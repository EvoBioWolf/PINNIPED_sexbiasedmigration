import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import rowcol
from pyproj import Transformer
from scipy import ndimage
from skimage.graph import MCP_Geometric
from pyproj import Geod

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

def plot_coastal_path_with_map(
    geoloc_file,
    bathy_tif,
    area,
    popA,
    popB,
    depth_thresholds=[-800.0, -1500.0, -3000.0],
    depth_costs=[1.0, 1.1, 1.2],
    additional_filename='',
    output_csv='',
    figsize=(5,5)
):
    """
    Get depth‐weighted coastal path between popA and popB;
    plot it on a Cartopy map (with coastlines).
    """

    df = pd.read_csv(geoloc_file)
    lonA, latA = df.loc[df.Population == popA, ["Lon","Lat"]].iloc[0]
    lonB, latB = df.loc[df.Population == popB, ["Lon","Lat"]].iloc[0]

    with rasterio.open(bathy_tif) as src:
        bathy = src.read(1)
        mask  = src.read_masks(1)
        transform = src.transform
        crs       = src.crs

    cost = np.full(bathy.shape, np.inf, dtype=float)
    shallow   = (bathy < 0) & (bathy >= depth_thresholds[0])
    deep      = (bathy < depth_thresholds[0]) & (bathy >= depth_thresholds[1])
    extradeep = (bathy < depth_thresholds[1]) & (bathy >= depth_thresholds[2])
    cost[shallow]   = depth_costs[0]
    cost[deep]      = depth_costs[1]
    cost[extradeep] = depth_costs[2]

    # nearest‐water snap indices
    land_mask = ~np.isfinite(cost)
    _, (ind_r, ind_c) = ndimage.distance_transform_edt(
        land_mask, return_distances=True, return_indices=True
    )

    dx = transform.a
    dy = -transform.e
    mcp = MCP_Geometric(cost, sampling=[dx, dy])

    # CRS transformers
    to_utm   = Transformer.from_crs("EPSG:4326", crs, always_xy=True)
    to_wgs84 = Transformer.from_crs(crs, "EPSG:4326", always_xy=True)

    def snap(lon, lat):
        x, y = to_utm.transform(lon, lat)
        r, c = rowcol(transform, x, y)
        if not np.isfinite(cost[r, c]):
            r, c = ind_r[r, c], ind_c[r, c]
        return (r, c)

    start = snap(lonA, latA)
    end   = snap(lonB, latB)

    # path
    _, _ = mcp.find_costs([start], [end])
    path = mcp.traceback(end)

    # path -> lon/lat
    rows, cols = zip(*path)
    xs, ys = rasterio.transform.xy(transform, rows, cols, offset="center")
    lons, lats = to_wgs84.transform(xs, ys)



    def path_length_km(lons, lats, ellps="WGS84"):
        """
        Get the total geodesic distance along the polyline in kilometers
        between the points in the list of longitudes and latitudes (in degrees).
        """
        geod = Geod(ellps=ellps)
        total_m = 0.0
        for lon1, lat1, lon2, lat2 in zip(lons[:-1], lats[:-1], lons[1:], lats[1:]):
            _, _, dist = geod.inv(lon1, lat1, lon2, lat2)
            total_m += dist
        return total_m / 1e3 # m -> km

    length_km = round(path_length_km(lons, lats),2)

    # # plot
    # fig = plt.figure(figsize=figsize)
    # ax = plt.axes(projection=ccrs.PlateCarree())

    # # coastlines
    # ax.add_feature(cfeature.COASTLINE.with_scale("10m"), linewidth=0.5)

    # # coastal path
    # ax.plot(lons, lats, color="red", linewidth=1, transform=ccrs.PlateCarree(),
    #         label=f"{popA}-{popB} coastal path ({length_km}km)")

    # # endpoints
    # ax.scatter([lonA, lonB], [latA, latB],
    #            color=("red","red"), s=20, zorder=5,
    #            transform=ccrs.PlateCarree())
    # ax.text(lonA, latA, popA, color="black", va="bottom",
    #         transform=ccrs.PlateCarree())
    # ax.text(lonB, latB, popB, color="black", va="bottom",
    #         transform=ccrs.PlateCarree())

    # ax.set_extent(area,
    #               crs=ccrs.PlateCarree())

    # ax.legend(loc="upper right")
    # plt.tight_layout()
    # plt.savefig(f"{output_csv}/{additional_filename}.{popA}_{popB}.{length_km}km.jpg")
    # plt.show()
    fig = plt.figure(figsize=figsize)
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))

    ax.add_feature(cfeature.COASTLINE.with_scale("10m"), linewidth=0.5)

    ax.plot(lons, lats, color="red", linewidth=1, transform=ccrs.PlateCarree(),
            label=f"{popA}-{popB} coastal path ({length_km}km)")

    ax.scatter([lonA, lonB], [latA, latB],
                color=("red","red"), s=20, zorder=5,
                transform=ccrs.PlateCarree())
    ax.text(lonA, latA, popA, color="black", va="bottom",
            transform=ccrs.PlateCarree())
    ax.text(lonB, latB, popB, color="black", va="bottom",
            transform=ccrs.PlateCarree())

    ax.set_extent(area, crs=ccrs.PlateCarree())

    ax.legend(loc="upper right")
    plt.tight_layout()
    plt.savefig(f"{output_csv}/w.{popA}_{popB}.{length_km}km.jpg")
    plt.show()

    return length_km

def count_distances(pops, geoloc_file, bathy_tif, area, additional_filename="", output_csv="",figsize=(5,5)):
    dist_mat = pd.DataFrame(np.zeros((len(pops), len(pops))), index=pops, columns=pops)
    for i, popA in enumerate(pops):
        for popB in pops[i+1:]:
            length_km = plot_coastal_path_with_map(
                geoloc_file=geoloc_file,
                bathy_tif=bathy_tif,
                area=area,
                popA=popA,
                popB=popB,
                additional_filename=additional_filename,
                output_csv=output_csv,
                figsize=figsize
            )
            dist_mat.loc[popA, popB] = length_km
            dist_mat.loc[popB, popA] = length_km
    dist_mat.to_csv(output_csv+f"{additional_filename}_pairwise_depth_weighted_distances_km.csv") 

    return dist_mat

def run_experiment(geoloc_path, bathy_tif, area, additional_filename, output_csv, figsize=(5,5)):
    geoloc_file=f"{geoloc_path}{additional_filename}_geolocations.csv"
    pops=pd.read_csv(geoloc_file, index_col=0).index.to_list()
    count_distances(pops, geoloc_file, bathy_tif, area, additional_filename=additional_filename, output_csv=output_csv,figsize=figsize)


basepath=""
area_csl=[-123, -108, 20, 38]# CSL
area_gsl=[-92, -89, -1.5, 0.9]# Galapagos
area_ssl=[125, 240, 35, 65]# SSL

bathy_tif_csl = f"{basepath}/Mantel_test/tifs/csl_bathy_utm.tif"
bathy_tif_gsl = f"{basepath}/Mantel_test/tifs/gsl_bathy_utm.tif"
bathy_tif_ssl = f"{basepath}/Mantel_test/tifs/ssl_bathy_utm.tif"

output_csv=f"{basepath}/Mantel_test/distances/"
geoloc_path=f"{basepath}/Mantel_test/coordinates/"


run_experiment(geoloc_path, bathy_tif_csl, area_csl, "csl_pops", output_csv)
run_experiment(geoloc_path, bathy_tif_gsl, area_gsl, "gsl_pops", output_csv)
run_experiment(geoloc_path, bathy_tif_ssl, area_ssl, "ssl_pops", output_csv, figsize=(7,2))

