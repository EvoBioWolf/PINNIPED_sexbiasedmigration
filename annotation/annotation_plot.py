import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, Rectangle, FancyBboxPatch

def _draw_feature_arrows(ax, feats, y=0.0, height=0.6, head_frac=0.15, edgecolor="black"):
    # short ones appear on top
    feats_sorted = sorted(feats, key=lambda d: (d["end"] - d["start"]), reverse=True)
    for i, f in enumerate(feats_sorted):
        start, end, strand, color = f["start"], f["end"], f["strand"], f["color"]
        span = max(end - start, 1)
        
        head_len = max(5, int(span * head_frac))
        head_len = min(head_len, int(span * 0.45))

        if strand == 0 or span <= 2 * head_len:
            ax.add_patch(Rectangle((start, y - height/2), span, height,
                                   facecolor=color, edgecolor=edgecolor, lw=0.7, zorder=10+i))
        else:
            if strand > 0:
                x, dx = start, span
            else:
                x, dx = end, -span
            ax.add_patch(FancyArrow(
                x, y - height/2, dx, 0,
                width=height, head_width=height, head_length=head_len,
                length_includes_head=True,
                facecolor=color, edgecolor=edgecolor, lw=0.7, zorder=10+i
            ))

def _label_features(ax, feats, labels, sides, y_above=1.0, y_below=-1.0,
                    rotation=90, fontsize=10, x_offsets=None):
    xs = np.array([(f["start"] + f["end"]) * 0.5 for f in feats])
    if x_offsets:
        for idx, dx in x_offsets.items():
            if 0 <= idx < len(xs): xs[idx] += dx

    for i, (x, lab, side) in enumerate(zip(xs, labels, sides)):
        is_above = str(side).strip().lower() in {"above","top","up","a"}
        y = y_above if is_above else y_below
        va = "bottom" if is_above else "top"
        ax.text(x, y, lab, rotation=rotation, ha="center", va=va,
                fontsize=fontsize, clip_on=False, zorder=100)

def plot_genomic_tracks(
    tracks,
    x_max=None,
    xlabel="Position (Mb)",
    figsize=None,
    coverage_height=2.0,
    annotation_height=10.0,
    dpi=300
):
    spans = []
    for t in tracks:
        if t["type"] == "ticks":
            a = np.asarray(t["data"]);  spans.append((np.min(a), np.max(a))) if a.size else None
        elif t["type"] == "intervals":
            df = pd.read_csv(t["data"], sep="\t", header=None)
            if not df.empty:
                starts, ends = df.iloc[:,1].to_numpy(), df.iloc[:,2].to_numpy()
                if starts.size: spans.append((starts.min(), ends.max()))
        elif t["type"] == "coverage_line":
            df = pd.read_csv(t["data"], sep="\t", header=None)
            if not df.empty:
                starts, ends = df.iloc[:,1].to_numpy(), df.iloc[:,2].to_numpy()
                if starts.size: spans.append((starts.min(), ends.max()))
        elif t["type"] == "features":
            ends = [f["end"] for f in t["features"]]
            starts = [f["start"] for f in t["features"]]
            if starts: spans.append((min(starts), max(ends)))
        else:
            raise ValueError(f"Unknown track type {t['type']}")
    if x_max is None:
        x_max = max((b for _, b in spans), default=0)

    n = len(tracks)
    if figsize is None:
        figsize = (12, max(3, 1.0 * n))
    height_ratios = []
    for t in tracks:
        if t["type"] == "coverage_line":
            height_ratios.append(coverage_height)
        elif t["type"] == "features":
            height_ratios.append(annotation_height)
        else:
            height_ratios.append(1.0)

    fig, axes = plt.subplots(
        n, 1, sharex=True, figsize=figsize,
        gridspec_kw={"height_ratios": height_ratios},
        dpi=dpi
    )
    if n == 1: axes = [axes]

    # DRAW TRACKS
    for ax, t in zip(axes, tracks):
        ttype = t["type"]
        color = t.get("color", "k")
        label = t.get("label", "")

        if ttype == "ticks":
            pos = np.asarray(t["data"])
            if pos.size:
                ax.eventplot([pos], orientation="horizontal",
                             lineoffsets=0, linelengths=0.8, colors=color)
            ax.set_ylim(-1, 1)
            ax.set_ylabel(label, rotation=0, ha="right", va="center", labelpad=10, fontsize=12)
            ax.tick_params(axis="y", left=False, labelleft=False)

        elif ttype == "intervals":
            df = pd.read_csv(t["data"], sep="\t", header=None)
            if not df.empty:
                starts = df.iloc[:,1].to_numpy()
                ends   = df.iloc[:,2].to_numpy()
                widths = (ends - starts)
                if starts.size:
                    ax.broken_barh(list(zip(starts, widths)),
                                   (0.1, 0.8),
                                   facecolors=color, edgecolors="none", linewidth=0)
            ax.set_ylim(0, 1)
            ax.set_ylabel(label, rotation=0, ha="right", va="center", labelpad=10, fontsize=12)
            ax.tick_params(axis="y", left=False, labelleft=False)

        elif ttype == "coverage_line":
            df = pd.read_csv(t["data"], sep="\t", header=None)
            if not df.empty:
                mids = ((df.iloc[:,1] + df.iloc[:,2]) / 2).to_numpy()
                vals = df.iloc[:,3].to_numpy()
                ax.plot(mids, vals, color=color, lw=0.8, zorder=3)
                if "ylim" in t:
                    ax.set_ylim(*t["ylim"])
                elif vals.size:
                    ymin, ymax = float(vals.min()), float(vals.max())
                    pad = 0.05 * (ymax - ymin if ymax > ymin else (ymax or 1))
                    ax.set_ylim(max(0, ymin - pad), ymax + pad)
            ax.yaxis.tick_right()
            ax.set_ylabel(label or "Coverage", rotation=0, ha="right", va="center", labelpad=10, fontsize=12)
            ax.tick_params(axis="y", right=True, labelright=True, left=False, labelleft=False, colors=color)

        elif ttype == "features":
            feats   = t["features"]
            labels_ = t.get("labels", [""]*len(feats))
            sides   = t.get("sides",  ["above"]*len(feats))
            height  = t.get("height_px", None)
            head_fr = t.get("head_frac", 0.15)
            base_y  = 0.0
            bar_h   = t.get("bar_height", 0.6)
            if t.get("baseline", True):
                ax.hlines(-0.28, 0, x_max, color="black", lw=1, zorder=1)
            _draw_feature_arrows(ax, feats, y=base_y, height=bar_h, head_frac=head_fr, edgecolor="black")
            _label_features(ax, feats, labels_, sides,
                            y_above=1.0, y_below=-1.0,
                            rotation=t.get("label_rotation", 90),
                            fontsize=t.get("label_size", 10),
                            x_offsets=t.get("x_offsets", None))
            ax.set_ylim(-3, 3)
            ax.set_ylabel(label, rotation=0, ha="right", va="center", labelpad=10, fontsize=12)
            ax.tick_params(left=False, labelleft=False)

        else:
            raise ValueError(f"Unknown track type: {ttype}")

        for side in ("top","right","left"):
            ax.spines[side].set_visible(False)
        ax.spines["bottom"].set_visible(True)
        ax.spines["bottom"].set_alpha(0.3)
        ax.tick_params(axis="x", direction="out")

    axes[-1].set_xlim(0, x_max)
    axes[-1].set_xlabel(xlabel)
    xt = axes[-1].get_xticks()
    axes[-1].set_xticklabels([f"{x/1e6:.1f}" for x in xt])
    for ax in axes[:-1]:
        ax.tick_params(axis="x", labelbottom=False)

    plt.tight_layout()
    plt.show()
