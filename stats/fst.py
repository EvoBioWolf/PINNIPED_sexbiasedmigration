import numpy as np
import pandas as pd
import allel
from pathlib import Path

def fst_permutation_test(
    vcf_path,
    samples_csv,
    pop_col="Population",
    sex_col="Sex",
    id_col="SampleID",
    method="hudson", # "hudson" or "wc84"
    biallelic_only=False,
    min_mac_all=1,
    clip_negative=False, # clip negatives to 0
    n_perm=1000,
    stratify_by_sex=True,
    random_state=42
):
    rng = np.random.default_rng(random_state)

    # read vcf
    callset = allel.read_vcf(
        vcf_path,
        fields=["samples", "calldata/GT", "variants/ALT"]
    )
    samples = callset["samples"].tolist()
    gt = allel.GenotypeArray(callset["calldata/GT"]) # (n_var, n_samp, ploidy)

    # biallelic + global MAC filters (to make sure there are none left)
    keep = np.ones(gt.shape[0], dtype=bool)
    if biallelic_only and "variants/ALT" in callset:
        keep &= ((callset["variants/ALT"] != '').sum(axis=1) == 1)

    ac_all = gt.count_alleles(max_allele=1) # (n_var, 2)
    mac = np.minimum(ac_all[:, 0], ac_all[:, 1])
    keep &= (mac >= min_mac_all)

    gt = gt.compress(keep, axis=0)
    if gt.n_variants == 0:
        raise ValueError("No variants left after biallelic/MAC filtering.")

    # intersect with metadata
    mapping = pd.read_csv(samples_csv, sep=None, engine="python", dtype=str)
    if id_col is None:
        id_col = "SampleID" if "SampleID" in mapping.columns else mapping.columns[0]
    mapping = mapping.set_index(id_col).reindex(samples)
    if mapping.isnull().any().any():
        missing = mapping.index[mapping.isnull().any(axis=1)].tolist()
        raise ValueError(f"Missing metadata for samples (first 10): {missing[:10]}")

    pops  = mapping[pop_col].to_numpy()
    sexes = mapping[sex_col].str.upper().str[0].to_numpy() if sex_col in mapping.columns else np.array(["U"]*len(pops))
    levels = list(pd.unique(pops))
    k = len(levels)

    # compute ratio-of-sums Fst for a given label assignment
    def _pair_fst_for_labels(labels):
        mat = np.full((k, k), np.nan, float)
        for a in range(k):
            for b in range(a+1, k):
                idx1 = np.where(labels == levels[a])[0]
                idx2 = np.where(labels == levels[b])[0]

                if method.lower() == "hudson":
                    ac1 = gt.count_alleles(subpop=idx1, max_allele=1)
                    ac2 = gt.count_alleles(subpop=idx2, max_allele=1)
                    # drop sites with zero callable chromosomes in either pop
                    ok = (ac1.sum(axis=1) > 0) & (ac2.sum(axis=1) > 0)
                    if not np.any(ok):
                        continue
                    num, den = allel.hudson_fst(ac1.compress(ok,0), ac2.compress(ok,0))
                    den_sum = np.nansum(den)
                    if den_sum == 0:
                        fst = np.nan
                    else:
                        fst = np.nansum(num) / den_sum
                        if clip_negative and np.isfinite(fst) and fst < 0:
                            fst = 0.0

                elif method.lower() == "wc84":
                    g = gt.to_n_alt()
                    g1, g2 = g[:, idx1], g[:, idx2]
                    ok = np.isfinite(g1).any(axis=1) & np.isfinite(g2).any(axis=1)
                    if not np.any(ok):
                        continue
                    a_, b_, c_ = allel.weir_cockerham_fst(g1[ok], g2[ok])
                    denom = a_ + b_ + c_
                    w = np.isfinite(denom) & (denom > 0)
                    if not np.any(w):
                        fst = np.nan
                    else:
                        fst = np.nansum(a_[w]) / np.nansum(denom[w])
                        if clip_negative and np.isfinite(fst) and fst < 0:
                            fst = 0.0
                else:
                    raise ValueError("method must be 'hudson' or 'wc84'.")

                mat[a, b] = mat[b, a] = fst
        np.fill_diagonal(mat, 0.0)
        return mat

    # Fst matrix
    obs_mat = _pair_fst_for_labels(pops.copy())

    # build permuted label generator, preserving pop counts
    def permute_labels():
        lab = pops.copy()
        if stratify_by_sex and np.unique(sexes).size > 1:
            for s in np.unique(sexes):
                idx = np.where(sexes == s)[0]
                lab[idx] = rng.permutation(lab[idx])
        else:
            lab = rng.permutation(lab)
        return lab

    # run permutations: for each pair return right-tailed p-value
    pval = np.full_like(obs_mat, np.nan, dtype=float)
    for a in range(k):
        for b in range(a+1, k):
            obs = obs_mat[a, b]
            if not np.isfinite(obs):
                continue
            ge = 0
            total = 0
            for _ in range(n_perm):
                lab = permute_labels()
                mat_perm = _pair_fst_for_labels(lab)
                val = mat_perm[a, b]
                if np.isfinite(val):
                    total += 1
                    if val >= obs:
                        ge += 1
            if total == 0:
                p = np.nan
            else:
                # add-one smoothing
                p = (ge + 1) / (total + 1)
            pval[a, b] = pval[b, a] = p

    return pd.DataFrame(obs_mat, index=levels, columns=levels), pd.DataFrame(pval, index=levels, columns=levels)



n_perm=1000
basepath=""
samples_csv=f"{basepath}/lists/sexlist_GSLSSL_nonrelatives.csv"

vcf_gsl_mt=f"{basepath}/VARCALL/vcf/GSL_MF_CM019821.1/GSL_MF_CM019821.1.MQ30.BQ30.snps.Q40_DP25000_MQ40_NA10_het25_bias25_PI_57pos.woEC01.vcf"
vcf_gsl_x=f"{basepath}/VARCALL/vcf/GSL_MF_CM019819.2/GSL_MF_CM019819.2.MQ30.BQ30.snps._hom._repOtZC._par.Q40_DP10000_MQ40_NA10_bias50_PI_30928pos.woEC01.vcf"
vcf_gsl_aut=f"{basepath}/VARCALL/vcf/GSL_MF_CM019810.2/GSL_MF_CM019810.2.MQ30.BQ30.snps._hom._repOtZC.Q40_DP10000_MQ40_NA10_bias50_PI_112730pos.woEC01.vcf"
vcf_gsl_y=f"{basepath}/VARCALL/vcf/GSL_M_CM019820.2/GSL_M_CM019820.2.MQ30.BQ30.snps.rep_Otariidae_ZalCal.fem90.Q60_DP1000mean25_MQ60_NA10_het50_bias50_PI_50pos.vcf"

vcf_ssl_mt=f"{basepath}/VARCALL/vcf/SSL_CM019821.1/SSL_CM019821.1.MQ30.BQ30.snps._hom._repOtZC.Q40_DP10000_NA0.1_bias0.25_het0.25_PI_50pos.vcf"
vcf_ssl_x=f"{basepath}/VARCALL/vcf/SSL_CM019819.2/SSL_CM019819.2.MQ30.BQ30.snps._hom._repOtZC._par.Q40_DP10000mean0.5_NA0.1_bias0.5_het0.5_PI_27459pos.vcf"
vcf_ssl_aut=f"{basepath}/VARCALL/vcf/SSL_CM019810.2/SSL_CM019810.2.MQ30.BQ30.snps._hom._repOtZC.Q40_DP10000mean0.5_NA0.1_bias0.5_het0.5_PI_116707pos.vcf"
vcf_ssl_y=f"{basepath}/VARCALL/vcf/SSL_CM019820.2/SSL_CM019820.2.MQ30.BQ30.snps._hom._repOtZC.Q60_DP300mean0.25_NA0.1_bias0.5_het0.5_PI_49pos.vcf"

for vcf in [vcf_ssl_aut]: #vcf_gsl_mt, vcf_gsl_x, vcf_gsl_aut, vcf_gsl_y, vcf_ssl_mt, vcf_ssl_x, vcf_ssl_aut, vcf_ssl_y
    pop_col="Population" #"Population" "Stock_K3" "Stock_K2"
    suffix=f"{Path(vcf).parent.name}_{pop_col}"
    print(suffix)
    obs, p = fst_permutation_test(
        vcf,
        samples_csv,
        pop_col=pop_col,
        sex_col="Sex",
        method="hudson",
        n_perm=n_perm,
        stratify_by_sex=False
    )
    outpath=f"{basepath}/Fst"
    obs.to_csv(f"{outpath}/fst_{suffix}.csv")
    p.to_csv(f"{outpath}/pval_fst_perm_{n_perm}_{suffix}.csv")
    print("Observed FST:\n", obs.round(4))
    print("P-values (right-tailed):\n", p.round(4))
