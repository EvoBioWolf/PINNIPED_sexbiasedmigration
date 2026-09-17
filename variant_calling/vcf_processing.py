import pandas as pd
import pysam
import tempfile
from pathlib import Path
import subprocess
import os

def open_vcf(vcf_instance):
    """ Opens a VCF file if given a path or returns it if already a pysam.VariantFile. """
    if isinstance(vcf_instance, pysam.VariantFile):
        return vcf_instance
    elif isinstance(vcf_instance, str) and Path(vcf_instance).exists():
        return pysam.VariantFile(vcf_instance)
    else:
        raise FileNotFoundError(f"VCF file not found or invalid instance: {vcf_instance}")

def report_saved_massarray_positions(vcf_records, massarray_positions):
    """
    Returns a set of positions that are missing from the massarray file.

    Args:
        vcf_records (list of VariantRecord): variant records.
        massarray_positions (list): The list with informative MassArray positions.

    Returns:
        set: A set of positions that are missing from the massarray file.
    """
    # massarray_positions = pd.read_csv(massarray_file, header=None, sep="\t", names=["Chromosome","Position"])
    # massarray_positions_set = set(massarray_positions["Position"])
    massarray_positions_set = set(massarray_positions)

    vcf_positions_set = set([record.pos for record in vcf_records])
    missing_positions = massarray_positions_set - vcf_positions_set

    return len(missing_positions)

# Filter positions based on numeric parameters (e.g., QUAL, DP, MQ)
def filter_positions_by_parameters(vcf_in, massarray_positions, **thresholds):
    """
    Filters positions for which the QUAL and INFO fields satisfy thresholds.
    Expects thresholds like QUAL=30, DP=10, MQ=20.
    Returns a list of filtered records and a stats dict.

    Args:
        vcf_in (list or pysam.VariantFile): A list of pysam.VariantRecord objects 
                                            or a path to a VCF file.
        massarray_positions (list): The list with informative MassArray positions.
        **thresholds (dict): A dictionary where keys are field names and values are 
                             threshold values to filter the VCF by.
                             For example: QUAL=30, DP=10

    Returns:
        list: A list of filtered pysam.VariantRecord objects.
        dict: A dictionary containing the number of total positions, 
              the number of positions kept after filtering, 
              the number of positions dropped, and the number of lost massarray positions.
    """
    filtered_records = []
    total = 0
    for record in vcf_in:
        total += 1
        include = True
        for field, thresh in thresholds.items():
            if field.upper() == "QUAL":
                if record.qual is None or record.qual < thresh:
                    include = False
                    break
            elif field.upper() == "DP":
                field_val = record.info.get(field, 0)
                if field_val is None or field_val > thresh:
                    include = False
                    break
            else:
                field_val = record.info.get(field, 0)
                if field_val is None or field_val < thresh:
                    include = False
                    break
        if include:
            filtered_records.append(record)
    stats = {"total": total, 
             "kept": len(filtered_records), 
             "dropped": total - len(filtered_records),
             "lost massarray positions": report_saved_massarray_positions(filtered_records, massarray_positions)}
    return filtered_records, stats

def filter_by_mean_depth(records, massarray_positions, mean_depth_distribution_file, mean_depth_threshold=0.5):
    """
    Filters positions based on mean depth.
    Returns a list of filtered records and a stats dict.

    Args:
        records (list): A list of pysam.VariantRecord objects.
        massarray_positions (list): The list with informative MassArray positions.
        mean_depth_distribution_file (str): The path to the file with mean depth distribution.
        mean_depth_threshold (float): The threshold for the proportion of high depth positions.

    Returns:
        list: A list of filtered pysam.VariantRecord objects.
        dict: A dictionary containing the number of total positions, 
              the number of positions kept after filtering, 
              the number of positions dropped, and the number of lost massarray positions.
    """

    # Read mean depth distribution file and create a dictionary
    def get_mean_depth(mean_depth_distribution_file):
        df = pd.read_csv(mean_depth_distribution_file,
                         header=0,
                         sep="\t",
                         usecols=["Sample", "85%"]).rename(columns={"Sample": "sample", "85%": "threshold"})
        return df.set_index("sample")["threshold"].to_dict()
    
    filtered_records = []
    total = 0
    kept = 0
    mean_dp = get_mean_depth(mean_depth_distribution_file)

    for record in records:
        total += 1
        dp_cnt = 0
        hdp_cnt = 0
        for sample in record.samples:
            dp = record.samples[sample].get("DP")
            if dp is not None and dp > 0:
                dp_cnt += 1
                if dp > mean_dp[sample]:
                    hdp_cnt += 1
        if hdp_cnt/dp_cnt < mean_depth_threshold:
            filtered_records.append(record)
            kept += 1
    stats = {"total": total, 
             "kept": kept, 
             "dropped": total - kept,
             "lost massarray positions": report_saved_massarray_positions(filtered_records, massarray_positions)}
    return filtered_records, stats


# Remove positions if a high proportion of genotypes are missing
def remove_high_missing_genotypes(records, massarray_positions, missing_threshold=0.5):
    """
    Removes positions (records) in which the proportion of missing genotypes is >= missing_threshold.
    Assumes that for a record, missing genotype is defined as GT containing None.
    Returns the filtered list and statistics.

    Args:
        records (list): A list of pysam.VariantRecord objects.
        massarray_positions (list): The list file with informative MassArray positions.    
        missing_threshold (float): The threshold for the proportion of missing genotypes.

    Returns:
        list: A list of filtered pysam.VariantRecord objects.
        dict: A dictionary containing the number of total positions, 
              the number of positions kept after filtering, 
              the number of positions dropped, and the number of lost massarray positions.
    """
    filtered_records = []
    total = 0
    kept = 0
    for record in records:
        total += 1
        samples = record.samples.keys()
        missing_count = sum(1 for sample in samples if record.samples[sample].get("GT") is None or None in record.samples[sample]["GT"])
        proportion_missing = missing_count / len(samples) if samples else 0
        if proportion_missing < missing_threshold:
            filtered_records.append(record)
            kept += 1
    stats = {"total": total, 
             "kept": kept, 
             "dropped": total - kept,
             "lost massarray positions": report_saved_massarray_positions(filtered_records, massarray_positions)}
    return filtered_records, stats

# Remove heterozygous positions
def remove_heterozygous_positions(records, massarray_positions, het_threshold=0.5):
    """
    Removes positions that are heterozygous in any sample.
    Returns filtered list and stats.

    Args:
        records (list): A list of pysam.VariantRecord objects.
        massarray_positions (list): The list with informative MassArray positions.    
        het_threshold (float): The threshold for the proportion of heterozygous genotypes.

    Returns:
        list: A list of filtered pysam.VariantRecord objects.
        dict: A dictionary containing the number of total positions, 
              the number of positions kept after filtering, 
              the number of positions dropped, and the number of lost massarray positions.
    """
    filtered_records = []
    total = 0
    kept = 0
    for record in records:
        total += 1
        het_cnt = 0
        for sample in record.samples:
            gt = record.samples[sample].get("GT")
            ad = record.samples[sample].get("AD")
            if gt is not None and gt[0] == 2 or ad is not None and ad[0] > 0 and ad[1] > 0:
                het_cnt += 1
        if het_cnt/len(record.samples) < het_threshold and len(record.alts) == 1:
            filtered_records.append(record)
            kept += 1
    stats = {"total": total, 
             "kept": kept, 
             "dropped": total - kept,
             "lost massarray positions": report_saved_massarray_positions(filtered_records, massarray_positions)}
    return filtered_records, stats

# Remove positions with high strand bias
def remove_strand_bias(records, massarray_positions, bias_threshold=0.1):
    """
    Removes positions with strand bias above the specified threshold.
    This example assumes that the INFO field 'SP' (strand bias) is numeric.
    Returns filtered records and statistics.

    Args:
        records (list): A list of pysam.VariantRecord objects.
        massarray_positions (list): The list with informative MassArray positions.    
        bias_threshold (float): The threshold for strand bias.

    Returns:
        list: A list of filtered pysam.VariantRecord objects.
        dict: A dictionary containing the number of total positions, 
              the number of positions kept after filtering, 
              the number of positions dropped, and the number of lost massarray positions.
    """

    # Returns True if for EITHER the reference allele (index 0)
    # or the alternative allele (index 1) the counts show strong strand bias:
    # i.e. one strand has 0 and the other has >1.
    def has_strand_bias(ADF, ADR):
        ref_bias = (ADF[0] == 0 and ADR[0] > 0) or (ADR[0] == 0 and ADF[0] > 0)
        if len(ADF) == 1 and len(ADR) == 1: # in case of invariant sites
            return ref_bias
        alt_bias = (ADF[1] == 0 and ADR[1] > 0) or (ADR[1] == 0 and ADF[1] > 0)
        return ref_bias or alt_bias

    filtered_records = []
    total = 0
    kept = 0
    for record in records:
        total += 1
        bias_rate = 0
        for sample in record.samples:
            adf = record.samples[sample].get("ADF")
            adr = record.samples[sample].get("ADR")
            if adf is not None and adr is not None and has_strand_bias(adf, adr):
                bias_rate += 1
        if bias_rate/len(record.samples) < bias_threshold:
            filtered_records.append(record)
            kept += 1
    stats = {"total": total, 
             "kept": kept, 
             "dropped": total - kept,
             "lost massarray positions": report_saved_massarray_positions(filtered_records, massarray_positions)}
    return filtered_records, stats


def impute_missing_genotypes(records,
                             populations_file,
                             genetic_identity,
                             impute_by,
                             impute_by_KNN=5):
    """
    Impute missing (haploid) GTs in-place for a list of VariantRecord objects.
    1) Try each subgroup in order, for example list Population, Group, Species.
    2) If still missing, try closest samples by KNN.
    Missing = GT is None, or GT[0] is None or ".".
    
    Returns:
        records  – list of VariantRecord (modified in-place)
        stats    – dict of counts: total, kept, dropped
    """
    # Load once
    pop_df = pd.read_csv(populations_file, sep="\t").set_index("SampleID") if os.path.isfile(populations_file) else None
    identity_df = pd.read_csv(genetic_identity, index_col=0, sep=",") if os.path.isfile(genetic_identity) else None
    
    # Helpers
    def is_missing(gt):
        return gt is None or gt[0] is None or (isinstance(gt[0], str) and gt[0] == ".")
    
    def most_common(gt_list):
        if len(set(gt_list)) == 1:
            return gt_list[0]
        return None
    
    def subgroup_samples(sample, level):
        val = pop_df.at[sample, level]
        # only male, not self
        sel = pop_df[
            (pop_df[level] == val) &
            (pop_df.Sex == "m")
        ].index
        return [s for s in sel if s != sample]

    s = {
        "total_missing": 0,
        "imputed_by_subgroup": 0,
        "imputed_by_identity": 0,
        "still_missing": 0
    }
    
    for record in records:
        for sample in record.samples:
            gt = record.samples[sample]["GT"]
            if not is_missing(gt):
                continue
            
            s["total_missing"] += 1
            imputed = None
            
            # 1) by subgroup(s)
            if pop_df is not None and impute_by != []:
                for lvl in impute_by:
                    gts = [
                        record.samples[s]["GT"]
                        for s in subgroup_samples(sample, lvl)
                        if not is_missing(record.samples[s]["GT"])
                    ]
                    imputed = most_common(gts)
                    if imputed is not None:
                        s["imputed_by_subgroup"] += 1
                        break
            
            # 2) by genetic distance
            if imputed is None and identity_df is not None:
                row = identity_df.loc[sample].drop(sample)
                candidates = row[row >= 85]
                neighbours = candidates.nlargest(impute_by_KNN).index
                gts = [
                    record.samples[s]["GT"]
                    for s in neighbours
                    if not is_missing(record.samples[s]["GT"])
                ]
                imputed = most_common(gts)
                if imputed is not None:
                    s["imputed_by_identity"] += 1
            
            # 3) assign or leave missing
            if imputed is not None:
                record.samples[sample]["GT"] = imputed
            else:
                s["still_missing"] += 1
        stats = {
            "total": f"total missing: {s['total_missing']}",
            "kept": f"still missing: {s['still_missing']}",
            "dropped": "; ".join([f"imputed by {impute_by}: {s["imputed_by_subgroup"]}" if impute_by != [] else "",
                                  f"imputed by KNN: {s["imputed_by_identity"]}" if identity_df is not None else ""]),
            "lost massarray positions": ""}
    
    return records, stats
                

# Remove non-parsimony informative positions
def remove_non_parsimony_informative_positions(records, massarray_positions):
    """
    Removes positions that are not parsimony informative.
    For example, for a diploid VCF, a site is parsimony informative 
    if at least two different alleles appear in at least two individuals each.
    Returns filtered records and stats.

    Args:
        records (list): A list of pysam.VariantRecord objects.
        massarray_positions (list): The list with informative MassArray positions.    

    Returns:
        list: A list of filtered pysam.VariantRecord objects.
        dict: A dictionary containing the number of total positions, 
              the number of positions kept after filtering, 
              the number of positions dropped, and the number of lost massarray positions.
    """
    filtered_records = []
    total = 0
    kept = 0
    for record in records:
        total += 1
        allele_counts = {}
        for sample in record.samples:
            gt = record.samples[sample].get("GT")
            # Convert None genotype to placeholder
            if gt is not None:
                for allele in gt:
                    allele_counts[allele] = allele_counts.get(allele, 0) + 1
        # Remove missing allele counts (None)
        allele_counts.pop(None, None)
        # Check parsimony informativeness:
        # A variant is parsimony informative if at least two alleles appear at least twice.
        num_informative = sum(1 for count in allele_counts.values() if count >= 2)
        if num_informative >= 2:
            filtered_records.append(record)
            kept += 1
    stats = {"total": total, 
             "kept": kept, 
             "dropped": total - kept,
             "lost massarray positions": report_saved_massarray_positions(filtered_records, massarray_positions)}
    return filtered_records, stats


def remove_close_positions(records, massarray_positions, distance_threshold=100):
    """
    Filters out positions that are too close to each other by grouping consecutive records 
    (based on a distance threshold) and keeping only the record with the highest QUAL from each group.
    
    Missing QUAL (None) is treated as 0.
    
    Args:
        records (list): List of pysam.VariantRecord objects.
        massarray_positions (list): The list with informative MassArray positions.
        distance_threshold (int): Maximum allowed gap (in bases) between consecutive records in the same group.
                                  If the gap between the current record and the previous record is less than this 
                                  threshold, they are considered part of the same group.
    
    Returns:
        list: A list of filtered pysam.VariantRecord objects.
        dict: A dictionary containing the number of total positions, 
              the number of positions kept after filtering, 
              the number of positions dropped, and the number of lost massarray positions.
    """
    if not records:
        return []
    # Sort records by position.
    sorted_records = sorted(records, key=lambda r: r.pos)
    
    # Initialize the best record of the current group with the first record.
    filtered_records = []
    total = 0
    best_record = sorted_records[0]
    prev_record = sorted_records[0]
    
    for record in sorted_records[1:]:
        total += 1
        # If the current record is close to the previous one, it is part of the same group.
        if record.pos - prev_record.pos < distance_threshold:
            rec_qual = record.qual if record.qual is not None else 0
            best_qual = best_record.qual if best_record.qual is not None else 0
            if rec_qual > best_qual:
                best_record = record
        else:
            filtered_records.append(best_record)
            best_record = record
        prev_record = record
    
    # Append the best record from the last group.
    filtered_records.append(best_record)    
    kept = len(filtered_records)
    stats = {"total": total, 
            "kept": kept, 
            "dropped": total - kept,
            "lost massarray positions": report_saved_massarray_positions(filtered_records, massarray_positions)}
    return filtered_records, stats

def remove_positions_vcf(records, massarray_positions, positions_to_remove):
    pos_set = set(positions_to_remove)
    filtered_records = [rec for rec in records if rec.pos not in pos_set]
    total = len(records)
    kept = len(filtered_records)
    stats = {"total": total, 
            "kept": kept, 
            "dropped": total - kept,
            "lost massarray positions": report_saved_massarray_positions(filtered_records, massarray_positions)}
    return filtered_records, stats
    

def remove_non_biallelic(records, massarray_positions):
    filtered_records = [rec for rec in records if len(rec.alts) == 1]
    total = len(records)
    kept = len(filtered_records)
    stats = {"total": total, 
            "kept": kept, 
            "dropped": total - kept,
            "lost massarray positions": report_saved_massarray_positions(filtered_records, massarray_positions)}
    return filtered_records, stats

# Pipeline function to run steps based on user selection
def run_pipeline(vcf_input, steps_config, massarray_positions, mean_depth_distribution_file="", 
                 populations_file="", positions_to_remove=[], genetic_identity="", **kwargs):
    """
    Runs a series of filtering steps on a VCF.
    
    Args:
      vcf_input: Path to VCF file or pysam.VariantFile.
      steps_config (dict): Keys are step names, and the value is True/False indicating whether to run.
          Expected keys: "filter_parameters", "remove_high_missing", "remove_het", "remove_strand_bias", "remove_non_parsimony".
      massarray_positions (list): The list with informative MassArray positions.
      kwargs: Additional thresholds for each step (e.g., QUAL=30, DP=10, missing_threshold, bias_threshold).
    
    Returns:
      filtered_records: Final filtered list of records.
      stats: A dictionary with statistics for each step.
    """
    vcf_obj = open_vcf(vcf_input)
    records = list(vcf_obj)
    
    stats = pd.DataFrame(columns=["thresholds", "total", "kept", "dropped", "lost massarray positions"])

    # Remove positions from the list
    if steps_config.get("remove_positions", False):
        records, s = remove_positions_vcf(records, massarray_positions, positions_to_remove)
        s["thresholds"] = ""
        stats.loc["Remove positions from the list"] = s
    
    # Remove non-biallelic positions
    if steps_config.get("biallelic_only", False):
        records, s = remove_non_biallelic(records, massarray_positions)
        s["thresholds"] = ""
        stats.loc["Biallelic filter"] = s
    
    # Filter positions by quality thresholds
    if steps_config.get("filter_parameters", False):
        records, s = filter_positions_by_parameters(records, massarray_positions,
                                                    QUAL=kwargs.get("QUAL", 30), 
                                                    DP=kwargs.get("DP", 5000), 
                                                    MQ=kwargs.get("MQ", 20))
        s["thresholds"] = f"QUAL>{kwargs.get('QUAL', 30)}, DP<{kwargs.get('DP', 5000)}, MQ>{kwargs.get('MQ', 20)}"
        stats.loc["Quality filter"] = s

    # Filter positions by mean depth threshold for each sample
    if steps_config.get("filter_by_mean_depth", False):
        if mean_depth_distribution_file != '':
            records, s = filter_by_mean_depth(records, 
                                              massarray_positions,
                                              mean_depth_distribution_file, 
                                              mean_depth_threshold=kwargs.get("mean_depth_threshold", 0.5))
            s["thresholds"] = kwargs.get('mean_depth_threshold', 0.5)
            stats.loc["Remove positions with > 85% percentile sample depth"] = s

    # Remove positions with high frequency of missing genotypes
    if steps_config.get("remove_high_missing", False):
        records, s = remove_high_missing_genotypes(records, 
                                                   massarray_positions, 
                                                   missing_threshold=kwargs.get("missing_threshold", 0.5))
        s["thresholds"] = kwargs.get('missing_threshold', 0.5)
        stats.loc["Remove high missing positions"] = s
    
    # Remove heterozygous positions
    if steps_config.get("remove_het", False):
        records, s = remove_heterozygous_positions(records, 
                                                   massarray_positions, 
                                                   het_threshold=kwargs.get("het_threshold", 0.5))
        s["thresholds"] = kwargs.get('het_threshold', 0.5)
        stats.loc["Remove heterozygous positions"] = s
    
    # Remove strand bias
    if steps_config.get("remove_strand_bias", False):
        records, s = remove_strand_bias(records, 
                                        massarray_positions, 
                                        bias_threshold=kwargs.get("bias_threshold", 0.1))
        s["thresholds"] = kwargs.get('bias_threshold', 0.1)
        stats.loc["Remove positions with strand bias"] = s

    # Impute missing genotypes
    if steps_config.get("impute", False): 
        impute_by=kwargs.get("impute_by", [])
        impute_by_KNN=kwargs.get("impute_by_KNN", 5)
        records, s = impute_missing_genotypes(records, 
                                              populations_file, 
                                              genetic_identity, 
                                              impute_by, 
                                              impute_by_KNN)
        s["thresholds"] = "; ".join([f"impute by {impute_by}" if impute_by != [] else "",
                                     f"KNN={impute_by_KNN}" if genetic_identity != "" else ""])
        stats.loc["Impute missing genotypes"] = s

    # Remove non-parsimony informative positions
    if steps_config.get("remove_non_parsimony", False):
        records, s = remove_non_parsimony_informative_positions(records, massarray_positions)
        s["thresholds"] = ""
        stats.loc["Remove non-parsimony informative positions"] = s

    # Remove close positions
    if steps_config.get("remove_close_variants", False):
        records, s = remove_close_positions(records, 
                                            massarray_positions, 
                                            distance_threshold=kwargs.get("distance_threshold", 100))
        s["thresholds"] = kwargs.get('distance_threshold', 100)
        stats.loc["Remove close positions"] = s

    return records, stats


def write_variant_records(variant_records, header_vcf_path, output_vcf_path):
    """
    Write a list of pysam.VariantRecord objects to a VCF file using the header
    from an existing VCF file.

    Args:
        variant_records (list): List of pysam.VariantRecord objects.
        header_vcf_path (str): Path to a VCF file from which to copy the header.
        output_vcf_path (str): Path where the output VCF file will be saved.
    """
    # Open the header VCF file to get its header.
    with pysam.VariantFile(header_vcf_path, "r") as header_vcf:
        header = header_vcf.header

    # Create output VCF file using this header.
    with pysam.VariantFile(output_vcf_path, "w", header=header) as vcf_out:
        for record in variant_records:
            vcf_out.write(record)


def variant_records_to_fasta(records, output_fasta):
    """
    Converts a list of pysam.VariantRecord objects to FASTA format.
    
    Each record describes several samples. For each sample, we construct a sequence 
    by concatenating alleles from each variant position. If a sample's genotype is 
    missing (None or "."), "N" is inserted.

    For haploid data:
      - If GT is missing (None or (None, ) or (".",)), "N" is used.
      - If GT[0] is 0, then the reference allele is used.
      - If GT[0] > 0, then the corresponding alternate allele is used 
        (i.e. record.alts[allele-1]).
        
    Args:
        records (list): List of pysam.VariantRecord objects.
        output_fasta (str): Path to the output FASTA file.
    """
    if not records:
        raise ValueError("No records provided.")

    # Get sample names from the header of the first record.
    sample_names = list(records[0].header.samples)
    
    # Dictionary to collect sequences per sample.
    sample_seqs = {sample: [] for sample in sample_names}
    
    for record in records:
        for sample in sample_names:
            # Extract the genotype tuple from the record for this sample.
            gt = record.samples[sample].get("GT")
            
            # For haploid data we expect GT to be a tuple of length 1.
            if gt is None or gt[0] is None or (isinstance(gt[0], str) and gt[0] == "."):
                # Missing genotype -> use "N"
                sample_seqs[sample].append("N")
            else:
                allele = gt[0]
                if allele == 0:
                    sample_seqs[sample].append(record.ref)
                elif allele > 0:
                    # Try to get the alternate allele; if it fails, fallback to "N"
                    try:
                        sample_seqs[sample].append(record.alts[allele - 1])
                    except (IndexError, TypeError):
                        sample_seqs[sample].append("N")
                else:
                    sample_seqs[sample].append("N")
    
    # Write the sequences to a FASTA file.
    with open(output_fasta, "w") as fasta:
        for sample, seq_list in sample_seqs.items():
            sequence = "".join(seq_list)
            fasta.write(f">{sample}\n{sequence}\n")


def get_vcf_positions(vcf_path):
    """
    Returns a list of all POS values from the VCF.
    """
    vcf = pysam.VariantFile(vcf_path)
    return [record.pos for record in vcf]

