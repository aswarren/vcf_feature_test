import pysam  # PyVCF library
import pandas as pd
from collections import defaultdict
from scipy.stats import fisher_exact
from statsmodels.sandbox.stats.multicomp import multipletests
import argparse
import os


#python vcf_enrichment.py --vcf_dir vcf_data_dir --output_prefix my_analysis --feature_names variant_type genes

def extract_features_from_vcf(vcf_file, vcf_type, group_id, min_qual=None, max_fs=None, min_dp=None):
    """
    Extracts position-ALT mutation features from a VCF file using pysam,
    with filtering based on quality metrics.

    Args:
        vcf_file (str): Path to the VCF file.
        vcf_type (str): Type of VCF (e.g., 'human', 'virus1', 'tumor').
        group_id (str): Group identifier for the samples in this VCF.
        min_qual (float, optional): Minimum QUAL score to pass filter. Defaults to None (no QUAL filter).
        max_fs (float, optional): Maximum FS (Fisher Strand bias) to pass filter. Defaults to None (no FS filter).
        min_dp (int, optional): Minimum DP (Depth) to pass filter. Defaults to None (no DP filter).

    Returns:
        pandas.DataFrame: DataFrame with position-ALT mutation features.
    """
    try:
        vcf_reader = pysam.VariantFile(vcf_file)
    except pysam.SamtoolsError as e:
        print(f"Error opening VCF file {vcf_file} with pysam: {e}")
        return pd.DataFrame()

    features = []
    for record in vcf_reader.fetch():
        sample_id = "NA"
        if record.samples:
            sample_id = list(record.samples.keys())[0] if record.samples.keys() else "NA"

        # --- Filtering based on quality metrics ---
        if min_qual is not None and (record.qual is None or record.qual < min_qual):
            continue # Skip variant if QUAL is below threshold
        if max_fs is not None and 'FS' in record.info and record.info['FS'] > max_fs:
            continue # Skip if FS is above threshold
        if min_dp is not None and 'DP' in record.info and record.info['DP'] < min_dp:
            continue # Skip if DP is below threshold
        if record.filter.keys(): # Skip if variant is FILTERed (i.e., not PASS)
            continue

        # --- Create position-ALT mutation feature ---
        chrom = record.chrom
        pos = record.pos
        ref = record.ref
        alt_alleles = ",".join([str(alt) for alt in record.alts]) if record.alts else "." # Handle no ALT allele (e.g., for reference records, though unlikely in your VCF snippet)
        mutation_feature = f"{chrom}_{pos}_{ref}->{alt_alleles}"

        features.append({
            'sample_id': sample_id,
            'group': group_id,
            'vcf_type': vcf_type,
            'mutation': mutation_feature, # Position-ALT mutation as the feature
            # You can still add other INFO field values if you want to analyze them alongside mutations
            'qual': record.qual,
            'depth': record.info.get('DP', "NA"),
            'fs': record.info.get('FS', "NA"),
            'position': f"{chrom}_{pos}"
        })
    vcf_reader.close()
    return pd.DataFrame(features)


def analyze_group_enrichment(feature_df, feature_name='mutation', group_name='group'): # Feature name is now fixed to 'mutation'
    """
    Analyzes enrichment of position-ALT mutations across groups.

    Args:
        feature_df (pandas.DataFrame): DataFrame containing mutation features (output of extract_features_from_vcf_pysam_position_alt).
        feature_name (str):  Fixed to 'mutation' in this version.

    Returns:
        pandas.DataFrame: DataFrame with enrichment analysis results for mutations.
    """
    enrichment_results = []
    groups = sorted(feature_df[group_name].unique())
    all_vcf_types = sorted(feature_df['vcf_type'].unique())

    #make sure feature_name column is set to string
    feature_df[feature_name] = feature_df[feature_name].astype(str)

    for vcf_type in all_vcf_types:
        vcf_type_df = feature_df[feature_df['vcf_type'] == vcf_type]
        if vcf_type_df.empty:
            continue

        mutation_values = sorted(vcf_type_df[feature_name].unique()) # Unique mutations
        if not mutation_values:
            continue

        for mutation in mutation_values:
            if not isinstance(mutation, str):
                continue

            for group in groups:
                group_counts = vcf_type_df[vcf_type_df[group_name] == group][feature_name].value_counts()
                mutation_count_in_group = group_counts.get(mutation, 0)
                total_in_group = len(vcf_type_df[vcf_type_df[group_name] == group])

                other_groups_df = vcf_type_df[vcf_type_df[group_name] != group]
                other_groups_counts = other_groups_df[feature_name].value_counts()
                mutation_count_in_others = other_groups_counts.get(mutation, 0)
                total_in_others = len(other_groups_df)

                observed = [[mutation_count_in_group, total_in_group - mutation_count_in_group], #count this mutation in this group, count all other mutations in this group
                            [mutation_count_in_others, total_in_others - mutation_count_in_others]] #count this mutation in other groups, count all other mutations in other groups

                if mutation_count_in_group == 0:
                    continue

                oddsratio, pvalue = fisher_exact(observed, alternative='greater')

                enrichment_results.append({
                    'vcf_type': vcf_type,
                    'mutation': mutation, # Mutation feature
                    'group': group,
                    'mutation_count_in_group': mutation_count_in_group,
                    'total_in_group': total_in_group,
                    'mutation_count_in_others': mutation_count_in_others,
                    'total_in_others': total_in_others,
                    'odds_ratio': oddsratio,
                    'p_value': pvalue
                })

    enrichment_df = pd.DataFrame(enrichment_results)

    if not enrichment_df.empty:
        p_values = enrichment_df['p_value']
        reject, pvals_corrected, _, _ = multipletests(p_values, method='fdr_bh', is_sorted=False)
        enrichment_df['corrected_p_value'] = pvals_corrected
        enrichment_df['significant'] = reject

    return enrichment_df


def main():
    parser = argparse.ArgumentParser(description="Analyze VCF data for group-level feature enrichment.")
    parser.add_argument("--vcf_dir", type=str, required=True, help="Directory containing VCF files, organized by group.")
    parser.add_argument("--output_prefix", type=str, default="enrichment_analysis", help="Prefix for output files.")
    parser.add_argument("--feature_names", nargs='+', default=['mutation'], help="Features to analyze (columns in feature DataFrame).") # Analyze variant_type and genes by default
    #binary parameter group_mode either True or False
    parser.add_argument("--progression_mode", type=str, choices=['progression1', 'progression2', 'group'], default='group', help="Mode of analysis: 'progression1', 'progression2', or 'group'.")

    args = parser.parse_args()

    all_features_df = pd.DataFrame() # Initialize empty DataFrame to store all features

    group_dirs = [d for d in os.listdir(args.vcf_dir) if os.path.isdir(os.path.join(args.vcf_dir, d)) and d.startswith("Group")] # Assuming group directories start with "Group"

    for group_dir_name in group_dirs:
        group_id = group_dir_name
        group_path = os.path.join(args.vcf_dir, group_dir_name)

        virus_paths={"EBV1": os.path.join(group_path, "EBV1"),
        "EBV2" : os.path.join(group_path, "EBV2"),
        "GK18" : os.path.join(group_path, "GK18")}
        vcf_files = {}
        for virus, virus_path in virus_paths.items():
            #get all vcf files in the virus path
            vcf_files[virus]=[os.path.join(virus_path, f) for f in os.listdir(virus_path) if f.endswith('.vcf')]

        for vcf_type, vcf_file_paths in vcf_files.items():
            for vcf_file_path in vcf_file_paths:
                if os.path.exists(vcf_file_path): # Check if file exists before processing
                    print(f"Processing {vcf_type} VCF in {group_id} from: {vcf_file_path}")
                    features_df = extract_features_from_vcf(vcf_file_path, vcf_type, group_id)
                    if not features_df.empty:
                        all_features_df = pd.concat([all_features_df, features_df], ignore_index=True)
                    else:
                        print(f"No features extracted from {vcf_type} VCF in {group_id}.")
                else:
                    print(f"Warning: {vcf_type} VCF file not found in {group_id}: {vcf_file_path}")
    #construct progression column Group1 and Group2 = 'early', Group3 and Group4 = 'late'
    if args.progression_mode == 'progression1':
        all_features_df['progression'] = all_features_df['group'].apply(lambda x: 'early' if x in ['Group_01', 'Group_02'] else 'late' if x in ['Group_03', 'Group_04'] else 'other')
    elif args.progression_mode == 'progression2':
        all_features_df['progression'] = all_features_df['group'].apply(lambda x: 'early' if x in ['Group_01', 'Group_02'] else 'late' if x in ['Group_03', 'Group_04','Group_05'] else 'other')
        group_name = 'progression'
    else:
        group_name = 'group'
    if not all_features_df.empty:
        print("\nPerforming Group Enrichment Analysis...")
        for feature_name in args.feature_names:
            print(f"Analyzing feature: {feature_name}")
            enrichment_results_df = analyze_group_enrichment(all_features_df, feature_name, group_name)
            if not enrichment_results_df.empty:
                output_file = f"{args.output_prefix}_{feature_name}_{args.progression_mode}_enrichment.tsv"
                enrichment_results_df.to_csv(output_file, sep='\t', index=False)
                print(f"Enrichment results for {feature_name} saved to: {output_file}")
            else:
                print(f"No enrichment results for feature: {feature_name}.")

    else:
        print("No features extracted from any VCF files. Analysis stopped.")


if __name__ == "__main__":
    main()