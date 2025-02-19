A script for processing VCF files grouped by folder to look for enrichment of features using fischer's exact test.  


Currently the folder groupings and contrasts are hard coded but could easily be made a parameter.  


python vcf_enrichment.py --vcf_dir ./vcf_folder/ --output_prefix mut_posv2 --feature_names mutation position  

python vcf_enrichment.py --vcf_dir /Users/anwarren/Documents/ndssl/saber/ST_Tumor_viruses/ --output_prefix mut_posv2 --progression_mode progression2 --feature_names mutation position  

