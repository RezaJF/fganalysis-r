

duckdb -c "COPY (select * from read_csv('finngen_R13_service_sector_detailed_longitudinal_1.0.txt.gz',delim='\t') ) to 'finngen_R13_service_sector_detailed_longitudinal_1.0.parquet' (FORMAT parquet,  PARTITION_BY (SOURCE))"

duckdb -c 'COPY (select * from read_parquet("finngen_R13_kanta_lab_1.0.parquet")) to "finngen_R13_kanta_lab_1.0.parquet" (format PARQUET)'


duckdb -c "COPY (select * from read_csv('finngen_R13_minimum_extended_1.0.txt.gz',delim='\t') ) to 'finngen_R13_minimum_extended_1.0.parquet' (FORMAT parquet)"


duckdb -c "COPY (select * from read_csv('R13_COV_PHENO_V0.txt.gz',delim='\t', nullstr="NA") ) to 'R13_COV_PHENO_V0.parquet' (FORMAT parquet)"

scripts/process_data.R
duckdb -c "COPY (select * from read_csv('finngen_R13_hilmo_avohilmo_extended_1.0_dedup_fixbp.txt.gz',delim='\t', nullstr="NA") ) to 'finngen_R13_hilmo_avohilmo_extended_1.0_dedup_fixbp.parquet' (FORMAT parquet)"


scripts/process_drugs.py finngen_R13_drug_events_2.0.csv finngen_R13_drug_events_2.0.simple.csv

duckdb -c "COPY (select * from read_csv('finngen_R13_drug_events_2.0.simple.csv.gz',delim='\t' columns = {'VNR'='VARCHAR'}) order by ATC) to 'finngen_R13_drug_events_2.0.simple.parquet' (FORMAT parquet)"