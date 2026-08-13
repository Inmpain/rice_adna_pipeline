The Rice 3000 Genomes Project Data.
===================================

Files:

1 - consortium_list.csv
contains the list of consortium members at the time of releasing the sequence data

2 - rice_line_metadata_20140521.tsv
The list of all details relating to the 3000 rice lines sequenced, in a tab seperated list

3 - rice_line_metadata_20140521.xlsx
the 3000 rice lines metadata in excel spreadsheet format, split into 2 sheets, one containing rice lines donated by IRRI, the other sheet for CAAS donated lines

4 - rice3K_sequence_data_at_EBI
The 3000 rice genomes sequence data is no longer available for download from GigaDB, it is hosted by the Sequence Read Archives (SRA) at EBI, DDBJ and NCBI. Please choose the appropriate geographical location:

Europe:
http://www.ebi.ac.uk/ena/data/view/PRJEB6180

USA:
http://www.ncbi.nlm.nih.gov/sra/?term=PRJEB6180

Asia:
http://trace.ddbj.nig.ac.jp/DRASearch/study?acc=ERP005654

5 - seq_file_mapping_to_SRA.txt
Contains all the file names (as originally stored in GigaDB) with the mappings to the SRA accessions for run, experiment, sample and project, as well as the specific file locations on the EBI FTP server. Its in tab delimited (TSV) format with the following column headings:

Project_Acc - 3000 rice lines project accession (PRJEB6180)
Sample_alias - rice line sample name used by the rice consortium
Sample_Acc - rice line Sample ID provided by SRA
Sample_taxid - NCBI taxonomic ID (4530)
Experiment_Acc - SRA accession given to each individual illumina run
Experiment_alias - arbitrary "experiment" name given during submission
Run_acc - SRA accession given to each individual illumina run
Run_alias - name used by consortium for each illumina run
read file - individual forward or reverse sequence file names
submission accession - SRA submission accession (submission done in batches)
files as submitted(from ENA) - FTP server location at EBI for the raw sequence file as submitted
ENA processed FastQ file - FTP server location at EBI for the processed FastQ file

