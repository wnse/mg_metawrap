
echo "BMTAGGER_DB=/mnt/data/metawrap_db/BMTAGGER_INDEX/" >>/root/metaWRAP/bin/config-metawrap
echo "KRAKEN2_DB=/mnt/data/metawrap_db/KRAKEN2_DB/KRAKEN2_DB/" >>/root/metaWRAP/bin/config-metawrap
echo "BLASTDB=/mnt/data/metawrap_db/NCBI_nt" >>/root/metaWRAP/bin/config-metawrap
echo "TAXDUMP=/mnt/data/metawrap_db/NCBI_tax" >>/root/metaWRAP/bin/config-metawrap


mkdir /mnt/data/mg_metawrap/test_out/test_samples/READ_QC
metawrap read_qc -1 /mnt/data/mg_metawrap/test_data/test3_1.fastq -2 /mnt/data/mg_metawrap/test_data/test3_2.fastq -t 24 -o /mnt/data/mg_metawrap/test_out/test_samples/READ_QC/test3
metawrap read_qc -1 /mnt/data/mg_metawrap/test_data/test1_1.fastq -2 /mnt/data/mg_metawrap/test_data/test1_2.fastq -t 24 -o /mnt/data/mg_metawrap/test_out/test_samples/READ_QC/test1
metawrap read_qc -1 /mnt/data/mg_metawrap/test_data/test2_1.fastq -2 /mnt/data/mg_metawrap/test_data/test2_2.fastq -t 24 -o /mnt/data/mg_metawrap/test_out/test_samples/READ_QC/test2

mkdir /mnt/data/mg_metawrap/test_out/test_samples/CLEAN_READS
mv /mnt/data/mg_metawrap/test_out/test_samples/READ_QC/test3/final_pure_reads_1.fastq /mnt/data/mg_metawrap/test_out/test_samples/CLEAN_READS/test3_1.fastq
mv /mnt/data/mg_metawrap/test_out/test_samples/READ_QC/test3/final_pure_reads_2.fastq /mnt/data/mg_metawrap/test_out/test_samples/CLEAN_READS/test3_2.fastq
mv /mnt/data/mg_metawrap/test_out/test_samples/READ_QC/test1/final_pure_reads_1.fastq /mnt/data/mg_metawrap/test_out/test_samples/CLEAN_READS/test1_1.fastq
mv /mnt/data/mg_metawrap/test_out/test_samples/READ_QC/test1/final_pure_reads_2.fastq /mnt/data/mg_metawrap/test_out/test_samples/CLEAN_READS/test1_2.fastq
mv /mnt/data/mg_metawrap/test_out/test_samples/READ_QC/test2/final_pure_reads_1.fastq /mnt/data/mg_metawrap/test_out/test_samples/CLEAN_READS/test2_1.fastq
mv /mnt/data/mg_metawrap/test_out/test_samples/READ_QC/test2/final_pure_reads_2.fastq /mnt/data/mg_metawrap/test_out/test_samples/CLEAN_READS/test2_2.fastq

mkdir /mnt/data/mg_metawrap/test_out/test_samples/CLEAN_READS_T
cat /mnt/data/mg_metawrap/test_out/test_samples/CLEAN_READS/*_1.fastq > /mnt/data/mg_metawrap/test_out/test_samples/CLEAN_READS_T/ALL_READS_1.fastq
cat /mnt/data/mg_metawrap/test_out/test_samples/CLEAN_READS/*_2.fastq > /mnt/data/mg_metawrap/test_out/test_samples/CLEAN_READS_T/ALL_READS_2.fastq

metawrap assembly -m 100 -t 42 -1 /mnt/data/mg_metawrap/test_out/test_samples/CLEAN_READS_T/ALL_READS_1.fastq -2 /mnt/data/mg_metawrap/test_out/test_samples/CLEAN_READS_T/ALL_READS_2.fastq -o /mnt/data/mg_metawrap/test_out/test_samples/ASSEMBLY 
metawrap kraken2 -t 42 -s 1000000 -o /mnt/data/mg_metawrap/test_out/test_samples/KRAKEN /mnt/data/mg_metawrap/test_out/test_samples/CLEAN_READS/*fastq /mnt/data/mg_metawrap/test_out/test_samples/ASSEMBLY/final_assembly.fasta
metawrap binning -m 100 -t 42 --metabat2 --maxbin2 --concoct -o /mnt/data/mg_metawrap/test_out/test_samples/INITIAL_BINNING -a /mnt/data/mg_metawrap/test_out/test_samples/ASSEMBLY/final_assembly.fasta /mnt/data/mg_metawrap/test_out/test_samples/CLEAN_READS/*fastq
metawrap bin_refinement -c 50 -x 10 -t 42 -o /mnt/data/mg_metawrap/test_out/test_samples/BIN_REFINEMENT -A /mnt/data/mg_metawrap/test_out/test_samples/INITIAL_BINNING/metabat2_bins/ -B /mnt/data/mg_metawrap/test_out/test_samples/INITIAL_BINNING/maxbin2_bins/ -C /mnt/data/mg_metawrap/test_out/test_samples/INITIAL_BINNING/concoct_bins/
metawrap blobology -t 42 -a /mnt/data/mg_metawrap/test_out/test_samples/ASSEMBLY/final_assembly.fasta -o /mnt/data/mg_metawrap/test_out/test_samples/BLOBOLOGY/ --bins /mnt/data/mg_metawrap/test_out/test_samples/BIN_REFINEMENT/metawrap_50_10_bins /mnt/data/mg_metawrap/test_out/test_samples/CLEAN_READS/*fastq
metawrap quant_bins -b /mnt/data/mg_metawrap/test_out/test_samples/BIN_REFINEMENT/metawrap_50_10_bins -o /mnt/data/mg_metawrap/test_out/test_samples/QUANT_BINS/ -a /mnt/data/mg_metawrap/test_out/test_samples/ASSEMBLY/final_assembly.fasta /mnt/data/mg_metawrap/test_out/test_samples/CLEAN_READS/*fastq
metawrap reassemble_bins -t 42 -m 100 -c 50 -x 10 -o /mnt/data/mg_metawrap/test_out/test_samples/BIN_REASSEMBLY -1 /mnt/data/mg_metawrap/test_out/test_samples/CLEAN_READS_T/ALL_READS_1.fastq -2 /mnt/data/mg_metawrap/test_out/test_samples/CLEAN_READS_T/ALL_READS_2.fastq -b /mnt/data/mg_metawrap/test_out/test_samples/BIN_REFINEMENT/metawrap_50_10_bins
metaWRAP annotate_bins -t 42 -o /mnt/data/mg_metawrap/test_out/test_samples/FUNCT_ANNOT/ -b /mnt/data/mg_metawrap/test_out/test_samples/BIN_REASSEMBLY/reassembled_bins
metawrap classify_bins -t 42 -o /mnt/data/mg_metawrap/test_out/test_samples/BIN_CLASSIFICATION -b /mnt/data/mg_metawrap/test_out/test_samples/BIN_REASSEMBLY/reassembled_bins

