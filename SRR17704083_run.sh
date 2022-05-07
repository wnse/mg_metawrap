
echo "SRR17704083 "
echo "BMTAGGER_DB=/mnt/data/metawrap_db/BMTAGGER_INDEX/" >>/root/metaWRAP/bin/config-metawrap
echo "KRAKEN2_DB=/mnt/data/metawrap_db/KRAKEN2_DB/KRAKEN2_DB/" >>/root/metaWRAP/bin/config-metawrap
echo "BLASTDB=/mnt/data/metawrap_db/NCBI_nt" >>/root/metaWRAP/bin/config-metawrap
echo "TAXDUMP=/mnt/data/metawrap_db/NCBI_tax" >>/root/metaWRAP/bin/config-metawrap


metawrap read_qc -t 42 -1 /mnt/data/mg_metawrap/test_data/SRR17704083_1.fastq -2 /mnt/data/mg_metawrap/test_data/SRR17704083_2.fastq -o /mnt/data/mg_metawrap/test_out/SRR17704083/READ_QC 

metawrap assembly -m 100 -t 42 -1 /mnt/data/mg_metawrap/test_out/SRR17704083/READ_QC/final_pure_reads_1.fastq -2 /mnt/data/mg_metawrap/test_out/SRR17704083/READ_QC/final_pure_reads_2.fastq -o /mnt/data/mg_metawrap/test_out/SRR17704083/ASSEMBLY 

metawrap kraken2 -o /mnt/data/mg_metawrap/test_out/SRR17704083/KRAKEN -t 42 -s 1000000 /mnt/data/mg_metawrap/test_out/SRR17704083/READ_QC/final_pure_reads_1.fastq /mnt/data/mg_metawrap/test_out/SRR17704083/READ_QC/final_pure_reads_2.fastq /mnt/data/mg_metawrap/test_out/SRR17704083/ASSEMBLY/final_assembly.fasta 

metawrap binning -m 100 -t 42 --metabat2 --maxbin2 --concoct -o /mnt/data/mg_metawrap/test_out/SRR17704083/INITIAL_BINNING -a /mnt/data/mg_metawrap/test_out/SRR17704083/ASSEMBLY/final_assembly.fasta /mnt/data/mg_metawrap/test_out/SRR17704083/READ_QC/final_pure_reads_1.fastq /mnt/data/mg_metawrap/test_out/SRR17704083/READ_QC/final_pure_reads_2.fastq

metawrap bin_refinement -c 50 -x 10 -t 42 -o /mnt/data/mg_metawrap/test_out/SRR17704083/BIN_REFINEMENT -A /mnt/data/mg_metawrap/test_out/SRR17704083/INITIAL_BINNING/metabat2_bins/ -B /mnt/data/mg_metawrap/test_out/SRR17704083/INITIAL_BINNING/maxbin2_bins/ -C /mnt/data/mg_metawrap/test_out/SRR17704083/INITIAL_BINNING/concoct_bins/

metawrap reassemble_bins -o /mnt/data/mg_metawrap/test_out/SRR17704083/BIN_REASSEMBLY -1 /mnt/data/mg_metawrap/test_out/SRR17704083/READ_QC/final_pure_reads_1.fastq -2 /mnt/data/mg_metawrap/test_out/SRR17704083/READ_QC/final_pure_reads_2.fastq -t 42 -m 100 -c 50 -x 10 -b /mnt/data/mg_metawrap/test_out/SRR17704083/BIN_REFINEMENT/metawrap_50_10_bins

metaWRAP annotate_bins -o /mnt/data/mg_metawrap/test_out/SRR17704083/FUNCT_ANNOT/ -t 42 -b /mnt/data/mg_metawrap/test_out/SRR17704083/BIN_REASSEMBLY/reassembled_bins/

metawrap blobology -a /mnt/data/mg_metawrap/test_out/SRR17704083/ASSEMBLY/final_assembly.fasta -t 42 -o /mnt/data/mg_metawrap/test_out/SRR17704083/BLOBOLOGY/ --bins /mnt/data/mg_metawrap/test_out/SRR17704083/BIN_REFINEMENT/metawrap_50_10_bins /mnt/data/mg_metawrap/test_out/SRR17704083/READ_QC/final_pure_reads_1.fastq /mnt/data/mg_metawrap/test_out/SRR17704083/READ_QC/final_pure_reads_2.fastq

metawrap classify_bins -b /mnt/data/mg_metawrap/test_out/SRR17704083/BIN_REASSEMBLY/reassembled_bins -o /mnt/data/mg_metawrap/test_out/SRR17704083/BIN_CLASSIFICATION -t 42

metawrap quant_bins -b /mnt/data/mg_metawrap/test_out/SRR17704083/BIN_REFINEMENT/metawrap_50_10_bins -o /mnt/data/mg_metawrap/test_out/SRR17704083/QUANT_BINS/ -a /mnt/data/mg_metawrap/test_out/SRR17704083/ASSEMBLY/final_assembly.fasta /mnt/data/mg_metawrap/test_out/SRR17704083/READ_QC/final_pure_reads_1.fastq /mnt/data/mg_metawrap/test_out/SRR17704083/READ_QC/final_pure_reads_2.fastq
