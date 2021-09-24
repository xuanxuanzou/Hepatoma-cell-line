PISA parse -t 4 -f -q 4 -dropN -config BGI_droplet_scRNA_readStructureV2_T1.json -cbdis /barcode_counts_raw.txt -report sequencing_report.csv DP8400011449TR_L01_1_1.fq.gz DP8400011449TR_L01_1_2.fq.gz -1 /reads.fq

STAR --outSAMmultNmax 1 --outStd SAM --outSAMunmapped Within --runThreadN 4 --genomeDir star_index --readFilesIn reads.fq --outFileNamePrefix /temp/ 1> aln.sam && \
PISA sam2bam -t 4 -k -o /temp/aln.bam -report alignment_report.csv aln.sam && rm -f aln.sam

PISA attrcnt -cb CB -tags UB,GN -@ 4 -dedup -o /temp/cell_stat.txt -list /temp/barcode_raw_list.txt final.bam

