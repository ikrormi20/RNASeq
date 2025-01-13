#This RNASeq pipeline maps RNASeq to Hg19 and qauntifies gene expression level
#To generate genome index
#/home/ngs/Downloads/STAR-2.7.0f/bin/Linux_x86_64/STAR --runThreadN 12 --runMode genomeGenerate --genomeDir hg19_index/ --genomeFastaFiles hg19_fasta/Homo_sapiens.GRCh37.dna.primary_assembly.fa --sjdbGTFfile hg19_annotationfile/Homo_sapiens.GRCh37.87.gtf --sjdbOverhang 99;
#To perform mapping with STAR
/home/ngs/Downloads/STAR-2.7.0f/bin/Linux_x86_64/STAR --runThreadN 12 --genomeDir hg19_index/ --readFilesCommand gunzip -c --readFilesIn /mnt/Diag/WGS-ikrormi/RNASeq/Prav_AD_S14_L001_R1_001.fastq.gz /mnt/Diag/WGS-ikrormi/RNASeq/Prav_AD_S14_L001_R2_001.fastq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix Prav;

#samtools view -bS  Aligned.out.sam > Aligned.bam;
#samtools index Aligned.bam;
#samtools sort Aligned.bam > sorted_Aligned.bam;
#samtools index sorted_Aligned.bam;

#*******************************************************************************************************************
#Command for RSEM
#To build reference
rsem-prepare-reference --gtf /mnt/Diag/RNASeq/hg19_annotationfile/Homo_sapiens.GRCh37.87.gtf --star --num-threads 12 --star-path /home/ngs/Downloads/STAR-2.7.0f/bin/Linux_x86_64  /mnt/Diag/RNASeq/hg19_fasta/Homo_sapiens.GRCh37.dna.primary_assembly.fa --star-sjdboverhang 99 ref/human_ensembl;

# To calculate expression levels

rsem-calculate-expression -p 8 --paired-end --star --num-threads 12 --star-path /home/ngs/Downloads/STAR-2.7.0f/bin/Linux_x86_64 --estimate-rspd --append-names --output-genome-bam /mnt/Diag/WGS-ikrormi/RNASeq_GATK/Prav_S13_L001_R1_001.fastq /mnt/Diag/WGS-ikrormi/RNASeq_GATK/Prav_S13_L001_R2_001.fastq ref/human_ensembl out/Prav;

rsem-calculate-expression -p 8 --paired-end --bam --estimate-rspd --append-names --output-genome-bam out/Prav.bam ref/human_ensembl out/Prav;

# To generate a diagnostic pdf
rsem-plot-model out/Prav Prav_diagnostic.pdf;

# R command to generate first 10 highly expressed genes

data = read.table("Prav.genes.results", header=T, stringsAsFactors=F)
#For decreasing order
idx_high = order(data["TPM"], decreasing=T) 
high = data[idx_high[1:10], c("gene_id", "expected_count", "TPM")]
write.csv(high, file = "high.csv")
#For increasing order
idx_low = order(data["TPM"], decreasing=F)
low = data[idx_low[1:10], c("gene_id", "expected_count", "TPM")]
write.csv(low, file = "low.csv")
