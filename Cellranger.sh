# D8-CO
cat D8-CO_S1_L007_R1_001.fastq.gz D8-CO_S2_L008_R1_001.fastq.gz D8-CO_S3_L004_R1_001.fastq.gz D8-CO_S4_L001_R1_001.fastq.gz > D8-CO_S1_L005_R1_001.fastq.gz
cat D8-CO_S1_L007_R2_001.fastq.gz D8-CO_S2_L008_R2_001.fastq.gz D8-CO_S3_L004_R2_001.fastq.gz D8-CO_S4_L001_R2_001.fastq.gz > D8-CO_S1_L005_R2_001.fastq.gz

# D20-CO
cat D20-CO_S1_L007_R1_001.fastq.gz D20-CO_S2_L008_R1_001.fastq.gz D20-CO_S3_L004_R1_001.fastq.gz D20-CO_S4_L001_R1_001.fastq.gz > D20-CO_S1_L005_R1_001.fastq.gz
cat D20-CO_S1_L007_R2_001.fastq.gz D20-CO_S2_L008_R2_001.fastq.gz D20-CO_S3_L004_R2_001.fastq.gz D20-CO_S4_L001_R2_001.fastq.gz > D20-CO_S1_L005_R2_001.fastq.gz

# D20-VCO
cat D20-VCO_S1_L007_R1_001.fastq.gz D20-VCO_S2_L002_R1_001.fastq.gz D20-VCO_S3_L008_R1_001.fastq.gz D20-VCO_S4_L007_R1_001.fastq.gz > D20-VCO_S1_L005_R1_001.fastq.gz
cat D20-VCO_S1_L007_R2_001.fastq.gz D20-VCO_S2_L002_R2_001.fastq.gz D20-VCO_S3_L008_R2_001.fastq.gz D20-VCO_S4_L007_R2_001.fastq.gz > D20-VCO_S1_L005_R2_001.fastq.gz


module load cellranger
SAMPLE_NAME=D8-CO
ID=${SAMPLE_NAME}_outs
FASTQS=./
SAMPLE=$SAMPLE_NAME
REFERENCE=~/references/refdata-gex-GRCh38-2020-A
cellranger count \
--id=$ID \
--fastqs=$FASTQS \
--sample=$SAMPLE \
--transcriptome=$REFERENCE \
--localmem=60 \
--localcores=10