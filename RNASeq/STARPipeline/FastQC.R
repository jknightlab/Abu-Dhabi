# Fastqc

fastqc data/M120_*R1*.fastq

fastqc -o ${OUTPUT_DIRECTORY_FASTQC} ${SAMPLE_PATH}_1.fastq.gz ${SAMPLE_PATH}_2.fastq.gz===========-

fastqc -o ${OUTPUT_DIRECTORY_FASTQC} ${SAMPLE_PATH}"_1.fastq.gz" ${SAMPLE_PATH}"_2.fastq.gz" > ${OUTPUT_DIRECTORY_LOGS}"/fastQC_"${SAMPLE_NAME}".log"

fastqc stdin

for	F	in	`ls	*.fastq`
do
fastqc	-o	$FASTQ	$F
done

module load java/1.8.0_latest

java -Xmx13G -jar /apps/well/nameofapp
