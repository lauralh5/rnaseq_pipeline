#$ -S /bin/bash
#$ -cwd
#$ -j yes
#$ -V

## Sample alignment and assambling script

## Defining parameters

i=$1
echo $i >> ../log/parameters
INS_DIR=$2
echo $INS_DIR >> ../log/parameters
WD=$3
echo $WD >> ../log/parameters
MF=$4
echo $MF >> ../log/parameters
SAMPLE=$5
echo $SAMPLE >> ../log/parameters
ACC_SAMPLES=$6
echo $ACC_SAMPLES >> ../log/parameters
ANNOTATION=$7
echo $ANNOTATION >> ../log/parameters
NUM_SAMPLES=$8
echo $NUM_SAMPLES >> ../log/parameters
INDEX=$9
echo $INDEX >> ../log/parameters

## Access to each sample directory and download the GEO .fastq files of that sample

cd $WD/$MF/samples/sample$i
echo "Access to sample directory" >> ../../log/samples

cp $ACC_SAMPLES .
gunzip ${SAMPLE}.gz

echo "Sample copied" >> ../../log/samples

## Map reads to the reference genome depending on single or paired end

##if [ -f ${SAMPLE} ]
##then
##	fastqc ${SAMPLE}
##	fastqc ${SAMPLE}_
##	hisat2 --dta -x $INDEX -1 ${SAMPLE} -2 ${SAMPLE} -S sample_$i.sam

##else
fastqc ${SAMPLE}
hisat2 --dta -x $INDEX -U ${SAMPLE} -S sample$i.sam

##fi

## sort and create  bam file

samtools sort -o sample$i.bam sample$i.sam

## Transcript assembly of the sample

stringtie -G $ANNOTATION -o sample$i.gtf -l sample$i sample$i.bam


## Synchronization before next task. If all samples are already processed, next script or task can be uploaded.

## Write on the blackboard
echo sample$i done >> ../../log/blackboard 

## Read from the blackboard
SAMPLES_DONE=$(wc -l ../../log/blackboard | awk '{ print $1 }')

## Synchronization point

if [ ${SAMPLES_DONE} -eq ${NUM_SAMPLES} ]
then
	qsub -N transcriptome_merging -o transcriptome_merging $INS_DIR/merge.sh $WD $MF $NUM_SAMPLES $INS_DIR

fi



