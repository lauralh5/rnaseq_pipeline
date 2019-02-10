#$ -S /bin/bash
#$ -cwd
#$ -j yes
#$ -V

i=$1
WD=$2
MF=$3
INS_DIR=$4

# Change to the sample directory

cd $WD/$MF/samples/sample$i

# Cuantificación de los niveles de expresión
stringtie -e -B -G ../../annotation/annotation.gtf -o sample$i.gtf sample$i.bam

# Eliminar ficheros sam y fastq
rm *.sam
rm *.bam
rm *.fastq

echo $i DONE >> ../../log/quantification
SAMPLE_DONE=$(wc -l ../../log/quantification | awk '{ print $1 }')

if [ ${SAMPLE_DONE} -eq 4 ]
then
	cp -a $INS_DIR/examen_ballgown.R $WD/$MF/samples
	Rscript --vanilla $WD/$MF/samples/examen_ballgown.R $INS_DIR/pheno_data.csv
fi

