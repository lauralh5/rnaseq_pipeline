#$ -S /bin/bash
#$ -cwd
#$ -j yes
#$ -V

WD=$1
MF=$2
NUM_SAMPLES=$3
INS_DIR=$4

## Change directory results

cd $WD/$MF/results

# Generar el fichero mergelist.txt con las rutas a los transcripts.gtf

i=1

while [ $i -le $NUM_SAMPLES ]
do
	echo $WD/$MF/samples/sample$i/sample$i.gtf >> mergelist.txt
	((i++))
done


# Combinar los trasncriptomas parciales de cada muestra
stringtie --merge -G ../annotation/annotation.gtf -o stringtie_merged.gtf mergelist.txt

# Comparar transcriptoma completo con la referencia
gffcompare -r ../annotation/annotation.gtf -G -o merged mergelist.txt


## Parallel quantification of each sample

i=1

while [ $i -le ${NUM_SAMPLES} ]
do
        qsub -N quant_$i -o quant_$i $INS_DIR/sample_quantification.sh $i $WD $MF $INS_DIR
	((i++))
done
