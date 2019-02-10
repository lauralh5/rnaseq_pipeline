###############################################################
## Tecnolog??as ??micas y Bioinform??tica 2017/2018             ##
## M??ster en Gen??tica Molecular y Biotecnolog??a              ##         
## Universidad de Sevilla                                    ##
##                                                           ##
## An??lisis transcript??micos masivos basados en RNA-seq.     ##  
##                                                           ##
## Prof. Francisco J. Romero-Campero fran@us.es              ##
###############################################################

## El paquete de bioconductor ballgown proporciona las funciones necesarias para 
## realizar un an??lisis de expresi??n g??nica diferencial y visualizar los resultados
## a partir de procesamiento de los datos brutos de secuenciaci??n realizados con 
## hisat2 y stringtie. 

## Para ejecutar con ??xito este script es necesario descarga la carpeta samples
## completa a tu ordenador, mover este script a la carpeta samples y fijar el 
## Working Directory To Source File Location. 

## Instalaci??n y carga de los paquetes necesarios. S??lo es necesario instalar los
## paquetes la primera vez que se ejecuta este script en un ordenador el resto de las
## veces bastar?? cargar los paquetes simplemente. 

library(ballgown)
library(genefilter)
library(plyranges)
library(CNVtools)
library(affy)
library(simpleaffy)
library(affyPLM)
library(limma)
library(annaffy)
library(ath1121501cdf)

## Para cargar los datos es necesario crear previamente un fichero tabular
## que contenga como primera columna los nombres de las carpetas donde se guarda
## cada muesra t??picamente sample1, sample2, etc ... El resto de columnas
## haran referencia al genotipo, tratamiento y dem??s caracteri??ticas de cada muestra. 
pheno.data <- read.csv("pheno_data.csv")
pheno.data

## La funci??n ballgown se usa para cargar o leer los datos. Es necesario especificar
## el directorio donde se encuentran las muestras. En nuestro caso especificamos .
## para indicar que se encuentran en el actual directorio. 
bg.data <- ballgown(dataDir = ".", samplePattern = "sample", pData=pheno.data)
bg.data  #9288 transcritos
sampleNames(bg.data)

## La funci??n gexpr extrae los niveles de expresi??n g??nicos en cada muestra
## medidos como FPKM (Fragments Per Kilobase of exon and Million of mapped reads)
gene.expression <- gexpr(bg.data)
head(gene.expression)

## Nombramos las columnas con los nombres de nuestras muestras. 
colnames(gene.expression) <- c("abc_1","abc_2","col0_1","col0_2")
head(gene.expression)


## Construimos un boxplot para comprobar que las distribuciones globales de las
## muestras son similares y comparables.
boxplot(log2(gene.expression + 1),col=rainbow(ncol(gene.expression)),ylab="log2(FPKM + 1)",cex.lab=1.5)

## Calculamos la matrix de expresi??n media. 
col0 <- (gene.expression[,"col0_1"] + gene.expression[,"col0_2"])/2
abc <- (gene.expression[,"abc_1"] + gene.expression[,"abc_2"])/2

mean.expression <- matrix(c(col0,abc),ncol=2)
colnames(mean.expression) <- c("col0","abc")
rownames(mean.expression) <- rownames(gene.expression)
head(mean.expression)

## Previsualizamos el efecto de la mutaci??n en un scatterplot.
plot(log2(col0+1),log2(abc+1),pch=19,cex=0.7,xlab="col0",ylab=substitute(italic("abc")),cex.lab=1.5)

## Pasamos a realizar un an??lisis de expresi??n g??nica diferencial. 
## Primero debido a la alta sensibilidad del RNA-seq eliminamos genes con muy 
## bajos niveles de expresi??n o que var??an muy poco. 
bg.data.filt <- subset(bg.data,"rowSums(texpr(bg.data)) >1 & rowVars(texpr(bg.data)) >1",genomesubset=TRUE)
bg.data.filt  #6169 transcritos
results.genes <- stattest(bg.data.filt,feature="gene",covariate="genotype",getFC=TRUE, meas="FPKM")
class(results.genes)
head(results.genes)

experimental.design <- model.matrix(~ -1+factor(c(1,1,2,2)))
experimental.design
colnames(experimental.design) <- c("col0_1", "col0_2", "abc_1", "abc_2")
experimental.design

linear.fit <- lmFit(gene.expression, experimental.design)

contrast.matrix <- makeContrasts(col0-abc,levels=c("col0","abc"))

contrast.linear.fit <- contrasts.fit(linear.fit, contrast.matrix)
contrast.results <- eBayes(contrast.linear.fit)


abc <- topTable(contrast.results, number=9288,coef=1,sort.by="logFC")
head(abc)

fold.change <- abc[["logFC"]]
q.value<- abc [["adj.P.Val"]]
genes.ids <- rownames(abc)


activated.genes <- genes.ids[fold.change > 1 & q.value < 0.05]
repressed.genes <- genes.ids[fold.change < - 1 & q.value < 0.05]

length(activated.genes)  #703
length(repressed.genes)  #788



## Resaltar activados y reprimidos
plot(log2(wt+1),log2(phyA+1),pch=19,cex=0.7,xlab="WT",ylab=substitute(italic("phyA")),cex.lab=1.5)
points(log2(wt[activated.genes]+1),log2(phyA[activated.genes]+1),pch=19,cex=0.7,col="red",cex.lab=1.5)
points(log2(wt[repressed.genes]+1),log2(phyA[repressed.genes]+1),pch=19,cex=0.7,col="blue",cex.lab=1.5)

plot(log2(wt),log2(phyA),pch=19,cex=0.7,xlab="WT",ylab=substitute(italic("phyA")),cex.lab=1.5,xlim=c(0,15),ylim=c(0,15))
points(log2(wt[activated.genes]),log2(phyA[activated.genes]),pch=19,cex=0.7,col="red",cex.lab=1.5)
points(log2(wt[repressed.genes]),log2(phyA[repressed.genes]),pch=19,cex=0.7,col="blue",cex.lab=1.5)

#Listas de genes activados y reprimidos
write.table(x = activated.genes,file = "activated_genes.txt",quote = FALSE,row.names = FALSE,sep = "\t")
write.table(x = repressed.genes,file = "repressed_genes.txt",quote = FALSE,row.names = FALSE,sep = "\t")

## De forma similar podemos realizar un an??lisis de expresi??n diferencial de transcritos
## en lugar de genes. 
transcript.expression <- texpr(bg.data, 'FPKM')
class(transcript.expression)
head(transcript.expression)

i <- 1
plotMeans(activated.genes[i], bg.data,groupvar="genotype",legend=FALSE)
i <- i + 1

i <- 16
plotMeans(activated.genes[i], bg.data,groupvar="genotype",legend=FALSE)

i <- 51
plotMeans(activated.genes[i], bg.data,groupvar="genotype",legend=FALSE)

i <- i + 1
gene.name <- activated.genes[i]
control.samples <- c("col0_1","col0_2")
treatment.samples <- c("abc_1","abc_2")
condition.names <- c("Col-0","AtBMI1abc")

head(gene.expression)
control.expression <- gene.expression[gene.name,control.samples]
treatment.expression <- gene.expression[gene.name,treatment.samples]

control.mean <- mean(control.expression)
treatment.mean <- mean(treatment.expression)

control.sd <- sd(control.expression)
treatment.sd <- sd(treatment.expression)

means <- c(control.mean,treatment.mean)
sds <- c(control.sd,treatment.sd)

ymax <- means + sds

pos <- barplot(means,col=rainbow(length(condition.names)),names.arg = condition.names,ylim=c(0,max(ymax)*1.1),main=gene.name,cex.main=2)
for(i in 1:length(condition.names))
{
  arrows(pos[i],means[i]-sds[i],pos[i],means[i]+sds[i],code=3,angle=90,length=0.1)  
}


