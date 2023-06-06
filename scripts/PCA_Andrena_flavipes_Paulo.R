library("SNPRelate")
library("RColorBrewer")
library(randomcoloR)
library(adegenet)
library(ggplot2)

#input
vcf.fn<-"Aflavipes_Paulo_Final.recode.vcf" #load the data
#snpgdsVCF2GDS(vcf.fn,"test_pca1",method="copy.num.of.ref") # to extract and store dosage (0, 1, 2) of the reference allele for all variant sites, including bi-allelic SNPs, multi-allelic SNPs, indels and structural variants.
snpgdsVCF2GDS(vcf.fn,"Aflavipes",method="biallelic.only") # to exact bi-allelic and polymorhpic SNP data (excluding monomorphic variants)
snpgdsSummary("test_pca2")


#open the translated vcf imported
genfile_1 <-snpgdsOpen("Aflavipes")

#PCA
pca_1<-snpgdsPCA(genfile_1, autosome.only = FALSE, eigen.method=c("DSPEV"), bayesian = TRUE, need.genmat = TRUE) # Bayesian Normalization
#pca_2<-snpgdsPCA(genfile_1, autosome.only = FALSE, eigen.method=c("DSPEV"), bayesian = FALSE) # No Bayesian Normalization
View(pca_1$genmat) # Genetic Covariance: Weight of each Principal Component in each sample
SNPLoad<-snpgdsPCASNPLoading(pca_1, genfile_1, verbose=TRUE)
#SNPLoad2<- snpgdsPCASNPLoading(pca_2, genfile_1, verbose = TRUE)
plot(SNPLoad$snploading[,1], type="h", ylab="PC 1")
pca_1$sample.id

#import population names
pop_code <- read.table("Andrena_flavipes_Paulo.txt",header = FALSE, sep = "") ## File 1
pop_code
#popmap<- scan("Andrena agilissima.txt", what=character(), ) ## File 1
#popmap
table(pop_code$V1) # Number of individuals per specie 
head(pop_code)
View(pop_code)
str(pop_code)

# variance proportion (%)
pca_1$varprop
pc.percent<- pca_1$varprop * 100
print(round(pc.percent, 2))
lbls <- paste("PC", 1:5, "\n", format(pc.percent[1:5], digits=2), "%", sep="")
pairs(pca_1$eigenvect[,1:5], col=colors, labels=lbls)

#get sample ids
sample.id<-read.gdsn(index.gdsn(genfile_1,"sample.id"))
sample.id

#data frame of eigenvectors
tab <- data.frame(sample.id = pca_1$sample.id,
                  pop = pop_code$V1[match(pca_1$sample.id, sample.id)],
                  EV1 = pca_1$eigenvect[,1],    # the first eigenvector
                  EV2 = pca_1$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
tab$pop

tab$color<- c(1, 2, 3, 4,4, 5, 6, 7,7, 8, 9, 10,10,10, 11, 12,12,12, 9, 13,13, 10, 11,11,11, 8, 5,5, 3,3, 1, 2,2, 4, 5, 6, 7, 8, 9, 13, 4, 1,1,1,1, 2, 6, 7, 8, 9,9,9, 13, 14, 12,12,12)
tab$color<- c("gray50", "yellow1", "darkblue", "green3", "green3", "coral", "tan4", "darkorchid", "darkorchid", "black", "red1", "springgreen", "springgreen", "springgreen", "deepskyblue", "olivedrab", "olivedrab", "olivedrab", "red1", "darkgoldenrod1", "darkgoldenrod1", "springgreen", "deepskyblue", "deepskyblue", "deepskyblue","black", "coral", "coral", "darkblue", "darkblue", "gray50", "yellow1", "yellow1", "green3", "coral", "tan4", "darkorchid", "black", "red1", "darkgoldenrod1", "green3", "gray50", "gray50", "gray50", "gray50", "yellow1", "tan4", "darkorchid", "black", "red1", "red1", "red1", "darkgoldenrod1", "magenta", "olivedrab", "olivedrab", "olivedrab")
label_legend_flavipes<- c("AND", "BEJ", "BET", "CAT", "CLM", "FVF", "GAL", "JAC", "MAL", "ODE", "SIN","SMA", "MGR", "SDA")
color_legend<- c("gray50", "yellow1", "darkblue", "green3", "coral", "tan4", "darkorchid", "black", "red1", "springgreen", "deepskyblue", "olivedrab", "darkgoldenrod1", "magenta")
color_legend

#tab_PC3 <- data.frame(sample.id = pca_1$sample.id,
                      pop = pop_code$V1[match(pca_1$sample.id, sample.id)],
                      EV1 = pca_1$eigenvect[,1],    # the first eigenvector
                      EV3 = pca_1$eigenvect[,3],    # the second eigenvector
                      stringsAsFactors = FALSE)

#Plot PC1-PC2 and PC1-PC3

par(mfrow = c(1, 2))
par(xpd = T, mar = par()$mar + c(0,0,0,7))
ev1=paste("PC1", round(pc.percent[1], digits = 2), "%")
ev2=paste("PC2", round(pc.percent[2], digits = 2), "%")
#ev3=paste("PC3", round(pc.percent[3], digits = 2), "%")
plot(tab$EV1, tab$EV2, col=tab$color, xlab=ev1 ,ylab=ev2, main="PCA for Andrena flavipes", pch=16)
legend("topright", legend=label_legend_flavipes, col=color_legend, cex = 0.6, pch = 16)
text(tab$EV1, tab$EV2, row.names(tab))

# Amostras que se separam: Todas as de JAC

plot(tab_PC3$EV1, tab_PC3$EV3, col=color_vector_altitude, xlab=ev1 ,ylab=ev3 , pch=myPch_1, main="PCA for Bombus terrestris")
legend(2.8, 0,6, legend=pop_code$V1, col=pop_code$V2, cex = 0.5) 
dev.off()














