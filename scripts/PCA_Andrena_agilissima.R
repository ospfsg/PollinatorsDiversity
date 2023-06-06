library("SNPRelate")
library("RColorBrewer")
library(randomcoloR)
library(adegenet)
library(ggplot2)

#input
vcf.fn<-"AagilissimaFinal.recode.vcf" #load the data
#snpgdsVCF2GDS(vcf.fn,"test_pca1",method="copy.num.of.ref") # to extract and store dosage (0, 1, 2) of the reference allele for all variant sites, including bi-allelic SNPs, multi-allelic SNPs, indels and structural variants.
snpgdsVCF2GDS(vcf.fn,"Aagilissima_final",method="biallelic.only") # to exact bi-allelic and polymorhpic SNP data (excluding monomorphic variants)
snpgdsSummary("test_pca2")


#open the translated vcf imported
genfile_1 <-snpgdsOpen("Aagilissima_final")

#import population names
pop_code <- read.table("Andrena agilissima.txt",header = FALSE, sep = "") ## File 1
pop_code
#popmap<- scan("Andrena agilissima.txt", what=character(), ) ## File 1
#popmap
table(pop_code$V1) # Number of individuals per specie 
head(pop_code)
View(pop_code)
str(pop_code)

#PCA
pca_1<-snpgdsPCA(genfile_1, autosome.only = FALSE, eigen.method=c("DSPEV"), bayesian = TRUE, need.genmat = TRUE) # Bayesian Normalization
pca_2<-snpgdsPCA(genfile_1, autosome.only = FALSE, eigen.method=c("DSPEV"), bayesian = FALSE) # No Bayesian Normalization
View(pca_1$genmat) # Genetic Covariance: Weight of each Principal Component in each sample
SNPLoad<-snpgdsPCASNPLoading(pca_1, genfile_1, verbose=TRUE)
SNPLoad2<- snpgdsPCASNPLoading(pca_2, genfile_1, verbose = TRUE)
plot(SNPLoad$snploading[,1], type="h", ylab="PC 1")
pca_1$sample.id

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
View(tab)

tab$color<- c("darkblue", "darkblue", "darkblue", "darkblue", "green3", "darkorange3", "darkorchid", "darkorchid", "darkorchid", "darkorchid", "black", "black")
label_legend_agilissima<- c("BET", "CAT", "FZZ", "GAL", "JAC")
color_legend<- c("darkblue", "green3", "darkorange3", "darkorchid", "black")
color_legend

tab_PC3 <- data.frame(sample.id = pca_1$sample.id,
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
plot(tab$EV1, tab$EV2, col=tab$color, xlab=ev1 ,ylab=ev2, main="PCA for Andrena agilissima", pch=16)
legend("topleft", legend=label_legend_agilissima, col=color_legend, cex = 0.8, pch = 16)
text(tab$EV1, tab$EV2, row.names(tab))

# Amostras que se separam: And_Bet_007_Poll123, And_Bet_013_Poll123

plot(tab_PC3$EV1, tab_PC3$EV3, col=color_vector_altitude, xlab=ev1 ,ylab=ev3 , pch=myPch_1, main="PCA for Bombus terrestris")
legend(2.8, 0,6, legend=pop_code$V1, col=pop_code$V2, cex = 0.5) 
dev.off()














