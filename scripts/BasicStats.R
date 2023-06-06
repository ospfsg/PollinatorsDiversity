library(adegenet)
library(hierfstat)
library(poppr)
library(mmod)

#Fst v
#Gst v
#D Jost v 
#Fis v
#Ho v
#He v
#Ar v
#private Alleles v

path <- "/home/paulosousa/Working/Bees/Analysis/Aagilissima/BasicStats/Aagilissima_input.gen"
Genind_Aagilissima<- read.genepop(path, ncode = 3)
Aagilissima <- genind2hierfstat(Genind_Aagilissima, pop=NULL) #Converts genind objects from adegenet into a hierfstat data frame

print(levels(Aagilissima$pop)) #Change population names

levels(Aagilissima$pop) <- c(levels(Aagilissima$pop), "Baetic")
Aagilissima$pop[Aagilissima$pop == 'Bet_03'] <- 'Baetic'
levels(Aagilissima$pop) <- c(levels(Aagilissima$pop), "Ebro")
Aagilissima$pop[Aagilissima$pop == 'Jac_02'] <- 'Ebro'
levels(Aagilissima$pop) <- c(levels(Aagilissima$pop), "West")
Aagilissima$pop[Aagilissima$pop == 'Gal_04'] <- 'West'
Aagilissima$pop <- factor(Aagilissima$pop)
table(Aagilissima$pop)


path <- "/home/paulosousa/Working/Bees/Analysis/Aflavipes/BasicStats/Aflavipes_input.gen"
Genind_Aflavipes <- read.genepop(path, ncode= 3)
Aflavipes <- genind2hierfstat(Genind_Aflavipes, pop=NULL) #Converts genind objects from adegenet into a hierfstat data frame

print(levels(Aflavipes$pop)) #Change population names

levels(Aflavipes$pop) <- c(levels(Aflavipes$pop), "Iberia")
Aflavipes$pop[Aflavipes$pop == 'Sma_06'] <- 'Iberia'
levels(Aflavipes$pop) <- c(levels(Aflavipes$pop), "Pyrenees")
Aflavipes$pop[Aflavipes$pop == 'Jac_02'] <- 'Pyrenees'
Aflavipes$pop <- factor(Aflavipes$pop)
table(Aflavipes$pop)


path <- "/home/paulosousa/Working/Bees/Analysis/Lmalachurum/BasicStats/Lmalachurum_input.gen"
Genind_Lmalachurum <- read.genepop(path, ncode= 3)
Lmalachurum <- genind2hierfstat(Genind_Lmalachurum, pop=NULL) #Converts genind objects from adegenet into a hierfstat data frame

print(levels(Lmalachurum$pop)) #Change population names

levels(Lmalachurum$pop) <- c(levels(Lmalachurum$pop), "South")
Lmalachurum$pop[Lmalachurum$pop == 'Fzz_01'] <- 'South'
levels(Lmalachurum$pop) <- c(levels(Lmalachurum$pop), "Ebro")
Lmalachurum$pop[Lmalachurum$pop == 'Cat_04'] <- 'Ebro'
levels(Lmalachurum$pop) <- c(levels(Lmalachurum$pop), "Central")
Lmalachurum$pop[Lmalachurum$pop == 'Lgl_Gua_014_Poll119'] <- 'Central'
Lmalachurum$pop <- factor(Lmalachurum$pop)
table(Lmalachurum$pop)

dim(Lmalachurum)
dim(Aflavipes)
dim(Aagilissima)

#Ho, Hs and Fis (hierfstat package)
BS_Lmalachurum <- basic.stats(Lmalachurum, diploid = TRUE) 
BS_Aagilissima <- basic.stats(Aagilissima, diploid = TRUE)
BS_Aflavipes <- basic.stats(Aflavipes, diploid = TRUE)

#Aflavipes
mean(BS_Aflavipes$Ho[,1], na.rm = TRUE)
sd(BS_Aflavipes$Ho[,1], na.rm = TRUE)
mean(BS_Aflavipes$Ho[,2], na.rm = TRUE)
sd(BS_Aflavipes$Ho[,2], na.rm = TRUE)
mean(BS_Aflavipes$perloc[,1], na.rm= TRUE)
sd(BS_Aflavipes$perloc[,1], na.rm= TRUE)

mean(BS_Aflavipes$Hs[,1], na.rm = TRUE)
sd(BS_Aflavipes$Hs[,1], na.rm = TRUE)
mean(BS_Aflavipes$Hs[,2], na.rm = TRUE)
sd(BS_Aflavipes$Hs[,2], na.rm = TRUE)
mean(BS_Aflavipes$perloc[,2], na.rm= TRUE)
sd(BS_Aflavipes$perloc[,2], na.rm= TRUE)

mean(BS_Aflavipes$Fis[,1], na.rm = TRUE)
sd(BS_Aflavipes$Fis[,1], na.rm = TRUE)
mean(BS_Aflavipes$Fis[,2], na.rm = TRUE)
sd(BS_Aflavipes$Fis[,2], na.rm = TRUE)
mean(BS_Aflavipes$perloc[,9], na.rm= TRUE)
sd(BS_Aflavipes$perloc[,9], na.rm= TRUE)

#Aagilissima
View(BS_Aagilissima$Ho)
mean(BS_Aagilissima$Ho[,1], na.rm = TRUE)
sd(BS_Aagilissima$Ho[,1], na.rm = TRUE)
mean(BS_Aagilissima$Ho[,2], na.rm = TRUE)
sd(BS_Aagilissima$Ho[,2], na.rm = TRUE)
mean(BS_Aagilissima$Ho[,3], na.rm = TRUE)
sd(BS_Aagilissima$Ho[,3], na.rm = TRUE)
mean(BS_Aagilissima$perloc[,1], na.rm = TRUE)
sd(BS_Aagilissima$perloc[,1], na.rm = TRUE)

mean(BS_Aagilissima$Hs[,1], na.rm = TRUE)
sd(BS_Aagilissima$Hs[,1], na.rm = TRUE)
mean(BS_Aagilissima$Hs[,2], na.rm = TRUE)
sd(BS_Aagilissima$Hs[,2], na.rm = TRUE)
mean(BS_Aagilissima$Ho[,3], na.rm = TRUE)
sd(BS_Aagilissima$Ho[,3], na.rm = TRUE)
mean(BS_Aagilissima$perloc[,2], na.rm=TRUE)
sd(BS_Aagilissima$perloc[,2], na.rm = TRUE)

mean(BS_Aagilissima$Fis[,1], na.rm = TRUE)
sd(BS_Aagilissima$Fis[,1], na.rm = TRUE)
mean(BS_Aagilissima$Fis[,2], na.rm = TRUE)
sd(BS_Aagilissima$Fis[,2], na.rm = TRUE)
mean(BS_Aagilissima$Fis[,3], na.rm = TRUE)
sd(BS_Aagilissima$Fis[,3], na.rm = TRUE)
mean(BS_Aagilissima$perloc[,9], na.rm=TRUE)
sd(BS_Aagilissima$perloc[,9], na.rm = TRUE)

#Lmalachurum
View(BS_Lmalachurum$Ho)
mean(BS_Lmalachurum$Ho[,1], na.rm = TRUE)
sd(BS_Lmalachurum$Ho[,1], na.rm = TRUE)
mean(BS_Lmalachurum$Ho[,2], na.rm = TRUE)
sd(BS_Lmalachurum$Ho[,2], na.rm = TRUE)
mean(BS_Lmalachurum$perloc[,1], na.rm=TRUE)
sd(BS_Lmalachurum$perloc[,1], na.rm=TRUE)

mean(BS_Lmalachurum$Hs[,1], na.rm = TRUE)
sd(BS_Lmalachurum$Hs[,1], na.rm = TRUE)
mean(BS_Lmalachurum$Hs[,2], na.rm = TRUE)
sd(BS_Lmalachurum$Hs[,2], na.rm = TRUE)
mean(BS_Lmalachurum$perloc[,2], na.rm = TRUE)
sd(BS_Lmalachurum$perloc[,2], na.rm=TRUE)

mean(BS_Lmalachurum$Fis[,1], na.rm = TRUE)
sd(BS_Lmalachurum$Fis[,1], na.rm = TRUE)
mean(BS_Lmalachurum$Fis[,2], na.rm = TRUE)
sd(BS_Lmalachurum$Fis[,2], na.rm = TRUE)
mean(BS_Lmalachurum$perloc[,9], na.rm = TRUE)
sd(BS_Lmalachurum$perloc[,9], na.rm=TRUE)

#He (popr package)
Aflavipes_stats <- poppr(Genind_Aflavipes)
Aflavipes_stat$
Aflavipes_stats$Hexp

Aagilissima_stats <- poppr(Genind_Aagilissima)
Aagilissima_stats$Pop
Aagilissima_stats$Hexp

Lmalachurum_stats <- poppr(Genind_Lmalachurum)
Lmalachurum_stats$Pop
Lmalachurum_stats$Hexp


#Fst and Fis totals (hierfstat package)
Fst_Aflavipes <- wc(Aflavipes, diploid = TRUE)
Fst_Aagilissima <- wc(Aagilissima, diploid = TRUE)
Fst_Lmalachurum <- wc(Lmalachurum, diploid = TRUE)


#Pairwise Fst Weir and Cockerham (1984) (hierfstat package)
PairFst_Aflavipes <- pairwise.WCfst(Aflavipes, diploid = TRUE)
PairFst_Aagilissima <- pairwise.WCfst(Aagilissima, diploid = TRUE)
PairFst_Lmalachurum <- pairwise.WCfst(Lmalachurum, diploid = TRUE)

print(PairFst_Aflavipes)
print(PairFst_Aagilissima)
print(PairFst_Lmalachurum)

#Pairwise Gst Hedrick 2005 (mmod package)
Gst_Aflavipes <- pairwise_Gst_Hedrick(Genind_Aflavipes, linearized = FALSE)
Gst_Aagilissima <- pairwise_Gst_Hedrick(Genind_Aagilissima, linearized = FALSE)
Gst_Lmalachurum <- pairwise_Gst_Hedrick(Genind_Lmalachurum, linearized = FALSE)

#Gst Hedrick 2005 and Jost D 2008 (mmod package)
Diff_Aflavipes <- diff_stats(Genind_Aflavipes)
Diff_Aagilissima <- diff_stats(Genind_Aagilissima)
Diff_Lmalachurum <- diff_stats(Genind_Lmalachurum)


#Pairwise Jost's D 2008 (mmod package)
Jost_Aflavipes <- pairwise_D(Genind_Aflavipes)
Jost_Aagilissima <- pairwise_D(Genind_Aagilissima)
Jost_Lmalachurum <- pairwise_D(Genind_Lmalachurum)
