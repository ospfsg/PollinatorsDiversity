#Aagilissima
setwd("/home/paulosousa/Working/Bees/Analysis/Aagilissima/pixy/")

Aagilissima_pi <- read.csv("pixy_pi.txt", sep="\t")
mean(Aagilissima_pi$avg_pi)
sd(Aagilissima_pi$avg_pi)

pi_plot <- ggplot(Aagilissima_pi, aes(x=pop, y=avg_pi)) +
  geom_boxplot(fill='#A4A4A4', color="black", outlier.colour="red", outlier.shape=8, outlier.size=4)+
  labs(title="Nucleotide diversity of three A. agilissima populations",x="Population", y = "Average per site nucleotide diversity")+
  theme_classic()
pi_plot <- pi_plot  + theme(legend.position= c(0.8, 0.5), legend.background = element_rect(size=0.5, linetype = "solid", colour = "black"))
pi_plot

Bet <- subset(Aagilissima_pi, subset = c(pop == "SerrasBeticas"))
mean(Bet$avg_pi)
sd(Bet$avg_pi)
West <- subset(Aagilissima_pi, subset = c(pop == "West"))
mean(West$avg_pi)
sd(West$avg_pi)
Ebro <- subset(Aagilissima_pi, subset = c(pop == "Ebro"))
mean(Ebro$avg_pi)
sd(Ebro$avg_pi)

Aagilissima_dxy <- read.csv("pixy_dxy.txt", sep="\t")
mean(Aagilissima_dxy$avg_dxy)
sd(Aagilissima_dxy$avg_dxy)

library(ggplot2)
Bet_Ebro <- subset(Aagilissima_dxy, subset = c(pop1 =="SerrasBeticas" & pop2== "Ebro"))
mean(Bet_Ebro$avg_dxy)
sd(Bet_Ebro$avg_dxy)

dxy_plot_BetEbro <-qplot(Bet_Ebro$avg_dxy, geom = "histogram", 
                         main= "Nucleotide divergence between populations of A. agilissima (Serras Beticas and Ebro)", 
                         xlab= "Average per site nucleotide divergence", ylab = "Count",
                         fill=I("lightblue"),
                         col=I("black"))
dxy_plot_BetEbro <- dxy_plot_BetEbro + theme_gray(base_size = 14)
dxy_plot_BetEbro


Bet_West <- subset(Aagilissima_dxy, subset = c(pop1 =="SerrasBeticas" & pop2== "West"))
mean(Bet_West$avg_dxy)
sd(Bet_West$avg_dxy)

dxy_plot_BetWest <-qplot(Bet_West$avg_dxy, geom = "histogram", 
                         main= "Nucleotide divergence between populations of A. agilissima (Serras Beticas and West)", 
                         xlab= "Average per site nucleotide divergence", ylab = "Count",
                         fill=I("lightblue"),
                         col=I("black"))
dxy_plot_BetWest <- dxy_plot_BetWest + theme_gray(base_size = 14)
dxy_plot_BetWest

Ebro_West <- subset(Aagilissima_dxy, subset = c(pop1 =="West" & pop2== "Ebro"))
mean(Ebro_West$avg_dxy)
sd(Ebro_West$avg_dxy)

dxy_plot_WestEbro <-qplot(Ebro_West$avg_dxy, geom = "histogram", 
                          main= "Nucleotide divergence between populations of A. agilissima (West and Ebro)", 
                          xlab= "Average per site nucleotide divergence", ylab = "Count",
                          fill=I("lightblue"),
                          col=I("black"))
dxy_plot_WestEbro <- dxy_plot_WestEbro + theme_gray(base_size = 14)
dxy_plot_WestEbro
###########################################################################################################################
#Aflavipes
setwd("/home/paulosousa/Working/Bees/Analysis/Aflavipes/pixy/")

Aflavipes_pi <- read.csv("pixy_pi.txt", sep="\t")
mean(Aflavipes_pi$avg_pi, na.rm = T)
sd(Aflavipes_pi$avg_pi, na.rm = T)

pi_plot <- ggplot(Aflavipes_pi, aes(x=pop, y=avg_pi)) +
  geom_boxplot(fill='#A4A4A4', color="black", outlier.colour="red", outlier.shape=8, outlier.size=4)+
  labs(title="Nucleotide diversity of three A. agilissima populations",x="Population", y = "Average per site nucleotide diversity")+
  theme_classic()
pi_plot <- pi_plot  + theme(legend.position= c(0.8, 0.5), legend.background = element_rect(size=0.5, linetype = "solid", colour = "black"))
pi_plot

Iberia <- subset(Aflavipes_pi, subset = c(pop == "Iberia"))
mean(Iberia$avg_pi)
sd(Iberia$avg_pi)
Pyrenees <- subset(Aflavipes_pi, subset = c(pop == "Pyrenees"))
mean(Pyrenees$avg_pi, na.rm = T)
sd(Pyrenees$avg_pi, na.rm = T)


Aflavipes_dxy <- read.csv("pixy_dxy.txt", sep="\t")
mean(Aflavipes_dxy$avg_dxy, na.rm = T)
sd(Aflavipes_dxy$avg_dxy, na.rm = T)

dxy_plot <-qplot(Aflavipes_dxy$avg_dxy, geom = "histogram", 
                 main= "Nucleotide divergence between populations of A. flavipes", 
                 xlab= "Average per site nucleotide divergence", ylab = "Count",
                 fill=I("lightblue"),
                 col=I("black"))
dxy_plot <- dxy_plot + theme_gray(base_size = 14)
dxy_plot
#############################################################################################
#Lmalachurum

setwd("/home/paulosousa/Working/Bees/Analysis/Lmalachurum/pixy/")

Lmalachurum_pi <- read.csv("pixy_pi.txt", sep="\t")
mean(Lmalachurum_pi$avg_pi, na.rm = T)
sd(Lmalachurum_pi$avg_pi, na.rm = T)

pi_plot <- ggplot(Lmalachurum_pi, aes(x=pop, y=avg_pi)) +
  geom_boxplot(fill='#A4A4A4', color="black", outlier.colour="red", outlier.shape=8, outlier.size=4)+
  labs(title="Nucleotide diversity of three A. agilissima populations",x="Population", y = "Average per site nucleotide diversity")+
  theme_classic()
pi_plot <- pi_plot  + theme(legend.position= c(0.8, 0.5), legend.background = element_rect(size=0.5, linetype = "solid", colour = "black"))
pi_plot

South <- subset(Lmalachurum_pi, subset = c(pop == "South"))
mean(South$avg_pi)
sd(South$avg_pi)
Ebro <- subset(Lmalachurum_pi, subset = c(pop == "Ebro"))
mean(Ebro$avg_pi, na.rm = T)
sd(Ebro$avg_pi, na.rm = T)


Lmalachurum_dxy <- read.csv("pixy_dxy.txt", sep="\t")
mean(Lmalachurum_dxy$avg_dxy, na.rm = T)
sd(Lmalachurum_dxy$avg_dxy, na.rm = T)

dxy_plot <-qplot(Lmalachurum_dxy$avg_dxy, geom = "histogram", 
                 main= "Nucleotide divergence between populations of A. flavipes", 
                 xlab= "Average per site nucleotide divergence", ylab = "Count",
                 fill=I("lightblue"),
                 col=I("black"))
dxy_plot <- dxy_plot + theme_gray(base_size = 14)
dxy_plot