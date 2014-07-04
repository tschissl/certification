library(ggplot2)
# library(devtools) #needed for "order_by"
# install_github("plotflow", "trinker")
library(plotflow) #needed for "order_by" # install_github("plotflow", "trinker")
library(metRology) #for mandel k and h calculations
library(ape) # needed for varcomp (variance component) extraction
library(nlme) # needed for lme
require(plyr)

GAS <- read.csv("~/GitHub/certification/GAS.csv", sep=";")
OKUM <- read.csv("~/GitHub/certification/OKUM.csv", sep=";")
MUH <- read.csv("~/GitHub/certification/MUH.csv", sep=";")
OKUM.methods <- read.csv("~/GitHub/certification/OKUM methods.csv", sep=";")
GOM <- data.frame(OKUM, MUH[6:60], GAS[6:60])
GOM <- merge(GOM, OKUM.methods) # merging all data and all methods
#defining factors 
GOM$Lab <- as.factor(GOM$Lab)
GOM$Packet <- as.factor(GOM$Packet)
GOM$Prep <- as.factor(GOM$Prep)
GOM$day <- as.factor(GOM$day)
GOM$names <- as.factor(GOM$names)

GOMorig <- GOM  #keeping the original file

mytheme <- theme_grey() + theme(plot.title = element_text(colour = "black", size = rel(2))) + theme(axis.title.x = element_text(size = rel(1.8)))+ theme(axis.title.y = element_text(size = rel(1.8))) +theme(axis.text.x = element_text(size = rel(1.5))) + theme(axis.text.y = element_text(size = rel(1.5))) + theme(legend.title = element_text(size = rel(1))) +theme(legend.text = element_text(size = rel(0.8)))

measurand <- 'SiO2'

mean <- mean(tapply(GOM[[measurand]], GOM$Lab, mean, na.rm=TRUE), na.rm=TRUE)
median <- median(tapply(GOM[[measurand]], GOM$Lab, median, na.rm=TRUE), na.rm=TRUE)

(test <- tapply(GOM[[measurand]], GOM$Packet, mean, na.rm=TRUE)) #only gives the mean of each packet/all not per packet/lab

avgbymethod <- function(x) {
  element <- x
  '%p%' <- function(x, y) {as.character(paste (x, y, sep =""))}
  prep <-  'Prep.'
  method <- 'Method.'
  anal.prep <- prep %p% element
  anal.method <- method %p% element
  anal <- GOM[[element]]
  anal.prep <- GOM[[anal.prep]]
  anal.method <- GOM[[anal.method]]
  analyte <- data.frame(GOM$Lab, GOM$names, anal, anal.prep, anal.method )
  analyte <- na.omit(analyte)
  bymethod <- tapply(analyte$anal, analyte$anal.method, mean, na.rm=TRUE)
  bymethod.n <- ddply(analyte, c("anal.method"), summarise, N=length(anal))
  print(bymethod)
  avg <- mean(tapply(analyte$anal, analyte$GOM.Lab, mean, na.rm=TRUE), na.rm=TRUE)
  median <- median(tapply(analyte$anal, analyte$GOM.Lab, median, na.rm=TRUE), na.rm=TRUE)
  cat("mean of lab means: ","\n")
  cat(round(avg,3),"g/100g", "\n")
  cat("median of lab medians: ","\n")
  cat(round(median,3),"g/100g", "\n")
}

###################### metRology #########################################################################
############ k plot #########
#barplot.mandel.k example for SiO2
k <- with(GOM, mandel.kh(GOM[[measurand]], g=Lab, type="k"))
# k <- with(SBC, mandel.kh(SBC[6:9], g=names, type="k", method="robust")) # robust version       
barplot(k, las=2, col=1:4)
#legend("topleft", legend=names(k), fill=1:8, cex=0.7, bg="white")
#as boxplot
# boxplot(k, las=2, col=1:4)

############ h plot #########
#barplot.mandel.h example for SiO2, TiO2, Al2O3, Fe2O3T
h <- with(GOM, mandel.kh(GOM[6:9], g=names, type="h"))
# h <- with(SBC, mandel.kh(SBC[6:9], g=names, type="h", method="robust")) # robust version
barplot(h, las=2, col=1:4) #without colours barplot(h, las=2)
legend("bottomleft", legend=names(h), fill=1:8, cex=0.7, bg="white")


##### summary statistic with package "pastec" #######
sd <- stat.desc(SBC$SiO2, basic = TRUE, norm = TRUE)
sd <- round(sd, 6)


#### Removing ouliers lab #16 based on Youden plot

outlier <- c(16, 12, 7, 30)
leng <- length(outlier)
for(i in 1:leng)
{
  GOM[[measurand]] <- ifelse(GOM$Lab==outlier[i], NA, GOM[[measurand]])
  print(summary(GOM[[measurand]], na.rm=TRUE, digits=4))
  #message(outlier[i])
}


avgbymethod(measurand)

## creating database with median for each lab, analyte and RM
## example works but statistically not correct as three packages
library(plyr) 
medianLab <- function(x) median(x, na.rm=TRUE)
medianGOM.lab <- ddply(GOM, .(Lab), numcolwise(medianLab))
### for single analyte only
medianSiO2.lab <- ddply(GOM, c("Lab"), function(df)c(median(df$SiO2, na.rm=TRUE))

## median over packets within lab fits perfectly with EXCEL calculation
medianGOM.packet <- ddply(GOM, c("Lab", "Packet"), numcolwise(medianLab))
### for single analyte only
medianSiO2.packet <- ddply(GOM, c("Lab", "Packet"), function(df) median(df$SiO2, na.rm=TRUE))

## median over median of packets within lab fits perfectly with EXCEL calculation
medianGOM.packet.lab <- ddply(medianGOM.packet, c("Lab"), numcolwise(medianLab))

## median over all lab medians, does not work properly
medianGOM.all <- data.frame(apply(medianGOM.packet.lab, 2, medianLab))
median.measurand.all <- median(medianGOM.packet[[measurand]], na.rm=TRUE)

## making mandel plots according to suggestion of Steve Ellison
k <- with(medianGOM.packet, mandel.kh(medianGOM.packet[[measurand]], g=Lab, type="k"))
barplot(k, las=2, col=1:4)

measurand <- "SiO2"
library(ape) # needed for varcomp (variance component) extraction
library(nlme) # needed for lme
x <- 'La'
VC.lme <- function(x){ # variance compenent extraction of one level only
  measurand <- x
  medianGOM.packet$Lab <- as.factor(medianGOM.packet$Lab) # using only the median of the 3 packages per lab
  anal <- medianGOM.packet[[measurand]]
  DF.lme <- data.frame(medianGOM.packet$Lab, medianGOM.packet[[measurand]])
  DF.lme <- na.omit(DF.lme)
  names(DF.lme) <- c("Lab", "measurand")
  GOM.lme <- lme(measurand ~ 1, random = ~ 1|Lab, data=DF.lme) # linear model with random effects
  sL2 <- varcomp(GOM.lme, FALSE, FALSE)[[1]] # between-laboratory variance
  sr2 <- varcomp(GOM.lme, FALSE, FALSE)[[2]] # repeatability standard deviation
  n.p <- dim(DF.lme)[1] # number of obervations
  p <- length(unique(DF.lme$Lab)) # haven't found a better way how to extract the number of labs (number of groups)
  u1 <- sqrt(sL2/p+sr2/n.p) # calculating the standard uncertainty of characterization
  u2 <- attr(GOM.lme$fixDF,"varFixFact") # gives the same results as u, amazing!
  cat("sL^2 VC Lab", "\n"); print(sL2);print(n.p);print(p)
  cat("sr^2 VC repeatability", "\n"); print(sr2)
  cat("standard uncertainty1", "standard uncertainty2", "\n")
  print(round(u1,5))
  print(round(u2,5))
  plot(DF.lme)
}

VC.lme2 <- function(x){ # variance compenent extraction of two levels (lab and packet)
  measurand <- x
  medianGOM.packet$Lab <- as.factor(medianGOM.packet$Lab) # using only the median of the 3 packages per lab
  medianGOM.packet$Packet <- as.factor(medianGOM.packet$Packet)
  anal <- medianGOM.packet[[measurand]]
  DF.lme <- data.frame(medianGOM.packet$Lab, medianGOM.packet$Packet, medianGOM.packet[[measurand]])
  DF.lme <- na.omit(DF.lme)
  names(DF.lme) <- c("Lab", "Packet", "measurand")
  GOM.lme <- lme(measurand ~ 1, random = ~ 1|Lab/Packet, data=DF.lme) # linear model with random effects
  sL2 <- varcomp(GOM.lme, FALSE, FALSE)[[1]] # between-laboratory variance
  sbb2 <- varcomp(GOM.lme, FALSE, FALSE)[[2]] # between bottle standard deviation
  sr2 <- varcomp(GOM.lme, FALSE, FALSE)[[3]] # repeatability standard deviation
  n.p <- dim(DF.lme)[1] # number of obervations
  p <- length(unique(DF.lme$Lab)) # haven't found a better way how to extract the number of labs (number of groups)
  r <- length(unique(DF.lme$Packet))
  u1 <- sqrt(sL2/p+sbb2/p/r+sr2/p/r/4) # calculating the standard uncertainty of characterization
  u2 <- attr(GOM.lme$fixDF,"varFixFact") # gives the same results as u, amazing!
  cat("sL^2 VC Lab", "\n"); print(sL2);print(n.p);print(p)
  cat("sr^2 VC repeatability", "\n"); print(sr2)
  cat("standard uncertainty1", "standard uncertainty2", "\n")
  print(round(u1,5))
  print(round(u2,5))
  plot(DF.lme)
}

# just a test to compare Statgraphics results with lme results
OKUM.SiO2.7 <- read.csv("C:/Daten/projects/IAG/Certification Committee/Reference materials/OKUM and MUH-1/OKUM calulcations/OKUM.SiO2.7.csv", sep=";")
VC.lme2 <- function(x){ # variance compenent extraction of two levels (lab and packet)
  measurand <- x
  anal <- OKUM.SiO2.7[[measurand]]
  DF.lme <- data.frame(OKUM.SiO2.7$Lab, OKUM.SiO2.7$Packet, OKUM.SiO2.7$Prep, OKUM.SiO2.7[[measurand]])
  DF.lme <- na.omit(DF.lme)
  names(DF.lme) <- c("Lab", "Packet", "Prep", "measurand")
  GOM.lme <- lme(measurand ~ 1, random = ~ 1|Lab/Packet, data=DF.lme) # linear model with random effects
  sL2 <- varcomp(GOM.lme, FALSE, FALSE)[[1]] # between-laboratory variance
  sbb2 <- varcomp(GOM.lme, FALSE, FALSE)[[2]] # between bottle standard deviation
  sr2 <- varcomp(GOM.lme, FALSE, FALSE)[[3]] # repeatability standard deviation
  n.p <- dim(DF.lme)[1] # number of obervations
  p <- length(unique(DF.lme$Lab)) # haven't found a better way how to extract the number of labs (number of groups)
  q <- length(unique(DF.lme$Packet))
  r <- length(unique(DF.lme$Prep))
  u1 <- sqrt(sL2/p+sbb2/p/q+sr2/p/3/4)
  #u1 <- sqrt(sL2/p+sbb2/3/p+sr2/p/3/4) # calculating the standard uncertainty of characterization
  u2 <- attr(GOM.lme$fixDF,"varFixFact") # gives the same results as u, amazing!
  cat("sL^2 VC Lab", "\n"); print(sL2);print(n.p);print(p)
  cat("sr^2 VC repeatability", "\n"); print(sr2)
  cat("standard uncertainty1", "standard uncertainty2", "\n")
  print(round(u1,5))
  print(round(u2,5))
  plot(DF.lme)
}

qqnorm(medianGOM.packet.after$SiO2)

require(nortest)
shapiro.test(GOM.median.after[[measurand]])

## means over packets within lab 
meanGOM <- function(x) mean(x, na.rm=TRUE)
meanGOM.packet <- ddply(GOM, c("Lab", "Packet"), numcolwise(meanGOM))
## median over median of packets within lab
GOM.mean <- ddply(meanGOM.packet, c("Lab"), numcolwise(meanGOM))
## median over packets within lab 
require(plyr)
medianGOM <- function(x) median(x, na.rm=TRUE)
sdGOM <- function(x) sd(x, na.rm=TRUE)
medianGOM.packet <- ddply(GOM, c("Lab", "Packet"), numcolwise(medianGOM))
## median over median of packets within lab
GOM.median <- ddply(medianGOM.packet, c("Lab"), numcolwise(medianGOM))
GOM.sd <- ddply(medianGOM.packet, c("Lab"), numcolwise(sdGOM))

plot_youd <- function(a, y, z) { # a = measurand, y and z ref
  '%p%' <- function(a, b) {as.character(paste (a, b, sep =""))}
  switch(
    y,
    GAS = rm1 <- 2,
    MUH = rm1 <- 1,
    OKUM = rm1 <- 0
  )
  switch(
    z,
    GAS = rm2 <- 2,
    MUH = rm2 <- 1,
    OKUM = rm2 <- 0
  )
  if(rm1 > 0) 
  {RM1 <- a %p% '.' %p% rm1
  } else 
  {
    RM1 <- a
  }
  if(rm2 > 0) 
  {RM2 <- a %p% '.' %p% rm2
  } else 
  {
    RM2 <- a
  }
  
  RM1.s <- GOM.sd[[RM1]]/sd(GOM.median[[RM1]], na.rm = TRUE)
  RM1 <- GOM.median[[RM1]]
  RM1 <- (RM1-median(RM1, na.rm = TRUE))/sd(RM1, na.rm = TRUE) #calculating z-scores
  RM2.s <- GOM.sd[[RM2]]/sd(GOM.median[[RM2]], na.rm = TRUE)
  RM2 <- GOM.median[[RM2]]
  RM2 <- (RM2-median(RM2, na.rm = TRUE))/sd(RM2, na.rm = TRUE) #calculating z-scores
  RM <- data.frame(GOM.median$Lab, RM1, RM2, RM1.s, RM2.s)
  RM <- na.omit(RM)
  p <- ggplot(RM, aes(RM1, RM2, label = GOM.median.Lab))  + xlim(-5, 5) + ylim(-5,5) + geom_point(aes(colour=factor(GOM.median.Lab)), size = 4) + geom_errorbar(aes(ymin=RM2-RM2.s, ymax=RM2+RM2.s)) + geom_errorbarh(aes(xmin=RM1-RM1.s, xmax=RM1+RM1.s))
  p <- p + xlab(y) + ylab(z) + labs(title = a) + labs(colour = "GOM.median.Lab") + mytheme + 
    geom_abline(intercept = 0, slope = 1) +  geom_abline(intercept = 2.8284, slope = 1) +  
    geom_abline(intercept = - 2.8284, slope = 1) +  geom_abline(intercept = 2.8284, slope = - 1) +  
    geom_abline(intercept = -2.8284, slope = - 1) + geom_text(data = NULL, x = -0.5, y = 2.85, label = "z = 2") + 
    geom_text(aes(colour=factor(GOM.median.Lab)), hjust=1, vjust=0) 
  p + theme(legend.position="none") #+ stat_density2d(aes(fill = ..level..), geom="polygon") 
}

##test
y <- 'OKUM'
z <- 'MUH'
a <- 'SiO2'
RM1 <- GOM.median$SiO2.1
RM1.s <- GOM.sd$SiO2.1/sd(GOM.median$SiO2.1, na.rm = TRUE)
RM1 <- (RM1-median(RM1, na.rm = TRUE))/sd(RM1, na.rm = TRUE) #calculating z-scores
RM2 <- GOM.median$SiO2
RM2.s <- GOM.sd$SiO2/sd(GOM.median$SiO2, na.rm = TRUE)
RM2 <- (RM2-median(RM2, na.rm = TRUE))/sd(RM2, na.rm = TRUE) #calculating z-scores
RM <- data.frame(GOM.median$Lab, RM1, RM2, RM1.s, RM2.s)
RM <- na.omit(RM)
p <- ggplot(RM, aes(RM1, RM2, label = GOM.median.Lab))  + xlim(-5, 5) + ylim(-5,5) + geom_point(aes(colour=factor(GOM.median.Lab)), size = 4)
p <- p + xlab(y) + ylab(z) + labs(title = a) + labs(colour = "GOM.median.Lab") + mytheme + geom_abline(intercept = 0, slope = 1) +  geom_abline(intercept = 2.8284, slope = 1) +  geom_abline(intercept = - 2.8284, slope = 1) +  geom_abline(intercept = 2.8284, slope = - 1) +  geom_abline(intercept = -2.8284, slope = - 1) + geom_text(data = NULL, x = -0.5, y = 2.85, label = "z = 2") + geom_text(aes(colour=factor(GOM.median.Lab)), hjust=1, vjust=0)
p + theme(legend.position="none") #+ stat_density2d(aes(fill = ..level..), geom="polygon") 
p + geom_errorbar(aes(ymin=RM2 - RM2.s, ymax=RM2 + RM2.s)) + geom_errorbarh(aes(xmin=RM1-RM1.s, xmax=RM1+RM1.s))

plot_lab <- function(x, type) {
  element <- x
  '%p%' <- function(x, y) {as.character(paste (x, y, sep =""))}
  prep <-  'Prep.'
  method <- 'Method.'
  anal.prep <- prep %p% measur
  anal.method <- method %p% measur
  anal <- GOM[[element]]
  anal.prep <- GOM[[anal.prep]]
  anal.method <- GOM[[anal.method]]
  analyte <- data.frame(GOM$Lab, GOM$names, anal, anal.prep, anal.method )
  analyte <- na.omit(analyte)
  avg <- median(tapply(analyte$anal, analyte$GOM.Lab, mean, na.rm=TRUE), na.rm=TRUE)
  switch(type, 
         M = hor <- 1*(0.01*avg)^0.8495, 
         T = hor <- 10000*1*(0.01*avg/10000)^0.8495)
  # hor <- 1*(0.01*avg)^0.8495
  u.hor <- avg + hor
  u.2hor <- avg + 2*hor
  l.hor <- avg - hor
  l.2hor <- avg - 2*hor
  p <- ggplot(reorder_by(GOM.Lab, ~anal, analyte, mean), aes(GOM.Lab, anal)) + geom_abline(intercept = avg, slope = 0) + geom_abline(intercept = u.2hor, slope = 0, linetype ="dotted") + geom_abline(intercept = l.2hor, slope = 0, linetype = "dotted") + geom_abline(intercept = u.hor, slope = 0, linetype ="dashed") + geom_abline(intercept = l.hor, slope = 0, linetype = "dashed")
  p <- p + geom_boxplot(aes(fill=anal.method)) + geom_point(size=4) + geom_point(aes(colour=anal.prep), size = 3.5) 
  p + xlab("lab") + ylab(unit) + labs(title = measur) + labs(colour = "Prep") + geom_smooth(aes(group=1)) + mytheme
}

sequence <- seq(from = 1, to = length(names(OKUM.outlier)), by = 2)
# col <- OKUM.outlier[,c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,57,59,61,63,65,67,69,71,73,75,77,79,81,83,85,87,79,91,93,95,97,99,101,103,105,107,109)]
col <- OKUM.outlier[,c(sequence)]
col.names <- colnames(col)
for (m in col.names) print(col.names)

