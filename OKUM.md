Analyzing submitted data for OKUM
========================================================


```r
measurand <- "TiO2"
```


## measurand selected TiO2


#### Importing the data and assigning factors


```r
library(ggplot2)
# library(devtools) #needed for 'order_by' install_github('plotflow',
# 'trinker')
library(plotflow)  #needed for 'order_by' # install_github('plotflow', 'trinker')
```

```
## Loading required package: gridExtra
## Loading required package: grid
```

```r
library(metRology)  #for mandel k and h calculations
```

```
## 
## Attaching package: 'metRology'
## 
## Die folgenden Objekte sind maskiert from 'package:base':
## 
##     cbind, rbind
```

```r
library(ape)  # needed for varcomp (variance component) extraction
library(nlme)  # needed for lme
require(plyr)
```

```
## Loading required package: plyr
```

```r

GAS <- read.csv("~/GitHub/certification/GAS.csv", sep = ";")
OKUM <- read.csv("~/GitHub/certification/OKUM.csv", sep = ";")
MUH <- read.csv("~/GitHub/certification/MUH.csv", sep = ";")
OKUM.methods <- read.csv("~/GitHub/certification/OKUM methods.csv", sep = ";")
GOM <- data.frame(OKUM, MUH[6:60], GAS[6:60])
GOM <- merge(GOM, OKUM.methods)  # merging all data and all methods
# defining factors
GOM$Lab <- as.factor(GOM$Lab)
GOM$Packet <- as.factor(GOM$Packet)
GOM$Prep <- as.factor(GOM$Prep)
GOM$day <- as.factor(GOM$day)
GOM$names <- as.factor(GOM$names)
GOMorig <- GOM  #keeping the original file
```


#### defining the plotting theme


```r
mytheme <- theme_grey() + theme(plot.title = element_text(colour = "black", 
    size = rel(2))) + theme(axis.title.x = element_text(size = rel(1.8))) + 
    theme(axis.title.y = element_text(size = rel(1.8))) + theme(axis.text.x = element_text(size = rel(1.5))) + 
    theme(axis.text.y = element_text(size = rel(1.5))) + theme(legend.title = element_text(size = rel(1))) + 
    theme(legend.text = element_text(size = rel(0.8)))
```

#### initial calculations with complete data set


```r
## means over packets within lab
meanGOM <- function(x) mean(x, na.rm = TRUE)
sdGOM <- function(x) sd(x, na.rm = TRUE)
meanGOM.packet <- ddply(GOM, c("Lab", "Packet"), numcolwise(meanGOM))
## median over median of packets within lab
GOM.mean <- ddply(meanGOM.packet, c("Lab"), numcolwise(meanGOM))
```


```r
## median over packets within lab
require(plyr)
medianGOM <- function(x) median(x, na.rm = TRUE)
medianGOM.packet <- ddply(GOM, c("Lab", "Packet"), numcolwise(medianGOM))
GOM.sd <- ddply(medianGOM.packet, c("Lab"), numcolwise(sdGOM))
## median over median of packets within lab
GOM.median <- ddply(medianGOM.packet, c("Lab"), numcolwise(medianGOM))
```

#### defining the function for plotting methods vs. measurand mass fraction. Sample preparation methods are also marked in the plot.

```r
### function plot_method defined here. Automated plot of methods boxplots
### sorted by increasing method mean. Enter with plot_method=('xx') #######
plot_method <- function(x) {
    element <- x
    "%p%" <- function(x, y) {
        as.character(paste(x, y, sep = ""))
    }
    prep <- "Prep."
    method <- "Method."
    anal.prep <- prep %p% element
    anal.method <- method %p% element
    anal <- GOM[[element]]
    anal.prep <- GOM[[anal.prep]]
    anal.method <- GOM[[anal.method]]
    analyte <- data.frame(GOM$Lab, GOM$names, anal, anal.prep, anal.method)
    analyte <- na.omit(analyte)
    p <- ggplot(reorder_by(anal.method, ~anal, analyte, mean), aes(anal.method, 
        anal)) + geom_boxplot(outlier.shape = NA) + geom_jitter(aes(colour = anal.prep), 
        size = 4.5, position = position_jitter(width = 0.2))
    p + xlab("method") + ylab("mg/kg or g/100g") + labs(title = element) + labs(colour = "Prep") + 
        mytheme
}
```


```r
###### function plot_lab #### defined here. Automated plot of lab boxplots sorted
###### by increasing lab mean and Horwitz function. Enter with plot_lab=('xx')
plot_lab <- function(x, type) {
    element <- x
    "%p%" <- function(x, y) {
        as.character(paste(x, y, sep = ""))
    }
    prep <- "Prep."
    method <- "Method."
    anal.prep <- prep %p% element
    anal.method <- method %p% element
    anal <- GOM[[element]]
    anal.prep <- GOM[[anal.prep]]
    anal.method <- GOM[[anal.method]]
    analyte <- data.frame(GOM$Lab, GOM$names, anal, anal.prep, anal.method)
    analyte <- na.omit(analyte)
    avg <- median(tapply(analyte$anal, analyte$GOM.Lab, mean, na.rm = TRUE), 
        na.rm = TRUE)
    switch(type, M = hor <- 1 * (0.01 * avg)^0.8495, T = hor <- 10000 * 1 * 
        (0.01 * avg/10000)^0.8495)
    # hor <- 1*(0.01*avg)^0.8495
    u.hor <- avg + hor
    u.2hor <- avg + 2 * hor
    l.hor <- avg - hor
    l.2hor <- avg - 2 * hor
    p <- ggplot(reorder_by(GOM.Lab, ~anal, analyte, mean), aes(GOM.Lab, anal)) + 
        geom_abline(intercept = avg, slope = 0) + geom_abline(intercept = u.2hor, 
        slope = 0, linetype = "dotted") + geom_abline(intercept = l.2hor, slope = 0, 
        linetype = "dotted") + geom_abline(intercept = u.hor, slope = 0, linetype = "dashed") + 
        geom_abline(intercept = l.hor, slope = 0, linetype = "dashed")
    p <- p + geom_boxplot(aes(fill = anal.method)) + geom_point(size = 4) + 
        geom_point(aes(colour = anal.prep), size = 3.5)
    p + xlab("lab") + ylab("mg/kg or g/100g") + labs(title = element) + labs(colour = "Prep") + 
        geom_smooth(aes(group = 1)) + mytheme
}
```


#### defining the function of Youden plots


```r
plot_youd <- function(a, y, z) {
    # a = measurand, y and z ref
    "%p%" <- function(a, b) {
        as.character(paste(a, b, sep = ""))
    }
    switch(y, GAS = rm1 <- 2, MUH = rm1 <- 1, OKUM = rm1 <- 0)
    switch(z, GAS = rm2 <- 2, MUH = rm2 <- 1, OKUM = rm2 <- 0)
    if (rm1 > 0) {
        RM1 <- a %p% "." %p% rm1
    } else {
        RM1 <- a
    }
    if (rm2 > 0) {
        RM2 <- a %p% "." %p% rm2
    } else {
        RM2 <- a
    }
    
    RM1.s <- GOM.sd[[RM1]]/sd(GOM.median[[RM1]], na.rm = TRUE)  # calculating the normalised standard deviations
    RM1 <- GOM.median[[RM1]]
    RM1 <- (RM1 - median(RM1, na.rm = TRUE))/sd(RM1, na.rm = TRUE)  #calculating z-scores
    RM2.s <- GOM.sd[[RM2]]/sd(GOM.median[[RM2]], na.rm = TRUE)  # calculating the normalised standard deviations
    RM2 <- GOM.median[[RM2]]
    RM2 <- (RM2 - median(RM2, na.rm = TRUE))/sd(RM2, na.rm = TRUE)  #calculating z-scores
    RM <- data.frame(GOM.median$Lab, RM1, RM2, RM1.s, RM2.s)  # creating a data frame for the measurand
    RM <- na.omit(RM)  # removing all 'NA'
    p <- ggplot(RM, aes(RM1, RM2, label = GOM.median.Lab)) + xlim(-5, 5) + ylim(-5, 
        5) + geom_point(aes(colour = factor(GOM.median.Lab)), size = 4) + geom_errorbar(aes(ymin = RM2 - 
        RM2.s, ymax = RM2 + RM2.s)) + geom_errorbarh(aes(xmin = RM1 - RM1.s, 
        xmax = RM1 + RM1.s))
    p <- p + xlab(y) + ylab(z) + labs(title = a) + labs(colour = "GOM.median.Lab") + 
        mytheme + geom_abline(intercept = 0, slope = 1) + geom_abline(intercept = 2.8284, 
        slope = 1) + geom_abline(intercept = -2.8284, slope = 1) + geom_abline(intercept = 2.8284, 
        slope = -1) + geom_abline(intercept = -2.8284, slope = -1) + geom_text(data = NULL, 
        x = -0.5, y = 2.85, label = "z = 2") + geom_text(aes(colour = factor(GOM.median.Lab)), 
        hjust = 1, vjust = 0)
    p + theme(legend.position = "none")  #+ stat_density2d(aes(fill = ..level..), geom='polygon') 
}
```



```r
summary(GOM[[measurand]], na.rm = TRUE, digits = 4)  # with values without outlier removal
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
##    0.35    0.37    0.38    0.39    0.39    0.46     126
```

```r
mean.before <- mean(GOM.mean[[measurand]], na.rm = TRUE)
median.before <- median(GOM.median[[measurand]], na.rm = TRUE)
```

The TiO2 mean of the Lab means is 0.3858 g/100 g  
The TiO2 median of the Lab+Package medians is 0.38 g/100 g  



```r
plot_method(measurand)
```

```
## Warning: Removed 12 rows containing missing values (geom_point).
```

![plot of chunk methods and lab plot as is](figure/methods_and_lab_plot_as_is1.png) 

```r
plot_lab(measurand, "M")  ## M for majors, T for traces
```

```
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
```

![plot of chunk methods and lab plot as is](figure/methods_and_lab_plot_as_is2.png) 


```r
plot_youd(measurand, "GAS", "OKUM")
```

```
## Error: 'data' muss vom Typ vector sein, war 'NULL'
```

```r
plot_youd(measurand, "MUH", "OKUM")
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3.png) 

#### Mandel k barplot displays the within lab performance relative to all participating labs using the median over packages

```r
k <- with(medianGOM.packet, mandel.kh(medianGOM.packet[[measurand]], g = Lab, 
    type = "k"))
barplot(k, las = 2, col = 1:4)
```

![plot of chunk mandel k](figure/mandel_k.png) 


#### Removing ouliers lab based on Youden plot and Mandel's k (lab performance)


```r
outlier <- c(10, 26, 31) ## defining the outlying lab here with lab#
leng <- length(outlier) ## counting the number of outliers for loop
for(i in seq(leng)) ##  looping
{
  GOM[[measurand]] <- ifelse(GOM$Lab==outlier[i], NA, GOM[[measurand]]) ## replacing values of outlying lab with "NA" and defining new GOM
  message("Lab ", outlier[i], " was removed")
  print(summary(GOM[[measurand]], na.rm=TRUE, digits=4))
  }
```

```
## Lab 10 was removed
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
##    0.35    0.37    0.38    0.39    0.39    0.46     138
```

```
## Lab 26 was removed
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
##    0.35    0.37    0.38    0.39    0.39    0.46     150
```

```
## Lab 31 was removed
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
##    0.35    0.37    0.38    0.39    0.39    0.46     162
```

```r
plot_lab(measurand, 'T') # replotting without outlier, M for majors and T for traces
```

```
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
```

![plot of chunk outlier removal](figure/outlier_removal.png) 

## Calculating mean and median of property value


```r
summary(GOM[[measurand]], na.rm = TRUE, digits = 4)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
##    0.35    0.37    0.38    0.39    0.39    0.46     162
```

```r
mean <- mean(tapply(GOM[[measurand]], GOM$Lab, mean, na.rm = TRUE), na.rm = TRUE)  # mean w/o outlier removal
median.before <- median(GOM.median[[measurand]], na.rm = TRUE)  # median of measurand w/o outlier removal
medianGOM.packet.after <- ddply(GOM, c("Lab", "Packet"), numcolwise(medianGOM))  # median Lab and Packets after outlier removal
## median over median of packets within lab
GOM.median.after <- ddply(medianGOM.packet.after, c("Lab"), numcolwise(medianGOM))  # median of labs after outlier removal
GOM.median.after <- merge(GOM.median.after, OKUM.methods, by = "Lab")
median.after <- median(GOM.median.after[[measurand]], na.rm = TRUE)  # median of measurand after outlier removal
## calculating method parameters
"%p%" <- function(x, y) {
    as.character(paste(x, y, sep = ""))
}
prep <- "Prep."
method <- "Method."
anal.prep <- prep %p% measurand
anal.method <- method %p% measurand
anal <- GOM.median.after[[measurand]]
anal.prep <- GOM.median.after[[anal.prep]]
anal.method <- GOM.median.after[[anal.method]]
analyte <- data.frame(GOM.median.after$Lab, GOM.median.after$names, anal, anal.prep, 
    anal.method)
analyte <- na.omit(analyte)
bymethod.n <- ddply(analyte, c("anal.method"), summarise, N = length(anal), 
    mean = round(mean(anal), 3), median = round(median(anal), 3), sd = round(sd(anal), 
        3), se = round(sd/sqrt(N), 3))
```


## Comparisons of property value calculations

The TiO2 median of the Lab+packet medians without outlier removal is 0.38 g/100 g  
The TiO2 mean of the Lab means after outlier removal of lab # 10, 26, 31 is 0.3857 g/100 g  
The TiO2 median of the Lab+packet medians after outlier removal is 0.3799 g/100 g  

## Comparisons of method parameters based on median of the Lab+packet medians after outlier removal

```r
print(bymethod.n)
```

```
##   anal.method  N  mean median    sd    se
## 1         AAS  1 0.370  0.370    NA    NA
## 2     ICP-AES  2 0.382  0.382 0.018 0.013
## 3      ICP-MS  3 0.391  0.373 0.040 0.023
## 4         XRF 17 0.386  0.380 0.020 0.005
```



## Measurement uncertainty estimations of property value


```r
medianGOM.packet$Lab <- as.factor(medianGOM.packet$Lab)  # using only the median of the 3 packages per lab
medianGOM.packet$Packet <- as.factor(medianGOM.packet$Packet)
anal <- medianGOM.packet[[measurand]]
DF.lme <- data.frame(medianGOM.packet$Lab, medianGOM.packet$Packet, medianGOM.packet[[measurand]])
DF.lme <- na.omit(DF.lme)
names(DF.lme) <- c("Lab", "Packet", "measurand")
GOM.lme <- lme(measurand ~ 1, random = ~1 | Lab/Packet, data = DF.lme)  # linear model with random effects
sL2 <- varcomp(GOM.lme, FALSE, FALSE)[[1]]  # between-laboratory variance
sbb2 <- varcomp(GOM.lme, FALSE, FALSE)[[2]]  # between bottle standard deviation
sr2 <- varcomp(GOM.lme, FALSE, FALSE)[[3]]  # repeatability standard deviation
n.p <- dim(DF.lme)[1]  # number of obervations
p <- length(unique(DF.lme$Lab))  # haven't found a better way how to extract the number of labs (number of groups)
r <- length(unique(DF.lme$Packet))
u1 <- sqrt(sL2/p + sbb2/p/r + sr2/p/r/4)  # calculating the standard uncertainty of characterization
u2 <- attr(GOM.lme$fixDF, "varFixFact")  # gives the same results as u, amazing!
plot(DF.lme)
```

![plot of chunk measurement uncertainty of data before outlier rejection](figure/measurement_uncertainty_of_data_before_outlier_rejection.png) 


```r
medianGOM.packet.after$Lab <- as.factor(medianGOM.packet$Lab)  # using only the median of the 3 packages per lab
medianGOM.packet.after$Packet <- as.factor(medianGOM.packet$Packet)
anal <- medianGOM.packet.after[[measurand]]
DF.lme <- data.frame(medianGOM.packet.after$Lab, medianGOM.packet.after$Packet, 
    medianGOM.packet.after[[measurand]])
DF.lme <- na.omit(DF.lme)
names(DF.lme) <- c("Lab", "Packet", "measurand")
GOM.lme <- lme(measurand ~ 1, random = ~1 | Lab/Packet, data = DF.lme)  # linear model with random effects
sL2.a <- varcomp(GOM.lme, FALSE, FALSE)[[1]]  # between-laboratory variance
sbb2.a <- varcomp(GOM.lme, FALSE, FALSE)[[2]]  # between bottle standard deviation
sr2.a <- varcomp(GOM.lme, FALSE, FALSE)[[3]]  # repeatability standard deviation
n.p <- dim(DF.lme)[1]  # number of obervations
p <- length(unique(DF.lme$Lab))  # haven't found a better way how to extract the number of labs (number of groups)
r <- length(unique(DF.lme$Packet))
t.value <- qt(0.975, df = p - 1)
u1.a <- sqrt(sL2.a/p + sbb2.a/p/r + sr2.a/p/r/4)  # calculating the standard uncertainty of characterization
u2.a <- attr(GOM.lme$fixDF, "varFixFact")  # gives the same results as u1, amazing!
plot(DF.lme)
```

![plot of chunk measurement uncertainty of data after outlier rejection](figure/measurement_uncertainty_of_data_after_outlier_rejection.png) 

### before outlier rejection

The between-laboratory variance for TiO2 is 4.0248 &times; 10<sup>-4</sup>   
The between-bottle variance for TiO2 is 6.3813 &times; 10<sup>-6</sup>   
The repeatability variance for TiO2 is 4.4695 &times; 10<sup>-6</sup>    
The standard uncertainty for the assigned value of TiO2 is 0.004  

### after outlier rejection
The between-laboratory variance for TiO2 is 4.4874 &times; 10<sup>-4</sup>   
The between-bottle variance for TiO2 is 7.0475 &times; 10<sup>-6</sup>   
The repeatability variance for TiO2 is 4.9493 &times; 10<sup>-6</sup>    
The standard uncertainty for the property value of TiO2 is 0.0044  
The standard uncertainty for the property value of TiO2 is 0.0044  (check with different way of calculation)
 
### tests for normal distribution

```r
qqnorm(GOM.median.after[[measurand]])
qqline(GOM.median.after[[measurand]])
```

![plot of chunk QQ plot](figure/QQ_plot.png) 


### final result based on median for property value and measurement uncertainty based on variance components
The TiO2 median of the Lab+packet medians after outlier removal is 0.3799 g/100 g   
The expanded standard uncertainty for the assigned value of TiO2 is 0.0092 
exluded labs for TiO2 is/are 10, 26, 31  
labs remaining for calculations 23  
#### comments
labs #10, #26 and #31 removed based on Youden plot
