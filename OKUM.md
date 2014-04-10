Analyzing submitted data for OKUM
========================================================


```r
measurand <- "CaO"
```


## measurand selected CaO


#### Importing the data and assigning factors


```r
library(ggplot2)
library(devtools)  #needed for 'order_by'
```

```
## WARNING: Rtools is required to build R packages, but is not currently installed.
## 
## Please download and install Rtools 3.1 from http://cran.r-project.org/bin/windows/Rtools/ and then run find_rtools().
```

```r
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



#### defining the function for plotting methods vs. measurand mass fraction. Sample preparation methods are also marked in the plot.

```r
##### function plot_method defined here. Automated plot of methods boxplots
##### sorted by increasing method mean. Enter with plot_method=('xx') #######
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
    p <- ggplot(order_by(anal.method, ~anal, analyte, mean), aes(anal.method, 
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
    p <- ggplot(order_by(GOM.Lab, ~anal, analyte, mean), aes(GOM.Lab, anal)) + 
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
    
    RM1 <- GOM[[RM1]]
    RM1 <- (RM1 - median(RM1, na.rm = TRUE))/sd(RM1, na.rm = TRUE)  #calculating z-scores
    RM2 <- GOM[[RM2]]
    RM2 <- (RM2 - median(RM2, na.rm = TRUE))/sd(RM2, na.rm = TRUE)  #calculating z-scores
    RM <- data.frame(GOM$Lab, GOM$names, RM1, RM2)
    RM <- na.omit(RM)
    p <- ggplot(RM, aes(RM1, RM2, label = GOM.Lab)) + xlim(-5, 5) + ylim(-5, 
        5) + geom_point(aes(colour = factor(GOM.Lab)), size = 4)
    p <- p + xlab(y) + ylab(z) + labs(title = a) + labs(colour = "Lab") + mytheme + 
        geom_abline(intercept = 0, slope = 1) + geom_abline(intercept = 2.8284, 
        slope = 1) + geom_abline(intercept = -2.8284, slope = 1) + geom_abline(intercept = 2.8284, 
        slope = -1) + geom_abline(intercept = -2.8284, slope = -1) + geom_text(data = NULL, 
        x = -0.5, y = 2.85, label = "z = 2") + geom_text(aes(colour = factor(GOM.Lab)), 
        hjust = 1, vjust = 0)
    p + theme(legend.position = "none")  #+ stat_density2d(aes(fill = ..level..), geom='polygon') 
}
```


```r
avgbymethod <- function(x) {
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
    bymethod <- tapply(analyte$anal, analyte$anal.method, mean, na.rm = TRUE)
    cat("method means", "\n")
    print(bymethod)
    avg <- mean(tapply(analyte$anal, analyte$GOM.Lab, mean, na.rm = TRUE), na.rm = TRUE)
    median <- median(tapply(analyte$anal, analyte$GOM.Lab, median, na.rm = TRUE), 
        na.rm = TRUE)
    # cat('mean of lab means: ','\n') cat(round(avg,3),'g/100g', '\n')
    # cat('median of lab medians: ','\n') cat(round(median,3),'g/100g', '\n')
}
```


```r
summary(GOM[[measurand]], na.rm = TRUE, digits = 4)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
##    7.21    7.78    7.86    7.88    7.95    8.83     162
```

```r
mean <- mean(tapply(GOM[[measurand]], GOM$Lab, mean, na.rm = TRUE), na.rm = TRUE)
median <- median(tapply(GOM[[measurand]], GOM$Lab, median, na.rm = TRUE), na.rm = TRUE)
```

The CaO mean of the Lab means is 7.8834 g/100 g  
The CaO median of the Lab medians is 7.855 g/100 g  


```r
avgbymethod(measurand)
```

```
## method means 
##               ICP-AES titrimetry        XRF 
##         NA      7.706      8.164      7.888
```




```r
plot_method(measurand)
```

```
## Warning: Removed 1 rows containing missing values (geom_point).
## Warning: Removed 12 rows containing missing values (geom_point).
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-101.png) 

```r
plot_lab(measurand, "M")
```

```
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-102.png) 


```r
plot_youd(measurand, "GAS", "OKUM")
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-111.png) 

```r
plot_youd(measurand, "MUH", "OKUM")
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-112.png) 

#### Mandel k barplot displays the within lab performance relative to all participating labs

```r
k <- with(GOM, mandel.kh(GOM[[measurand]], g = Lab, type = "k"))
barplot(k, las = 2, col = 1:4)
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12.png) 



#### Removing ouliers lab based on Youden plot


```r
outlier <- c(30, 12, 7, 16)
leng <- length(outlier)
for (i in 1:leng) {
    GOM[[measurand]] <- ifelse(GOM$Lab == outlier[i], NA, GOM[[measurand]])
    message("Lab ", outlier[i], " was removed")
    print(summary(GOM[[measurand]], na.rm = TRUE, digits = 4))
}
```

```
## Lab 30 was removed
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
##    7.21    7.78    7.86    7.88    7.93    8.83     174
```

```
## Lab 12 was removed
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
##    7.21    7.76    7.86    7.84    7.92    8.39     186
```

```
## Lab 7 was removed
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
##    7.21    7.74    7.85    7.83    7.91    8.39     198
```

```
## Lab 16 was removed
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
##    7.38    7.78    7.86    7.83    7.91    8.02     210
```

```r
plot_lab(measurand, "M")
```

```
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
```

![plot of chunk outlier removal](figure/outlier_removal.png) 

```r
# GOM[[measurand]] <- ifelse(GOM$Lab==remove[2], NA, GOM[[measurand]])
# summary(GOM[[measurand]], na.rm=TRUE, digits=4) plot_lab(measurand, 'M')
```


```r
avgbymethod(measurand)
```

```
## method means 
##               ICP-AES titrimetry        XRF 
##         NA      7.701         NA      7.840
```


```r
summary(GOM[[measurand]], na.rm = TRUE, digits = 4)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
##    7.38    7.78    7.86    7.83    7.91    8.02     210
```

```r
mean <- mean(tapply(GOM[[measurand]], GOM$Lab, mean, na.rm = TRUE), na.rm = TRUE)
median <- median(tapply(GOM[[measurand]], GOM$Lab, median, na.rm = TRUE), na.rm = TRUE)
```

The CaO mean of the Lab means after outlier removal is 7.8325 g/100 g  
The CaO median of the Lab medians after outlier removal is 7.846 g/100 g 
