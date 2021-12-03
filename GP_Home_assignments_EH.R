##########################################################
##### ASSIGMENT PART 1: PRINCIPLE COMPONENT ANALYSIS #####
##########################################################

# Package loading 
library(tidyverse)
library(gsheet)
library(psych)	
library(gridExtra)
library(dplyr)
library(ggplot2)
library(psych)

# Data collection

data_PAQ <- PAQ_Emma
view(data_PAQ)

# Reshaping of data from long to wide 

data_PAQ_long <- data_PAQ

library(reshape2) 
library(plyr)

duplicated(data_PAQ_long[, c("id", "var", "value")])
sum(duplicated(data_PAQ_long[, c("id", "var", "value")]))  

data_PAQ_wide <- reshape(data_PAQ_long, v.names = c("value"), timevar = "var", direction = "wide", sep = "")
view(data_PAQ_wide)

# Renaming of variables 
data_PAQ_wide <- rename(data_PAQ_wide, c("valuesex" = "sex"))
data_PAQ_wide <- rename(data_PAQ_wide, c("valueage" = "age"))
data_PAQ_wide <- rename(data_PAQ_wide, c("valueQ1_cry" = "Q1_cry"))
data_PAQ_wide <- rename(data_PAQ_wide, c("valueQ2_help" = "Q2_help"))
data_PAQ_wide <- rename(data_PAQ_wide, c("valueQ3_breathe" = "Q3_breathe"))
data_PAQ_wide <- rename(data_PAQ_wide, c("valueQ4_freeze" = "Q4_freeze"))
data_PAQ_wide <- rename(data_PAQ_wide, c("valueQ5_alien" = "Q5_alien"))
data_PAQ_wide <- rename(data_PAQ_wide, c("valueQ6_inferior" = "Q6_inferior"))
data_PAQ_wide <- rename(data_PAQ_wide, c("valueQ7_weep" = "Q7_weep"))
data_PAQ_wide <- rename(data_PAQ_wide, c("valueQ8_Support" = "Q8_support"))
data_PAQ_wide <- rename(data_PAQ_wide, c("valueQ9_Nerd" = "Q9_nerd"))

view(data_PAQ_wide)

data_anx <- data_PAQ_wide
view(data_anx)

# Removing rows with NA

data_anx <- na.omit(data_anx)
view(data_anx)

# Removing ID-column 
data_anx$id <- NULL
view(data_anx)

# Descriptive information about the dataset
describe(data_anx)
dim(data_anx)
str(data_anx)
names(data_anx)

row.names(data_anx)
nrow(data_anx)
ncol(data_anx)
head(data_anx)
tail(data_anx)

# Checking data for coding errors and outliers 
ggplot(data_anx, aes(y=age))+stat_boxplot(geom="errorbar")+geom_boxplot()
ggplot(data_anx, aes(y=sex))+stat_boxplot(geom="errorbar")+geom_boxplot()
ggplot(data_anx, aes(y=Q1_cry))+stat_boxplot(geom="errorbar")+geom_boxplot()
ggplot(data_anx, aes(y=Q2_help))+stat_boxplot(geom="errorbar")+geom_boxplot()
ggplot(data_anx, aes(y=Q3_breathe))+stat_boxplot(geom="errorbar")+geom_boxplot()
ggplot(data_anx, aes(y=Q4_freeze))+stat_boxplot(geom="errorbar")+geom_boxplot()
ggplot(data_anx, aes(y=Q5_alien))+stat_boxplot(geom="errorbar")+geom_boxplot()
ggplot(data_anx, aes(y=Q6_inferior))+stat_boxplot(geom="errorbar")+geom_boxplot()
ggplot(data_anx, aes(y=Q7_weep))+stat_boxplot(geom="errorbar")+geom_boxplot()
ggplot(data_anx, aes(y=Q8_support))+stat_boxplot(geom="errorbar")+geom_boxplot()
ggplot(data_anx, aes(y=Q9_nerd))+stat_boxplot(geom="errorbar")+geom_boxplot()

# Creating covariance and correlation matrix 

cov_matrix <- cov(x=data_anx)
cov_matrix

cor_matrix <- cor(x=data_anx)
cor_matrix

# Selecting questionnarie values for the PCA (excluding sex and age)
data_anx$sex <- NULL
view(data_anx)

data_anx$age <- NULL
view(data_anx)

#Principle component analysis

data_pcacov=princomp(cov_matrix,cor=FALSE)
summary(data_pcacov, loadings=TRUE)

data_pcacor=princomp(cor_matrix,cor=TRUE)
summary(data_pcacor, loadings=TRUE)

# Plotting covariance and correlation for variance and eigenvalue

plot(data_pcacov$sdev^2, xlab = "Component number", ylab = "Component variance", type = "l", main = "Scree diagram")
plot(log(data_pcacov$sdev^2), xlab = "Component number", ylab = "log(Component variance)", type="l", main = "Log(eigenvalue) diagram")

plot(data_pcacor$sdev^2, xlab = "Component number", ylab = "Component variance", type = "l", main = "Scree diagram")
plot(log(data_pcacor$sdev^2), xlab = "Component number", ylab = "log(Component variance)", type="l", main = "Log(eigenvalue) diagram")

# Based on cumulative proportion of variance, and the screeplot of eigenvalue, 
# component 1-3 will be used 

data_anx_Q3 <- data_anx[, c("Q1_cry", "Q2_help", "Q3_breathe")]
colMeans(data_anx)

cov(data_anx_Q3)

anx_pca <- princomp(x = data_anx_Q3)
anx_pca

print(summary(anx_pca), loadings = TRUE)

xlim <- range(anx_pca$scores[,1])
plot(anx_pca$scores, xlim = xlim, ylim = xlim)

# Biplot
biplot(anx_pca, col = c("grey", "black"))


######################################################
##### ASSIGMENT PART 2: MULTIDIMENSIONAL SCALING #####
######################################################

?cmdscale

library(smacof)
library(MASS)

view(Nations)

# Checking data for coding errors (numbers outside the scale) 
Nations[Nations>9]
Nations[Nations<1]

# Nations distance and distance matrix 
require(smacof)
Nations.distance = sim2diss(Nations, method = 1)
Nations.distance

Nations.distance <- dist(Nations.distance)
Nations.distance

# CMD scaling 
cmdscale(Nations.distance, k = 8, eig = TRUE)

max(abs(dist(Nations) - dist(cmdscale(Nations.distance, k = 8))))

# Multidimensional scaling of Nations
require(MASS)
Nations_mds = isoMDS(Nations.distance)
Nations_mds

# Non-metrics MDS and plotting 
x <- Nations_mds$points[,1]
y <- Nations_mds$points[,2]
plot(x, y, xlab = "Coordinate 1", ylab = "Coordinate 2", xlim = range(Nations_mds$points[,1])*1.2, type = "n")
text(x, y, labels = colnames(Nations), cex = 0.6)
Nations_sh <- Shepard(Nations[lower.tri(Nations)], Nations_mds$points)
Nations_sh

# Two dimensional solution 
(Nations_mds <- isoMDS(Nations.distance))

# Dissimilarities plot 
plot(Nations_sh, pch = ".", xlab = "Dissimilarity", ylab = "Distance", xlim = range(Nations_sh$x), ylim = range(Nations_sh$x)) +
lines(Nations_sh$x, Nations_sh$yf, type = "S")

