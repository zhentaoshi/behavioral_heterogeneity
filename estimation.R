library(optimx)
library(numDeriv)
library(sandwich)
library(zoo)
library(timeDate)
library(tseries)
library(urca)
library(MASS)

cat("###################################\n")
cat("\n time period = ", period, "\n")


source("func5.R")
load(paste("period", period,".Rdata", sep=""))

# tuning parameters
blocksize = floor(1.14*T^(1/3))
Del = p.lag1 - p.lag2  
Del = Del [abs(Del) < 0.40] # remove outlines
bw = 1.06*sd(Del)*(len^(-1/5));


source("CUE.R")

source("EL.R")




