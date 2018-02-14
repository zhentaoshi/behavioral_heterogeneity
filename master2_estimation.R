library(plyr)
library(optimx)
library(numDeriv)
library(sandwich)
library(zoo)
library(timeDate)
library(tseries)
library(urca)
library(MASS)


estimate = function(j){
  period <<- j # global variable
  
  pts0 = Sys.time()
  source("estimation")
  cat("=====", Sys.time()-pts0, "===\n\n\n\n\n\n")
  
  RES = matrix(0, 4, 6)

  
  RES[,1] = theta.gmm
  RES[,2] = theta.gmm.lower
  RES[,3] = theta.gmm.upper

  

  
  RES[,4] = theta.EL
  RES[,5] = theta.EL.lower
  RES[,6] = theta.EL.upper
  
  
  
  return(RES)
}


Estimat = ldply(.data = 1:3, .fun = estimate)
print(Estimat)

Estimat3 = format(Estimat, nsmall = 3, digits = 3, justify = "left")
write.csv(Estimat3, file = "Estimate.csv")
