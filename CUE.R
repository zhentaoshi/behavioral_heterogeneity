cat("################## GMM ################\n")

mm = 8 # total number of moments
multi_start = 100 
VAL = rep(0, multi_start)

set.seed(period+1) # random seed

lower = rep(0.001, 4)
upper = c(3,3,3,6)

init = matrix( runif(4*multi_start, min = 0.001, 
                     max = rep(upper, each = multi_start)), ncol = 4)


init[, 1] = sd(diff.mu)



for (i in 1:multi_start){
  print(i)
  val = optimx( init[i,], gmm.model1, method = "L-BFGS-B",  
                lower = lower, upper = upper )$value
  VAL[i] = val
}
print(VAL)

chosen.init = init[which.min(VAL), ]
opt <- optimx( chosen.init, gmm.model1, method = "L-BFGS-B",  
               lower = lower, upper = upper  ) 


theta.hat <- as.vector( as.numeric(opt[1:length(chosen.init)]) )



####################################################
# calculate the covariance matrix
gg = function(theta) colMeans( moments.components(theta)$moments )

gmm.jacobian <- jacobian(gg, theta.hat)
gmm.cov <- moments.components(theta.hat)$W

####################################################
theta.gmm <- theta.hat

variance =  ginv( t( gmm.jacobian ) %*% gmm.cov %*%  gmm.jacobian  ) 
Asd.theta.gmm = sqrt( diag(variance)/T )
theta.gmm.lower = theta.gmm - 1.96*Asd.theta.gmm
theta.gmm.upper = theta.gmm + 1.96*Asd.theta.gmm

print(rbind(theta.gmm, theta.gmm.lower, theta.gmm.upper))

# 
J = opt$value
cat( "\n J = ", J, "p-value = ", 1-pchisq(J, mm - 4), "\n" )
