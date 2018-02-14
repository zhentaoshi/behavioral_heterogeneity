library(Rmosek)
library(nloptr)
library(MASS)
library(numDeriv)

cat("################## EL ################\n")

eff.T = T
theta0 = theta.gmm


# this script apply the idea of empirical likelihood 
# through the R-MOSEk in the inner optimization

n = floor( T/blocksize )

eff.T = n * blocksize
p.lag1 = p.lag1[1:eff.T]
p.lag2 = p.lag2[1:eff.T]
mu_vec = mu_vec[1:eff.T]
mu_vec_lag = mu_vec_lag[1:eff.T]
diff.mu = diff.mu[1:eff.T]
R = R[1:eff.T]
p = p[1:eff.T]


################## functions 

ggblock = function(theta){
  # the parameters
  sigma.mu = theta[1]
  eta = theta[2]
  tau = theta[3]
  alpha = theta[4]
  
  
  
  # generate the uniform moments
  lD = post.model1(theta)  # the demand from the formula
  demean.lD = demean(lD)
  
  
  gg = rbind(  ( diff.mu^2 - theta[1]^2)*100,
               (R - lD)*0.5, # do not demean lD. only demean those of higher moments
               ( R^2 - demean.lD^2)*10,
               ( R^3 - demean.lD^3)*50,
               ( R^4 - demean.lD^4)*100  )
  
  
  # moments from thin-set
  dist1 = p.lag1 - p.lag2
  weight1 = kern(dist1)
  
  z1 = mu_vec_lag - p.lag1
  z2 = mu_vec - p.lag1
  ld.f = eta/(1+alpha) * z1 + eta * alpha/(1+alpha) * z2
  
  thin.moment1 = R - ld.f
  thin.moment2 = R^2 - demean(ld.f)^2

  dist2 = (1+alpha) * p.lag1 - mu_vec_lag - alpha * mu_vec
  weight2 = kern(dist2)
  
  z3 = tau * ( 2 * pnorm( sqrt(tau/eta) * abs( dist1 )/( alpha * sigma.mu)  ) - 1 ) * dist1
  thin.moment3 = abs(R) - abs(z3)
  
  gg = rbind( gg,
              weight1 * thin.moment1*10,
              weight1 * thin.moment2*10,
              weight2 * thin.moment3*10
  )
  
  
  num.moments = dim(gg)[1]
  gg.block = array( t(gg), dim = c(blocksize, n, num.moments ) )
  gg.block = apply(gg.block, c(2,3), mean)
  return(gg.block)
}




inner = function(theta, probb = FALSE){
  # the parameters
  sigma.mu = theta[1]
  eta = theta[2]
  tau = theta[3]
  alpha = theta[4]
  
  # print(theta)
  
  if (any(theta< 0.0001)) { J = 1000000} # restrict the parameters
  else  {    
    gg.block = ggblock(theta)
    # all moments
    
    prob$A = Matrix( rbind( 1,  t( gg.block ) ) ) # combine the sum of prob into the linear constraint
    
    prob$scopt <- list ( opro = opro  )
    
    
    r = mosek ( prob, opts = list(soldetail = 2, verbose = 0) )
    if (r$response$msg == "MSK_RES_ERR_USER_NLO_FUNC: The user-defined nonlinear function reported an error."){
      J = 1000}
    else{
      J = -r$sol$itr$pobjval # primal objective value
}
  }

  if (probb == FALSE){
    return( J ) 
  } else if (probb == TRUE){
    return( r$sol$itr$xx )
  }
}




############### execution 

mm = 8 # number of moments


prob = list(sense = "max")

prob$c = rep(0,n) 


prob$bc = rbind( blc = c( 1, rep(0,mm) ), buc = c( 1,  rep(0,mm) ) ) 
prob$bx = rbind( blx = rep(0,n), bux = rep(1,n) )

# the decision variables are the probablity 'p'
opro <- matrix ( list (), nrow =5, ncol = n )
rownames ( opro ) <- c(" type ","j","f","g","h")
for (i in 1:n){ opro[,i] = list("LOG", i, 1, 1, 0 ) }




opts = list("algorithm"="NLOPT_LN_NELDERMEAD",
            "xtol_rel"=1.0e-9, 
            "maxeval" = 1000, 
            "xtol_abs" = 1e-9)




res = nloptr( x0= theta0, 
                   eval_f=inner, probb = FALSE,
                   opts=opts)


theta.EL = res$solution

probs = inner(theta.EL, probb = TRUE)

print( res )


probs = probs/sum(probs)

# LR stat
Tb =    (T/blocksize  )
LR = 2* (-sum(  log(probs) ) - ( Tb ) *log( Tb ) )


cat( "LR = ", LR, "p-value = ", 1-pchisq(LR, mm - 4), "\n" )
cat("\n sum probs", sum(probs), "\n")

### inference #### 
g_cov = moment.p( theta.EL, TRUE)


el.jacobian <- jacobian(moment.p, method = "simple",  theta.EL)
Omega <-  t(el.jacobian) %*% ginv(g_cov) %*% el.jacobian

variance = ginv( Omega  ) 
Asd.theta.EL = sqrt( diag(variance/T) )

theta.EL.lower = theta.EL - 1.96*Asd.theta.EL
theta.EL.upper = theta.EL + 1.96*Asd.theta.EL


print(rbind(theta.EL, theta.EL.lower, theta.EL.upper))

