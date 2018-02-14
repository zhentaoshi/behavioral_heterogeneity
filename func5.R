
demean = function(x) { x - mean(x)  } # remove the sample mean
kern = function(x, w = bw)  exp( - 0.5 * (x/w)^2 ) # w is the bandwith


post.model1 <- function(theta, fraction = FALSE){
  # the function that inputs "theta" and returns a vector "D(t)"
  
  sigma.mu = theta[1]  
  eta = theta[2] 
  tau = theta[3]
  alpha = theta[4]
  
  
  A1 = sqrt(tau/eta) * abs(p.lag1 - p.lag2) / (alpha * sigma.mu)
  delta_t = ( (1+alpha) * p.lag1 - (mu_vec_lag + alpha * mu_vec) )/ (alpha * sigma.mu)
  
  yM = delta_t  + A1;
  ym = delta_t  - A1;
  
  m_vec = pnorm( yM ) - pnorm( ym );
  
  B1 = eta/(1+alpha) * sqrt(alpha) * sigma.mu * ( dnorm(yM) - dnorm(ym) - (1-m_vec) * delta_t  )
  B3 = tau * m_vec * (p.lag1 - p.lag2)
  
  lambda.D =  B1 + B3
  
  if (fraction == TRUE){
    return(m_vec)
  } else {
    return(lambda.D)
  }
}



moments.components <- function(theta, typ = "CUE"){
  
  sigma.mu = theta[1]  
  eta = theta[2] 
  tau = theta[3]
  alpha = theta[4]
  
  demean = function(x) { x - mean(x)  }
  
  
  if (typ == "CUE") {   
    lD = post.model1(theta)    
    demean.lD = demean(lD)
    
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
    
    z3 = tau * ( 2 * pnorm( sqrt(tau/eta) * abs( dist1 )/(alpha * sigma.mu)  ) - 1 ) * dist1
    thin.moment3 = abs(R) - abs(z3)
    
    #######
    
    moments = cbind(
      (diff.mu^2 - theta[1]^2)*100 ,
      weight1 * thin.moment1*10,
      weight1 * thin.moment2*10,
      weight2 * thin.moment3*10,
      R - lD, # do not demean lD. only demean those of higher moments
      ( R^2 - demean.lD^2),
      ( R^3 - demean.lD^3),
      ( R^4 - demean.lD^4)  
    )  
    
    # HAC
    # "lag_order" and "bartlett" are only useful in CUE
    lag_order = ceiling(T^(1/3)) 
    bartlett = sort((1:lag_order)/lag_order, decreasing = TRUE)
    ACOV = acf(moments, lag =lag_order - 1, type = "covariance", plot = F, demean = T)$acf
    LRvar = ACOV * bartlett
    LRvar = 2 * apply(LRvar, c(2,3), FUN = sum) - ACOV[1, ,] # the first ACOV is the variance

    W = solve( LRvar )
    
  }   
  
  return(list(moments = moments, W = W))
}


###########################################

gmm.model1 <- function(theta, typ = "CUE" ){
  
  if  ( any(theta < 0.001) ) { J = 1000000 }
  else {
    mW <- moments.components(theta, typ )
    
    mm <- colMeans(mW$moments)
    J = T * as.vector( t(mm) %*% mW$W %*% mm)
  }
  
  return( J )
}

### used in EL

moment.p = function(theta, COV = FALSE){
  gg.block = ggblock(theta)
  num.moments = dim(gg.block)[2]
  
  Prob.matrix = matrix(  probs, nrow = n, ncol = num.moments ) 
  Prob.array  = array( probs, dim = c( n, num.moments, num.moments) )
  
  g_each = gg.block * Prob.matrix 
  g_sum  = colSums(g_each)
  
  g_cov1 = array( t(gg.block), dim = c( n, num.moments, num.moments  ) )
  g_cov2 = aperm( g_cov1, perm = c(1, 3, 2) )
  g_COV  = g_cov1 * g_cov2 * Prob.array
  g_COV  = apply( g_COV, c(2,3), sum )
  
  if ( COV == TRUE){
    return( g_COV)
  } else if (COV == FALSE){
    return(  g_sum )
  }
}

