# this script generates the data of the three periods
# and save them to data files for analysis


## prepare

### the detrend function

detrend <- function(y){
  # remove the linear time trend 
  
  time_trend = 1:length(y)
  y_new = residuals( lm( y ~ time_trend ))
  return(y_new)
}

#############
# start the data generating


###### import the data
d0 <- read.csv("raw_data.csv")
names(d0)[1] <- "date"



# star to handle the data

for (period in 1:3) {
  
  # because of the lagged terms, the starting point is about three years as 
  # specified in the paper, so that the effective data is exactly the same 
  # as in the paper.
  
  if (period == 1){
    end.date <- 2013.99
    start.date <- 1988.09; 
  } else if (period == 2){
    end.date <- 1990.99
    start.date <- 1958.09;  
  } else if (period == 3){
    end.date <- 1960.99
    start.date <- 1908.09; 
  } 
  
  
  cat("\n start =", floor(start.date)+3, "Jan; end = ", floor(end.date) , "Dec \n")
  
  len = floor( end.date - start.date - 3) * 12
  
  number_of_years = end.date - start.date;
  
  
  d1 <- d0[ d0$date >= start.date & d0$date <= end.date, ]
  
  div <- d1$D
  p   <- d1$p
  mu_vec = log(div)
  p <- log(p)
  
  
  k = 1
  # generate the lag terms
  p.lag1 <- c(rep(0,k),p)
  p.lag1 <- p.lag1[ - ( (length(p.lag1)-k):length(p.lag1) ) ] # 1-month lag
  p.lag2 <- c(rep(0,k),p.lag1)
  p.lag2 <- p.lag2[ - ( (length(p.lag2)-k):length(p.lag2) ) ] # 2-month lag
  
  mu_vec_lag <- c(rep(0,k),mu_vec)
  mu_vec_lag <- mu_vec_lag[ - ( (length(mu_vec_lag)-k):length(mu_vec_lag) )]
  
  # align the time series by removing the first k observations
  
  p<- p[-(1:k)]
  p.lag1<- p.lag1[-(1:k)]
  p.lag2<- p.lag2[-(1:k)]
  mu_vec  <- mu_vec[-(1:k)]
  mu_vec_lag <- mu_vec_lag[-(1:k)]
  
  
  p.new <- numeric()
  p.new1 <- numeric()
  p.new2 <- numeric()
  m.new <- numeric()
  mu_lag_new <- numeric()
  
  
  # generate the moving average time series
  k = 12
  max.lag = 36
  
  for (t in 12:length(p.lag2)){
    
    p.new[t]  <-   p[ t  ]
    p.new1[t] <-   p.lag1[t] 
    p.new2[t] <-   mean( p.lag2[( t-12+1): t] )
  }
  
  
  
  # remove those first 60 observations
  p = p.new[-(1:(max.lag-1))]
  p.lag1  <- p.new1[-(1:(max.lag-1))]
  p.lag2  <- p.new2[-(1:(max.lag-1))]
  
  mu_vec  <- mu_vec[-(1:(max.lag-1))]
  mu_vec_lag <- mu_vec_lag[-(1:(max.lag-1))]
  mu_vec <- mu_vec[-(1:2)]
  mu_vec_lag <- mu_vec_lag[-1]
  
  
  
  ### detrend 
  p = detrend(p)
  p.lag1 = detrend(p.lag1)
  p.lag2 = detrend(p.lag2)
  
  
  mu_vec = detrend(mu_vec)
  mu_vec_lag = detrend(mu_vec_lag)

  ##### 
  diff.mu = mu_vec - mu_vec_lag;
  
  R = p - p.lag1
  
  T = length(R)
  save.image(file = paste( "period", period, ".Rdata", sep = ""))
}