z = rnorm(100, 0, 0.3)
x = 1
x_store=rep(NA, 10000)
for(i in 1:length(x_store)){
  x_prop = x + rnorm(1,0,1)
  mh = exp(dcauchy(x_prop,0,1,log=TRUE) + sum(dnorm(z,0,abs(x_prop),log=TRUE)) - dcauchy(x,0,1,log=TRUE) - sum(dnorm(z,0,abs(x),log=TRUE)))
  if(runif(1,0,1)<=mh){x = x_prop}
  x_store[i] = abs(x)
}

sigma = cauchy_test(z)
