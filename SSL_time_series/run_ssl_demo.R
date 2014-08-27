###
# Steller sea lion forecast extrapolation
###

library(agTrend)
library(ggplot2)
library(coda)
library(Rcpp)
Rcpp::sourceCpp('PCtimeMCMC.cpp')

data(wdpsNonpups)

mar_data = wdpsNonpups[wdpsNonpups$site=="MARMOT" & wdpsNonpups$year>1989,]

###
pred_window=20
lf_freq=5
hf_freq=5
###

knots_lf = seq(1990-lf_freq,2011+pred_window+lf_freq, by=lf_freq)
knots_hf = seq(1990-hf_freq,2011+pred_window+hf_freq, by=hf_freq)

K_lf = exp(-0.5*outer(1990:(2011+pred_window),knots_lf,"-")^2/lf_freq^2)
K_lf=sweep(K_lf, 2, colSums(K_lf), "/")

K_hf = exp(-0.5*outer(1990:(2011+pred_window),knots_hf,"-")^2/hf_freq^2)
K_hf=sweep(K_hf, 2, colSums(K_hf), "/")


post_smp=PCtimeMCMC(
  y=mar_data$count, 
  K_lf=K_lf, 
  K_hf=0.0*K_hf, 
  sampleIndex=as.numeric(1990:(2011+pred_window)%in%mar_data$year), 
  beta_mean=0, 
  beta_prec=0, 
  phi_lf_scale=100, 
  phi_hf_scale=100, 
  a_sigma = 0,
  b_sigma=0,
  block=1000,
  burn=10000, 
  iter=1000000,
  sample_sigma=TRUE,
  sample_hf=FALSE
  )

# xxx = apply(post_smp$jump_idx, 2, function(x){cumsum(x)/c(1:length(x))})
# apply(tail(post_smp$jump_idx, 10000), 2, mean)
# plot(xxx[,2], type='l')
# abline(h=0.234)
# plot(post_smp$tune[,1], type='l')
# 
# plot(mcmc(post_smp$beta))
# plot(mcmc(post_smp$alpha_f))
# 
# plot(mcmc(post_smp$phi_hf))


pred_data = data.frame(year=1990:(2011+pred_window), pred=apply(post_smp$pred, 2, median), HPDinterval(mcmc(post_smp$pred), prob=0.9))
pred_plot= ggplot() +
  geom_point(aes(x=year, y=count), data=mar_data) +
  geom_path(aes(x=year, y=pred), data=pred_data[pred_data$year<2012,]) + 
#   geom_ribbon(aes(ymin=lower, ymax=upper, x=year), data=pred_data, alpha=0.3)
  geom_ribbon(aes(ymin=lower, ymax=upper, x=year), data=pred_data[pred_data$year<2012,], alpha=0.3)
ggsave("ssl_pred_plot.pdf", pred_plot, width=6.5, height=4.5)

# ggplot() + geom_histogram(aes(x=post_smp$phi_lf))
# ggplot() + geom_histogram(aes(x=post_smp$phi_hf))
# ggplot() + geom_histogram(aes(x=post_smp$sigma))

Xaug = cbind(1, K_lf, 0.0*K_hf)
V = var(cbind(post_smp$beta, post_smp$alpha_lf, post_smp$alpha_hf))
givh = diag(Xaug%*%V%*%t(Xaug)) 
givh_data = data.frame(year=pred_data$year, sample=as.numeric(1990:(2011+pred_window)%in%mar_data$year), givh=givh)
givh_data$givh = givh_data$givh/max(givh_data$givh[givh_data$sample==1])

p = ggplot() + geom_path(aes(x=year, y=givh), data=givh_data, lwd=2) + 
  geom_point(aes(x=year, y=-5), data=givh_data[givh_data$sample==1,], cex=2)
ggsave("ssl_time_series.pdf", p, width=6.5, height=6.5)

