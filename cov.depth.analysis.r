library(countreg)
library(flexmix)
library(robustbase)
source('glmrob.nb.r') # main function, estimates coefficients and standard errors

cov.depth.analysis <- function(Y, method = "zinb", q.outliers = 0.999, plot.trunc = TRUE, plot = TRUE) {
  
  if (method == "pois") {
    # Outlier detection
    outlier.cutoff <- quantile(Y, probs = q.outliers)
    Y.new <- Y[Y<=outlier.cutoff]
    
    # Poisson regression
    mle.est.poi <- glm(Y.new~1, family = "poisson")
    print(mle.est.poi)
    
    poi.summary <- data.frame(outlier.cutoff,
                              pi = NA,
                              mu = exp(coef(mle.est.poi)),
                              sigma = exp(coef(mle.est.poi))
    )
    row.names(poi.summary) <- NULL
    
    
    if (plot.trunc) {
      hist <- hist(Y.new, freq = F, breaks = "FD", 
                   xlab = "Depth of coverage", ylab = "Density")
      lines(0:max(Y.new),
            dpois(0:max(Y.new), lambda = poi.summary$mu),
            col = "red", lwd = 2)
    } else {
      hist <- hist(Y, freq = F, breaks = "FD")
      lines(0:max(Y.new),
            dpois(0:max(Y.new), lambda = poi.summary$mu),
            col = "red", lwd = 2)
    }
    
    return(poi.summary)
  }
  
  if (method == "nb") {
    # Outlier detection
    outlier.cutoff <- quantile(Y, probs = q.outliers)
    Y.new <- Y[Y<=outlier.cutoff]
    
    # NB regression
    mle.est.nb <- glm.nb(Y.new~1)
    # print(mle.est.nb)
    
    nb.summary <- data.frame(outlier.cutoff,
                             pi = NA,
                             mu = exp(coef(mle.est.nb)),
                             sigma = mle.est.nb$theta
    )
    row.names(nb.summary) <- NULL
    
    if (plot.trunc & plot) {
      hist <- hist(Y.new, freq = F, breaks = "FD", 
                   xlab = "Depth of coverage", ylab = "Density")
      lines(0:max(Y.new),
            dnbinom(0:max(Y.new), mu = nb.summary$mu, size = nb.summary$sigma),
            col = "red", lwd = 2)
    } else if (plot) {
      hist <- hist(Y, freq = F, breaks = "FD")
      lines(0:max(Y.new),
            dnbinom(0:max(Y.new), mu = nb.summary$mu, size = nb.summary$sigma),
            col = "red", lwd = 2)
    }
    
    return(nb.summary )
  }
  
  if (method == "zinb") {
    # Y.gt0 <- Y[Y>0]
    # Initialized: Robust M-estimators
    # X=as.matrix(rep(1,1*length(Y.gt0)), ncol = 1)
    # rob.est.nb <- glmrob.nb(Y.gt0, X, 
    #                         bounding.func='T/T',
    #                         c.tukey.beta=6,c.tukey.sig=4,
    #                         maxit.sig=50,
    #                         weights.on.x='none')
    
    # Outlier detection
    # rob.est.nb.mu <- exp(rob.est.nb$coef[-1])
    # rob.est.nb.sigma <- rob.est.nb$coef[1]
    # outlier.cutoff <- qnbinom(q.outliers, mu = rob.est.nb.mu, size = 1 / rob.est.nb.sigma)
    outlier.cutoff <- quantile(Y, probs = q.outliers)
    Y.new <- Y[Y<=outlier.cutoff]

    # zero inflacted model
    # Y.new.df <- as.data.frame(Y.new)
    # Y.new.df$X <- 1
    rob.est.zero <- zeroinfl(Y.new ~ 1 | 1, dist = "negbin")
    # print(summary(rob.est.zero))
    # print(rob.est.zero)
    zinb.summary <- data.frame(
                               # mu.ini = rob.est.nb.mu,
                               # sigma.ini = rob.est.nb.sigma,
                               outlier.cutoff = outlier.cutoff,
                               pi = exp(coef(rob.est.zero)[2])/(1+exp(coef(rob.est.zero)[2])),
                               mu = exp(rob.est.zero$coefficients$count),
                               sigma = rob.est.zero$theta)
    row.names(zinb.summary) <- NULL
    
    if (plot.trunc & plot) {
      hist <- hist(Y.new, freq = F)
      lines(0:max(Y.new), 
            zinb.summary$pi * c(1, rep(0, max(Y.new))) + (1 - zinb.summary$pi) * dnbinom(0:max(Y.new), mu = zinb.summary$mu, size = zinb.summary$sigma), 
            col = "blue", lwd = 2)
    } else if (plot) {
      hist <- hist(Y, freq = F, breaks = "FD")
      lines(0:max(Y.new), 
            zinb.summary$pi * c(1, rep(0, max(Y.new))) + (1 - zinb.summary$pi) * dnbinom(0:max(Y.new), mu = zinb.summary$mu, size = zinb.summary$sigma), 
            col = "blue", lwd = 2)
    }
    
    return(zinb.summary)
  }
  
  if (method == "fmnb") {
    # Outlier detection
    outlier.cutoff <- quantile(Y, probs = q.outliers)
    Y.new <- Y[Y<=outlier.cutoff]
    
    # Flexible mixture model of nb
    fm.nb <- flexmix(Y.new ~ 1 | 1, k = 2, model = FLXMRnegbin()) 
    print(summary(fm.nb))
    print(parameters(fm.nb))
    fm.summary <- data.frame(outlier.cutoff = outlier.cutoff,
                             pi = tabulate(clusters(fm.nb))[1] / sum(tabulate(clusters(fm.nb))),
                             mu1 = exp(parameters(fm.nb)[1,1]),
                             sigma1 = parameters(fm.nb)[2,1],
                             mu2 = exp(parameters(fm.nb)[1,2]),
                             sigma2 = parameters(fm.nb)[2,2])
    row.names(fm.summary) <- NULL
    
    if (plot.trunc) {
      hist <- hist(Y.new, freq = F, breaks = "FD")
      lines(0:max(Y.new), 
            fm.summary$pi * dnbinom(0:max(Y.new), mu = fm.summary$mu1, size =  fm.summary$sigma1), 
            col = "blue", lwd = 2)
      lines(0:max(Y.new), 
            (1 - fm.summary$pi) * dnbinom(0:max(Y.new), mu = fm.summary$mu2, size = fm.summary$sigma2), 
            col = "blue", lwd = 2)
      lines(0:max(Y.new),
            fm.summary$pi * dnbinom(0:max(Y.new), mu = fm.summary$mu1, size = fm.summary$sigma1) + (1 - fm.summary$pi) * dnbinom(0:max(Y.new), mu = fm.summary$mu2, size = fm.summary$sigma2),
            col = "red", lwd = 2)
    } else {
      hist <- hist(Y, freq = F, breaks = "FD")
      lines(0:max(Y.new), 
            fm.summary$pi * dnbinom(0:max(Y.new), mu = fm.summary$mu1, size = 1 / fm.summary$sigma1), 
            col = "blue", lwd = 2)
      lines(0:max(Y.new), 
            (1 - fm.summary$pi) * dnbinom(0:max(Y.new), mu = fm.summary$mu2, size = 1 / fm.summary$sigma2), 
            col = "blue", lwd = 2)
      lines(0:max(Y.new),
            fm.summary$pi * dnbinom(0:max(Y.new), mu = fm.summary$mu1, size = 1 / fm.summary$sigma1) + (1 - fm.summary$pi) * dnbinom(0:max(Y.new), mu = fm.summary$mu2, size = 1 / fm.summary$sigma2),
            col = "red", lwd = 2)
    }

    return(fm.summary)
  }
  
  if (method == "rnor") {
    rob.est.nor <- lmrob(Y~1, method = "MM")
    print(summary(rob.est.nor))
    print(rob.est.nor)
    
    outlier.cutoff.min <- min(Y[which(rob.est.nor$rweights > 0)])
    outlier.cutoff.max <- max(Y[which(rob.est.nor$rweights > 0)])
    outlier.cutoff = paste0("min=", outlier.cutoff.min, "; max=", outlier.cutoff.max)
    rnor.summary <- data.frame(outlier.cutoff = outlier.cutoff,
                               pi = 0,
                               mu = rob.est.nor$coefficients,
                               sigma = rob.est.nor$coefficients ** 2 / (rob.est.nor$scale ** 2 - rob.est.nor$coefficients)
                              )
    
    Y.new <- Y[Y<=outlier.cutoff.max & Y>=outlier.cutoff.min]
    if (plot.trunc & plot) {
      hist <- hist(Y.new, freq = F, breaks = "FD")
      lines(min(Y.new):max(Y.new), dnbinom(min(Y.new):max(Y.new), mu = rnor.summary$mu, size = rnor.summary$sigma), col = "red", lwd = 2)
    } else if (plot) {
      hist <- hist(Y, freq = F, breaks = "FD")
      lines(min(Y.new):max(Y.new), dnbinom(min(Y.new):max(Y.new), mu = rnor.summary$mu, size = rnor.summary$sigma), col = "red", lwd = 2)
    }
    
    return(rnor.summary)
  }
  
  if (method=="nbnor") {
    # Outlier detection
    outlier.cutoff <- quantile(Y, probs = q.outliers)
    Y.new <- Y[Y<=outlier.cutoff]  
    
    # Mixture of NB and normal
    mix.nbnor <- flexmix(Y.new ~ 1 | 1, k = 2, model = list(FLXMRnegbin(theta = 3), FLXMRglm(family = "gaussian"))) 
    print(summary(mix.nbnor))
    print(parameters(mix.nbnor))
    nbnor.summary <- data.frame(outlier.cutoff = outlier.cutoff,
                                pi = tabulate(clusters(mix.nbnor))[1] / length(Y),
                                mu1 = parameters(mix.nbnor, component = 1)[[2]][1,1],
                                sigma1 = parameters(mix.nbnor, component = 1)[[2]][1,1] ** 2 / (parameters(mix.nbnor, component = 1)[[2]][2,1] ** 2 - parameters(mix.nbnor, component = 1)[[2]][1,1]),
                                mu2 = parameters(mix.nbnor, component = 2)[[2]][1,1],
                                sigma2 = parameters(mix.nbnor, component = 2)[[2]][1,1] **2 / (parameters(mix.nbnor, component = 2)[[2]][2,1] ** 2 - parameters(mix.nbnor, component = 2)[[2]][1,1])
    )
    
    if (plot.trunc) {
      hist <- hist(Y.new, freq = F, breaks = "FD")
    } else {
      hist <- hist(Y, freq = F, breaks = "FD")
    }
    lines(0:max(Y.new), nbnor.summary$pi * dnbinom(0:max(Y.new), mu = nbnor.summary$mu1, size = nbnor.summary$sigma1), col = "blue", lwd = 2)
    lines(0:max(Y.new), (1 - nbnor.summary$pi) * dnbinom(0:max(Y.new), mu = nbnor.summary$mu2, size = nbnor.summary$sigma2), col = "blue", lwd = 2)
    lines(0:max(Y.new), 
          nbnor.summary$pi * dnbinom(0:max(Y.new), mu = nbnor.summary$mu1, size = nbnor.summary$sigma1) + (1 - nbnor.summary$pi) * dnbinom(0:max(Y.new), mu = nbnor.summary$mu2, size = nbnor.summary$sigma2),
          col = "red",
          lwd = 2)
    
    return(nbnor.summary)
  }
}
