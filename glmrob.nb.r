glmrob.nb <- function(y,X,bounding.func='T/T',c.tukey.beta=5,c.tukey.sig=3,c.by.beta=4,weights.on.x='none',
                      a.hampel.beta=3,b.hampel.beta=4,c.hampel.beta=5,a.hampel.sig=2,b.hampel.sig=3,c.hampel.sig=4,
                      minsig=1e-2,maxsig=50,minmu=1e-10,maxmu=1e120,maxit=30,maxit.sig=50,sig.prec=1e-8,tol=1e-4,...){
  ### Written by William H. Aeberhard, February 2014
  ## Disclaimer: Users of these routines are cautioned that, while due care has been taken and they are
  ## believed accurate, they have not been rigorously tested and their use and results are
  ## solely the responsibilities of the user.
  #-------------------------------------------------------------------
  # General set up
  #-------------------------------------------------------------------
  n <- length(y)
  if (dim(X)[1]!=n){stop('length(y) does not match dim(X)[1].')}
  onevec <- rep(1,n)
  if (identical(X[,1],onevec)){X <- X[,-1]}
  res <- list()
  #-------------------------------------------------------------------
  # initial estimates: MLEs for beta and sigma
  #-------------------------------------------------------------------
  loglkhd <- function(sig,y,mu){
    sum(lgamma(y+1/sig)-lgamma(1/sig)-lgamma(y+1)-(1/sig)*log(sig*mu+1)+y*log(sig*mu/(sig*mu+1)))
  }
  score.sig.ML <- function(sig,y,mu){
    sum(digamma(y+1/sig)-digamma(1/sig)-log(sig*mu+1)-sig*(y-mu)/(sig*mu+1))
  }
  info.sig.ML <- function(sig,y,mu){
    (-1/sig^2)*(sum(trigamma(y+1/sig))-length(y)*trigamma(1/sig))-sum((sig*mu^2+y)/(sig*mu+1)^2)
  }
  invlink <- function(eta){exp(eta)}
  derivlink <- function(mu){1/mu}
  varfunc <- function(mu,sig){mu+sig*mu^2}
  # initial mu computed through Poisson GLM
  # glm.ini <- glm(y~X,family='poisson',...)
  if (dim(X)[2]==0) {
    glm.ini <- glm(y~1,family='poisson')
  } else {
    glm.ini <- glm(y~X,family='poisson',...)
  }
  eta <- glm.ini$lin
  mu <- invlink(eta)
  mu[which(mu>maxmu)] <- maxmu
  mu[which(mu<minmu)] <- minmu
  # sig MLE based on initial mu, with starting value = moment based
  sig <- sum((y/mu-1)^2)/length(y) 
  loglkhd.sig1 <- 0
  loglkhd.sig0 <- loglkhd.sig1+tol+1
  it.sig <- 0
  while (abs(loglkhd.sig1-loglkhd.sig0)>tol & it.sig<maxit){
    loglkhd.sig0 <- loglkhd(sig=sig,y=y,mu=mu)
    sig <- sig-score.sig.ML(sig=sig,y=y,mu=mu)/info.sig.ML(sig=sig,y=y,mu=mu)
    sig <- abs(sig)
    if (sig>maxsig){sig <- maxsig}
    loglkhd.sig1 <- loglkhd(sig=sig,y=y,mu=mu)
    it.sig <- it.sig+1
  }
  # ML estimations of mu and sig
  full.loglkhd1 <- 0
  full.loglkhd0 <- full.loglkhd1+tol+1
  it <- 0
  while(abs(full.loglkhd1-full.loglkhd0)>tol & it<maxit){
    full.loglkhd0 <- full.loglkhd1
    # mu
    loglkhd.mu1 <- 0
    loglkhd.mu0 <- loglkhd.mu1+tol+1
    it.mu <- 0
    while(abs(loglkhd.mu1-loglkhd.mu0)>tol & it.mu<maxit){
      loglkhd.mu0 <- loglkhd.mu1
      w.mat <- 1/((derivlink(mu))^2*varfunc(mu,sig))
      z <- eta+(y-mu)*derivlink(mu)
      # wls <- lm(z~X,weights=w.mat,...)
      wls <- lm(z~1,weights=w.mat,...)
      eta <- fitted(wls)
      mu <- invlink(eta)
      mu[which(mu>maxmu)] <- maxmu
      mu[which(mu<minmu)] <- minmu
      loglkhd.mu1 <- loglkhd(y=y,sig=sig,mu=mu)
      it.mu <- it.mu+1
    }
    # sigma
    loglkhd.sig1 <- 0
    loglkhd.sig0 <- loglkhd.sig1+tol+1
    it.sig <- 0
    while (abs(loglkhd.sig1-loglkhd.sig0)>tol & it.sig<maxit){
      loglkhd.sig0 <- loglkhd(sig=sig,y=y,mu=mu)
      sig <- sig-score.sig.ML(sig=sig,y=y,mu=mu)/info.sig.ML(sig=sig,y=y,mu=mu)
      sig <- abs(sig)
      if (sig>maxsig){sig <- maxsig}
      loglkhd.sig1 <- loglkhd(sig=sig,y=y,mu=mu)
      it.sig <- it.sig+1
    }
    full.loglkhd1 <- loglkhd(y=y,sig=sig,mu=mu)
    it <- it+1
  }
  #-------------------------------------------------------------------
  # Robust estimations
  #-------------------------------------------------------------------
  derivinvlink <- function(etai){exp(etai)}
  if (weights.on.x=='none'){
    weights.x <- onevec
  } else if (weights.on.x=='hard'){
    require(MASS) # for cov.rob
    Xrc <- cov.rob(X,quantile.used=floor(0.8*n))
    D2 <- mahalanobis(X,center=Xrc$center,cov=Xrc$cov) # copied from robustbase:::wts_RobDist
    qchi2 <- qchisq(p=0.95,df=dim(X)[2])
    weights.x <- ifelse(D2<=qchi2,1,0)
  } else {stop('Only "hard" and "none" are implemented for weights.on.x.')}
  derivinvlink <- function(eta){exp(eta)}
  psi.sig.ML <- function(r,mu,sig){
    digamma(r*sqrt(mu*(sig*mu+1))+mu+1/sig)-sig*r*sqrt(mu/(sig*mu+1))-digamma(1/sig)-log(sig*mu+1)
  }
  if (bounding.func=='T/T'){
    ### estimations
    tukeypsi <- function(r,c.tukey){
      ifelse(abs(r)>c.tukey,0,((r/c.tukey)^2-1)^2*r)
    }
    E.tukeypsi.1 <- function(mui,sig,c.tukey){
      sqrtVmui <- sqrt(varfunc(mui,sig))
      j1 <- max(c(ceiling(mui-c.tukey*sqrtVmui),0))
      j2 <- floor(mui+c.tukey*sqrtVmui)
      if (j1>j2){0}
      else {
        j12 <- j1:j2
        sum((((j12-mui)/(c.tukey*sqrtVmui))^2-1)^2*(j12-mui)*dnbinom(j12,mu=mui,size=1/sig))/sqrtVmui
      }
    }
    E.tukeypsi.2 <- function(mui,sig,c.tukey){
      sqrtVmui <- sqrt(varfunc(mui,sig))
      j1 <- max(c(ceiling(mui-c.tukey*sqrtVmui),0))
      j2 <- floor(mui+c.tukey*sqrtVmui)
      if (j1>j2){0}
      else {
        j12 <- j1:j2
        sum((((j12-mui)/(c.tukey*sqrtVmui))^2-1)^2*(j12-mui)^2*dnbinom(j12,mu=mui,size=1/sig))/sqrtVmui
      }
    }
    ai.sig.tukey <- function(mui,sig,c.tukey){
      psi.sig.ML.mod <- function(j,mui,invsig){
        digamma(j+invsig)-digamma(invsig)-log(mui/invsig+1)-(j-mui)/(mui+invsig)
      }
      sqrtVmui <- sqrt(mui*(sig*mui+1))
      invsig <- 1/sig
      j1 <- max(c(ceiling(mui-c.tukey*sqrtVmui),0))
      j2 <- floor(mui+c.tukey*sqrtVmui)
      if (j1>j2){0}
      else {
        j12 <- j1:j2
        sum((((j12-mui)/(c.tukey*sqrtVmui))^2-1)^2*psi.sig.ML.mod(j=j12,mui=mui,invsig=invsig)*dnbinom(x=j12,mu=mui,size=invsig))
      }
    }
    sig.rob.tukey <- function(sig,y,mu,c.tukey){
      r <- (y-mu)/sqrt(varfunc(mu,sig))
      wi <- tukeypsi(r=r,c.tukey=c.tukey)/r
      sum(wi*psi.sig.ML(r=r,mu=mu,sig=sig)-sapply(X=mu,FUN=ai.sig.tukey,sig=sig,c.tukey=c.tukey))
    }
    sig0 <- sig+1+tol
    beta11 <- 0
    beta00 <- beta11+tol+1
    it <- 0
    while(abs(sig-sig0)>tol | abs(max(beta11-beta00))>tol & it<maxit){
      sig0 <- sig
      beta00 <- beta11
      # estimate sigma given mu
      sig <- uniroot(f=sig.rob.tukey,interval=c(minsig,maxsig),tol=sig.prec,maxiter=maxit.sig,mu=mu,y=y,c.tukey=c.tukey.sig)$root
      if (is.na(sig) | sig>maxsig){sig <- maxsig}
      # estimate mu given sigma
      beta1 <- 0
      beta0 <- beta1+tol+1
      it.mu <- 0
      while(abs(max(beta1-beta0))>tol & it.mu<maxit){
        beta0 <- beta1
        bi <- sapply(X=mu,FUN=E.tukeypsi.2,sig=sig,c.tukey=c.tukey.beta)*varfunc(mu,sig)^(-3/2)*derivinvlink(eta)^2
        ei <- (tukeypsi(r=(y-mu)/sqrt(varfunc(mu,sig)),c.tukey=c.tukey.beta)-sapply(X=mu,FUN=E.tukeypsi.1,sig=sig,c.tukey=c.tukey.beta))/
          sapply(X=mu,FUN=E.tukeypsi.2,sig=sig,c.tukey=c.tukey.beta)*varfunc(mu,sig)*derivlink(mu)
        zi <- eta + ei
        # wls <- lm(zi~X,weights=bi,...)
        wls <- lm(zi~1,weights=bi,...)
        beta1 <- coef(wls)
        eta <- fitted(wls)
        mu <- invlink(eta)
        mu[which(mu>maxmu)] <- maxmu
        mu[which(mu==0)] <- minmu
        it.mu <- it.mu+1
      }
      beta11 <- beta1
      it <- it+1
    }
    res$coef <- c(sig,beta1)
    ### standard deviations
    fullscore.sig <- function(y,mui,sigma){
      (digamma(y+1/sigma)-digamma(1/sigma)-log(sigma*mui+1)-sigma*(y-mui)/(sigma*mui+1))/(-sigma^2)
    }
    all.expectations.tukey <- function(mui,sigma,c.tukey.beta,c.tukey.sigma){
      expec <- list()
      sqrtVmui <- sqrt(varfunc(mui,sigma))
      j1.beta <- max(c(ceiling(mui-c.tukey.beta*sqrtVmui),0))
      j2.beta <- floor(mui+c.tukey.beta*sqrtVmui)
      if (j1.beta>j2.beta){
        expec$tukeypsi2 <- 0
        expec$psibetascoresig.beta <- 0
        expec$tukeypsi13 <- 0
        expec$psibetaminuspsisig <- 0
        j1.sigma <- max(c(ceiling(mui-c.tukey.sigma*sqrtVmui),0))
        j2.sigma <- floor(mui+c.tukey.sigma*sqrtVmui)
        if (j1.sigma>j2.sigma){
          expec$psibetascoresig.sigma <- 0
          expec$psiscoresig2 <- 0
          expec$psiscoresig13 <- 0
        } else {
          j12.sigma <- j1.sigma:j2.sigma
          probNB.sigma <- dnbinom(j12.sigma,mu=mui,size=1/sigma)
          resi.sigma <- (j12.sigma-mui)/sqrtVmui
          tukeyresi.sigma <- tukeypsi(r=resi.sigma,c.tukey=c.tukey.sigma)
          fullscoresig.sigma <- fullscore.sig(y=j12.sigma,mui=mui,sigma=sigma)
          ## M21, M21=t(M12) if c.tukey.beta=c.tukey.sigma
          expec$psibetascoresig.sigma <- sum(tukeyresi.sigma*fullscoresig.sigma*probNB.sigma)/sqrtVmui # varfunc^(1/2) is from the rest of M21
          ## M22
          expec$psiscoresig2 <- sum(tukeyresi.sigma/resi.sigma*fullscoresig.sigma^2*probNB.sigma)
          ## Q22
          expec$psiscoresig13 <- sum((tukeyresi.sigma/resi.sigma*fullscoresig.sigma)^2*probNB.sigma)+
            -sum(tukeyresi.sigma/resi.sigma*fullscoresig.sigma*probNB.sigma)^2
        }
      } else {
        j12.beta <- j1.beta:j2.beta
        probNB.beta <- dnbinom(j12.beta,mu=mui,size=1/sigma)
        resi.beta <- (j12.beta-mui)/sqrtVmui
        tukeyresi.beta <- tukeypsi(r=resi.beta,c.tukey=c.tukey.beta)
        fullscoresig.beta <- fullscore.sig(y=j12.beta,mui=mui,sigma=sigma)
        ## M11
        expec$tukeypsi2 <- sum(tukeyresi.beta*(j12.beta-mui)*probNB.beta)/sqrtVmui^3 # varfunc^(-3/2) is from the rest of M11
        ## M12
        expec$psibetascoresig.beta <- sum(tukeyresi.beta*fullscoresig.beta*probNB.beta)/sqrtVmui # varfunc^(1/2) is from the rest of M12
        ## Q11
        expec$tukeypsi13 <- (sum(tukeyresi.beta^2*probNB.beta)-sum(tukeyresi.beta*probNB.beta)^2)/sqrtVmui^2 # varfunc^(-1) is from the rest of Q11
        j1.sigma <- max(c(ceiling(mui-c.tukey.sigma*sqrtVmui),0))
        j2.sigma <- floor(mui+c.tukey.sigma*sqrtVmui)
        if (j1.sigma>j2.sigma){
          expec$psibetascoresig.sigma <- 0
          expec$psiscoresig2 <- 0
          expec$psibetaminuspsisig <- 0
          expec$psiscoresig13 <- 0
        } else {
          j12.sigma <- j1.sigma:j2.sigma
          probNB.sigma <- dnbinom(j12.sigma,mu=mui,size=1/sigma)
          resi.sigma <- (j12.sigma-mui)/sqrtVmui
          tukeyresi.sigma <- tukeypsi(r=resi.sigma,c.tukey=c.tukey.sigma)
          fullscoresig.sigma <- fullscore.sig(y=j12.sigma,mui=mui,sigma=sigma)
          ## M21, M21=t(M12) if c.tukey.beta=c.tukey.sigma
          expec$psibetascoresig.sigma <- sum(tukeyresi.sigma*fullscoresig.sigma*probNB.sigma)/sqrtVmui # varfunc^(1/2) is from the rest of M21
          ## M22
          expec$psiscoresig2 <- sum(tukeyresi.sigma/resi.sigma*fullscoresig.sigma^2*probNB.sigma)
          ## Q12
          if (j2.beta<j2.sigma){ # we assume j1.beta=j1.sigma=0
            expec$psibetaminuspsisig <- (sum(tukeyresi.beta*tukeyresi.sigma[1:length(j12.beta)]/
                                               resi.sigma[1:length(j12.beta)]*fullscoresig.sigma[1:length(j12.beta)]*probNB.beta)+
                                           -sum(tukeyresi.beta*probNB.beta)*sum(tukeyresi.sigma/resi.sigma*fullscoresig.sigma*probNB.sigma))/sqrtVmui # varfunc^(1/2) is from the rest of Q12
          } else {
            expec$psibetaminuspsisig <- (sum(tukeyresi.beta[1:length(j12.sigma)]*tukeyresi.sigma/resi.sigma*fullscoresig.sigma*probNB.sigma)+
                                           -sum(tukeyresi.beta*probNB.beta)*sum(tukeyresi.sigma/resi.sigma*fullscoresig.sigma*probNB.sigma))/sqrtVmui # varfunc^(1/2) is from the rest of Q12
          }
          ## Q22
          expec$psiscoresig13 <- sum((tukeyresi.sigma/resi.sigma*fullscoresig.sigma)^2*probNB.sigma)+
            -sum(tukeyresi.sigma/resi.sigma*fullscoresig.sigma*probNB.sigma)^2
        }
      }
      return(expec)
    }
    X <- cbind(onevec,X)
    eta <- X%*%beta1
    mu <- invlink(eta)
    ## evaluate all expectations at each mui
    expect <- sapply(X=mu,FUN=all.expectations.tukey,sigma=sig,c.tukey.beta=c.tukey.beta,c.tukey.sigma=c.tukey.sig)
    ## compute all blocks of M and Q
    M11 <- t(X)%*%diag(as.numeric(unlist(expect['tukeypsi2',])*derivinvlink(eta)^2*weights.x))%*%X/n
    M12 <- as.numeric(t(X)%*%(unlist(expect['psibetascoresig.beta',])*derivinvlink(eta)*weights.x))/n
    M21 <- t(as.numeric(t(X)%*%(unlist(expect['psibetascoresig.sigma',])*derivinvlink(eta)*weights.x)))/n
    M22 <- t(weights.x)%*%unlist(expect['psiscoresig2',])/n
    Q11 <- t(X)%*%diag(as.numeric(unlist(expect['tukeypsi13',])*(derivinvlink(eta)*weights.x)^2))%*%X/n
    Q12 <- as.numeric(t(X)%*%(unlist(expect['psibetaminuspsisig',])*derivinvlink(eta)*weights.x^2))/n
    # Q21 <- t(Q12) # Q is always symmetric
    Q22 <- t(weights.x^2)%*%unlist(expect['psiscoresig13',])/n
    ## stdev from diag of full sandwich M^(-1)*Q*M^(-T)
    fullM <- rbind(cbind(M11,M12),t(c(M21,M22)))
    fullQ <- rbind(cbind(Q11,Q12),t(c(Q12,Q22)))
    stdev <- sqrt(diag(solve(fullM)%*%fullQ%*%solve(fullM))/n)
    res$stdev <- stdev[-length(stdev)]
    ### weights on response
    resid <- (y-mu)/sqrt(varfunc(mu,sig))
    res$weights.y <- as.numeric(ifelse(resid==0,1,tukeypsi(r=resid,c.tukey=c.tukey.beta)/resid))
    ### weights on design
    res$weights.x <- weights.x
  } else if (bounding.func=='BY/T'){
    deviance <- function(y,mu,sig){y*(log(y/mu)+log(sig*mu+1))-log(sig*y+1)*(y+1/sig)+log(sig*mu+1)/sig}
    findj0 <- function(mui,c.by,sigma,maxj=1e5){
      deviance.c <- function(j,mui,c.by,sig){j*(log(j/mui)+log(sig*mui+1))-log(sig*j+1)*(j+1/sig)+log(sig*mui+1)/sig-c.by}
      return(floor(uniroot(f=deviance.c,interval=c(mui,maxj),mui=mui,c.by=c.by,sig=sigma)$root))
    }
    ### estimations
    psi.by <- function(muyi,cc,sig){ # muyi = c(mui,yi)
      if (muyi[2]==0){
        devi <- log(sig*muyi[1]+1)/sig
        if (devi<=cc){
          1-devi/cc
        } else {0}
      } else {
        devi <- deviance(y=muyi[2],mu=muyi[1],sig=sig)
        if (devi<=cc){
          1-devi/cc
        } else {0}
      }
    }
    E.psi.by.1 <- function(mui,cc,sig){
      j0 <- findj0(mui=mui,c.by=cc,sigma=sig)
      if (j0>0){
        jrange <- 1:j0
        devi0 <- log(sig*mui+1)/sig
        if (devi0<=cc){
          term0 <- -(1-devi0/cc)*mui/(sig*mui+1)^(1/sig)
        } else {term0 <- 0}
        term0+sum((1-deviance(y=jrange,mu=mui,sig=sig)/cc)*(jrange-mui)*dnbinom(jrange,mu=mui,size=1/sig))
      } else {
        devi0 <- log(sig*mui+1)/sig
        if (devi0<=cc){
          -(1-devi0/cc)*mui/(sig*mui+1)^(1/sig)
        } else {0}
      }
    }
    E.psi.by.2 <- function(mui,cc,sig){
      j0 <- findj0(mui=mui,c.by=cc,sigma=sig)
      if (j0>0) {
        jrange <- 1:j0
        devi0 <- log(sig*mui+1)/sig
        if (devi0<=cc){
          term0 <- (1-devi0/cc)*mui^2/(sig*mui+1)^(1/sig)
        } else {term0 <- 0}
        term0+sum((1-deviance(y=jrange,mu=mui,sig=sig)/cc)*(jrange-mui)^2*dnbinom(jrange,mu=mui,size=1/sig))
      } else {
        devi0 <- log(sig*mui+1)/sig
        if (devi0<=cc){
          (1-devi0/cc)*mui^2/(sig*mui+1)^(1/sig)
        } else {0}
      }
    }
    tukeypsi <- function(r,c.tukey){
      ifelse(abs(r)>c.tukey,0,((r/c.tukey)^2-1)^2*r)
    }
    ai.sig.tukey <- function(mui,sig,c.tukey){
      psi.sig.ML.mod <- function(j,mui,invsig){
        digamma(j+invsig)-digamma(invsig)-log(mui/invsig+1)-(j-mui)/(mui+invsig)
      }
      sqrtVmui <- sqrt(mui*(sig*mui+1))
      invsig <- 1/sig
      j1 <- max(c(ceiling(mui-c.tukey*sqrtVmui),0))
      j2 <- floor(mui+c.tukey*sqrtVmui)
      if (j1>j2){0}
      else {
        j12 <- j1:j2
        sum((((j12-mui)/(c.tukey*sqrtVmui))^2-1)^2*psi.sig.ML.mod(j=j12,mui=mui,invsig=invsig)*dnbinom(x=j12,mu=mui,size=invsig))
      }
    }
    sig.rob.tukey <- function(sig,y,mu,c.tukey){
      r <- (y-mu)/sqrt(varfunc(mu,sig))
      wi <- tukeypsi(r=r,c.tukey=c.tukey)/r
      sum(wi*psi.sig.ML(r=r,mu=mu,sig=sig)-sapply(X=mu,FUN=ai.sig.tukey,sig=sig,c.tukey=c.tukey))
    }
    sig0 <- sig+1+tol
    beta11 <- 0
    beta00 <- beta11+tol+1
    it <- 0
    while(abs(sig-sig0)>tol | abs(max(beta11-beta00))>tol & it<maxit){
      sig0 <- sig
      beta00 <- beta11
      # estimate sigma given mu
      sig <- uniroot(f=sig.rob.tukey,interval=c(minsig,maxsig),tol=sig.prec,maxiter=maxit.sig,mu=mu,y=y,c.tukey=c.tukey.sig)$root
      if (is.na(sig) | sig>maxsig){sig <- maxsig}
      # estimate mu given sigma
      beta1 <- 0
      beta0 <- beta1+tol+1
      it.mu <- 0
      while(abs(max(beta1-beta0))>tol & it.mu<maxit){
        beta0 <- beta1
        bi <- sapply(X=mu,FUN=E.psi.by.2,cc=c.by.beta,sig=sig)/varfunc(mu,sig)^2*derivinvlink(eta)^2
        ei <- (apply(X=cbind(mu,y),MARGIN=1,FUN=psi.by,cc=c.by.beta,sig=sig)*(y-mu)-sapply(X=mu,FUN=E.psi.by.1,cc=c.by.beta,sig=sig))/sapply(X=mu,FUN=E.psi.by.2,cc=c.by.beta,sig=sig)*varfunc(mu,sig)*derivlink(mu)
        zi <- eta+ei
        wls <- lm(zi~X,weights=bi,...)
        beta1 <- coef(wls)
        eta <- fitted(wls)
        mu <- invlink(eta)
        mu[which(mu>maxmu)] <- maxmu
        mu[which(mu==0)] <- minmu
        it.mu <- it.mu+1
      }
      beta11 <- beta1
      it <- it+1
    }
    res$coef <- c(sig,beta1)
    ### standard deviations
    fullscore.sig <- function(y,mui,sigma){
      (digamma(y+1/sigma)-digamma(1/sigma)-log(sigma*mui+1)-sigma*(y-mui)/(sigma*mui+1))/(-sigma^2)
    }
    all.expectations.BY <- function(mui,sigma,c.by.b,c.tukey.s){
      expec <- list()
      sqrtVmui <- sqrt(varfunc(mui,sigma))
      invsig <- 1/sigma
      j0.b <- findj0(mui=mui,c.by=c.by.b,sigma=sigma)
      if (j0.b>0){
        jrange.b <- 0:j0.b
        pnb.b <- dnbinom(jrange.b,mu=mui,size=invsig)
        by.dev.b <- 1-c(log(sigma*mui+1)/sigma,deviance(y=jrange.b[-1],mu=mui,sig=sigma))/c.by.b
      } else if (j0.b==0){
        jrange.b <- 0
        pnb.b <- dnbinom(jrange.b,mu=mui,size=invsig)
        by.dev.b <- 1-(log(sigma*mui+1)/sigma)/c.by.b
      } else {
        jrange.b <- -1
        pnb.b <- 0
        by.dev.b <- 0
      }
      fullscoresig.b <- fullscore.sig(y=jrange.b,mui=mui,sigma=sigma)
      j1.s <- max(c(ceiling(mui-c.tukey.s*sqrtVmui),0))
      j2.s <- floor(mui+c.tukey.s*sqrtVmui)
      if (j1.s>j2.s){
        expec$psibetascoresig.sigma <- 0
        expec$psiscoresig2 <- 0
        expec$psiscoresig13 <- 0
      } else {
        j12.s <- j1.s:j2.s
        pnb.s <- dnbinom(j12.s,mu=mui,size=invsig)
        res.s <- (j12.s-mui)/sqrtVmui
        tukey.res.s <- tukeypsi(r=res.s,c.tukey=c.tukey.s)
        fullscoresig.s <- fullscore.sig(y=j12.s,mui=mui,sigma=sigma)
        ## M11
        expec$bypsi2 <- sum(by.dev.b*(jrange.b-mui)^2*pnb.b)/sqrtVmui^4 # varfunc^(-2) is from the rest of M11
        ## M12
        expec$psibetascoresig.beta <- sum(by.dev.b*(jrange.b-mui)*fullscoresig.b*pnb.b)/sqrtVmui^2 # varfunc^(-1) is from the rest of M12
        ## M21
        expec$psibetascoresig.sigma <- sum(tukey.res.s*fullscoresig.s*pnb.s)/sqrtVmui # varfunc^(-1/2) is from the rest of M21
        # M22
        expec$psiscoresig2 <- sum(tukey.res.s/res.s*fullscoresig.s^2*pnb.s)
        ## Q11
        expec$bypsi13 <- (sum(by.dev.b^2*(jrange.b-mui)^2*pnb.b)-sum(by.dev.b*(jrange.b-mui)*pnb.b)^2)/sqrtVmui^4 # varfunc^(-2) is from the rest of Q11
        ## Q12
        if (j0.b<=j2.s){ # we assume j1.s=0 and j0.b>0
          expec$psibetaminuspsisig <- (sum(by.dev.b*tukey.res.s[which(j12.s==0):which(j12.s==j0.b)]/
                                             res.s[which(j12.s==0):which(j12.s==j0.b)]*fullscoresig.s[which(j12.s==0):which(j12.s==j0.b)]*pnb.b)+
                                         -sum(by.dev.b*(jrange.b-mui)*pnb.b)*sum(tukey.res.s/res.s*fullscoresig.s*pnb.s))/sqrtVmui^2 # varfunc^(-1) is from the rest of Q12
        } else {
          expec$psibetaminuspsisig <- (sum(by.dev.b[1:(j2.s+1)]*tukey.res.s/res.s*fullscoresig.s*pnb.s)+
                                         -sum(by.dev.b*(jrange.b-mui)*pnb.b)*sum(tukey.res.s/res.s*fullscoresig.s*pnb.s))/sqrtVmui^2 # varfunc^(-1) is from the rest of Q12
        }
        ## Q22
        expec$psiscoresig13 <- sum((tukey.res.s/res.s*fullscoresig.s)^2*pnb.s)+
          -sum(tukey.res.s/res.s*fullscoresig.s*pnb.s)^2
      }
      return(expec)
    }
    X <- cbind(onevec,X)
    eta <- X%*%beta1
    mu <- invlink(eta)
    expect <- sapply(X=mu,FUN=all.expectations.BY,sigma=sig,c.by.b=c.by.beta,c.tukey.s=c.tukey.sig)
    ## compute all blocks of M and Q
    M11 <- t(X)%*%diag(as.numeric(unlist(expect['bypsi2',])*derivinvlink(eta)^2*weights.x))%*%X/n
    M12 <- as.numeric(t(X)%*%(unlist(expect['psibetascoresig.beta',])*derivinvlink(eta)*weights.x))/n
    M21 <- t(as.numeric(t(X)%*%(unlist(expect['psibetascoresig.sigma',])*derivinvlink(eta)*weights.x)))/n
    M22 <- t(weights.x)%*%unlist(expect['psiscoresig2',])/n
    Q11 <- t(X)%*%diag(as.numeric(unlist(expect['bypsi13',])*(derivinvlink(eta)*weights.x)^2))%*%X/n
    Q12 <- as.numeric(t(X)%*%(unlist(expect['psibetaminuspsisig',])*derivinvlink(eta)*weights.x^2))/n
    # Q21 <- t(Q12) # Q is always symmetric
    Q22 <- t(weights.x^2)%*%unlist(expect['psiscoresig13',])/n
    ## stdev from diag of full sandwich M^(-1)*Q*M^(-T)
    fullM <- rbind(cbind(M11,M12),t(c(M21,M22)))
    fullQ <- rbind(cbind(Q11,Q12),t(c(Q12,Q22)))
    stdev <- sqrt(diag(solve(fullM)%*%fullQ%*%solve(fullM))/n)
    res$stdev <- stdev[-length(stdev)]
    ### weights on response
    res$weights.y <- as.numeric(apply(X=cbind(mu,y),MARGIN=1,FUN=psi.by,cc=c.by.beta,sig=sig))
    ### weights on design
    res$weights.x <- weights.x  
  } else if (bounding.func=='HA/HA'){
    ### estimations
    hampelpsi <- function(r,a.hampel,b.hampel,c.hampel){
      res <- r
      absr <- abs(r)
      between.b.c <- which(absr<=c.hampel & absr>b.hampel)
      between.a.b <- which(absr<=b.hampel & absr>a.hampel)
      res[which(absr>c.hampel)] <- 0
      res[between.b.c] <- a.hampel*(c.hampel-absr[between.b.c])/(c.hampel-b.hampel)*sign(r[between.b.c])
      res[between.a.b] <- sign(r[between.a.b])*a.hampel
      res
    }
    E.hampelpsi.1 <- function(mui,sig,a.hampel,b.hampel,c.hampel){
      sqrtVmui <- sqrt(varfunc(mui,sig))
      invsig <- 1/sig
      j1 <- ceiling(mui-c.hampel*sqrtVmui)
      j2 <- ceiling(mui-b.hampel*sqrtVmui)
      j3 <- ceiling(mui-a.hampel*sqrtVmui)
      j4 <- floor(mui+a.hampel*sqrtVmui)
      j5 <- floor(mui+b.hampel*sqrtVmui)
      j6 <- floor(mui+c.hampel*sqrtVmui)
      sum.jminusmu.1 <- function(u,v,mui,sig,invsig){
        mui*(dnbinom(u-1,mu=mui,size=invsig)*(sig*u-sig+1)-dnbinom(v,mu=mui,size=invsig)*(sig*v+1))
      }
      a.hampel/(b.hampel-c.hampel)*(c.hampel*(pnbinom(j2-1,mu=mui,size=invsig)-pnbinom(j1-1,mu=mui,size=invsig)+
                                                -pnbinom(j6,mu=mui,size=invsig)+pnbinom(j5,mu=mui,size=invsig))+
                                      +(sum.jminusmu.1(u=j1,v=j2-1,mui=mui,sig=sig,invsig=invsig)+
                                          +sum.jminusmu.1(u=j5+1,v=j6,mui=mui,sig=sig,invsig=invsig))/sqrtVmui)+
        +a.hampel*(pnbinom(j5,mu=mui,size=invsig)-pnbinom(j4,mu=mui,size=invsig)+
                     -pnbinom(j3-1,mu=mui,size=invsig)+pnbinom(j2-1,mu=mui,size=invsig))+
        +sum.jminusmu.1(u=j3,v=j4,mui=mui,sig=sig,invsig=invsig)/sqrtVmui
    }
    E.hampelpsi.2 <- function(mui,sig,a.hampel,b.hampel,c.hampel){
      sqrtVmui <- sqrt(varfunc(mui,sig))
      invsig <- 1/sig
      j1 <- ceiling(mui-c.hampel*sqrtVmui)
      j2 <- ceiling(mui-b.hampel*sqrtVmui)
      j3 <- ceiling(mui-a.hampel*sqrtVmui)
      j4 <- floor(mui+a.hampel*sqrtVmui)
      j5 <- floor(mui+b.hampel*sqrtVmui)
      j6 <- floor(mui+c.hampel*sqrtVmui)
      sum.jminusmu.1 <- function(u,v,mui,sig,invsig){
        mui*(dnbinom(u-1,mu=mui,size=invsig)*(sig*u-sig+1)-dnbinom(v,mu=mui,size=invsig)*(sig*v+1))
      }
      sum.jminusmu.2 <- function(u,v,mui,sig,invsig){
        (sig*mui^2+mui)*(pnbinom(v-2,mu=mui,size=invsig)-pnbinom(u-3,mu=mui,size=invsig))+
          -dnbinom(v,mu=mui,size=invsig)*(sig*mui*v^2-2*sig*mui^2*v-mui^2)+
          +dnbinom(u-1,mu=mui,size=invsig)*(sig*mui*(u-1)^2-2*sig*mui^2*(u-1)-mui^2)+
          -dnbinom(v-1,mu=mui,size=invsig)*(sig*mui^2*(sig+1)*(v-1)-mui+mui^2)+
          +dnbinom(u-2,mu=mui,size=invsig)*(sig*mui^2*(sig+1)*(u-2)-mui+mui^2)
      }
      a.hampel/(b.hampel-c.hampel)*(c.hampel*(sum.jminusmu.1(u=j1,v=j2-1,mui=mui,sig=sig,invsig=invsig)+
                                                -sum.jminusmu.1(u=j5+1,v=j6,mui=mui,sig=sig,invsig=invsig))+
                                      +(sum.jminusmu.2(u=j1,v=j2-1,mui=mui,sig=sig,invsig=invsig)+
                                          +sum.jminusmu.2(u=j5+1,v=j6,mui=mui,sig=sig,invsig=invsig))/sqrtVmui)+
        +a.hampel*(sum.jminusmu.1(u=j4+1,v=j5,mui=mui,sig=sig,invsig=invsig)+
                     -sum.jminusmu.1(u=j2,v=j3-1,mui=mui,sig=sig,invsig=invsig))+
        +sum.jminusmu.2(u=j3,v=j4,mui=mui,sig=sig,invsig=invsig)/sqrtVmui
    }
    ai.sig.hampel <- function(mui,sig,a.hampel,b.hampel,c.hampel){
      psi.sig.ML.mod <- function(j,mui,invsig){
        digamma(j+invsig)-digamma(invsig)-log(mui/invsig+1)-(j-mui)/(mui+invsig)
      }
      sqrtVmui <- sqrt(mui*(sig*mui+1))
      invsig <- 1/sig
      j1 <- max(c(ceiling(mui-c.hampel*sqrtVmui),0))
      j2 <- max(c(ceiling(mui-b.hampel*sqrtVmui),0))
      j3 <- max(c(ceiling(mui-a.hampel*sqrtVmui),0))
      j4 <- floor(mui+a.hampel*sqrtVmui)
      j5 <- floor(mui+b.hampel*sqrtVmui)
      j6 <- floor(mui+c.hampel*sqrtVmui)
      if (j1==j2){
        if (j2==j3){
          if (j4==j5){
            if (j5==j6){
              j34 <- j3:j4
              sum(psi.sig.ML.mod(j=j34,mui=mui,invsig=invsig)*dnbinom(x=j34,mu=mui,size=invsig))
            } else {
              j34 <- j3:j4
              j56 <- (j5+1):j6
              sum(psi.sig.ML.mod(j=j34,mui=mui,invsig=invsig)*dnbinom(x=j34,mu=mui,size=invsig))+
                +sum((c.hampel*sqrtVmui/(j56-mui)-1)*a.hampel/(c.hampel-b.hampel)*psi.sig.ML.mod(j=j56,mui=mui,invsig=invsig)*dnbinom(x=j56,mu=mui,size=invsig))
            }
          } else {
            if (j5==j6){
              j34 <- j3:j4
              j45 <- (j4+1):j5
              a.hampel*sqrtVmui*sum(psi.sig.ML.mod(j=j45,mui=mui,invsig=invsig)/(j45-mui)*dnbinom(x=j45,mu=mui,size=invsig))+
                +sum(psi.sig.ML.mod(j=j34,mui=mui,invsig=invsig)*dnbinom(x=j34,mu=mui,size=invsig))
            } else {
              j34 <- j3:j4
              j45 <- (j4+1):j5
              j56 <- (j5+1):j6
              a.hampel*sqrtVmui*sum(psi.sig.ML.mod(j=j45,mui=mui,invsig=invsig)/(j45-mui)*dnbinom(x=j45,mu=mui,size=invsig))+
                +sum(psi.sig.ML.mod(j=j34,mui=mui,invsig=invsig)*dnbinom(x=j34,mu=mui,size=invsig))+
                +sum((c.hampel*sqrtVmui/(j56-mui)-1)*a.hampel/(c.hampel-b.hampel)*psi.sig.ML.mod(j=j56,mui=mui,invsig=invsig)*dnbinom(x=j56,mu=mui,size=invsig))
            }
          }
        } else {
          if (j4==j5){
            if (j5==j6){
              j23 <- j2:(j3-1)
              j34 <- j3:j4
              -a.hampel*sqrtVmui*sum(psi.sig.ML.mod(j=j23,mui=mui,invsig=invsig)/(j23-mui)*dnbinom(x=j23,mu=mui,size=invsig))+
                +sum(psi.sig.ML.mod(j=j34,mui=mui,invsig=invsig)*dnbinom(x=j34,mu=mui,size=invsig))
            } else {
              j23 <- j2:(j3-1)
              j34 <- j3:j4
              j56 <- (j5+1):j6
              -a.hampel*sqrtVmui*sum(psi.sig.ML.mod(j=j23,mui=mui,invsig=invsig)/(j23-mui)*dnbinom(x=j23,mu=mui,size=invsig))+
                +sum(psi.sig.ML.mod(j=j34,mui=mui,invsig=invsig)*dnbinom(x=j34,mu=mui,size=invsig))+
                +sum((c.hampel*sqrtVmui/(j56-mui)-1)*a.hampel/(c.hampel-b.hampel)*psi.sig.ML.mod(j=j56,mui=mui,invsig=invsig)*dnbinom(x=j56,mu=mui,size=invsig))
            }
          } else {
            if (j5==j6){
              j23 <- j2:(j3-1)
              j34 <- j3:j4
              j45 <- (j4+1):j5
              a.hampel*sqrtVmui*(sum(psi.sig.ML.mod(j=j45,mui=mui,invsig=invsig)/(j45-mui)*dnbinom(x=j45,mu=mui,size=invsig))+
                                   -sum(psi.sig.ML.mod(j=j23,mui=mui,invsig=invsig)/(j23-mui)*dnbinom(x=j23,mu=mui,size=invsig)))+
                +sum(psi.sig.ML.mod(j=j34,mui=mui,invsig=invsig)*dnbinom(x=j34,mu=mui,size=invsig))
            } else {
              j23 <- j2:(j3-1)
              j34 <- j3:j4
              j45 <- (j4+1):j5
              j56 <- (j5+1):j6
              a.hampel*sqrtVmui*(sum(psi.sig.ML.mod(j=j45,mui=mui,invsig=invsig)/(j45-mui)*dnbinom(x=j45,mu=mui,size=invsig))+
                                   -sum(psi.sig.ML.mod(j=j23,mui=mui,invsig=invsig)/(j23-mui)*dnbinom(x=j23,mu=mui,size=invsig)))+
                +sum(psi.sig.ML.mod(j=j34,mui=mui,invsig=invsig)*dnbinom(x=j34,mu=mui,size=invsig))+
                +sum((c.hampel*sqrtVmui/(j56-mui)-1)*a.hampel/(c.hampel-b.hampel)*psi.sig.ML.mod(j=j56,mui=mui,invsig=invsig)*dnbinom(x=j56,mu=mui,size=invsig))
            }
          }
        }
      } else {
        if (j2==j3){
          if (j4==j5){
            if (j5==j6){
              j12 <- j1:(j2-1)
              j34 <- j3:j4
              sum((c.hampel*sqrtVmui/(j12-mui)+1)*a.hampel/(b.hampel-c.hampel)*psi.sig.ML.mod(j=j12,mui=mui,invsig=invsig)*dnbinom(x=j12,mu=mui,size=invsig))+
                +sum(psi.sig.ML.mod(j=j34,mui=mui,invsig=invsig)*dnbinom(x=j34,mu=mui,size=invsig))
            } else {
              j12 <- j1:(j2-1)
              j34 <- j3:j4
              j56 <- (j5+1):j6
              sum((c.hampel*sqrtVmui/(j12-mui)+1)*a.hampel/(b.hampel-c.hampel)*psi.sig.ML.mod(j=j12,mui=mui,invsig=invsig)*dnbinom(x=j12,mu=mui,size=invsig))+
                +sum(psi.sig.ML.mod(j=j34,mui=mui,invsig=invsig)*dnbinom(x=j34,mu=mui,size=invsig))+
                +sum((c.hampel*sqrtVmui/(j56-mui)-1)*a.hampel/(c.hampel-b.hampel)*psi.sig.ML.mod(j=j56,mui=mui,invsig=invsig)*dnbinom(x=j56,mu=mui,size=invsig))
            }
          } else {
            if (j5==j6){
              j12 <- j1:(j2-1)
              j34 <- j3:j4
              j45 <- (j4+1):j5
              sum((c.hampel*sqrtVmui/(j12-mui)+1)*a.hampel/(b.hampel-c.hampel)*psi.sig.ML.mod(j=j12,mui=mui,invsig=invsig)*dnbinom(x=j12,mu=mui,size=invsig))+
                +a.hampel*sqrtVmui*sum(psi.sig.ML.mod(j=j45,mui=mui,invsig=invsig)/(j45-mui)*dnbinom(x=j45,mu=mui,size=invsig))+
                +sum(psi.sig.ML.mod(j=j34,mui=mui,invsig=invsig)*dnbinom(x=j34,mu=mui,size=invsig))
            } else {
              j12 <- j1:(j2-1)
              j34 <- j3:j4
              j45 <- (j4+1):j5
              j56 <- (j5+1):j6
              sum((c.hampel*sqrtVmui/(j12-mui)+1)*a.hampel/(b.hampel-c.hampel)*psi.sig.ML.mod(j=j12,mui=mui,invsig=invsig)*dnbinom(x=j12,mu=mui,size=invsig))+
                +a.hampel*sqrtVmui*sum(psi.sig.ML.mod(j=j45,mui=mui,invsig=invsig)/(j45-mui)*dnbinom(x=j45,mu=mui,size=invsig))+
                +sum(psi.sig.ML.mod(j=j34,mui=mui,invsig=invsig)*dnbinom(x=j34,mu=mui,size=invsig))+
                +sum((c.hampel*sqrtVmui/(j56-mui)-1)*a.hampel/(c.hampel-b.hampel)*psi.sig.ML.mod(j=j56,mui=mui,invsig=invsig)*dnbinom(x=j56,mu=mui,size=invsig))
            }
          }
        } else {
          if (j4==j5){
            if (j5==j6){
              j12 <- j1:(j2-1)
              j23 <- j2:(j3-1)
              j34 <- j3:j4
              sum((c.hampel*sqrtVmui/(j12-mui)+1)*a.hampel/(b.hampel-c.hampel)*psi.sig.ML.mod(j=j12,mui=mui,invsig=invsig)*dnbinom(x=j12,mu=mui,size=invsig))+
                -a.hampel*sqrtVmui*sum(psi.sig.ML.mod(j=j23,mui=mui,invsig=invsig)/(j23-mui)*dnbinom(x=j23,mu=mui,size=invsig))+
                +sum(psi.sig.ML.mod(j=j34,mui=mui,invsig=invsig)*dnbinom(x=j34,mu=mui,size=invsig))
            } else {
              j12 <- j1:(j2-1)
              j23 <- j2:(j3-1)
              j34 <- j3:j4
              j56 <- (j5+1):j6
              sum((c.hampel*sqrtVmui/(j12-mui)+1)*a.hampel/(b.hampel-c.hampel)*psi.sig.ML.mod(j=j12,mui=mui,invsig=invsig)*dnbinom(x=j12,mu=mui,size=invsig))+
                -a.hampel*sqrtVmui*sum(psi.sig.ML.mod(j=j23,mui=mui,invsig=invsig)/(j23-mui)*dnbinom(x=j23,mu=mui,size=invsig))+
                +sum(psi.sig.ML.mod(j=j34,mui=mui,invsig=invsig)*dnbinom(x=j34,mu=mui,size=invsig))+
                +sum((c.hampel*sqrtVmui/(j56-mui)-1)*a.hampel/(c.hampel-b.hampel)*psi.sig.ML.mod(j=j56,mui=mui,invsig=invsig)*dnbinom(x=j56,mu=mui,size=invsig))
            }
          } else {
            if (j5==j6){
              j12 <- j1:(j2-1)
              j23 <- j2:(j3-1)
              j34 <- j3:j4
              j45 <- (j4+1):j5
              sum((c.hampel*sqrtVmui/(j12-mui)+1)*a.hampel/(b.hampel-c.hampel)*psi.sig.ML.mod(j=j12,mui=mui,invsig=invsig)*dnbinom(x=j12,mu=mui,size=invsig))+
                +a.hampel*sqrtVmui*(sum(psi.sig.ML.mod(j=j45,mui=mui,invsig=invsig)/(j45-mui)*dnbinom(x=j45,mu=mui,size=invsig))+
                                      -sum(psi.sig.ML.mod(j=j23,mui=mui,invsig=invsig)/(j23-mui)*dnbinom(x=j23,mu=mui,size=invsig)))+
                +sum(psi.sig.ML.mod(j=j34,mui=mui,invsig=invsig)*dnbinom(x=j34,mu=mui,size=invsig))
            } else {
              j12 <- j1:(j2-1)
              j23 <- j2:(j3-1)
              j34 <- j3:j4
              j45 <- (j4+1):j5
              j56 <- (j5+1):j6
              sum((c.hampel*sqrtVmui/(j12-mui)+1)*a.hampel/(b.hampel-c.hampel)*psi.sig.ML.mod(j=j12,mui=mui,invsig=invsig)*dnbinom(x=j12,mu=mui,size=invsig))+
                +a.hampel*sqrtVmui*(sum(psi.sig.ML.mod(j=j45,mui=mui,invsig=invsig)/(j45-mui)*dnbinom(x=j45,mu=mui,size=invsig))+
                                      -sum(psi.sig.ML.mod(j=j23,mui=mui,invsig=invsig)/(j23-mui)*dnbinom(x=j23,mu=mui,size=invsig)))+
                +sum(psi.sig.ML.mod(j=j34,mui=mui,invsig=invsig)*dnbinom(x=j34,mu=mui,size=invsig))+
                +sum((c.hampel*sqrtVmui/(j56-mui)-1)*a.hampel/(c.hampel-b.hampel)*psi.sig.ML.mod(j=j56,mui=mui,invsig=invsig)*dnbinom(x=j56,mu=mui,size=invsig))
            }
          }
        }
      }
    }
    sig.rob.hampel <- function(sig,y,mu,a.hampel,b.hampel,c.hampel){
      r <- (y-mu)/sqrt(varfunc(mu,sig))
      wi <- hampelpsi(r=r,a.hampel=a.hampel,b.hampel=b.hampel,c.hampel=c.hampel)/r
      sum(wi*psi.sig.ML(r=r,mu=mu,sig=sig)-sapply(X=mu,FUN=ai.sig.hampel,sig=sig,a.hampel=a.hampel,b.hampel=b.hampel,c.hampel=c.hampel))
    }
    sig0 <- sig+1+tol
    beta11 <- 0
    beta00 <- beta11+tol+1
    it <- 0
    while(abs(sig-sig0)>tol | abs(max(beta11-beta00))>tol & it<maxit){
      sig0 <- sig
      beta00 <- beta11
      # estimate sigma given mu
      sig <- uniroot(f=sig.rob.hampel,interval=c(minsig,maxsig),tol=sig.prec,maxiter=maxit.sig,mu=mu,y=y,a.hampel=a.hampel.sig,b.hampel=b.hampel.sig,c.hampel=c.hampel.sig)$root
      if (is.na(sig) | sig>maxsig){sig <- maxsig}
      # estimate mu given sigma
      beta1 <- 0
      beta0 <- beta1+tol+1
      it.mu <- 0
      while(abs(max(beta1-beta0))>tol & it.mu<maxit){
        beta0 <- beta1
        bi <- sapply(X=mu,FUN=E.hampelpsi.2,sig=sig,a.hampel=a.hampel.beta,b.hampel=b.hampel.beta,c.hampel=c.hampel.beta)*varfunc(mu,sig)^(-3/2)*derivinvlink(eta)^2
        ei <- (hampelpsi(r=(y-mu)/sqrt(mu*(sig*mu+1)),a.hampel=a.hampel.beta,b.hampel=b.hampel.beta,c.hampel=c.hampel.beta)-sapply(X=mu,FUN=E.hampelpsi.1,sig=sig,a.hampel=a.hampel.beta,b.hampel=b.hampel.beta,c.hampel=c.hampel.beta))/
          sapply(X=mu,FUN=E.hampelpsi.2,sig=sig,a.hampel=a.hampel.beta,b.hampel=b.hampel.beta,c.hampel=c.hampel.beta)*varfunc(mu,sig)*derivlink(mu)
        zi <- eta+ei
        wls <- lm(zi~X,weights=bi,...)
        beta1 <- coef(wls)
        eta <- fitted(wls)
        mu <- invlink(eta)
        mu[which(mu>maxmu)] <- maxmu
        mu[which(mu==0)] <- minmu
        it.mu <- it.mu+1
      }
      beta11 <- beta1
      it <- it+1
    }
    res$coef <- c(sig,beta1)
    ### standard deviations
    fullscore.sig <- function(y,mui,sigma){
      (digamma(y+1/sigma)-digamma(1/sigma)-log(sigma*mui+1)-sigma*(y-mui)/(sigma*mui+1))/(-sigma^2)
    }
    hampelpsi <- function(r,tuning){
      res <- r
      absr <- abs(r)
      between.b.c <- which(absr<=tuning[3] & absr>tuning[2])
      between.a.b <- which(absr<=tuning[2] & absr>tuning[1])
      res[which(absr>tuning[3])] <- 0
      res[between.b.c] <- tuning[1]*(tuning[3]-absr[between.b.c])/(tuning[3]-tuning[2])*sign(r[between.b.c])
      res[between.a.b] <- sign(r[between.a.b])*tuning[1]
      return(res)
    }
    all.expectations.hampel <- function(mui,sigma,tuning.b,tuning.s){
      sqrtVmui <- sqrt(varfunc(mui,sigma))
      invsig <- 1/sigma
      j1.b <- max(c(ceiling(mui-tuning.b[3]*sqrtVmui),0))
      j6.b <- floor(mui+tuning.b[3]*sqrtVmui)
      j1.s <- max(c(ceiling(mui-tuning.s[3]*sqrtVmui),0))
      j6.s <- floor(mui+tuning.s[3]*sqrtVmui)
      # setup for beta
      j16.b <- j1.b:j6.b
      res.b <- (j16.b-mui)/sqrtVmui
      pnb16.b <- dnbinom(j16.b,mu=mui,size=invsig)
      hampel.res.b <- hampelpsi(r=res.b,tuning=tuning.b)
      fullscoresig.b <- fullscore.sig(y=j16.b,mui=mui,sigma=sigma)
      # setup for sigma
      j16.s <- j1.s:j6.s
      res.s <- (j16.s-mui)/sqrtVmui
      pnb16.s <- dnbinom(j16.s,mu=mui,size=invsig)
      hampel.res.s <- hampelpsi(r=res.s,tuning=tuning.s)
      fullscoresig.s <- fullscore.sig(y=j16.s,mui=mui,sigma=sigma)
      expec <- list()
      ## M11
      expec$hampelpsi2 <- sum(hampel.res.b*(j16.b-mui)*pnb16.b)/sqrtVmui^3 # varfunc^(-3/2) is from the rest of M11
      ## M12
      expec$psibetascoresig.beta <- sum(hampel.res.b*fullscoresig.b*pnb16.b)/sqrtVmui # varfunc^(1/2) is from the rest of M12
      ## M21, M21=t(M12) if tuning.b=tuning.s
      expec$psibetascoresig.sigma <- sum(hampel.res.s*fullscoresig.s*pnb16.s)/sqrtVmui # varfunc^(1/2) is from the rest of M21
      ## M22
      expec$psiscoresig2 <- sum(hampel.res.s*sqrtVmui/(j16.s-mui)*fullscoresig.s^2*pnb16.s)
      ## Q11
      expec$hampelpsi13 <- (sum(hampel.res.b^2*pnb16.b)-sum(hampel.res.b*pnb16.b)^2)/sqrtVmui^2 # varfunc^(-1) is from the rest of Q11
      ## Q12
      if (j6.b<=j6.s){ # we fill the remaining elements of hampel.res.b with zeros, to match length(hampel.res.s)
        if (j1.b<=j1.s){
          expec$psibetaminuspsisig <- (sum(c(rep(0,j1.s-j1.b),hampel.res.b,rep(0,j6.s-j6.b))*hampel.res.s/res.s*fullscoresig.s*pnb16.s)+
                                         -sum(hampel.res.b*pnb16.b)*sum(hampel.res.s/res.s*fullscoresig.s*pnb16.s))/sqrtVmui
        } else {
          expec$psibetaminuspsisig <- (sum(c(hampel.res.b,rep(0,j6.s-j6.b))*c(rep(0,j1.b-j1.s),hampel.res.s)/c(rep(1,j1.b-j1.s),res.s)*
                                             c(rep(0,j1.b-j1.s),fullscoresig.s)*c(rep(0,j1.b-j1.s),pnb16.s))+
                                         -sum(hampel.res.b*pnb16.b)*sum(hampel.res.s/res.s*fullscoresig.s*pnb16.s))/sqrtVmui
        }
      } else {
        if (j1.b<=j1.s){
          expec$psibetaminuspsisig <- (sum(c(rep(0,j1.s-j1.b),hampel.res.b)*c(hampel.res.s,rep(0,j6.b-j6.s))/c(res.s,rep(1,j6.b-j6.s))*
                                             c(fullscoresig.s,rep(0,j6.b-j6.s))*c(pnb16.s,rep(0,j6.b-j6.s)))+
                                         -sum(hampel.res.b*pnb16.b)*sum(hampel.res.s/res.s*fullscoresig.s*pnb16.s))/sqrtVmui
        } else {
          expec$psibetaminuspsisig <- (sum(hampel.res.b*c(rep(0,j1.b-j1.s),hampel.res.s,rep(0,j6.b-j6.s))/
                                             c(rep(1,j1.b-j1.s),res.s,rep(1,j6.b-j6.s))*
                                             c(rep(0,j1.b-j1.s),fullscoresig.s,rep(0,j6.b-j6.s))*pnb16.b)+
                                         -sum(hampel.res.b*pnb16.b)*sum(hampel.res.s/res.s*fullscoresig.s*pnb16.s))/sqrtVmui
        }
      }
      ## Q22
      expec$psiscoresig13 <- sum((hampel.res.s/res.s*fullscoresig.s)^2*pnb16.s)-sum(hampel.res.s/res.s*fullscoresig.s*pnb16.s)^2
      return(expec)
    }
    X <- cbind(onevec,X)
    eta <- X%*%beta1
    mu <- invlink(eta)
    is.whole <- function(x,tol=1e-8){abs(x-round(x))<tol}
    mu <- ifelse(is.whole(mu),mu+1e-8,mu) # we don't allow exact integers for mu, otherwise may divide by 0
    expect <- sapply(X=mu,FUN=all.expectations.hampel,sigma=sig,tuning.b=c(a.hampel.beta,b.hampel.beta,c.hampel.beta),
                     tuning.s=c(a.hampel.sig,b.hampel.sig,c.hampel.sig))
    ## compute all blocks of M and Q
    M11 <- t(X)%*%diag(as.numeric(unlist(expect['hampelpsi2',])*derivinvlink(eta)^2*weights.x))%*%X/n
    M12 <- as.numeric(t(X)%*%(unlist(expect['psibetascoresig.beta',])*derivinvlink(eta)*weights.x))/n
    M21 <- t(as.numeric(t(X)%*%(unlist(expect['psibetascoresig.sigma',])*derivinvlink(eta)*weights.x)))/n
    M22 <- t(weights.x)%*%unlist(expect['psiscoresig2',])/n
    Q11 <- t(X)%*%diag(as.numeric(unlist(expect['hampelpsi13',])*(derivinvlink(eta)*weights.x)^2))%*%X/n
    Q12 <- as.numeric(t(X)%*%(unlist(expect['psibetaminuspsisig',])*derivinvlink(eta)*weights.x^2))/n
    # Q21 <- t(Q12) # Q is always symmetric
    Q22 <- t(weights.x^2)%*%unlist(expect['psiscoresig13',])/n
    ## stdev from diag of full sandwich M^(-1)*Q*M^(-T)
    fullM <- rbind(cbind(M11,M12),t(c(M21,M22)))
    fullQ <- rbind(cbind(Q11,Q12),t(c(Q12,Q22)))
    # return(list('M'=fullM,'Q'=fullQ))
    stdev <- sqrt(diag(solve(fullM)%*%fullQ%*%solve(fullM))/n)
    res$stdev <- stdev[-length(stdev)]
    ### weights on response
    resid <- (y-mu)/sqrt(varfunc(mu,sig))
    res$weights.y <- as.numeric(ifelse(resid==0,1,hampelpsi(r=resid,tuning=c(a.hampel.beta,b.hampel.beta,c.hampel.beta))/resid))
    ### weights on design
    res$weights.x <- weights.x
  } else {stop('Available bounding.func are "T/T", "BY/T" and "HA/HA".')}
  return(res)
}
#--- END glmrob.nb
