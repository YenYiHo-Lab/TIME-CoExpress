
library(copula)
library(gamlss.dist)
library(GJRM)
library(pracma)


dgamma_musig <- function(x, mu, sigma){
  
  shape_use <- 1/sigma^2
  
  scale_use <- mu*sigma^2
  
  dgamma(x=x, shape=shape_use, scale=scale_use)
  
}


pgamma_musig <- function(q, mu, sigma){
  
  shape_use <- 1/sigma^2
  
  scale_use <- mu*sigma^2
  
  pgamma(q=q, shape=shape_use, scale=scale_use)
  
}


qgamma_musig <- function(p, mu, sigma){
  
  shape_use <- 1/sigma^2
  
  scale_use <- mu*sigma^2
  
  qgamma(p=p, shape=shape_use, scale=scale_use)
  
}


rgamma_musig <- function(n, mu, sigma){
  
  shape_use <- 1/sigma^2
  
  scale_use <- mu*sigma^2
  
  rgamma(n=n, shape=shape_use, scale=scale_use)
  
}



zinfgammacop.loglik <- function(beta_mu1,
                       beta_mu2,
                       beta_sig1,
                       beta_sig2,
                       beta_rho,
                       p1,
                       p2,
                       xmat_mu1,
                       xmat_mu2,
                       xmat_sig1,
                       xmat_sig2,
                       xmat_rho,
                       y1,
                       y2
){
  
  neitherzero <- which((y1 > 0) & (y2 > 0))
  justy1zero <- which((y1 == 0) & (y2 > 0))
  justy2zero <- which((y1 > 0) & (y2 == 0))
  bothzero <- which((y1 == 0) & (y2 == 0))
  
  mu1 <- exp(xmat_mu1%*%beta_mu1)
  mu2 <- exp(xmat_mu2%*%beta_mu2)
  sig1 <- exp(xmat_sig1%*%beta_sig1)
  sig2 <- exp(xmat_sig2%*%beta_sig2)
  rho <- tanh(xmat_rho%*%beta_rho)
  
  q1 <- qnorm(pgamma_musig(q=y1, mu=mu1, sigma=sig1))
  q2 <- qnorm(pgamma_musig(q=y2, mu=mu2, sigma=sig2))
  
  dmarg1 <- dgamma_musig(x=y1, mu=mu1, sigma=sig1)
  dmarg2 <- dgamma_musig(x=y2, mu=mu2, sigma=sig2)
  
  frontpiece <- 1/(2*pi*sqrt(1-rho^2))
  mainpiece <- exp((-1/(2*(1-rho^2)))*(q1^2 + q2^2 - 2*rho*q1*q2))
  
  coppiece <- frontpiece*mainpiece*(dmarg1/dnorm(q1))*(dmarg2/dnorm(q2))
  
  resvec <- rep(0, length(y1))
  
  resvec[neitherzero] <- log((1-p1)*(1-p2)*coppiece)[neitherzero]
  resvec[justy2zero] <- log((1-p1)*p2*dmarg1)[justy2zero]
  resvec[justy1zero] <- log(p1*(1-p2)*dmarg2)[justy1zero]
  resvec[bothzero] <- log(p1*p2)[bothzero]
  
  sum(resvec)
  
}




zinfgammacop.gradhess <- function(beta_mu1,
                            beta_mu2,
                            beta_sig1,
                            beta_sig2,
                            beta_rho,
                            p1,
                            p2,
                            xmat_mu1,
                            xmat_mu2,
                            xmat_sig1,
                            xmat_sig2,
                            xmat_rho,
                            y1,
                            y2){
  
  mu1_idcs <- 1:length(beta_mu1)
  mu2_idcs <- (mu1_idcs[2]+1):(mu1_idcs[2] + length(beta_mu2))
  sig1_idcs <- (mu2_idcs[2]+1):(mu2_idcs[2] + length(beta_sig1))
  sig2_idcs <- (sig1_idcs[2]+1):(sig1_idcs[2] + length(beta_sig2))
  rho_idcs <- (sig2_idcs[2]+1):(sig2_idcs[2] + length(beta_rho))
  
  
  loglik_temp <- function(betavec){
    
    zinfgammacop.loglik(beta_mu1=betavec[mu1_idcs],
                  beta_mu2=betavec[mu2_idcs],
                  beta_sig1=betavec[sig1_idcs],
                  beta_sig2=betavec[sig2_idcs],
                  beta_rho=betavec[rho_idcs],
                  p1=p1,
                  p2=p2,
                  xmat_mu1=xmat_mu1,
                  xmat_mu2=xmat_mu2,
                  xmat_sig1=xmat_sig1,
                  xmat_sig2=xmat_sig2,
                  xmat_rho=xmat_rho,
                  y1=y1,
                  y2=y2)
    
  }
  
  retlist <- list()
  
  retlist[[1]] <- loglik_temp(betavec=c(beta_mu1, beta_mu2, beta_sig1, beta_sig2, beta_rho))
  
  retlist[[2]] <- pracma::grad(f=loglik_temp, x0=c(beta_mu1, beta_mu2, beta_sig1, beta_sig2, beta_rho))

  retlist[[3]] <- pracma::hessian(loglik_temp, x0=c(beta_mu1, beta_mu2, beta_sig1, beta_sig2, beta_rho))
  
  retlist
    
}




zinfgammacop.sim <- function(mu1, mu2, sig1, sig2, rho, p1, p2){
  
  paramlens <- c(length(mu1), length(mu2), length(sig1), length(sig2), length(rho), length(p1), length(p2))
  
  if (max(paramlens) != min(paramlens)){
    stop("one of mu1, mu2, sig1, sig2, rho, p1, p2 has a different length than the rest")
  } else {
    n <- max(paramlens)
  }
  
  oursamp <- matrix(0, nrow=n, ncol=2)
  
  for (i in 1:n){
    
    copsetup <- normalCopula(param = rho[i])
    
    speclist1 <- list( mu = mu1[i], sigma = sig1[i] )
    speclist2 <- list( mu = mu2[i], sigma = sig2[i] )
    
    
    spec <- mvdc(copula = copsetup, c("GA", "GA"), list(speclist1, speclist2) )
    
    tempsamp <- rMvdc(1, spec)*rbinom(n=2, size=1, prob=1-c(p1[i], p2[i]))
    
    oursamp[i, ] <- tempsamp
    
  }

  oursamp
  
}



zinfgammacop.override <- function (params, respvec, VC, ps, AT = FALSE) 
{
  p1 <- p2 <- pdf1 <- pdf2 <- c.copula.be2 <- c.copula.be1 <- c.copula2.be1be2 <- NA
  eta1 <- VC$X1 %*% params[1:VC$X1.d2]
  eta2 <- VC$X2 %*% params[(VC$X1.d2 + 1):(VC$X1.d2 + VC$X2.d2)]
  etad <- etas1 <- etas2 <- l.ln <- NULL
  if (is.null(VC$X3)) {
    sigma21.st <- etas1 <- params[(VC$X1.d2 + VC$X2.d2 + 
                                     1)]
    sigma22.st <- etas2 <- params[(VC$X1.d2 + VC$X2.d2 + 
                                     2)]
    teta.st <- etad <- params[(VC$X1.d2 + VC$X2.d2 + 3)]
  }
  if (!is.null(VC$X3)) {
    sigma21.st <- etas1 <- VC$X3 %*% params[(VC$X1.d2 + 
                                               VC$X2.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2)]
    sigma22.st <- etas2 <- VC$X4 %*% params[(VC$X1.d2 + 
                                               VC$X2.d2 + VC$X3.d2 + 1):(VC$X1.d2 + VC$X2.d2 + 
                                                                           VC$X3.d2 + VC$X4.d2)]
    teta.st <- etad <- VC$X5 %*% params[(VC$X1.d2 + VC$X2.d2 + 
                                           VC$X3.d2 + VC$X4.d2 + 1):(VC$X1.d2 + VC$X2.d2 + 
                                                                       VC$X3.d2 + VC$X4.d2 + VC$X5.d2)]
  }
  sstr1 <- esp.tr(sigma21.st, VC$margins[1])
  sstr2 <- esp.tr(sigma22.st, VC$margins[2])
  sigma21.st <- sstr1$vrb.st
  sigma22.st <- sstr2$vrb.st
  sigma21 <- sstr1$vrb
  sigma22 <- sstr2$vrb
  eta1 <- eta.tr(eta1, VC$margins[1])
  eta2 <- eta.tr(eta2, VC$margins[2])
  resT <- teta.tr(VC, teta.st)
  teta.st1 <- teta.st2 <- teta.st <- resT$teta.st
  teta1 <- teta2 <- teta <- resT$teta
  Cop1 <- Cop2 <- VC$BivD
  teta.ind1 <- as.logical(c(1, 0, round(runif(VC$n - 2))))
  teta.ind2 <- teta.ind1 == FALSE
  if (!(VC$BivD %in% VC$BivD2) && length(teta.st) > 1) {
    teta.st1 <- teta.st[teta.ind1]
    teta.st2 <- teta.st[teta.ind2]
    teta1 <- teta[teta.ind1]
    teta2 <- teta[teta.ind2]
  }
  if (VC$BivD %in% VC$BivD2) {
    if (VC$BivD %in% VC$BivD2[c(1:4, 13:16)]) 
      teta.ind1 <- ifelse(VC$my.env$signind * teta > exp(VC$zerov), 
                          TRUE, FALSE)
    if (VC$BivD %in% VC$BivD2[5:12]) 
      teta.ind1 <- ifelse(VC$my.env$signind * teta > exp(VC$zerov) + 
                            1, TRUE, FALSE)
    teta.ind2 <- teta.ind1 == FALSE
    VC$my.env$signind <- ifelse(teta.ind1 == TRUE, 1, -1)
    teta1 <- teta[teta.ind1]
    teta2 <- -teta[teta.ind2]
    teta.st1 <- teta.st[teta.ind1]
    teta.st2 <- teta.st[teta.ind2]
    if (length(teta) == 1) 
      teta.ind2 <- teta.ind1 <- rep(TRUE, VC$n)
    Cop1Cop2R <- Cop1Cop2(VC$BivD)
    Cop1 <- Cop1Cop2R$Cop1
    Cop2 <- Cop1Cop2R$Cop2
  }
  dHs1 <- distrHs(respvec$y1, eta1, sigma21, sigma21.st, nu = 1, 
                  nu.st = 1, margin2 = VC$margins[1], naive = FALSE, min.dn = VC$min.dn, 
                  min.pr = VC$min.pr, max.pr = VC$max.pr)
  dHs2 <- distrHs(respvec$y2, eta2, sigma22, sigma22.st, nu = 1, 
                  nu.st = 1, margin2 = VC$margins[2], naive = FALSE, min.dn = VC$min.dn, 
                  min.pr = VC$min.pr, max.pr = VC$max.pr)
  pdf1 <- dHs1$pdf2
  pdf2 <- dHs2$pdf2
  p1 <- dHs1$p2
  p2 <- dHs2$p2
  if (length(teta1) != 0) 
    dH1 <- copgHsAT(p1[teta.ind1], p2[teta.ind1], teta1, 
                    Cop1, Ln = TRUE, min.dn = VC$min.dn, min.pr = VC$min.pr, 
                    max.pr = VC$max.pr)
  if (length(teta2) != 0) 
    dH2 <- copgHsAT(p1[teta.ind2], p2[teta.ind2], teta2, 
                    Cop2, Ln = TRUE, min.dn = VC$min.dn, min.pr = VC$min.pr, 
                    max.pr = VC$max.pr)
  c.copula2.be1be2 <- NA
  if (length(teta1) != 0) 
    c.copula2.be1be2[teta.ind1] <- dH1$c.copula2.be1be2
  if (length(teta2) != 0) 
    c.copula2.be1be2[teta.ind2] <- dH2$c.copula2.be1be2
  l.par <- VC$weights * (log(pdf1) + log(pdf2) + log(c.copula2.be1be2))
  
  # replacing old loglik with zinf loglik

  neitherzero <- 1*((respvec$y1 > 0) & (respvec$y2 > 0))
  justy1zero <- 1*((respvec$y1 == 0) & (respvec$y2 > 0))
  justy2zero <- 1*((respvec$y1 > 0) & (respvec$y2 == 0))
  bothzero <- 1*((respvec$y1 == 0) & (respvec$y2 == 0))
  
  # l.par.new <- log(
  #   (1-pr1)*(1-pr2)*c.copula2.be1be2*pdf1*pdf2*neitherzero +
  #     (1-pr1)*pr2*pdf1*justy2zero +
  #     pr1*(1-pr2)*pdf2*justy1zero + 
  #     pr1*pr2*bothzero
  # )
  
  l.par.new <- log(
    c.copula2.be1be2*pdf1*pdf2*neitherzero +
      pdf1*justy2zero +
      pdf2*justy1zero + 
      bothzero
  )
  
  ourlik <- exp(l.par.new)
  theirlik <- exp(l.par)
  
  # end of modification

  derpdf1.dereta1 <- dHs1$derpdf2.dereta2
  derpdf1.dersigma21.st <- dHs1$derpdf2.dersigma2.st
  derpdf2.dereta2 <- dHs2$derpdf2.dereta2
  derpdf2.dersigma22.st <- dHs2$derpdf2.dersigma2.st
  derp1.dereta1 <- dHs1$derp2.dereta2
  derp1.dersigma21.st <- dHs1$derp2.dersigma.st
  derp2.dereta2 <- dHs2$derp2.dereta2
  derp2.dersigma22.st <- dHs2$derp2.dersigma.st
  if (length(teta1) != 0) 
    BITS1 <- copgHsCont(p1[teta.ind1], p2[teta.ind1], teta1, 
                        teta.st1, Cop1, Cont = TRUE)
  if (length(teta2) != 0) 
    BITS2 <- copgHsCont(p1[teta.ind2], p2[teta.ind2], teta2, 
                        teta.st2, Cop2, Cont = TRUE)
  der2h.derp1p1 <- NA
  if (length(teta1) != 0) 
    der2h.derp1p1[teta.ind1] <- BITS1$der2h.derp1p1
  if (length(teta2) != 0) 
    der2h.derp1p1[teta.ind2] <- BITS2$der2h.derp1p1
  derc.dereta1 <- der2h.derp1p1 * derp1.dereta1
  derc.dersigma21.st <- der2h.derp1p1 * derp1.dersigma21.st
  der2h.derp1p2 <- NA
  if (length(teta1) != 0) 
    der2h.derp1p2[teta.ind1] <- BITS1$der2h.derp1p2
  if (length(teta2) != 0) 
    der2h.derp1p2[teta.ind2] <- BITS2$der2h.derp1p2
  derc.dereta2 <- der2h.derp1p2 * derp2.dereta2
  derc.dersigma22.st <- der2h.derp1p2 * derp2.dersigma22.st
  der2h.derp1teta <- NA
  derteta.derteta.st <- NA
  if (length(teta1) != 0) 
    der2h.derp1teta[teta.ind1] <- BITS1$der2h.derp1teta
  if (length(teta2) != 0) 
    der2h.derp1teta[teta.ind2] <- BITS2$der2h.derp1teta
  if (length(teta1) != 0) 
    derteta.derteta.st[teta.ind1] <- BITS1$derteta.derteta.st
  if (length(teta2) != 0) 
    derteta.derteta.st[teta.ind2] <- BITS2$derteta.derteta.st
  der2h.derp1teta.st <- der2h.derp1teta * derteta.derteta.st
  
  # here we have dl/d(XB) where l is log(c(F1(y1),F2(y2))*f1(y1)f2(y2))

  dl.dbe1 <- VC$weights * (derpdf1.dereta1/pdf1 + derc.dereta1/c.copula2.be1be2)
  dl.dbe2 <- VC$weights * (derpdf2.dereta2/pdf2 + derc.dereta2/c.copula2.be1be2)
  dl.dsigma21.st <- VC$weights * (derpdf1.dersigma21.st/pdf1 + 
                                    derc.dersigma21.st/c.copula2.be1be2)
  dl.dsigma22.st <- VC$weights * (derpdf2.dersigma22.st/pdf2 + 
                                    derc.dersigma22.st/c.copula2.be1be2)
  dl.dteta.st <- VC$weights * (der2h.derp1teta.st/c.copula2.be1be2)
  
  # seeing as our loglik is 
  # log((1-p1)(1-p2)c(F1(y1),F2(y2))*f1(y1)f2(y2)*(y1 > 0 & y2 > 0) 
  #     + (1-p1)p2*f1(y1)*(y1 > 0 & y2 = 0)
  #     + p1(1-p2)*f2(y2)*(y1 = 0 & y2 > 0)
  #     + p1p2(y1 = 0 & y2 = 0)
  # we need to find ways to get d/d(BX) c(F1(y1),F2(y2))*f1(y1)f2(y2)
  # and d/d(BX) f1(y1) etc. from existing components
  # well dl/d(XB) = (1/L) * dL/d(XB) so dL/d(XB) = L * dl/d(XB), so
  # d/d(XB) c(F1(y1),F2(y2))*f1(y1)f2(y2) = dl.dbe1 * pdf1*pdf2*c.copula2.be1be2
  # d/d(XB) f1(y1) = derpdf1.dereta1
  
  # dl.dbe1.new <- VC$weights * (1/ourlik)*(
  #   (1-pr1)*(1-pr2)*theirlik*dl.dbe1*neitherzero +
  #     (1-pr1)*pr2*(derpdf1.dereta1)*justy2zero +
  #     0 + 
  #     0)
  
  dl.dbe1.new <- VC$weights * (1/ourlik)*(
    theirlik*dl.dbe1*neitherzero +
    (derpdf1.dereta1)*justy2zero +
      0 +
      0)
  
  # dl.dbe2.new <- VC$weights * (1/ourlik)*(
  #   (1-pr1)*(1-pr2)*theirlik*dl.dbe2*neitherzero +
  #     0 +
  #     pr1*(1-pr2)*(derpdf2.dereta2)*justy1zero +
  #     0)
  
  dl.dbe2.new <- VC$weights * (1/ourlik)*(
    theirlik*dl.dbe2*neitherzero +
      0 +
    (derpdf2.dereta2)*justy1zero +
      0)
  
  # dl.dsigma21.st.new <- VC$weights * (1/ourlik)*(
  #   (1-pr1)*(1-pr2)*theirlik*dl.dsigma21.st*neitherzero +
  #     (1-pr1)*pr2*(derpdf1.dersigma21.st)*justy2zero +
  #     0 +
  #     0)
  
  dl.dsigma21.st.new <- VC$weights * (1/ourlik)*(
    theirlik*dl.dsigma21.st*neitherzero +
    (derpdf1.dersigma21.st)*justy2zero +
      0 +
      0)
  
  # dl.dsigma22.st.new <- VC$weights * (1/ourlik)*(
  #   (1-pr1)*(1-pr2)*theirlik*dl.dsigma22.st*neitherzero +
  #     0 +
  #     pr1*(1-pr2)*(derpdf2.dersigma22.st)*justy1zero +
  #     0)
  
  dl.dsigma22.st.new <- VC$weights * (1/ourlik)*(
    theirlik*dl.dsigma22.st*neitherzero +
      0 +
    (derpdf2.dersigma22.st)*justy1zero +
      0)
  
  # dl.dteta.st.new <- VC$weights * (1/ourlik)*(
  #   (1-pr1)*(1-pr2)*theirlik*dl.dteta.st*neitherzero +
  #     0 +
  #     0 +
  #     0)
  
  dl.dteta.st.new <- VC$weights * (1/ourlik)*(
    theirlik*dl.dteta.st*neitherzero +
      0 +
      0 +
      0)
  
  # end of modification 
  
  # check if this matches up with numerical derivatives:
  
  # whatres <- zinfgammacop.gradhess(beta_mu1=params[1:VC$X1.d2],
  #                            beta_mu2=params[(VC$X1.d2 + 1):(VC$X1.d2 + VC$X2.d2)],
  #                            beta_sig1=params[(VC$X1.d2 + 
  #                                                VC$X2.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2)],
  #                            beta_sig2=params[(VC$X1.d2 + 
  #                                                VC$X2.d2 + VC$X3.d2 + 1):(VC$X1.d2 + VC$X2.d2 + 
  #                                                                            VC$X3.d2 + VC$X4.d2)],
  #                            beta_rho=params[(VC$X1.d2 + VC$X2.d2 + 
  #                                               VC$X3.d2 + VC$X4.d2 + 1):(VC$X1.d2 + VC$X2.d2 + 
  #                                                                           VC$X3.d2 + VC$X4.d2 + VC$X5.d2)],
  #                            p1=pr1,
  #                            p2=pr2,
  #                            xmat_mu1=VC$X1,
  #                            xmat_mu2=VC$X2,
  #                            xmat_sig1=VC$X3,
  #                            xmat_sig2=VC$X4,
  #                            xmat_rho=VC$X5,
  #                            y1=respvec$y1,
  #                            y2=respvec$y2)
  
  
  der2c.derrho.derrho <- NA
  der2c.derp1.derp1 <- NA
  der2c.derp2.derp2 <- NA
  der2c.derp1.derp2 <- NA
  der2c.derp1.derrho <- NA
  der2c.derp2.derrho <- NA
  der2teta.derteta.stteta.st <- NA
  if (length(teta1) != 0) {
    der2c.derrho.derrho[teta.ind1] <- BITS1$der2c.derrho.derrho
    der2c.derp1.derp1[teta.ind1] <- BITS1$der2c.derp1.derp1
    der2c.derp2.derp2[teta.ind1] <- BITS1$der2c.derp2.derp2
    der2c.derp1.derp2[teta.ind1] <- BITS1$der2c.derp1.derp2
    der2c.derp1.derrho[teta.ind1] <- BITS1$der2c.derp1.derrho
    der2c.derp2.derrho[teta.ind1] <- BITS1$der2c.derp2.derrho
  }
  if (length(teta2) != 0) {
    der2c.derrho.derrho[teta.ind2] <- BITS2$der2c.derrho.derrho
    der2c.derp1.derp1[teta.ind2] <- BITS2$der2c.derp1.derp1
    der2c.derp2.derp2[teta.ind2] <- BITS2$der2c.derp2.derp2
    der2c.derp1.derp2[teta.ind2] <- BITS2$der2c.derp1.derp2
    der2c.derp1.derrho[teta.ind2] <- BITS2$der2c.derp1.derrho
    der2c.derp2.derrho[teta.ind2] <- BITS2$der2c.derp2.derrho
  }
  if (length(teta1) != 0) 
    der2teta.derteta.stteta.st[teta.ind1] <- BITS1$der2teta.derteta.stteta.st
  if (length(teta2) != 0) 
    der2teta.derteta.stteta.st[teta.ind2] <- BITS2$der2teta.derteta.stteta.st
  der2pdf1.dereta1 <- dHs1$der2pdf2.dereta2
  der2pdf2.dereta2 <- dHs2$der2pdf2.dereta2
  der2pdf1.dersigma21.st2 <- dHs1$der2pdf2.dersigma2.st2
  der2pdf2.dersigma22.st2 <- dHs2$der2pdf2.dersigma2.st2
  der2p1.dereta1eta1 <- dHs1$der2p2.dereta2eta2
  der2p2.dereta2eta2 <- dHs2$der2p2.dereta2eta2
  der2p1.dersigma21.st2 <- dHs1$der2p2.dersigma2.st2
  der2p2.dersigma22.st2 <- dHs2$der2p2.dersigma2.st2
  der2pdf1.dereta1dersigma21.st <- dHs1$der2pdf2.dereta2dersigma2.st
  der2pdf2.dereta2dersigma22.st <- dHs2$der2pdf2.dereta2dersigma2.st
  der2p1.dereta1dersigma21.st <- dHs1$der2p2.dereta2dersigma2.st
  der2p2.dereta2dersigma22.st <- dHs2$der2p2.dereta2dersigma2.st
  d2l.be1.be1 <- -VC$weights * ((der2pdf1.dereta1 * pdf1 - 
                                   derpdf1.dereta1^2)/pdf1^2 + ((der2c.derp1.derp1 * derp1.dereta1^2 + 
                                                                   der2h.derp1p1 * der2p1.dereta1eta1) * c.copula2.be1be2 - 
                                                                  derc.dereta1^2)/c.copula2.be1be2^2)
  d2l.be2.be2 <- -VC$weights * ((der2pdf2.dereta2 * pdf2 - 
                                   derpdf2.dereta2^2)/pdf2^2 + ((der2c.derp2.derp2 * derp2.dereta2^2 + 
                                                                   der2h.derp1p2 * der2p2.dereta2eta2) * c.copula2.be1be2 - 
                                                                  derc.dereta2^2)/c.copula2.be1be2^2)
  d2l.rho.rho <- -VC$weights * (((der2c.derrho.derrho * derteta.derteta.st^2 + 
                                    der2h.derp1teta * der2teta.derteta.stteta.st) * c.copula2.be1be2 - 
                                   der2h.derp1teta.st^2)/c.copula2.be1be2^2)
  d2l.sigma21.sigma21 <- -VC$weights * ((der2pdf1.dersigma21.st2 * 
                                           pdf1 - derpdf1.dersigma21.st^2)/pdf1^2 + ((der2c.derp1.derp1 * 
                                                                                        derp1.dersigma21.st^2 + der2h.derp1p1 * der2p1.dersigma21.st2) * 
                                                                                       c.copula2.be1be2 - derc.dersigma21.st^2)/c.copula2.be1be2^2)
  d2l.sigma22.sigma22 <- -VC$weights * ((der2pdf2.dersigma22.st2 * 
                                           pdf2 - derpdf2.dersigma22.st^2)/pdf2^2 + ((der2c.derp2.derp2 * 
                                                                                        derp2.dersigma22.st^2 + der2h.derp1p2 * der2p2.dersigma22.st2) * 
                                                                                       c.copula2.be1be2 - derc.dersigma22.st^2)/c.copula2.be1be2^2)
  d2l.be1.be2 <- -VC$weights * ((der2c.derp1.derp2 * derp1.dereta1 * 
                                   derp2.dereta2 * c.copula2.be1be2 - derc.dereta1 * derc.dereta2)/c.copula2.be1be2^2)
  d2l.be1.rho <- -VC$weights * ((der2c.derp1.derrho * derp1.dereta1 * 
                                   derteta.derteta.st * c.copula2.be1be2 - derc.dereta1 * 
                                   der2h.derp1teta * derteta.derteta.st)/c.copula2.be1be2^2)
  d2l.be2.rho <- -VC$weights * ((der2c.derp2.derrho * derp2.dereta2 * 
                                   derteta.derteta.st * c.copula2.be1be2 - derc.dereta2 * 
                                   der2h.derp1teta * derteta.derteta.st)/c.copula2.be1be2^2)
  d2l.be1.sigma21 <- -VC$weights * ((der2pdf1.dereta1dersigma21.st * 
                                       pdf1 - derpdf1.dereta1 * derpdf1.dersigma21.st)/pdf1^2 + 
                                      ((der2c.derp1.derp1 * derp1.dereta1 * derp1.dersigma21.st + 
                                          der2h.derp1p1 * der2p1.dereta1dersigma21.st) * c.copula2.be1be2 - 
                                         derc.dereta1 * derc.dersigma21.st)/c.copula2.be1be2^2)
  d2l.be2.sigma22 <- -VC$weights * ((der2pdf2.dereta2dersigma22.st * 
                                       pdf2 - derpdf2.dereta2 * derpdf2.dersigma22.st)/pdf2^2 + 
                                      ((der2c.derp2.derp2 * derp2.dereta2 * derp2.dersigma22.st + 
                                          der2h.derp1p2 * der2p2.dereta2dersigma22.st) * c.copula2.be1be2 - 
                                         derc.dereta2 * derc.dersigma22.st)/c.copula2.be1be2^2)
  d2l.be2.sigma21 <- -VC$weights * ((der2c.derp1.derp2 * derp2.dereta2 * 
                                       derp1.dersigma21.st * c.copula2.be1be2 - derc.dereta2 * 
                                       derc.dersigma21.st)/c.copula2.be1be2^2)
  d2l.be1.sigma22 <- -VC$weights * ((der2c.derp1.derp2 * derp1.dereta1 * 
                                       derp2.dersigma22.st * c.copula2.be1be2 - derc.dereta1 * 
                                       derc.dersigma22.st)/c.copula2.be1be2^2)
  d2l.rho.sigma21 <- -VC$weights * ((der2c.derp1.derrho * 
                                       derp1.dersigma21.st * derteta.derteta.st * c.copula2.be1be2 - 
                                       derc.dersigma21.st * der2h.derp1teta * derteta.derteta.st)/c.copula2.be1be2^2)
  d2l.rho.sigma22 <- -VC$weights * ((der2c.derp2.derrho * 
                                       derp2.dersigma22.st * derteta.derteta.st * c.copula2.be1be2 - 
                                       derc.dersigma22.st * der2h.derp1teta * derteta.derteta.st)/c.copula2.be1be2^2)
  d2l.sigma21.sigma22 <- -VC$weights * ((der2c.derp1.derp2 * 
                                           derp1.dersigma21.st * derp2.dersigma22.st * c.copula2.be1be2 - 
                                           derc.dersigma21.st * derc.dersigma22.st)/c.copula2.be1be2^2)
  
  
  # now we will calculate the second derivatives for OUR zinf log likelihood. The negative sign on the d2l components is to reverse the -VC$weights factor
  
  # d2l.be1.be1.new <- -VC$weights * (-(dl.dbe1.new^2) + (1/ourlik)*((1-pr1)*(1-pr2)*(theirlik)*(-d2l.be1.be1 + dl.dbe1^2)*neitherzero + (1-pr1)*pr2*der2pdf1.dereta1*justy2zero))
  # 
  # d2l.be2.be2.new <- -VC$weights * (-(dl.dbe2.new^2) + (1/ourlik)*((1-pr1)*(1-pr2)*(theirlik)*(-d2l.be2.be2 + dl.dbe2^2)*neitherzero + pr1*(1-pr2)*der2pdf2.dereta2*justy1zero))
  # 
  # d2l.rho.rho.new <- -VC$weights * (-(dl.dteta.st.new^2) + (1/ourlik)*((1-pr1)*(1-pr2)*(theirlik)*(-d2l.rho.rho + dl.dteta.st^2)*neitherzero))
  # 
  # d2l.sigma21.sigma21.new <- -VC$weights * (-(dl.dsigma21.st.new^2) + (1/ourlik)*((1-pr1)*(1-pr2)*(theirlik)*(-d2l.sigma21.sigma21 + dl.dsigma21.st^2)*neitherzero + (1-pr1)*pr2*der2pdf1.dersigma21.st2*justy2zero))
  # 
  # d2l.sigma22.sigma22.new <- -VC$weights * (-(dl.dsigma22.st.new^2) + (1/ourlik)*((1-pr1)*(1-pr2)*(theirlik)*(-d2l.sigma22.sigma22 + dl.dsigma22.st^2)*neitherzero + pr1*(1-pr2)*der2pdf2.dersigma22.st2*justy1zero))
  # 
  # d2l.be1.be2.new <- -VC$weights * (-(dl.dbe1.new*dl.dbe2.new) + (1/ourlik)*((1-pr1)*(1-pr2)*(theirlik)*(-d2l.be1.be2 + dl.dbe1*dl.dbe2)*neitherzero))
  # 
  # d2l.be1.rho.new <- -VC$weights * (-(dl.dbe1.new*dl.dteta.st.new) + (1/ourlik)*((1-pr1)*(1-pr2)*(theirlik)*(-d2l.be1.rho + dl.dbe1*dl.dteta.st)*neitherzero))
  # 
  # d2l.be2.rho.new <- -VC$weights * (-(dl.dbe2.new*dl.dteta.st.new) + (1/ourlik)*((1-pr1)*(1-pr2)*(theirlik)*(-d2l.be2.rho + dl.dbe2*dl.dteta.st)*neitherzero))
  # 
  # d2l.be1.sigma21.new <- -VC$weights * (-(dl.dbe1.new*dl.dsigma21.st.new) + (1/ourlik)*((1-pr1)*(1-pr2)*(theirlik)*(-d2l.be1.sigma21 + dl.dbe1*dl.dsigma21.st)*neitherzero + (1-pr1)*pr2*der2pdf1.dereta1dersigma21.st*justy2zero))
  # 
  # d2l.be2.sigma22.new <- -VC$weights * (-(dl.dbe2.new*dl.dsigma22.st.new) + (1/ourlik)*((1-pr1)*(1-pr2)*(theirlik)*(-d2l.be2.sigma22 + dl.dbe2*dl.dsigma22.st)*neitherzero + pr1*(1-pr2)*der2pdf2.dereta2dersigma22.st*justy1zero))
  # 
  # d2l.be2.sigma21.new <- -VC$weights * (-(dl.dbe2.new*dl.dsigma21.st.new) + (1/ourlik)*((1-pr1)*(1-pr2)*(theirlik)*(-d2l.be2.sigma21 + dl.dbe2*dl.dsigma21.st)*neitherzero))
  # 
  # d2l.be1.sigma22.new <- -VC$weights * (-(dl.dbe1.new*dl.dsigma22.st.new) + (1/ourlik)*((1-pr1)*(1-pr2)*(theirlik)*(-d2l.be1.sigma22 + dl.dbe1*dl.dsigma22.st)*neitherzero))
  # 
  # d2l.rho.sigma21.new <- -VC$weights * (-(dl.dsigma21.st.new*dl.dteta.st.new) + (1/ourlik)*((1-pr1)*(1-pr2)*(theirlik)*(-d2l.rho.sigma21 + dl.dsigma21.st*dl.dteta.st)*neitherzero))
  # 
  # d2l.rho.sigma22.new <- -VC$weights * (-(dl.dsigma22.st.new*dl.dteta.st.new) + (1/ourlik)*((1-pr1)*(1-pr2)*(theirlik)*(-d2l.rho.sigma22 + dl.dsigma22.st*dl.dteta.st)*neitherzero))
  # 
  # d2l.sigma21.sigma22.new <- -VC$weights * (-(dl.dsigma21.st.new*dl.dsigma22.st.new) + (1/ourlik)*((1-pr1)*(1-pr2)*(theirlik)*(-d2l.sigma21.sigma22 + dl.dsigma21.st*dl.dsigma22.st)*neitherzero))

  
  d2l.be1.be1.new <- -VC$weights * (-(dl.dbe1.new^2) + (1/ourlik)*((theirlik)*(-d2l.be1.be1 + dl.dbe1^2)*neitherzero + der2pdf1.dereta1*justy2zero))
  
  d2l.be2.be2.new <- -VC$weights * (-(dl.dbe2.new^2) + (1/ourlik)*((theirlik)*(-d2l.be2.be2 + dl.dbe2^2)*neitherzero + der2pdf2.dereta2*justy1zero))
  
  d2l.rho.rho.new <- -VC$weights * (-(dl.dteta.st.new^2) + (1/ourlik)*((theirlik)*(-d2l.rho.rho + dl.dteta.st^2)*neitherzero))
  
  d2l.sigma21.sigma21.new <- -VC$weights * (-(dl.dsigma21.st.new^2) + (1/ourlik)*((theirlik)*(-d2l.sigma21.sigma21 + dl.dsigma21.st^2)*neitherzero + der2pdf1.dersigma21.st2*justy2zero))
  
  d2l.sigma22.sigma22.new <- -VC$weights * (-(dl.dsigma22.st.new^2) + (1/ourlik)*((theirlik)*(-d2l.sigma22.sigma22 + dl.dsigma22.st^2)*neitherzero + der2pdf2.dersigma22.st2*justy1zero))
  
  d2l.be1.be2.new <- -VC$weights * (-(dl.dbe1.new*dl.dbe2.new) + (1/ourlik)*((theirlik)*(-d2l.be1.be2 + dl.dbe1*dl.dbe2)*neitherzero))
  
  d2l.be1.rho.new <- -VC$weights * (-(dl.dbe1.new*dl.dteta.st.new) + (1/ourlik)*((theirlik)*(-d2l.be1.rho + dl.dbe1*dl.dteta.st)*neitherzero))
  
  d2l.be2.rho.new <- -VC$weights * (-(dl.dbe2.new*dl.dteta.st.new) + (1/ourlik)*((theirlik)*(-d2l.be2.rho + dl.dbe2*dl.dteta.st)*neitherzero))
  
  d2l.be1.sigma21.new <- -VC$weights * (-(dl.dbe1.new*dl.dsigma21.st.new) + (1/ourlik)*((theirlik)*(-d2l.be1.sigma21 + dl.dbe1*dl.dsigma21.st)*neitherzero + der2pdf1.dereta1dersigma21.st*justy2zero))
  
  d2l.be2.sigma22.new <- -VC$weights * (-(dl.dbe2.new*dl.dsigma22.st.new) + (1/ourlik)*((theirlik)*(-d2l.be2.sigma22 + dl.dbe2*dl.dsigma22.st)*neitherzero + der2pdf2.dereta2dersigma22.st*justy1zero))
  
  d2l.be2.sigma21.new <- -VC$weights * (-(dl.dbe2.new*dl.dsigma21.st.new) + (1/ourlik)*((theirlik)*(-d2l.be2.sigma21 + dl.dbe2*dl.dsigma21.st)*neitherzero))
  
  d2l.be1.sigma22.new <- -VC$weights * (-(dl.dbe1.new*dl.dsigma22.st.new) + (1/ourlik)*((theirlik)*(-d2l.be1.sigma22 + dl.dbe1*dl.dsigma22.st)*neitherzero))
  
  d2l.rho.sigma21.new <- -VC$weights * (-(dl.dsigma21.st.new*dl.dteta.st.new) + (1/ourlik)*((theirlik)*(-d2l.rho.sigma21 + dl.dsigma21.st*dl.dteta.st)*neitherzero))
  
  d2l.rho.sigma22.new <- -VC$weights * (-(dl.dsigma22.st.new*dl.dteta.st.new) + (1/ourlik)*((theirlik)*(-d2l.rho.sigma22 + dl.dsigma22.st*dl.dteta.st)*neitherzero))
  
  d2l.sigma21.sigma22.new <- -VC$weights * (-(dl.dsigma21.st.new*dl.dsigma22.st.new) + (1/ourlik)*((theirlik)*(-d2l.sigma21.sigma22 + dl.dsigma21.st*dl.dsigma22.st)*neitherzero))
  
    
  
  l.par <- l.par.new
  
  dl.dbe1 <- dl.dbe1.new
  dl.dbe2 <- dl.dbe2.new
  dl.dsigma21.st <- dl.dsigma21.st.new
  dl.dsigma22.st <- dl.dsigma22.st.new
  dl.dteta.st <- dl.dteta.st.new
  
  d2l.be1.be1 <- d2l.be1.be1.new
  d2l.be2.be2 <- d2l.be2.be2.new
  d2l.rho.rho <- d2l.rho.rho.new
  d2l.sigma21.sigma21 <- d2l.sigma21.sigma21.new
  d2l.sigma22.sigma22 <- d2l.sigma22.sigma22.new
  d2l.be1.be2 <- d2l.be1.be2.new
  d2l.be1.rho <- d2l.be1.rho.new
  d2l.be2.rho <- d2l.be2.rho.new
  d2l.be1.sigma21 <- d2l.be1.sigma21.new
  d2l.be2.sigma22 <- d2l.be2.sigma22.new
  d2l.be2.sigma21 <- d2l.be2.sigma21.new
  d2l.be1.sigma22 <- d2l.be1.sigma22.new
  d2l.rho.sigma21 <- d2l.rho.sigma21.new
  d2l.rho.sigma22 <- d2l.rho.sigma22.new
  d2l.sigma21.sigma22 <- d2l.sigma21.sigma22.new

  # end of modification
  
  if (is.null(VC$X3)) {
    G <- -c(colSums(c(dl.dbe1) * VC$X1), colSums(c(dl.dbe2) * 
                                                   VC$X2), sum(dl.dsigma21.st), sum(dl.dsigma22.st), 
            sum(dl.dteta.st))
    be1.be1 <- crossprod(VC$X1 * c(d2l.be1.be1), VC$X1)
    be2.be2 <- crossprod(VC$X2 * c(d2l.be2.be2), VC$X2)
    be1.be2 <- crossprod(VC$X1 * c(d2l.be1.be2), VC$X2)
    be1.rho <- t(t(rowSums(t(VC$X1 * c(d2l.be1.rho)))))
    be1.sigma21 <- t(t(rowSums(t(VC$X1 * c(d2l.be1.sigma21)))))
    be1.sigma22 <- t(t(rowSums(t(VC$X1 * c(d2l.be1.sigma22)))))
    be2.rho <- t(t(rowSums(t(VC$X2 * c(d2l.be2.rho)))))
    be2.sigma21 <- t(t(rowSums(t(VC$X2 * c(d2l.be2.sigma21)))))
    be2.sigma22 <- t(t(rowSums(t(VC$X2 * c(d2l.be2.sigma22)))))
    H <- rbind(cbind(be1.be1, be1.be2, be1.sigma21, be1.sigma22, 
                     be1.rho), cbind(t(be1.be2), be2.be2, be2.sigma21, 
                                     be2.sigma22, be2.rho), cbind(t(be1.sigma21), t(be2.sigma21), 
                                                                  sum(d2l.sigma21.sigma21), sum(d2l.sigma21.sigma22), 
                                                                  sum(d2l.rho.sigma21)), cbind(t(be1.sigma22), t(be2.sigma22), 
                                                                                               sum(d2l.sigma21.sigma22), sum(d2l.sigma22.sigma22), 
                                                                                               sum(d2l.rho.sigma22)), cbind(t(be1.rho), t(be2.rho), 
                                                                                                                            sum(d2l.rho.sigma21), sum(d2l.rho.sigma22), sum(d2l.rho.rho)))
  }
  if (!is.null(VC$X3)) {
    G <- -c(colSums(c(dl.dbe1) * VC$X1), colSums(c(dl.dbe2) * 
                                                   VC$X2), colSums(c(dl.dsigma21.st) * VC$X3), colSums(c(dl.dsigma22.st) * 
                                                                                                         VC$X4), colSums(c(dl.dteta.st) * VC$X5))
    be1.be1 <- crossprod(VC$X1 * c(d2l.be1.be1), VC$X1)
    be2.be2 <- crossprod(VC$X2 * c(d2l.be2.be2), VC$X2)
    be1.be2 <- crossprod(VC$X1 * c(d2l.be1.be2), VC$X2)
    be1.rho <- crossprod(VC$X1 * c(d2l.be1.rho), VC$X5)
    be2.rho <- crossprod(VC$X2 * c(d2l.be2.rho), VC$X5)
    be1.sigma21 <- crossprod(VC$X1 * c(d2l.be1.sigma21), 
                             VC$X3)
    be1.sigma22 <- crossprod(VC$X1 * c(d2l.be1.sigma22), 
                             VC$X4)
    be2.sigma21 <- crossprod(VC$X2 * c(d2l.be2.sigma21), 
                             VC$X3)
    be2.sigma22 <- crossprod(VC$X2 * c(d2l.be2.sigma22), 
                             VC$X4)
    sigma21.sigma21 <- crossprod(VC$X3 * c(d2l.sigma21.sigma21), 
                                 VC$X3)
    sigma21.sigma22 <- crossprod(VC$X3 * c(d2l.sigma21.sigma22), 
                                 VC$X4)
    rho.sigma21 <- crossprod(VC$X3 * c(d2l.rho.sigma21), 
                             VC$X5)
    sigma22.sigma22 <- crossprod(VC$X4 * c(d2l.sigma22.sigma22), 
                                 VC$X4)
    rho.sigma22 <- crossprod(VC$X4 * c(d2l.rho.sigma22), 
                             VC$X5)
    rho.rho <- crossprod(VC$X5 * c(d2l.rho.rho), VC$X5)
    H <- rbind(cbind(be1.be1, be1.be2, be1.sigma21, be1.sigma22, 
                     be1.rho), cbind(t(be1.be2), be2.be2, be2.sigma21, 
                                     be2.sigma22, be2.rho), cbind(t(be1.sigma21), t(be2.sigma21), 
                                                                  sigma21.sigma21, sigma21.sigma22, rho.sigma21), 
               cbind(t(be1.sigma22), t(be2.sigma22), t(sigma21.sigma22), 
                     sigma22.sigma22, rho.sigma22), cbind(t(be1.rho), 
                                                          t(be2.rho), t(rho.sigma21), t(rho.sigma22), 
                                                          rho.rho))
  }
  res <- -sum(l.par)
  if (VC$extra.regI == "pC") 
    H <- regH(H, type = 1)
  S.h <- ps$S.h
  if (length(S.h) != 1) {
    S.h1 <- 0.5 * crossprod(params, S.h) %*% params
    S.h2 <- S.h %*% params
  }
  else S.h <- S.h1 <- S.h2 <- 0
  S.res <- res
  res <- S.res + S.h1
  G <- G + S.h2
  H <- H + S.h
  if (VC$extra.regI == "sED") 
    H <- regH(H, type = 2)
  if (VC$margins[1] == "LN" || VC$margins[2] == "LN") {
    if (VC$margins[1] == "LN") 
      dHs1 <- distrHsAT(exp(respvec$y1), eta1, sigma21, 
                        1, margin2 = VC$margins[1], min.dn = VC$min.dn, 
                        min.pr = VC$min.pr, max.pr = VC$max.pr)
    if (VC$margins[2] == "LN") 
      dHs2 <- distrHsAT(exp(respvec$y2), eta2, sigma22, 
                        1, margin2 = VC$margins[2], min.dn = VC$min.dn, 
                        min.pr = VC$min.pr, max.pr = VC$max.pr)
    pdf1 <- dHs1$pdf2
    pdf2 <- dHs2$pdf2
    p1 <- dHs1$p2
    p2 <- dHs2$p2
    if (length(teta1) != 0) 
      dH1 <- copgHsAT(p1[teta.ind1], p2[teta.ind1], teta1, 
                      Cop1, Ln = TRUE, min.dn = VC$min.dn, min.pr = VC$min.pr, 
                      max.pr = VC$max.pr)
    if (length(teta2) != 0) 
      dH2 <- copgHsAT(p1[teta.ind2], p2[teta.ind2], teta2, 
                      Cop2, Ln = TRUE, min.dn = VC$min.dn, min.pr = VC$min.pr, 
                      max.pr = VC$max.pr)
    c.copula2.be1be2 <- NA
    if (length(teta1) != 0) 
      c.copula2.be1be2[teta.ind1] <- dH1$c.copula2.be1be2
    if (length(teta2) != 0) 
      c.copula2.be1be2[teta.ind2] <- dH2$c.copula2.be1be2
    l.ln <- -sum(VC$weights * (log(pdf1) + log(pdf2) + log(c.copula2.be1be2)))
  }
  list(value = res, gradient = G, hessian = H, S.h = S.h, 
       S.h1 = S.h1, S.h2 = S.h2, l = S.res, l.ln = l.ln, l.par = l.par, 
       ps = ps, eta1 = eta1, eta2 = eta2, etad = etad, etas1 = etas1, 
       etas2 = etas2, BivD = VC$BivD, p1 = p1, p2 = p2, pdf1 = pdf1, 
       pdf2 = pdf2, c.copula.be2 = c.copula.be2, c.copula.be1 = c.copula.be1, 
       c.copula2.be1be2 = c.copula2.be1be2, dl.dbe1 = dl.dbe1, 
       dl.dbe2 = dl.dbe2, dl.dsigma21.st = dl.dsigma21.st, 
       dl.dsigma22.st = dl.dsigma22.st, dl.dteta.st = dl.dteta.st, 
       teta.ind2 = teta.ind2, teta.ind1 = teta.ind1, Cop1 = Cop1, 
       Cop2 = Cop2, teta1 = teta1, teta2 = teta2)
}












zinfgammacop <- 
function (formula, data = list(), weights = NULL, subset = NULL, 
          copula = "N", copula2 = "N", margins, model, dof = 3, dof2 = 3, 
          ordinal = FALSE, surv = FALSE, cens1 = NULL, cens2 = NULL, 
          cens3 = NULL, dep.cens = FALSE, upperBt1 = NULL, upperBt2 = NULL, 
          gamlssfit = FALSE, fp = FALSE, infl.fac = 1, rinit = 1, 
          rmax = 100, iterlimsp = 50, tolsp = 1e-07, gc.l = FALSE, 
          parscale, extra.regI = "t", k1.tvc = 0, k2.tvc = 0, knots = NULL, 
          penCor = "unpen", sp.penCor = 3, Chol = FALSE, gamma = 1, 
          w.alasso = NULL, drop.unused.levels = TRUE, min.dn = 1e-40, 
          min.pr = 1e-16, max.pr = 0.999999) 
{
  
  BivD <- copula
  BivD2 <- copula2
  Model <- model
  if (dep.cens == TRUE) 
    stop("The dependent censoring case is work in progress. \nGet in touch should you wish to know more.")
  if (model %in% c("TSS", "TESS")) 
    stop("This trivariate model case has not been made available yet. \nGet in touch should you wish to know more.")
  if (model == "BSS" && margins[2] == "TW") 
    stop("This model case has not been made available yet. \nGet in touch should you wish to know more.")
  if (surv == TRUE && !(margins[1] %in% c("PH", "PO", "probit"))) 
    stop("This model case has not been made available yet. \nGet in touch should you wish to know more.")
  if (margins[1] %in% c("PO", "ZTP", "NBI", "NBII", "PIG", 
                        "DGP", "DGPII", "DGP0") && margins[2] %in% c("N", "TW", 
                                                                     "LN", "GU", "rGU", "LO", "WEI", "iG", "GA", "DAGUM", 
                                                                     "SM", "BE", "FISK", "GP", "GPII", "GPo")) 
    stop("This model case has not been made available yet. \nGet in touch should you wish to know more.")
  if (missing(margins)) 
    stop("You must choose the margins' values.")
  if (missing(model)) 
    stop("You must choose a model type.")
  if (margins[1] == "PH" && surv == TRUE) 
    margins[1] <- "cloglog"
  if (margins[1] == "PO" && surv == TRUE) 
    margins[1] <- "logit"
  if (margins[2] == "PH" && surv == TRUE) 
    margins[2] <- "cloglog"
  if (margins[2] == "PO" && surv == TRUE) 
    margins[2] <- "logit"
  bl <- c("probit", "logit", "cloglog")
  v.rB1 <- upperBt1
  v.rB2 <- upperBt2
  if (Model == "ROY") {
    L <- eval(substitute(SemiParROY(formula, data, weights, 
                                    subset, BivD1 = BivD, BivD2, margins, dof1 = dof, 
                                    dof2, gamlssfit, fp, infl.fac, rinit, rmax, iterlimsp, 
                                    tolsp, gc.l, parscale, extra.regI, knots = knots, 
                                    drop.unused.levels = drop.unused.levels, min.dn = min.dn, 
                                    min.pr = min.pr, max.pr = max.pr), list(weights = weights)))
  }
  else {
    if (surv == FALSE && ordinal == FALSE) {
      if ((margins[1] %in% bl && margins[2] %in% bl && 
           is.na(margins[3])) || (margins[1] %in% bl && 
                                  !(margins[2] %in% bl) && Model == "B" && is.na(margins[3]))) {
        L <- eval(substitute(SemiParBIV(formula, data, 
                                        weights, subset, Model, BivD, margins, dof, 
                                        gamlssfit, fp, hess = TRUE, infl.fac, rinit, 
                                        rmax, iterlimsp, tolsp, gc.l, parscale, extra.regI, 
                                        intf = TRUE, theta.fx = NULL, knots = knots, 
                                        drop.unused.levels = drop.unused.levels, min.dn = min.dn, 
                                        min.pr = min.pr, max.pr = max.pr), list(weights = weights)))
      }
    }
    if (surv == FALSE && ordinal == TRUE) {
      if ((margins[1] %in% bl && margins[2] %in% bl && 
           is.na(margins[3])) || (margins[1] %in% bl && 
                                  !(margins[2] %in% bl) && is.na(margins[3]))) {
        L <- eval(substitute(CopulaCLM(formula, data, 
                                       weights, subset, Model, BivD, margins, dof, 
                                       gamlssfit, fp, hess = TRUE, infl.fac, rinit, 
                                       rmax, iterlimsp, tolsp, gc.l, parscale, extra.regI, 
                                       intf = TRUE, theta.fx = NULL, knots = knots, 
                                       drop.unused.levels = drop.unused.levels, min.dn = min.dn, 
                                       min.pr = min.pr, max.pr = max.pr), list(weights = weights)))
      }
      else {
        stop("The first margin must be ordinal and the second either ordinal or continuous.")
      }
    }
    if (margins[1] %in% bl && !(margins[2] %in% bl) && surv == 
        FALSE && is.na(margins[3]) && Model == "BSS" && 
        ordinal == FALSE) {
      L <- eval(substitute(copulaSampleSel(formula, data, 
                                           weights, subset, BivD, margins, dof, fp, infl.fac, 
                                           rinit, rmax, iterlimsp, tolsp, gc.l, parscale, 
                                           extra.regI, knots, drop.unused.levels = drop.unused.levels, 
                                           min.dn = min.dn, min.pr = min.pr, max.pr = max.pr), 
                           list(weights = weights)))
    }
    if (!is.na(margins[3])) {
      if (margins[1] %in% bl && margins[2] %in% bl && 
          margins[3] %in% bl && surv == FALSE && ordinal == 
          FALSE) {
        L <- eval(substitute(SemiParTRIV(formula, data, 
                                         weights, subset, Model, margins, penCor, sp.penCor, 
                                         approx = FALSE, Chol, infl.fac, gamma, w.alasso, 
                                         rinit, rmax, iterlimsp, tolsp, gc.l, parscale, 
                                         extra.regI, knots, drop.unused.levels = drop.unused.levels, 
                                         min.dn = min.dn, min.pr = min.pr, max.pr = max.pr), 
                             list(weights = weights)))
      }
      else {
        stop("The model currently support only binary outcomes.")
      }
    }
    if ((!(margins[1] %in% bl) || surv == TRUE) && ordinal == 
        FALSE) {
      robust <- FALSE
      t.c = 3
      sp <- qu.mag <- y1.y2 <- y1.cy2 <- cy1.y2 <- cy1.cy2 <- cy <- cy1 <- gamlss1 <- gamlss2 <- gam1 <- gam2 <- y1m <- y2m <- indexTeq1 <- indexTeq2 <- NULL
      i.rho <- log.sig2.2 <- log.nu.2 <- log.nu.1 <- log.sig2.1 <- dof.st <- NULL
      end <- X3.d2 <- X4.d2 <- X5.d2 <- X6.d2 <- X7.d2 <- X8.d2 <- l.sp3 <- l.sp4 <- l.sp5 <- l.sp6 <- l.sp7 <- l.sp8 <- l.sp9 <- 0
      sp1 <- sp2 <- NULL
      sp3 <- gp3 <- gam3 <- X3 <- sp4 <- gp4 <- gam4 <- X4 <- sp5 <- gp5 <- gam5 <- X5 <- gam9 <- NULL
      sp6 <- gp6 <- gam6 <- X6 <- sp7 <- gp7 <- gam7 <- X7 <- sp8 <- gp8 <- gam8 <- X8 <- sp9 <- NULL
      c11 <- c10 <- c01 <- c00 <- NA
      cens1Mix <- cens2Mix <- NULL
      Sl.sf <- NULL
      sp.method <- "perf"
      Xd1 <- Xd2 <- mono.sm.pos1 <- mono.sm.pos2 <- mono.sm.pos <- NULL
      surv.flex <- FALSE
      Deq1 <- pos.pbeq1 <- Deq2 <- pos.pbeq2 <- list()
      BivD2 <- c("C0C90", "C0C270", "C180C90", "C180C270", 
                 "J0J90", "J0J270", "J180J90", "J180J270", "G0G90", 
                 "G0G270", "G180G90", "G180G270", "GAL0GAL90", 
                 "GAL0GAL270", "GAL180GAL90", "GAL180GAL270")
      opc <- c("N", "C0", "C90", "C180", "C270", "J0", 
               "J90", "J180", "J270", "G0", "G90", "G180", 
               "G270", "F", "AMH", "FGM", "T", "PL", "HO", 
               "GAL0", "GAL90", "GAL180", "GAL270")
      scc <- c("C0", "C180", "GAL0", "GAL180", "J0", "J180", 
               "G0", "G180", BivD2)
      sccn <- c("C90", "C270", "GAL90", "GAL270", "J90", 
                "J270", "G90", "G270")
      m2 <- c("N", "GU", "rGU", "LO", "LN", "WEI", "iG", 
              "GA", "BE", "FISK", "GP", "GPII", "GPo")
      m3 <- c("DAGUM", "SM", "TW")
      m1d <- c("PO", "ZTP", "DGP0")
      m2d <- c("NBI", "NBII", "PIG", "DGP", "DGPII")
      m3d <- c("DEL", "SICHEL")
      ct <- data.frame(c(opc), c(1:14, 55, 56, 57, 60, 
                                 61, 62:65))
      cta <- data.frame(c(opc), c(1, 3, 23, 13, 33, 6, 
                                  26, 16, 36, 4, 24, 14, 34, 5, 55, 56, 2, 60, 
                                  61, 62:65))
      if (BivD %in% BivD2) {
        if (BivD %in% BivD2[1:4]) 
          BivDt <- "C0"
        if (BivD %in% BivD2[5:12]) 
          BivDt <- "J0"
        if (BivD %in% BivD2[13:16]) 
          BivDt <- "C0"
        nC <- ct[which(ct[, 1] == BivDt), 2]
        nCa <- cta[which(cta[, 1] == BivDt), 2]
      }
      if (!(BivD %in% BivD2)) {
        nC <- ct[which(ct[, 1] == BivD), 2]
        nCa <- cta[which(cta[, 1] == BivD), 2]
      }
      if (!is.list(formula)) 
        stop("You must specify a list of equations.")
      l.flist <- length(formula)
      form.check(formula, l.flist)
      cl <- match.call()
      mf <- match.call(expand.dots = FALSE)
      pred.varR <- pred.var(formula, l.flist)
      v1 <- pred.varR$v1
      v2 <- pred.varR$v2
      pred.n <- pred.varR$pred.n
      if (!is.null(v.rB1)) 
        pred.n <- c(pred.n, v.rB1)
      if (!is.null(v.rB2)) 
        pred.n <- c(pred.n, v.rB2)
      fake.formula <- paste(v1[1], "~", paste(pred.n, 
                                              collapse = " + "))
      environment(fake.formula) <- environment(formula[[1]])
      mf$formula <- fake.formula
      mf$upperBt1 <- mf$upperBt2 <- mf$min.dn <- mf$min.pr <- mf$max.pr <- mf$dep.cens <- mf$ordinal <- mf$Model <- mf$model <- mf$knots <- mf$k1.tvc <- mf$k2.tvc <- mf$surv <- mf$BivD <- mf$copula <- mf$copula2 <- mf$margins <- mf$fp <- mf$dof <- mf$infl.fac <- mf$rinit <- mf$rmax <- mf$iterlimsp <- mf$tolsp <- mf$gc.l <- mf$parscale <- mf$extra.regI <- mf$gamlssfit <- NULL
      mf$drop.unused.levels <- drop.unused.levels
      if (surv == TRUE) 
        mf$na.action <- na.pass
      mf[[1]] <- as.name("model.frame")
      data <- eval(mf, parent.frame())
      zinf_formula <- formula[6:7]
      formula <- formula[1:5]
      if (surv == TRUE) {
        if (!("(cens1)" %in% names(data)) && margins[1] %in% 
            bl) 
          stop("You must provide both censoring indicators.")
        if (!("(cens2)" %in% names(data)) && margins[2] %in% 
            bl) 
          stop("You must provide both censoring indicators.")
      }
      if (gc.l == TRUE) 
        gc()
      if (!("(weights)" %in% names(data))) {
        weights <- rep(1, dim(data)[1])
        data$weights <- weights
        names(data)[length(names(data))] <- "(weights)"
      }
      else weights <- data[, "(weights)"]
      if (!("(cens1)" %in% names(data))) {
        cens1 <- rep(0, dim(data)[1])
        data$cens1 <- cens1
        names(data)[length(names(data))] <- "(cens1)"
      }
      else cens1 <- data[, "(cens1)"]
      if (!("(cens2)" %in% names(data))) {
        cens2 <- rep(0, dim(data)[1])
        data$cens2 <- cens2
        names(data)[length(names(data))] <- "(cens2)"
      }
      else cens2 <- data[, "(cens2)"]
      if (!("(cens3)" %in% names(data))) {
        cens3 <- rep(0, dim(data)[1])
        data$cens3 <- cens3
        names(data)[length(names(data))] <- "(cens3)"
      }
      else cens3 <- data[, "(cens3)"]
      if (surv == TRUE) {
        if (is.factor(cens1) && !is.factor(cens2)) 
          stop("Both censoring indicators have to be factor variables for mixed censoring case.")
        if (!is.factor(cens1) && is.factor(cens2)) 
          stop("Both censoring indicators have to be factor variables for mixed censoring case.")
      }
      if (surv == TRUE && is.factor(cens1) && !is.null(v.rB1)) 
        data[!(cens1 == "I"), v.rB1] <- data[!(cens1 == 
                                                 "I"), v1[1]]
      if (surv == TRUE && is.factor(cens2) && !is.null(v.rB2)) 
        data[!(cens2 == "I"), v.rB2] <- data[!(cens2 == 
                                                 "I"), v2[1]]
      if (surv == TRUE) {
        if (any(is.na(data[, v1[1]]) | is.na(data[, 
                                                  v2[1]]))) 
          stop("Time to event with NA's. Please check your time covariates.")
        actual.NAs = as.numeric(which(apply(apply(data, 
                                                  1, is.na), 2, any)))
        data <- na.omit(data)
        if (length(actual.NAs) > 0) {
          cens1 <- cens1[-actual.NAs]
          cens2 <- cens2[-actual.NAs]
        }
      }
      n <- dim(data)[1]
      if (surv == TRUE && is.factor(cens1) && is.factor(cens2)) {
        cens1Mix <- cens1
        cens2Mix <- cens2
        cens1 <- cens2 <- rep(1, n)
      }
      M <- list(m1d = m1d, m2 = m2, m2d = m2d, m3 = m3, 
                m3d = m3d, BivD = BivD, bl = bl, robust = robust, 
                opc = opc, extra.regI = extra.regI, margins = margins, 
                BivD2 = BivD2, dof = dof, surv = surv, c1 = cens1, 
                c2 = cens2, c3 = cens3, dep.cens = dep.cens)
      M$K1 <- NULL
      M$type.cens1 <- M$type.cens2 <- "R"
      # pream.wm(formula, margins, M, l.flist)
      formula.eq1 <- formula[[1]]
      formula.eq2 <- formula[[2]]
      form.eq12R <- form.eq12(formula.eq1, data, v1, margins[1], 
                              m1d, m2d)
      formula.eq1 <- form.eq12R$formula.eq1
      formula.eq1r <- form.eq12R$formula.eq1r
      y1 <- form.eq12R$y1
      y1.test <- form.eq12R$y1.test
      y1m <- form.eq12R$y1m
      
      # making it so the starting values are found using only the nonzero portions,
      # otherwise it will confuse mgcv which is assuming the marginals are non-zero-inflated GA

      y1weights <- weights*(y1 > 0)

      # end of modification
      
      if (surv == FALSE) 
        # making weights as y1weights so that zero values are ignored
        gam1 <- eval(substitute(gam(formula.eq1, gamma = infl.fac, 
                                    weights = y1weights, data = data, knots = knots, 
                                    drop.unused.levels = drop.unused.levels), 
                                list(y1weights = y1weights)))
        # end of modification
      if (surv == TRUE && margins[1] %in% c(m2, m3) && 
          margins[2] %in% bl) 
        gam1 <- eval(substitute(gam(formula.eq1, gamma = infl.fac, 
                                    weights = weights, data = data, knots = knots, 
                                    drop.unused.levels = drop.unused.levels), 
                                list(weights = weights)))
      else {
        if (surv == TRUE && !(margins[1] %in% bl)) 
          gam1 <- eval(substitute(gam(formula.eq1, gamma = infl.fac, 
                                      weights = weights * cens1, data = data, 
                                      knots = knots, drop.unused.levels = drop.unused.levels), 
                                  list(weights = weights, cens1 = cens1)))
      }
      if (surv == TRUE && margins[1] %in% bl) {
        surv.flex <- TRUE
        f.eq1 <- form.eq12R$f.eq1
        data$urcfcphmwicu <- seq(-10, 10, length.out = dim(data)[1])
        tempb <- eval(substitute(gam(f.eq1, family = cox.ph(), 
                                     data = data, weights = cens1, drop.unused.levels = drop.unused.levels), 
                                 list(cens1 = cens1)))
        data$Sh <- as.vector(mm(predict(tempb, type = "response"), 
                                min.pr = min.pr, max.pr = max.pr))
        cens11 <- ifelse(cens1 == 0, 1e-07, cens1)
        gam1 <- eval(substitute(scam(formula.eq1, gamma = infl.fac, 
                                     weights = weights * cens11, data = data), 
                                list(weights = weights, cens11 = cens11)))
        lsgam1 <- length(gam1$smooth)
        if (lsgam1 == 0) 
          stop("You must use at least a monotonic smooth function of time in the first equation.")
        clsm <- ggr <- NA
        for (i in 1:lsgam1) {
          clsm[i] <- class(gam1$smooth[[i]])[1]
        }
        if (sum(as.numeric(clsm[1] %in% c("mpi.smooth"))) == 
            0) 
          stop("You must have a monotonic smooth of time and it has to be the first to be included.")
        l.sp1 <- length(gam1$sp)
        if (l.sp1 != 0) 
          sp1 <- gam1$sp
        sp1[1] <- 1
        gam.call <- gam1$call
        gam.call$sp <- sp1
        gam1 <- eval(gam.call)
        j <- 1
        for (i in 1:lsgam1) {
          if (max(as.numeric(grepl(v1[1], gam1$smooth[[i]]$term))) != 
              0 && clsm[i] == "mpi.smooth") 
            mono.sm.pos1 <- c(mono.sm.pos1, c(gam1$smooth[[i]]$first.para:gam1$smooth[[i]]$last.para))
        }
        X1 <- predict(gam1, type = "lpmatrix")
        if (!is.null(indexTeq1) && k1.tvc != 0) {
          if (range(X1[, indexTeq1])[1] < 0) 
            stop("Check design matrix for smooth(s) of tvc term(s) in eq. 1.")
        }
        Xd1 <- Xdpred(gam1, data, v1[1])
        gam1$y <- data[, v1[1]]
        st.v1 <- c(gam1$coefficients)
        if (!is.null(indexTeq1)) {
          st.v1[mono.sm.pos1] <- exp(st.v1[mono.sm.pos1])
          while (range(Xd1 %*% st.v1)[1] < 0) st.v1[indexTeq1] <- 0.999 * 
              st.v1[indexTeq1]
          gam1$coefficients <- gam1$coefficients.t <- st.v1
          gam1$coefficients.t[mono.sm.pos1] <- exp(gam1$coefficients.t[mono.sm.pos1])
        }
      }
      gam1$formula <- formula.eq1r
      lsgam1 <- length(gam1$smooth)
      y1 <- y1.test
      if (margins[1] %in% c("LN")) 
        y1 <- log(y1)
      attr(data, "terms") <- NULL
      if (!(surv == TRUE && margins[1] %in% bl)) {
        names(gam1$model)[1] <- as.character(formula.eq1r[2])
        X1 <- predict(gam1, type = "lpmatrix")
        l.sp1 <- length(gam1$sp)
        sp1 <- gam1$sp
      }
      gp1 <- gam1$nsdf
      X1.d2 <- dim(X1)[2]
      form.eq12R <- form.eq12(formula.eq2, data, v2, margins[2], 
                              m1d, m2d)
      formula.eq2 <- form.eq12R$formula.eq1
      formula.eq2r <- form.eq12R$formula.eq1r
      y2 <- form.eq12R$y1
      y2.test <- form.eq12R$y1.test
      y2m <- form.eq12R$y1m
      
      # making it so the starting values are found using only the nonzero portions,
      # otherwise it will confuse mgcv which is assuming the marginals are non-zero-inflated GA
      
      y2weights <- weights*(y2 > 0)
      
      # end of modification
      
      if (surv == FALSE) 
        # making weights as y2weights so that zero values are ignored
        gam2 <- eval(substitute(gam(formula.eq2, gamma = infl.fac, 
                                    weights = y2weights, data = data, knots = knots, 
                                    drop.unused.levels = drop.unused.levels), 
                                list(y2weights = y2weights)))
        # end of modification
      if (surv == TRUE && !(margins[2] %in% bl)) 
        gam2 <- eval(substitute(gam(formula.eq2, gamma = infl.fac, 
                                    weights = weights * cens2, data = data, knots = knots, 
                                    drop.unused.levels = drop.unused.levels), 
                                list(weights = weights, cens2 = cens2)))
      if (surv == TRUE && margins[2] %in% bl) {
        surv.flex <- TRUE
        f.eq2 <- form.eq12R$f.eq1
        data$urcfcphmwicu <- seq(-10, 10, length.out = dim(data)[1])
        tempb <- eval(substitute(gam(f.eq2, family = cox.ph(), 
                                     data = data, weights = cens2, drop.unused.levels = drop.unused.levels), 
                                 list(cens2 = cens2)))
        data$Sh <- as.vector(mm(predict(tempb, type = "response"), 
                                min.pr = min.pr, max.pr = max.pr))
        cens22 <- ifelse(cens2 == 0, 1e-07, cens2)
        gam2 <- eval(substitute(scam(formula.eq2, gamma = infl.fac, 
                                     weights = weights * cens22, data = data), 
                                list(weights = weights, cens22 = cens22)))
        lsgam2 <- length(gam2$smooth)
        if (lsgam2 == 0) 
          stop("You must use at least a monotonic smooth function of time in the second equation.")
        clsm <- ggr <- NA
        for (i in 1:lsgam2) {
          clsm[i] <- class(gam2$smooth[[i]])[1]
        }
        if (sum(as.numeric(clsm[1] %in% c("mpi.smooth"))) == 
            0) 
          stop("You must have a monotonic smooth of time and it has to be the first to be included.")
        l.sp2 <- length(gam2$sp)
        if (l.sp2 != 0) 
          sp2 <- gam2$sp
        sp2[1] <- 1
        gam.call <- gam2$call
        gam.call$sp <- sp2
        gam2 <- eval(gam.call)
        j <- 1
        for (i in 1:lsgam2) {
          if (max(as.numeric(grepl(v2[1], gam2$smooth[[i]]$term))) != 
              0 && clsm[i] == "mpi.smooth") 
            mono.sm.pos2 <- c(mono.sm.pos2, c(gam2$smooth[[i]]$first.para:gam2$smooth[[i]]$last.para))
        }
        X2 <- predict(gam2, type = "lpmatrix")
        if (!is.null(indexTeq2) && k2.tvc != 0) {
          if (range(X2[, indexTeq2])[1] < 0) 
            stop("Check design matrix for smooth(s) of tvc term(s) in eq. 2.")
        }
        Xd2 <- Xdpred(gam2, data, v2[1])
        gam2$y <- data[, v2[1]]
        st.v2 <- c(gam2$coefficients)
        if (!is.null(indexTeq2)) {
          st.v2[mono.sm.pos2] <- exp(st.v2[mono.sm.pos2])
          while (range(Xd2 %*% st.v2)[1] < 0) st.v2[indexTeq2] <- 0.999 * 
              st.v2[indexTeq2]
          gam2$coefficients <- gam2$coefficients.t <- st.v2
          gam2$coefficients.t[mono.sm.pos2] <- exp(gam2$coefficients.t[mono.sm.pos2])
        }
      }
      gam2$formula <- formula.eq2r
      lsgam2 <- length(gam2$smooth)
      y2 <- y2.test
      if (margins[2] %in% c("LN")) 
        y2 <- log(y2)
      attr(data, "terms") <- NULL
      if (!(surv == TRUE && margins[2] %in% bl)) {
        names(gam2$model)[1] <- as.character(formula.eq2r[2])
        X2 <- predict(gam2, type = "lpmatrix")
        l.sp2 <- length(gam2$sp)
        sp2 <- gam2$sp
      }
      
      gp2 <- gam2$nsdf
      X2.d2 <- dim(X2)[2]
      
      # calculate logistic regression using mgcv for p1, p2
      data$resp1zero <- 1*(y1 == 0)
      data$resp2zero <- 1*(y2 == 0)

      p1gam <- gam(as.formula(paste("resp1zero~", zinf_formula[[1]][2], sep = "")), data=data, family=binomial())
      p2gam <- gam(as.formula(paste("resp2zero~", zinf_formula[[2]][2], sep="")), data=data, family=binomial())
      # end of modification
      
      
      # the correlation should only be measured where y1 and y2 are both nonzero
      res1 <- residuals(gam1)[(y1 > 0) & (y2 > 0)]
      res2 <- residuals(gam2)[(y1 > 0) & (y2 > 0)]
      # end of modification
      ass.s <- cor(res1, res2, method = "kendall")
      ass.s <- sign(ass.s) * ifelse(abs(ass.s) > 0.9, 
                                    0.9, abs(ass.s))
      i.rho <- ass.dp(ass.s, BivD, scc, sccn, nCa)
      dof.st <- log(dof - 2)
      names(dof.st) <- "dof.star"
      if (!(margins[1] %in% c(m1d, bl))) {
        # this sigma estimate needs to only come from the positive y1 values
        # not the zero inflated portion
        start.snR <- startsn(margins[1], y1[y1 > 0])
        # end of modification
        log.sig2.1 <- start.snR$log.sig2.1
        names(log.sig2.1) <- "sigma1.star"
        if (margins[1] %in% c(m3)) {
          log.nu.1 <- start.snR$log.nu.1
          names(log.nu.1) <- "nu.1.star"
        }
      }
      if (!(margins[2] %in% c(m1d, bl))) {
        # this sigma estimate needs to only come from the positive y1 values
        # not the zero inflated portion
        start.snR <- startsn(margins[2], y2[y2 > 0])
        # end of modification
        log.sig2.2 <- start.snR$log.sig2.1
        names(log.sig2.2) <- "sigma2.star"
        if (margins[2] %in% c(m3)) {
          log.nu.2 <- start.snR$log.nu.1
          names(log.nu.2) <- "nu.2.star"
        }
      }
      vo <- list(gam1 = gam1, gam2 = gam2, i.rho = i.rho, 
                 log.sig2.2 = log.sig2.2, log.nu.2 = log.nu.2, 
                 log.nu.1 = log.nu.1, log.sig2.1 = log.sig2.1, 
                 dof.st = dof.st, n = n, drop.unused.levels = drop.unused.levels)
      start.v <- overall.sv(margins, M, vo)
      if (l.flist > 2) {
        overall.svGR <- overall.svG(formula, data, ngc = 2, 
                                    margins, M, vo, gam1, gam2, knots = knots)
        start.v = overall.svGR$start.v
        X3 = overall.svGR$X3
        X4 = overall.svGR$X4
        X5 = overall.svGR$X5
        X6 = overall.svGR$X6
        X7 = overall.svGR$X7
        X8 = overall.svGR$X8
        X3.d2 = overall.svGR$X3.d2
        X4.d2 = overall.svGR$X4.d2
        X5.d2 = overall.svGR$X5.d2
        X6.d2 = overall.svGR$X6.d2
        X7.d2 = overall.svGR$X7.d2
        X8.d2 = overall.svGR$X8.d2
        gp3 = overall.svGR$gp3
        gp4 = overall.svGR$gp4
        gp5 = overall.svGR$gp5
        gp6 = overall.svGR$gp6
        gp7 = overall.svGR$gp7
        gp8 = overall.svGR$gp8
        gam3 = overall.svGR$gam3
        gam4 = overall.svGR$gam4
        gam5 = overall.svGR$gam5
        gam6 = overall.svGR$gam6
        gam7 = overall.svGR$gam7
        gam8 = overall.svGR$gam8
        l.sp3 = overall.svGR$l.sp3
        l.sp4 = overall.svGR$l.sp4
        l.sp5 = overall.svGR$l.sp5
        l.sp6 = overall.svGR$l.sp6
        l.sp7 = overall.svGR$l.sp7
        l.sp8 = overall.svGR$l.sp8
        sp3 = overall.svGR$sp3
        sp4 = overall.svGR$sp4
        sp5 = overall.svGR$sp5
        sp6 = overall.svGR$sp6
        sp7 = overall.svGR$sp7
        sp8 = overall.svGR$sp8
      }
      GAM <- list(gam1 = gam1, gam2 = gam2, gam3 = gam3, 
                  gam4 = gam4, gam5 = gam5, gam6 = gam6, gam7 = gam7, 
                  gam8 = gam8, gam9 = gam9)
      if ((l.sp1 != 0 || l.sp2 != 0 || l.sp3 != 0 || l.sp4 != 
           0 || l.sp5 != 0 || l.sp6 != 0 || l.sp7 != 0 || 
           l.sp8 != 0) && fp == FALSE) {
        L.GAM <- list(l.gam1 = length(gam1$coefficients), 
                      l.gam2 = length(gam2$coefficients), l.gam3 = length(gam3$coefficients), 
                      l.gam4 = length(gam4$coefficients), l.gam5 = length(gam5$coefficients), 
                      l.gam6 = length(gam6$coefficients), l.gam7 = length(gam7$coefficients), 
                      l.gam8 = length(gam8$coefficients), l.gam9 = 0)
        L.SP <- list(l.sp1 = l.sp1, l.sp2 = l.sp2, l.sp3 = l.sp3, 
                     l.sp4 = l.sp4, l.sp5 = l.sp5, l.sp6 = l.sp6, 
                     l.sp7 = l.sp7, l.sp8 = l.sp8, l.sp9 = l.sp9)
        sp <- c(sp1, sp2, sp3, sp4, sp5, sp6, sp7, sp8, 
                sp9)
        qu.mag <- S.m(GAM, L.SP, L.GAM)
      }
      if (missing(parscale)) 
        parscale <- 1
      respvec <- respvec2 <- respvec3 <- list(y1 = y1, 
                                              y2 = y2, y1.y2 = NULL, y1.cy2 = NULL, cy1.y2 = NULL, 
                                              cy1.cy2 = NULL, cy1 = NULL, cy = NULL, univ = 0)
      my.env <- new.env()
      my.env$signind <- 1
      lsgam3 <- length(gam3$smooth)
      lsgam4 <- length(gam4$smooth)
      lsgam5 <- length(gam5$smooth)
      lsgam6 <- length(gam6$smooth)
      lsgam7 <- length(gam7$smooth)
      lsgam8 <- length(gam8$smooth)
      lsgam9 <- length(gam9$smooth)
      indUR <- indUL <- indUI <- indUU <- indRR <- indRL <- indRI <- indRU <- indLR <- indLL <- indLI <- indLU <- indIR <- indIL <- indII <- indIU <- rep(0, 
                                                                                                                                                          n)
      if (surv == TRUE && dep.cens == FALSE) {
        if ((surv == TRUE && margins[1] %in% bl && margins[2] %in% 
             bl && !is.factor(cens1) && !is.factor(cens2)) || 
            (surv == TRUE && margins[1] %in% m2 && margins[2] %in% 
             m2)) {
          c11 <- cens1 * cens2
          c10 <- cens1 * (1 - cens2)
          c01 <- (1 - cens1) * cens2
          c00 <- (1 - cens1) * (1 - cens2)
        }
        if (surv == TRUE && margins[1] %in% c(m2, m3) && 
            margins[2] %in% bl) {
          c11 <- cens2
          c10 <- 1 - cens2
          c01 <- NULL
          c00 <- NULL
        }
        if (!is.null(cens1Mix) && !is.null(cens2Mix)) {
          if (surv == TRUE && margins[1] %in% bl && 
              margins[2] %in% bl && is.factor(cens1Mix) && 
              is.factor(cens2Mix)) {
            gamlssfit <- TRUE
            indUR <- as.numeric(cens1Mix == "U") * as.numeric(cens2Mix == 
                                                                "R")
            indUL <- as.numeric(cens1Mix == "U") * as.numeric(cens2Mix == 
                                                                "L")
            indUI <- as.numeric(cens1Mix == "U") * as.numeric(cens2Mix == 
                                                                "I")
            indUU <- as.numeric(cens1Mix == "U") * as.numeric(cens2Mix == 
                                                                "U")
            indRR <- as.numeric(cens1Mix == "R") * as.numeric(cens2Mix == 
                                                                "R")
            indRL <- as.numeric(cens1Mix == "R") * as.numeric(cens2Mix == 
                                                                "L")
            indRI <- as.numeric(cens1Mix == "R") * as.numeric(cens2Mix == 
                                                                "I")
            indRU <- as.numeric(cens1Mix == "R") * as.numeric(cens2Mix == 
                                                                "U")
            indLR <- as.numeric(cens1Mix == "L") * as.numeric(cens2Mix == 
                                                                "R")
            indLL <- as.numeric(cens1Mix == "L") * as.numeric(cens2Mix == 
                                                                "L")
            indLI <- as.numeric(cens1Mix == "L") * as.numeric(cens2Mix == 
                                                                "I")
            indLU <- as.numeric(cens1Mix == "L") * as.numeric(cens2Mix == 
                                                                "U")
            indIR <- as.numeric(cens1Mix == "I") * as.numeric(cens2Mix == 
                                                                "R")
            indIL <- as.numeric(cens1Mix == "I") * as.numeric(cens2Mix == 
                                                                "L")
            indII <- as.numeric(cens1Mix == "I") * as.numeric(cens2Mix == 
                                                                "I")
            indIU <- as.numeric(cens1Mix == "I") * as.numeric(cens2Mix == 
                                                                "U")
          }
        }
      }
      if (surv == TRUE && dep.cens == TRUE) {
        c11 <- NULL
        c10 <- cens1
        c01 <- cens2
        c00 <- cens3
      }
      my.env$k1 <- k1.tvc
      my.env$k2 <- k2.tvc
      VC <- list(lsgam1 = lsgam1, indexTeq1 = indexTeq1, 
                 indexTeq2 = indexTeq2, lsgam2 = lsgam2, Deq1 = Deq1, 
                 pos.pbeq1 = pos.pbeq1, Deq2 = Deq2, pos.pbeq2 = pos.pbeq2, 
                 lsgam3 = lsgam3, robust = FALSE, sp.fixed = NULL, 
                 lsgam4 = lsgam4, Sl.sf = Sl.sf, sp.method = sp.method, 
                 lsgam5 = lsgam5, K1 = NULL, lsgam6 = lsgam6, 
                 lsgam7 = lsgam7, lsgam8 = lsgam8, lsgam9 = lsgam9, 
                 X1 = X1, X2 = X2, X3 = X3, X4 = X4, X5 = X5, 
                 X6 = X6, X7 = X7, X8 = X8, X1.d2 = X1.d2, X2.d2 = X2.d2, 
                 X3.d2 = X3.d2, X4.d2 = X4.d2, X5.d2 = X5.d2, 
                 X6.d2 = X6.d2, X7.d2 = X7.d2, X8.d2 = X8.d2, 
                 gp1 = gp1, gp2 = gp2, gp3 = gp3, gp4 = gp4, 
                 gp5 = gp5, gp6 = gp6, gp7 = gp7, gp8 = gp8, 
                 l.sp1 = l.sp1, l.sp2 = l.sp2, l.sp3 = l.sp3, 
                 l.sp4 = l.sp4, l.sp5 = l.sp5, l.sp6 = l.sp6, 
                 l.sp7 = l.sp7, l.sp8 = l.sp8, l.sp9 = 0, my.env = my.env, 
                 infl.fac = infl.fac, weights = weights, fp = fp, 
                 gamlssfit = gamlssfit, hess = NULL, Model = "CC", 
                 univ.gamls = FALSE, model = model, end = end, 
                 BivD = BivD, nCa = nCa, copula = copula, copula2 = copula2, 
                 nC = nC, gc.l = gc.l, n = n, extra.regI = extra.regI, 
                 parscale = parscale, margins = margins, Cont = "YES", 
                 ccss = "no", m2 = m2, m3 = m3, m1d = m1d, m2d = m2d, 
                 m3d = m3d, bl = bl, triv = FALSE, y1m = y1m, 
                 y2m = y2m, tc = t.c, i.rho = i.rho, dof = dof, 
                 dof.st = dof.st, BivD2 = BivD2, cta = cta, ct = ct, 
                 zerov = -10, c11 = c11, c10 = c10, c01 = c01, 
                 c00 = c00, indUR = indUR, indUL = indUL, indUI = indUI, 
                 indUU = indUU, indRR = indRR, indRL = indRL, 
                 indRI = indRI, indRU = indRU, indLR = indLR, 
                 indLL = indLL, indLI = indLI, indLU = indLU, 
                 indIR = indIR, indIL = indIL, indII = indII, 
                 indIU = indIU, surv = surv, Xd1 = Xd1, Xd2 = Xd2, 
                 mono.sm.pos1 = mono.sm.pos1, mono.sm.pos2 = mono.sm.pos2, 
                 surv.flex = surv.flex, mono.sm.pos = mono.sm.pos, 
                 gp2.inf = NULL, informative = "no", zero.tol = 0.01, 
                 min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)
      if (gc.l == TRUE) 
        gc()
      

      if (gamlssfit == TRUE) {
        type.cens1 <- type.cens2 <- "R"
        surv1 <- surv2 <- surv
        form.gamlR <- form.gaml(formula, l.flist, M)
        if (surv == TRUE && margins[1] %in% c(m2, m3) && 
            margins[2] %in% bl) 
          surv1 <- FALSE
        if (surv == TRUE && margins[1] %in% bl && margins[2] %in% 
            bl && is.factor(cens1Mix) && is.factor(cens2Mix)) {
          cens1 <- cens1Mix
          cens2 <- cens2Mix
          type.cens1 <- type.cens2 <- "mixed"
          M$type.cens1 = type.cens1
          M$type.cens2 = type.cens2
        }
        # changed weights to y1weights so that zeros would not be used in fit
        gamlss1 <- eval(substitute(gamlss(form.gamlR$formula.gamlss1, 
                                          data = data, weights = y1weights, subset = subset, 
                                          margin = margins[1], surv = surv1, cens = cens1, 
                                          type.cens = type.cens1, upperB = upperBt1, 
                                          infl.fac = infl.fac, rinit = rinit, rmax = rmax, 
                                          iterlimsp = iterlimsp, tolsp = tolsp, gc.l = gc.l, 
                                          parscale = 1, extra.regI = extra.regI, k.tvc = k1.tvc, 
                                          drop.unused.levels = drop.unused.levels), 
                                   list(y1weights = y1weights, cens1 = cens1)))
        # end of modification
        # changed weights to y2weights so that zeros would not be used in fit
        gamlss2 <- eval(substitute(gamlss(form.gamlR$formula.gamlss2, 
                                          data = data, weights = y2weights, subset = subset, 
                                          margin = margins[2], surv = surv2, cens = cens2, 
                                          type.cens = type.cens2, upperB = upperBt2, 
                                          infl.fac = infl.fac, rinit = rinit, rmax = rmax, 
                                          iterlimsp = iterlimsp, tolsp = tolsp, gc.l = gc.l, 
                                          parscale = 1, extra.regI = extra.regI, k.tvc = k2.tvc, 
                                          drop.unused.levels = drop.unused.levels), 
                                   list(y2weights = y2weights, cens2 = cens2)))
        # end of modification
        SP <- list(sp1 = sp1, sp2 = sp2, sp3 = sp3, 
                   sp4 = sp4, sp5 = sp5, sp6 = sp6, sp7 = sp7, 
                   sp8 = sp8)
        gamls.upsvR <- gamls.upsv(gamlss1, gamlss2, 
                                  margins, M, l.flist, nstv = names(start.v), 
                                  VC, GAM, SP)
        sp <- gamls.upsvR$sp
        start.v <- gamls.upsvR$start.v
        VC$X1 <- gamlss1$VC$X1
        VC$Xd1 <- gamlss1$VC$Xd1
        VC$X1.2 <- gamlss1$VC$X2
        VC$X2 <- gamlss2$VC$X1
        VC$Xd2 <- gamlss2$VC$Xd1
        VC$X2.2 <- gamlss2$VC$X2
        rangeSurv1 <- gamlss1$rangeSurv
        rangeSurv2 <- gamlss2$rangeSurv
      }
      
      # replacing func.opt with zinf version
      # func.opt <- func.OPT(margins, M)
      
      func.opt <- zinfgammacop.override
      
      # end of modification

      SemiParFit <- SemiParBIV.fit(func.opt = func.opt, 
                                   start.v = start.v, rinit = rinit, rmax = rmax, 
                                   iterlim = 100, iterlimsp = iterlimsp, tolsp = tolsp, 
                                   respvec = respvec, VC = VC, sp = sp, qu.mag = qu.mag)
      SemiParFit.p <- copulaReg.fit.post(SemiParFit = SemiParFit, 
                                         VC = VC, GAM)
      y1.m <- y1
      if (margins[1] == "LN") 
        y1.m <- exp(y1)
      y2.m <- y2
      if (margins[2] == "LN") 
        y2.m <- exp(y2)
      SemiParFit <- SemiParFit.p$SemiParFit
      if (gc.l == TRUE) 
        gc()
      cov.c(SemiParFit)
      formula[6:7] <- zinf_formula
      gam1$call$data <- gam2$call$data <- gam3$call$data <- gam4$call$data <- gam5$call$data <- gam6$call$data <- gam7$call$data <- gam8$call$data <- cl$data
      L <- list(fit = SemiParFit$fit, dataset = NULL, 
                n = n, gamlss1 = gamlss1, gamlss2 = gamlss2, 
                formula = formula, robust = FALSE, edf11 = SemiParFit.p$edf11, 
                surv = surv, gam1 = gam1, gam2 = gam2, gam3 = gam3, 
                gam4 = gam4, gam5 = gam5, gam6 = gam6, gam7 = gam7, 
                gam8 = gam8, coefficients = SemiParFit$fit$argument, 
                coef.t = SemiParFit.p$coef.t, iterlimsp = iterlimsp, 
                weights = weights, cens1 = cens1, cens2 = cens2, 
                cens3 = cens3, sp = SemiParFit.p$sp, iter.sp = SemiParFit$iter.sp, 
                l.sp1 = l.sp1, l.sp2 = l.sp2, l.sp3 = l.sp3, 
                l.sp4 = l.sp4, l.sp5 = l.sp5, l.sp6 = l.sp6, 
                l.sp7 = l.sp7, l.sp8 = l.sp8, bl = bl, l.sp9 = l.sp9, 
                gam9 = gam9, fp = fp, iter.if = SemiParFit$iter.if, 
                iter.inner = SemiParFit$iter.inner, theta = SemiParFit.p$theta, 
                theta.a = SemiParFit.p$theta.a, sigma21 = SemiParFit.p$sigma21, 
                sigma22 = SemiParFit.p$sigma22, sigma21.a = SemiParFit.p$sigma21.a, 
                sigma22.a = SemiParFit.p$sigma22.a, sigma1 = SemiParFit.p$sigma21, 
                sigma2 = SemiParFit.p$sigma22, sigma1.a = SemiParFit.p$sigma21.a, 
                sigma2.a = SemiParFit.p$sigma22.a, nu1 = SemiParFit.p$nu1, 
                nu2 = SemiParFit.p$nu2, nu1.a = SemiParFit.p$nu1.a, 
                nu2.a = SemiParFit.p$nu2.a, dof.a = SemiParFit.p$dof.a, 
                dof = SemiParFit.p$dof, X1 = X1, X2 = X2, X3 = X3, 
                X4 = X4, X5 = X5, X6 = X6, X7 = X7, X8 = X8, 
                X1.d2 = X1.d2, X2.d2 = X2.d2, X3.d2 = X3.d2, 
                X4.d2 = X4.d2, X5.d2 = X5.d2, X6.d2 = X6.d2, 
                X7.d2 = X7.d2, X8.d2 = X8.d2, He = SemiParFit.p$He, 
                HeSh = SemiParFit.p$HeSh, Vb = SemiParFit.p$Vb, 
                Ve = SemiParFit.p$Ve, F = SemiParFit.p$F, F1 = SemiParFit.p$F1, 
                t.edf = SemiParFit.p$t.edf, edf = SemiParFit.p$edf, 
                edf1 = SemiParFit.p$edf1, edf2 = SemiParFit.p$edf2, 
                edf3 = SemiParFit.p$edf3, edf4 = SemiParFit.p$edf4, 
                edf5 = SemiParFit.p$edf5, edf6 = SemiParFit.p$edf6, 
                edf7 = SemiParFit.p$edf7, edf8 = SemiParFit.p$edf8, 
                edf1.1 = SemiParFit.p$edf1.1, edf1.2 = SemiParFit.p$edf1.2, 
                edf1.3 = SemiParFit.p$edf1.3, edf1.4 = SemiParFit.p$edf1.4, 
                edf1.5 = SemiParFit.p$edf1.5, edf1.6 = SemiParFit.p$edf1.6, 
                edf1.7 = SemiParFit.p$edf1.7, edf1.8 = SemiParFit.p$edf1.8, 
                R = SemiParFit.p$R, bs.mgfit = SemiParFit$bs.mgfit, 
                conv.sp = SemiParFit$conv.sp, wor.c = SemiParFit$wor.c, 
                eta1 = SemiParFit$fit$eta1, eta2 = SemiParFit$fit$eta2, 
                etad = SemiParFit$fit$etad, etas1 = SemiParFit$fit$etas1, 
                etas2 = SemiParFit$fit$etas2, y1 = y1.m, y2 = y2.m, 
                BivD = BivD, margins = margins, copula = copula, 
                copula2 = copula2, logLik = SemiParFit.p$logLik, 
                nC = nC, respvec = respvec, hess = TRUE, qu.mag = qu.mag, 
                gp1 = gp1, gp2 = gp2, gp3 = gp3, gp4 = gp4, 
                gp5 = gp5, gp6 = gp6, gp7 = gp7, gp8 = gp8, 
                VC = VC, magpp = SemiParFit$magpp, gamlssfit = gamlssfit, 
                Cont = "YES", tau = SemiParFit.p$tau, tau.a = SemiParFit.p$tau.a, 
                l.flist = l.flist, v1 = v1, v2 = v2, triv = FALSE, 
                univar.gamlss = FALSE, BivD2 = BivD2, call = cl, 
                surv = surv, surv.flex = surv.flex, Vb.t = SemiParFit.p$Vb.t, 
                coef.t = SemiParFit.p$coef.t, Model = "CC", 
                model = model, p1gam=p1gam, p2gam=p2gam)
      if (BivD %in% BivD2) {
        L$teta1 <- SemiParFit$fit$teta1
        L$teta.ind1 <- SemiParFit$fit$teta.ind1
        L$teta2 <- SemiParFit$fit$teta2
        L$teta.ind2 <- SemiParFit$fit$teta.ind2
        L$Cop1 <- SemiParFit$fit$Cop1
        L$Cop2 <- SemiParFit$fit$Cop2
      }
      class(L) <- c("zinfgammacop")
    }
  }
  L
}





form.eq12 <- function (formula.eq1, data, v1, margins, m1d, m2d, copSS = FALSE, 
          inde = NULL) 
{
  y1m <- f.eq1 <- NULL
  formula.eq1r <- formula.eq1
  y1 <- y1.test <- data[, v1[1]]
  if (copSS == TRUE) 
    y1 <- y1.test <- y1[inde]
  if (v1[1] != as.character(formula.eq1r[2])) 
    y1.test <- try(data[, as.character(formula.eq1r[2])], 
                   silent = TRUE)
  if (inherits(y1.test, "try-error")) 
    stop("Please check the syntax of the equations' responses.")
  if (margins %in% c(m1d, m2d) && min(y1.test, na.rm = TRUE) < 
      0) 
    stop("The response of one or more margins must be positive.")
  if (margins %in% c(m1d, m2d)) {
    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - 
                                                                       round(x)) < tol
    if (sum(as.numeric(is.wholenumber(y1.test))) != length(y1.test)) 
      stop("The response of one or more margins must be discrete.")
  }
  if (margins %in% c("ZTP") && min(y1.test, na.rm = TRUE) < 
      1) 
    stop("The response of one or more margins must be greater than 0.")
  
  # modifying so that the zeros in the zero-inflated gammas don't raise an error
  
  # if (margins %in% c("probit", "logit", "cloglog", "LN", "WEI", 
  #                    "GO", "iG", "GA", "GAi", "DAGUM", "SM", "FISK") && min(y1.test, 
  #                                                                           na.rm = TRUE) <= 0) 
  #   stop("The response of one or more margins must be positive.")
  
  # end of modification
  
  if (margins %in% c("TW") && min(y1.test, na.rm = TRUE) < 
      0) 
    stop("The response of one or more margins must be >= 0.")
  if (margins %in% c("BE") && (min(y1.test, na.rm = TRUE) <= 
                               0 || max(y1.test, na.rm = TRUE) >= 1)) 
    stop("The response of one or more margins must be in the interval (0,1).")
  if (margins == "GEVlink" && length(table(y1.test)) != 2) 
    stop("The response must be binary.")
  if (margins %in% c("NBIa", "NBIIa", "NBI", "PO", "ZTP", 
                     "DGP", "DGPII", "DGP0")) {
    ly1 <- length(y1)
    y1m <- list()
    my1 <- max(y1)
    if (margins != "ZTP") 
      for (i in 1:ly1) {
        y1m[[i]] <- seq(0, y1[i])
        length(y1m[[i]]) <- my1 + 1
      }
    if (margins == "ZTP") 
      for (i in 1:ly1) {
        y1m[[i]] <- seq(1, y1[i])
        length(y1m[[i]]) <- my1
      }
    y1m <- do.call(rbind, y1m)
    if (max(y1) > 170 && margins %in% c("PO", "ZTP", "DGP0")) 
      y1m <- mpfr(y1m, pmax(53, getPrec(y1)))
  }
  if (margins %in% c("N", "LO", "GU", "rGU", "GAi")) 
    formula.eq1 <- update(formula.eq1, (. + mean(.))/2 ~ 
                            .)
  if (margins %in% c(m1d, m2d) && margins != "GEVlink" && 
      margins != "DGP") 
    formula.eq1 <- update(formula.eq1, log((. + mean(.))/2) ~ 
                            .)
  if (margins %in% c("LN")) 
    formula.eq1 <- update(formula.eq1, (log(.) + mean(log(.)))/2 ~ 
                            .)
  if (margins %in% c("iG", "GA", "GGA", "DAGUM", "SM", "FISK", 
                     "TW")) 
    formula.eq1 <- update(formula.eq1, log((. + mean(.))/2) ~ 
                            .)
  if (margins %in% c("WEI")) 
    formula.eq1 <- update(formula.eq1, log(exp(log(.) + 
                                                 0.5772/(1.283/sqrt(var(log(.)))))) ~ .)
  if (margins %in% c("BE")) 
    formula.eq1 <- update(formula.eq1, qlogis((. + mean(.))/2) ~ 
                            .)
  if (margins %in% c("GP", "GPII", "GPo", "DGP", "DGPII", 
                     "DGP0")) 
    formula.eq1 <- update(formula.eq1, estobXiGP ~ .)
  f.eq1LI <- temp.respV ~ urcfcphmwicu
  if (margins %in% c("probit")) {
    f.eq1 <- update(formula.eq1r, . ~ urcfcphmwicu)
    formula.eq1 <- update(formula.eq1, -qnorm(Sh) ~ .)
  }
  if (margins %in% c("logit")) {
    f.eq1 <- update(formula.eq1r, . ~ urcfcphmwicu)
    formula.eq1 <- update(formula.eq1, -qlogis(Sh) ~ .)
  }
  if (margins %in% c("cloglog")) {
    f.eq1 <- update(formula.eq1r, . ~ urcfcphmwicu)
    formula.eq1 <- update(formula.eq1, log(-log(Sh)) ~ .)
  }
  list(f.eq1LI = f.eq1LI, formula.eq1 = formula.eq1, formula.eq1r = formula.eq1r, 
       y1 = y1, y1.test = y1.test, y1m = y1m, f.eq1 = f.eq1)
}









post.check <- 
function (x, main = "Histogram and Density Estimate of Residuals", 
          main2 = "Histogram and Density Estimate of Residuals", xlab = "Quantile Residuals", 
          xlab2 = "Quantile Residuals", intervals = FALSE, n.sim = 100, 
          prob.lev = 0.05, ...) 
{
  y1m <- y2m <- y3m <- NA
  qr <- qr1 <- qr2 <- NULL
  if (x$Cont == "YES") {
    y1 <- x$y1
    y2 <- x$y2
    if (x$VC$margins[1] %in% c(x$VC$m2, x$VC$m3)) {
      p1 <- distrHsAT(x$y1, x$eta1, x$sigma21, x$nu1, 
                      x$margins[1], min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, 
                      max.pr = x$VC$max.pr)$p2
      y1zero <- which(x$y1 == 0)
      p1[y1zero] <- x$p1gam$fitted.values[y1zero]
      if (x$VC$margins[1] %in% c("TW")) {
        if (any(x$y1 == 0) == TRUE) 
          p1[x$y1 == 0] <- runif(sum(x$y1 == 0), 
                                 min = 0, max = p1[x$y1 == 0])
      }
    }
    if (x$VC$margins[2] %in% c(x$VC$m2, x$VC$m3)) {
      p2 <- distrHsAT(x$y2, x$eta2, x$sigma22, x$nu2, 
                      x$margins[2], min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, 
                      max.pr = x$VC$max.pr)$p2
      y2zero <- which(x$y2 == 0)
      p2[y2zero] <- x$p2gam$fitted.values[y2zero]
      if (x$VC$margins[2] %in% c("TW")) {
        if (any(x$y2 == 0) == TRUE) 
          p2[x$y2 == 0] <- runif(sum(x$y2 == 0), 
                                 min = 0, max = p2[x$y2 == 0])
      }
    }
    par(mfrow = c(2, 2))
    qr1 <- qnorm(p1)
    hist(qr1, freq = FALSE, main = main, xlab = xlab, 
         ylab = "Density", ...)
    lines(density(qr1, adjust = 2), lwd = 2)
    if (intervals == FALSE) {
      qqnorm(qr1)
      abline(0, 1, col = "red")
    }
    if (intervals == TRUE) 
      int.postcheck(x, x$VC$margins[1], n.rep = n.sim, 
                    prob.lev = prob.lev, y2m = y1m, eq = 1)
    qr2 <- qnorm(p2)
    hist(qr2, freq = FALSE, main = main2, xlab = xlab2, 
         ylab = "Density", ...)
    lines(density(qr2, adjust = 2), lwd = 2)
    if (intervals == FALSE) {
      qqnorm(qr2)
      abline(0, 1, col = "red")
    }
    if (intervals == TRUE) 
      int.postcheck(x, x$VC$margins[2], n.rep = n.sim, 
                    prob.lev = prob.lev, y2m = y2m, eq = 2)
  }
  L <- list(qr = qr, qr1 = qr1, qr2 = as.numeric(qr2))
  invisible(L)
  par(mfrow = c(1, 1))
}
















predict.zinfgammacop <- 
function (object, eq, ...) 
{
  if (missing(eq)) 
    stop("You must provide the equation number.")
  # if (eq > object$l.flist) 
  #   stop("The fitted model has a smaller number of equations.")
  if (eq <= 5){
    if (eq == 1) {
      ss.pred <- object$gam1
      ind <- 1:object$X1.d2
    }
    if (eq == 2) {
      ss.pred <- object$gam2
      ind <- (object$X1.d2 + 1):(object$X1.d2 + object$X2.d2)
    }
    if (eq == 3) {
      ss.pred <- object$gam3
      ind <- (object$X1.d2 + object$X2.d2 + 1):(object$X1.d2 + 
                                                  object$X2.d2 + object$X3.d2)
    }
    if (eq == 4) {
      ss.pred <- object$gam4
      ind <- (object$X1.d2 + object$X2.d2 + object$X3.d2 + 
                1):(object$X1.d2 + object$X2.d2 + object$X3.d2 + 
                      object$X4.d2)
    }
    if (eq == 5) {
      ss.pred <- object$gam5
      ind <- (object$X1.d2 + object$X2.d2 + object$X3.d2 + 
                object$X4.d2 + 1):(object$X1.d2 + object$X2.d2 + 
                                     object$X3.d2 + object$X4.d2 + object$X5.d2)
    }
    ss.pred$coefficients <- object$coefficients[ind]
    ss.pred$coefficients.t <- object$coef.t[ind]
    ss.pred$Vp <- object$Vb[ind, ind]
    ss.pred$Vp.t <- object$Vb.t[ind, ind]
    ss.pred$sig2 <- 1
    ss.pred$scale.estimated <- FALSE    
  } else {
    
    if (eq == 6) {
      ss.pred <- object$p1gam
    }
    if (eq == 7) {
      ss.pred <- object$p2gam
    }
    
  }

  predict(ss.pred, ...)
}




summary.zinfgammacop <- function (object, n.sim = 100, prob.lev = 0.05, ...) 
{
  bs <- SE <- Vb <- epds <- sigma2.st <- sigma2 <- nu.st <- nu <- est.RHOb <- XX <- Xt <- V <- 1
  cont1par <- c(object$VC$m1d, object$VC$bl)
  cont2par <- c(object$VC$m2, object$VC$m2d)
  cont3par <- c(object$VC$m3, object$VC$m3d)
  n <- object$n
  lf <- length(object$coefficients)
  Vb <- object$Vb
  SE <- sqrt(diag(Vb))
  bs <- rMVN(n.sim, mean = object$coefficients, sigma = Vb)
  susutsnR <- susutsn(object, bs, lf, cont1par, cont2par, 
                      cont3par, prob.lev)
  CIrs <- susutsnR$CIrs
  CIkt <- susutsnR$CIkt
  CIsig21 <- susutsnR$CIsig21
  CIsig22 <- susutsnR$CIsig22
  CInu1 <- susutsnR$CInu1
  CInu2 <- susutsnR$CInu2
  CIdof <- susutsnR$CIdof
  if (object$VC$gc.l == TRUE) 
    gc()
  susuR <- susu(object, SE, Vb)
  tableN <- susuR$tableN
  table <- susuR$table
  rm(bs, SE, Vb, XX, Xt, V)
  # incorporating the two zero inflation parameters
  p1summ <- summary(object$p1gam)
  p2summ <- summary(object$p2gam)
  # this includes changing tableP6, tableNP7, formula6, l.sp6 etc. etc.
  res <- list(tableP1 = table[[1]], tableP2 = table[[2]], 
              tableP3 = table[[3]], tableP4 = table[[4]], tableP5 = table[[5]], 
              tableP6 = p1summ$p.table, tableP7 = p2summ$p.table, tableP8 = table[[8]], 
              tableNP1 = tableN[[1]], tableNP2 = tableN[[2]], tableNP3 = tableN[[3]], 
              tableNP4 = tableN[[4]], tableNP5 = tableN[[5]], tableNP6 = p1summ$s.table, 
              tableNP7 = p2summ$s.table, tableNP8 = tableN[[8]], Model = object$Model, 
              n = n, theta = object$theta, theta.a = object$theta.a, 
              dof = object$dof, dof.a = object$dof.a, sigma21 = object$sigma21, 
              sigma22 = object$sigma22, sigma1 = object$sigma21, sigma2 = object$sigma22, 
              nu1 = object$nu1, nu2 = object$nu2, tau = object$tau, 
              sigma21.a = object$sigma21.a, sigma22.a = object$sigma22.a, 
              sigma1.a = object$sigma21.a, sigma2.a = object$sigma22.a, 
              nu1.a = object$nu1.a, nu2.a = object$nu2.a, tau.a = object$tau.a, 
              formula = object$formula, formula1 = object$gam1$formula, 
              formula2 = object$gam2$formula, formula3 = object$gam3$formula, 
              formula4 = object$gam4$formula, formula5 = object$gam5$formula, 
              formula6 = object$p1gam$formula, formula7 = object$p2gam$formula, 
              formula8 = object$gam8$formula, t.edf = object$t.edf, 
              CItheta = CIrs, CIsig1 = CIsig21, CIsig2 = CIsig22, 
              CInu1 = CInu1, CInu2 = CInu2, CItau = CIkt, CIdof = CIdof, 
              BivD = object$BivD, margins = object$margins, l.sp1 = object$l.sp1, 
              l.sp2 = object$l.sp2, l.sp3 = object$l.sp3, l.sp4 = object$l.sp4, 
              l.sp5 = object$l.sp5, l.sp6 = length(object$p1gam$sp), l.sp7 = length(object$p2gam$sp), 
              l.sp8 = object$l.sp8, X3.null = is.null(object$X3), 
              univar.gamlss = FALSE, m2 = object$VC$m2, m3 = object$VC$m3, 
              surv = object$surv, surv.flex = object$surv.flex, K1 = NULL)
  class(res) <- "summary.zinfgammacop"
  res
}





print.summary.zinfgammacop <- 
function (x, digits = max(3, getOption("digits") - 3), signif.stars = getOption("show.signif.stars"), 
          ...) 
{
  ppR <- pp(x)
  cont1par <- ppR$cont1par
  cont2par <- ppR$cont2par
  cont3par <- ppR$cont3par
  cop <- ppR$cop
  lind <- ppR$lind
  m1l <- ppR$m1l
  m2l <- ppR$m2l
  main.t <- "\nCOPULA:  "
  cp <- "  theta = "
  as.p <- x$theta.a
  dof <- x$dof.a
  ct <- "  tau = "
  kt.p <- x$tau.a
  s1 <- "sigma.1 = "
  s1.p <- x$sigma21.a
  s2 <- "sigma.2 = "
  s2.p <- x$sigma22.a
  n1 <- "nu.1 = "
  n1.p <- x$nu1.a
  n2 <- "nu.2 = "
  n2.p <- x$nu2.a
  cat(main.t, cop)
  pscr0(x)
  pscr(x, lind, m1l, m2l, cont1par, cont2par, cont3par, type = "copR", 
       digits, signif.stars, ...)
  CIrs <- colMeans(x$CItheta, na.rm = TRUE)
  CIkt <- colMeans(x$CItau, na.rm = TRUE)
  if (x$margins[1] %in% c(cont2par, cont3par)) 
    CIsig21 <- colMeans(x$CIsig1, na.rm = TRUE)
  if (x$margins[2] %in% c(cont2par, cont3par)) 
    CIsig22 <- colMeans(x$CIsig2, na.rm = TRUE)
  if (x$margins[1] %in% cont3par) 
    CInu1 <- colMeans(x$CInu1, na.rm = TRUE)
  if (x$margins[2] %in% cont3par) 
    CInu2 <- colMeans(x$CInu2, na.rm = TRUE)
  BivD <- x$BivD
  if (BivD == "T" && x$margins[1] %in% c(x$m2, x$m3) && x$margins[2] %in% 
      c(x$m2, x$m3)) 
    CIdof <- colMeans(x$CIdof, na.rm = TRUE)
  nodi <- 3
    if (x$margins[1] %in% cont2par && x$margins[2] %in% 
        cont2par) 
      cat(s1, format(s1.p, digits = nodi), "(", format(CIsig21[1], 
                                                       digits = nodi), ",", format(CIsig21[2], digits = nodi), 
          ")", "  ", s2, format(s2.p, digits = nodi), 
          "(", format(CIsig22[1], digits = nodi), ",", 
          format(CIsig22[2], digits = nodi), ")", "\ntheta = ", 
          format(as.p, digits = nodi), "(", format(CIrs[1], 
                                                   digits = nodi), ",", format(CIrs[2], digits = nodi), 
          ")", ct, format(kt.p, digits = nodi), "(", format(CIkt[1], 
                                                            digits = nodi), ",", format(CIkt[2], digits = nodi), 
          ")", "\nn = ", x$n, "  total edf = ", format(x$t.edf, 
                                                       digits = nodi), "\n\n", sep = "")
 invisible(x)
}








pscr <- 
function (x, lind, m1l, m2l, cont1par, cont2par, cont3par, type = "copR", 
          digits, signif.stars, m3l = NULL, lind2 = NULL, ...) 
{
      doff <- "log(· - 2)"
      cat("\n\nEQUATION 1")
      if (x$surv.flex == FALSE || (x$surv.flex == TRUE && 
                                   x$margins[1] %in% c(x$m2, x$m3))) 
        cat("\nLink function for mu.1:", m1l, "\n")
      if (x$surv.flex == TRUE && !(x$surv.flex == TRUE && 
                                   x$margins[1] %in% c(x$m2, x$m3))) 
        cat("\n")
      cat("Formula: ")
      print(x$formula[[1]])
      cat("\n")
      cat("Parametric coefficients:\n")
      printCoefmat(x$tableP1, digits = digits, signif.stars = signif.stars, 
                   na.print = "NA", ...)
      cat("\n")
      if (x$l.sp1 != 0) {
        cat("Smooth components' approximate significance:\n")
        printCoefmat(x$tableNP1, digits = digits, signif.stars = signif.stars, 
                     has.Pvalue = TRUE, na.print = "NA", cs.ind = 1, 
                     ...)
        cat("\n")
      }
      cat("\nEQUATION 2")
      if (x$surv.flex == FALSE) 
        cat("\nLink function for mu.2:", m2l, "\n")
      if (x$surv.flex == TRUE) 
        cat("\n")
      cat("Formula: ")
      print(x$formula[[2]])
      cat("\n")
      cat("Parametric coefficients:\n")
      printCoefmat(x$tableP2, digits = digits, signif.stars = signif.stars, 
                   na.print = "NA", ...)
      cat("\n")
      if (x$l.sp2 != 0) {
        cat("Smooth components' approximate significance:\n")
        printCoefmat(x$tableNP2, digits = digits, signif.stars = signif.stars, 
                     has.Pvalue = TRUE, na.print = "NA", cs.ind = 1, 
                     ...)
        cat("\n")
      }

            cat("\nEQUATION 3")
            if (x$margins[1] != "BE") 
              cat("\nLink function for sigma.1:", "log", 
                  "\n")
            else cat("\nLink function for sigma.1:", 
                     "qlogis", "\n")
            cat("Formula: ")
            print(x$formula[[3]])
            cat("\n")
            cat("Parametric coefficients:\n")
            printCoefmat(x$tableP3, digits = digits, 
                         signif.stars = signif.stars, na.print = "NA", 
                         ...)
            cat("\n")
            if (x$l.sp3 != 0) {
              cat("Smooth components' approximate significance:\n")
              printCoefmat(x$tableNP3, digits = digits, 
                           signif.stars = signif.stars, has.Pvalue = TRUE, 
                           na.print = "NA", cs.ind = 1, ...)
              cat("\n")
            }
            cat("\nEQUATION 4")
            if (x$margins[2] != "BE") 
              cat("\nLink function for sigma.2:", "log", 
                  "\n")
            else cat("\nLink function for sigma.2:", 
                     "qlogis", "\n")
            cat("Formula: ")
            print(x$formula[[4]])
            cat("\n")
            cat("Parametric coefficients:\n")
            printCoefmat(x$tableP4, digits = digits, 
                         signif.stars = signif.stars, na.print = "NA", 
                         ...)
            cat("\n")
            if (x$l.sp4 != 0) {
              cat("Smooth components' approximate significance:\n")
              printCoefmat(x$tableNP4, digits = digits, 
                           signif.stars = signif.stars, has.Pvalue = TRUE, 
                           na.print = "NA", cs.ind = 1, ...)
              cat("\n")
            }
            cat("\nEQUATION 5")
            cat("\nLink function for theta:", lind, 
                "\n")
            cat("Formula: ")
            print(x$formula[[5]])
            cat("\n")
            cat("Parametric coefficients:\n")
            printCoefmat(x$tableP5, digits = digits, 
                         signif.stars = signif.stars, na.print = "NA", 
                         ...)
            cat("\n")
            if (x$l.sp5 != 0) {
              cat("Smooth components' approximate significance:\n")
              printCoefmat(x$tableNP5, digits = digits, 
                           signif.stars = signif.stars, has.Pvalue = TRUE, 
                           na.print = "NA", cs.ind = 1, ...)
              cat("\n")
            }
            cat("\nEQUATION 6")
            cat("\nLink function for kappa: logit \n")
            cat("Formula: ")
            print(x$formula[[6]])
            cat("\n")
            cat("Parametric coefficients:\n")
            printCoefmat(x$tableP6, digits = digits, 
                         signif.stars = signif.stars, na.print = "NA", 
                         ...)
            cat("\n")
            if (x$l.sp6 != 0) {
              cat("Smooth components' approximate significance:\n")
              printCoefmat(x$tableNP6, digits = digits, 
                           signif.stars = signif.stars, has.Pvalue = TRUE, 
                           na.print = "NA", cs.ind = 1, ...)
              cat("\n")
            }
            cat("\nEQUATION 7")
            cat("\nLink function for kappa: logit \n")
            cat("Formula: ")
            print(x$formula[[7]])
            cat("\n")
            cat("Parametric coefficients:\n")
            printCoefmat(x$tableP7, digits = digits, 
                         signif.stars = signif.stars, na.print = "NA", 
                         ...)
            cat("\n")
            if (x$l.sp7 != 0) {
              cat("Smooth components' approximate significance:\n")
              printCoefmat(x$tableNP7, digits = digits, 
                           signif.stars = signif.stars, has.Pvalue = TRUE, 
                           na.print = "NA", cs.ind = 1, ...)
              cat("\n")
            }
}











tmpfun <- get("form.eq12", envir = asNamespace("GJRM"))

environment(form.eq12) <- environment(tmpfun)

assignInNamespace("form.eq12", form.eq12, ns = "GJRM")



tmpfun <- get("post.check", envir = asNamespace("GJRM"))

environment(post.check) <- environment(tmpfun)

assignInNamespace("post.check", post.check, ns = "GJRM")



tmpfun <- get("gjrm", envir = asNamespace("GJRM"))

environment(predict.zinfgammacop) <- environment(tmpfun)



tmpfun <- get("pscr", envir = asNamespace("GJRM"))

environment(pscr) <- environment(tmpfun)

assignInNamespace("pscr", pscr, ns = "GJRM")



tmpfun <- get("gjrm", envir = asNamespace("GJRM"))

environment(summary.zinfgammacop) <- environment(tmpfun)



tmpfun <- get("gjrm", envir = asNamespace("GJRM"))

environment(zinfgammacop) <- environment(tmpfun)


tmpfun <- get("gjrm", envir = asNamespace("GJRM"))

environment(zinfgammacop.override) <- environment(tmpfun)



