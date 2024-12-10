
source('zinfgammacop_functions.R')

n <- 5000
set.seed(1)
x1 <- runif(n,0,27) 
x2 <- rbinom(n=n, size=1, prob=0.55)

mu1 <- (1-x2)*(0.008*x1^2 -0.2*x1 +2.5) + x2*(exp(0.08*exp(0.1*x1)+0.2))
mu2 <- (1-x2)*(exp(sin(0.1*(x1-16)))) +  x2*(exp(-0.02*x1+0.6))
sigma1 <- 0.1*x2+0.1
sigma2 <- 0.1*x2+0.2
rho <- (1-x2)*(tanh(0.4*(0.2*x1-1))) +  x2*(tanh(-0.01*x1+0.2))
p1 <- (1-x2)*(sigmoid(-0.004*x1^2+0.037*x1-0.2)) +  x2*(sigmoid(0.01*x1-0.9))
p2 <- (1-x2)*(sigmoid(-0.05*x1+0.3)) +  x2*(sigmoid(0.01*x1^2-0.3*x1-0.1))

y1y2.sim <- zinfgammacop.sim(mu1=mu1,mu2=mu2,sig1=sigma1,sig2=sigma2,rho=rho,p1=p1,p2=p2)
x2 <- as.factor(x2)
simdata <- data.frame(y1=y1y2.sim[,1], y2=y1y2.sim[,2], x1=x1, x2=x2)

eq1 <- y1 ~ x2 + s(x1,by=x2) 
eq2 <- y2 ~ x2 + s(x1,by=x2) 
eq3 <- ~ x2
eq4 <- ~ x2
eq5 <- ~ x2 + s(x1,by=x2) 
eq6 <- ~ x2 + s(x1,by=x2)
eq7 <- ~ x2 + s(x1,by=x2)
eqlist <- list(eq1, eq2, eq3, eq4, eq5, eq6, eq7)

obj <- zinfgammacop(formula=eqlist, data=simdata, copula="N", model="B", margins=c("GA","GA"), gamlssfit=FALSE)
