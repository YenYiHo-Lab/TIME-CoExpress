
source('TIME-CoExpress.R')

eq1 <- y1 ~ group + s(t)+s(t, by=ogroup)
eq2 <- y2 ~ group + s(t)+s(t, by=ogroup)
eq3 <- ~ 1
eq4 <- ~ 1
eq5 <- ~ group + s(t)+s(t, by=ogroup)
eqlist <- list(eq1, eq2, eq3, eq4, eq5)

obj <- TIME_CoExpress(formula=eqlist, data=dat, copula="N", model="B", margins=c("GA","GA"), gamlssfit=TRUE, rinit=50, rmax=10000)

pvalue <- summary(obj)$tableNP5[8]