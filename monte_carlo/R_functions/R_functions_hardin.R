## Downloaded and adapted from:
## http://www.biomedcentral.com/content/supplementary/1471-2105-8-220-S1.r

#######################################################
#### R code for computing the Biweight Correlation ####
#### J. Hardin, A. Mitani, L. Hicks, B. VanKoten   ####
#### last updated 6/05/2007                        ####
#######################################################

# #################
# ##   EXAMPLE   ##
# #################
#
# library(rrcov)
# library(MASS)

# samp.data <- mvrnorm(30,mu=c(0,0),Sigma=matrix(c(1,.75,.75,1),ncol=2)) # you need MASS
#
# r<-0.2 # breakdown
# n1<-30 # number of samples
#
# c1<-rejpt.bw(p=2,r)[1]
# b0<-erho.bw(p=2,c1)[1]
#
# samp.mcd <- covMcd(samp.data) #you need rrcov
# samp.bw <- biwt.est(samp.data,n1,p=2,r,c1,b0,samp.mcd)
# samp.bw.corr <- samp.bw$biwt.sig[1,2] / sqrt(samp.bw$biwt.sig[1,1]*samp.bw$biwt.sig[2,2])
#
# ##############
# # to speed up the calculations, use the median/mad for the initialization:
# ##############
#
# samp.init <- list()
# 	samp.init$cov <- diag(apply(samp.data,2,mad,na.rm=T))
# 	samp.init$center <- apply(samp.data,2,median,na.rm=T)
# samp.bw <- biwt.est(samp.data,n1,p=2,r,c1,b0,samp.init)
# samp.bw.corr <- samp.bw$biwt.sig[1,2] / sqrt(samp.bw$biwt.sig[1,1]*samp.bw$biwt.sig[2,2])


##################
## The Function ##
##################

# this function calculates the biweight mean and covariance in any dimension
# given the appropriate c1 & b0
# we will always use p=2 when looking for a standard correlation

biwt.est <- function(x,n,p=2,r,c1,b0,med.init){

	NAid<-FALSE
	d <- sqrt(mahalanobis(x,med.init$center,med.init$cov))
	k <- ksolve(d,p,c1,b0)


	## !! Changed here !!
	## if is.na(k) solution in the equation (3) in the paper is not avilable
	## and the algorithm looks for another initial based on the MCD
	## I remove this because we never want to get into the MCD burden
# 	if(is.na(k)) {
# 		NAid<-TRUE
# 		med.init <- covMcd(x)
# 		d <- sqrt(mahalanobis(x,med.init$center,med.init$cov))
# 		k <- ksolve(d,p,c1,b0)}
#    # MCD is a more robust estimate of the center/shape
#    # than the median which is sometimes used


	eps <- 1e-5
	crit <- 100
	iter <- 1
	while (crit > eps & iter < 100) {
		d <- d/k
		biwt.mu <- apply(wtbw(d,c1)*x,2,sum,na.rm=T) / sum (wtbw(d,c1),na.rm=T)
		cent <- array(dim=c(n,p,p))
		for (i in 1:n){
		cent[i,,] <- (x[i,] - biwt.mu)%*%t(x[i,]-biwt.mu)}
		biwt.sig <- apply(cent*wtbw(d,c1),c(2,3),sum,na.rm=T)/
				sum(vbw(d,c1),na.rm=T)


		d2 <- sqrt(mahalanobis(x,biwt.mu,biwt.sig))
		k <- ksolve(d2,p,c1,b0)
		crit <- max(abs(d-(d2/k)),na.rm=T)
		d <- d2
		iter <-  iter+1}

return(list(biwt.mu=biwt.mu,biwt.sig=biwt.sig))}



##########################################################
#### Functions needed to run the biweight correlation ####
##########################################################

ksolve <- function(d,p,c1,b0){
    k <- 1
    iter <- 1
    crit <- 100
    eps <- 1e-10
    while ((crit > eps)&(iter<100)){
    k.old <- k
        fk <- mean(rhobw(d/k,c1),na.rm=T)-b0
        fkp <- -mean(psibw(d/k,c1)*d/k^2,na.rm=T)
		if (fkp==0) {k<-NA
		return(k)
		stop("no values close enough")}
        k <- k - fk/fkp
        if (k < 0)  k <- k.old/2
        crit <- abs(k-k.old)
        iter <- iter+1    }
    return(k) }

rhobw <- function(x,c1){
	ivec <- (abs(x)>c1)
	return((c1^2/6)*ivec +(1-ivec)*(x^2/2-x^4/(2*c1^2)+x^6/(6*c1^4)))}

psibw <- function(x,c1){
	ivec <- (abs(x)>c1)
    	return((1-ivec)*(x*(1-(x/c1)^2)^2))}

wtbw <- function(x,c1){
    	ivec <- (abs(x)>c1)
    	return((1-ivec)*(1-(x/c1)^2)^2)}

vbw <- function(x,c1) return(psibw(x,c1)*x)


erho.bw <- function(p,c1) # gives b0 = E(rho)
	return(chi.int(p,2,c1)/2-chi.int(p,4,c1)/(2*c1^2)+
     		chi.int(p,6,c1)/(6*c1^4)+c1^2*chi.int2(p,0,c1)/6)

erho.bw.p <- function(p,c1)
	return(chi.int.p(p,2,c1)/2-chi.int.p(p,4,c1)/(2*c1^2)+					2*chi.int(p,4,c1)/(2*c1^3)+chi.int.p(p,6,c1)/(6*c1^4)-
		4*chi.int(p,6,c1)/(6*c1^5)+c1^2*chi.int2.p(p,0,c1)/6
    		+2*c1*chi.int2(p,0,c1)/6)


chi.int <- function(p,a,c1)
	return( exp(lgamma((p+a)/2)-lgamma(p/2))*2^{a/2}*pchisq(c1^2,p+a) )
chi.int2 <- function(p,a,c1)
	return( exp(lgamma((p+a)/2)-lgamma(p/2))*2^{a/2}*(1-pchisq(c1^2,p+a)) )
chi.int.p <- function(p,a,c1)
	return( exp(lgamma((p+a)/2)-lgamma(p/2))*2^{a/2}*dchisq(c1^2,p+a)*2*c1 )
chi.int2.p <- function(p,a,c1)
	return( -exp(lgamma((p+a)/2)-lgamma(p/2))*2^{a/2}*dchisq(c1^2,p+a)*2*c1 )


rejpt.bw <- function(p,r){ # gives c1 = ARP
    c1 <- 2*p
    iter <- 1
    crit <- 100
    eps <- 1e-5
    while ((crit > eps)&(iter<100)){
	c1.old <- c1
        fc <- erho.bw(p,c1) - c1^2*r/6
        fcp <- erho.bw.p(p,c1) - c1*r/3
        c1 <- c1 - fc/fcp
        if (c1 < 0)  c1 <- c1.old/2
        crit <- abs(fc)
        iter <- iter+1    }
    return(c(c1,pchisq(c1^2,p),log10(1-pchisq(c1^2,p))))}





