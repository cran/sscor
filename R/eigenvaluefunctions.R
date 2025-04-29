shape <- function(shapema) {
p <- length(shapema[,1])
ablma <- diag(rep(1,p^2))
index <- seq(from=1,to=p^2,by=p+1)
for(i in 1:p^2) {
ablma[i,i] <- 1/shapema[floor((i-1)/p)+1,floor((i-1)/p)+1]^(1/2)/shapema[((i-1)%%p)+1,((i-1)%%p)+1]^(1/2)
}
for(i in 1:p^2) {
ablma[i,index[floor((i-1)/p)+1]] <- -1/2/shapema[floor((i-1)/p)+1,floor((i-1)/p)+1]^(3/2)/shapema[((i-1)%%p)+1,((i-1)%%p)+1]^(1/2)*shapema[floor((i-1)/p)+1,((i-1)%%p)+1]
}
for(i in 1:p^2) {
ablma[i,index[((i-1)%%p)+1]] <- -1/2/shapema[((i-1)%%p)+1,((i-1)%%p)+1]^(3/2)/shapema[floor((i-1)/p)+1,floor((i-1)/p)+1]^(1/2)*shapema[floor((i-1)/p)+1,((i-1)%%p)+1]
}
for (i in 1:p) {
ablma[index[i],index[i]] <- 0
}
return(ablma)
}

varicor <- function(sigma) {

varma <- varishape(sigma)

shapeab <- shape(sigma/eigen(sigma)$values[1])

return(shapeab%*%varma%*%t(shapeab))
}


varishape <- function(sigma) {
A <- eigenveccov(sigma)
B <- eigenvaluecov(sigma)
p1 <- length(A[1,])
p2 <- length(B[1,])
B <- B[2:p2,2:p2]
p2 <- p2-1
covi <- matrix(0,ncol=p1+p2,nrow=p1+p2)
covi[1:p1,1:p1] <- A
covi[(p1+1):(p1+p2),(p1+1):(p1+p2)] <- B
ableitung <- Ableitungalles(sigma)
return(ableitung%*%covi%*%t(ableitung))
}


Ableitungalles <- function(sigma) {
p <- length(sigma[,1])
evd <- eigen(sigma)
U <- evd$vectors
lamstan <- evd$values
lamstan <- lamstan/lamstan[1]
lamma <- diag(lamstan)
f1 <- t(kronecker(lamma%*%t(U),diag(rep(1,p))))
f2 <- cbind(diag(rep(1,p^2)),matrix(0,ncol=p-1,nrow=p^2))


S1 <- f1%*%f2
f1 <- kronecker(diag(rep(1,p)),U)
f21 <- t(kronecker(t(U),diag(rep(1,p))))
f22 <- matrix(0,nrow=p^2,ncol=p^2+p-1)
mama <- Ableitung(sigma)
index <- seq(from=1,to=p^2,by=p+1)
for (i in 2:p) {
f22[index[i],(p^2+1):(p^2+p-1)] <- mama[(i-1),]
}
S21 <- f21%*%f22
f31 <- kronecker(diag(rep(1,p)),lamma)
f32 <- matrix(0,ncol=p^2+p-1,nrow=p^2)
index <- as.numeric(t(outer((1:p),0:(p-1)*p,"+")))
for (i in  1:(p^2)) {
f32[i,index[i]] <- 1
}
S22 <- f31%*%f32
abl <- S1+f1%*%(S21+S22)
return(abl)
}



Ableitung <- function(sigma) {
p <- length(sigma[,1])
evd <- eigen(sigma)
werte <- evd$values/evd$values[1]
mama <- fuint4(werte)
mamainv <- solve(mama[2:p,2:p])
return(mamainv)
}

perm <- function(dime) {
index <- outer((1:dime),seq(from=0,to=(dime-1)*dime,by=dime),"+")
index <- as.numeric(t(index))
perma <- matrix(0,ncol=dime^2,nrow=dime^2)
for (i in 1:dime^2) {
perma[i,index[i]] <- 1
}
return(perma)
}

Ableitungalles2 <- function(sigma) {
p <- length(sigma[,1])
evd <- eigen(sigma)
U <- evd$vectors
lamstan <- evd$values
lamstan <- lamstan/lamstan[1]
lamma <- diag(lamstan)
f1 <- kronecker(U%*%lamma,diag(rep(1,p)))
f2 <- kronecker(diag(rep(1,p)),U%*%lamma)%*%perm(p)
ma1 <- f1+f2
index <- seq(from=1,to=p^2,by=p+1)
index <- index[-1]
ablma <- Ableitung(sigma)
abl <- matrix(0,nrow=p^2,ncol=p-1)
for (i in 1:p-1) {
abl[index[i],] <- ablma[i,]
}
ma2 <- kronecker(U,U)%*%abl
return(cbind(ma1,ma2))
}



Ableitung <- function(sigma) {
p <- length(sigma[,1])
evd <- eigen(sigma)
werte <- evd$values/evd$values[1]
mama <- fuint4(werte)
mamainv <- solve(mama[2:p,2:p])
return(mamainv)
}


eigenveccov <- function(sigma) {
p <- length(sigma[,1])
index <- seq(from=1,to=p^2,by=p+1)
evd <- eigen(sigma)
evsscm <- evShape2evSSCM(evd$values)
Wma <- asycov(diag(evd$values))
denominator <- outer(evsscm,evsscm,"-")
denominator <- as.numeric(denominator)%*%t(as.numeric(denominator))
Wma <- Wma/denominator
index <- seq(from=1,to=p^2,by=p+1)
Wma[index,] <- 0
Wma[,index] <- 0
UU <- kronecker(diag(1,p),evd$vectors)
covma <- UU%*%Wma%*%t(UU)
return(covma)
}


eigenvaluecov <- function(sigma) {
p <- length(sigma[,1])
index <- seq(from=1,to=p^2,by=p+1)
asycov(diag(eigen(sigma)$values))[index,index]
}


asycov <- function(sigma) {
p <- length(sigma[1,])
evd <- eigen(sigma)
Eta <- matrix(0,nrow=p^2,ncol=p^2)
etav <- eta(evd$values)
for(i in 1:p) {
    for(j in 1:p) {
        if (i != j) {
            Eta[(i-1)*p+j,(i-1)*p+j] <- etav[i,j]
            Eta[(i-1)*p+i,(j-1)*p+j] <- etav[i,j]
            Eta[(i-1)*p+j,(j-1)*p+i] <- etav[i,j]
        } else {
            Eta[(i-1)*p+i,(i-1)*p+i] <- etav[i,i]
        }
    } 

}

delta <- evShape2evSSCM(evd$values)
UU <- kronecker(evd$vectors,evd$vectors)
Etaz <- Eta - as.numeric(diag(delta))%*%t(as.numeric(diag(delta)))
Ws <- UU%*%Etaz%*%t(UU)
return(Ws)
}

fuwo <- function(sigma) {
eigen(sigma)

}






## Calculates eigenvalues of SSCM;
# input:
# evShape: eigenvalues of shape matrix
# output: eigenvalues of SSCM

eta <- function(evShape) {
  evShape <- evShape/sum(evShape)*length(evShape)
  eta <- fuint3(evShape)*(evShape%*%t(evShape))/4
  diag(eta) <- 3*diag(eta)
  return(eta)
}


## Calculates eigenvalues of SSCM;
# input:
# evShape: eigenvalues of shape matrix
# output: eigenvalues of SSCM

evShape2evSSCM <- function(evShape) {
  evShape <- evShape/sum(evShape)*length(evShape)
  delta <- fuint(evShape)*evShape/2
  return(delta)
}



## Calculates the eigenvalues of the shape matrix given the eigenvalues of the SSCM
# input: eigenvalues of the SSCM
# tol: absolute tolerance for the approximated eigenvalues
# itermax: maximal number of iterations of the approximation algorithm 
# output: list containing:
# eigenv: eigenvalues of the shape matrix
# iternumber: number of iterations of the approximation algorithm


evSSCM2evShape <- function(delta,tol=10^(-10),itermax=100) {
pm0 <- length(delta)
lambda <- numeric(pm0)
Index <- delta==0
p0 <- sum(Index)
delta <- delta[!Index]
p <- pm0-p0
if (p==2) {
lambda[1] <- 1
lambda[2] <- ((1-delta[1])/(delta[1]))^2
return(lambda)
}

pdata <- p
pdata <- min(pdata,200)
jacobiquad <- get(load(system.file("extdata", "jacobiquad", package = "sscor")))	
kno <- jacobiquad[[1]][pdata,]
wei <- jacobiquad[[2]][pdata,]

lambdaalt <- lambdaneu <- delta
lambdaalt[1] <- lambdaalt[1]+1

for(i in 1:itermax) {
lambdaalt <- lambdaneu
integr <- fuint(lambdaalt,wei=wei,kno=kno)
lambdaneu <- 2*delta/integr
lambdaneu <- lambdaneu/sum(lambdaneu)*p
deltaapprox <- lambdaalt*integr/2
if(sum(abs(deltaapprox-delta))<tol) break
if(i==itermax) warning("Maximal number of iterations reached. Approximated shape matrix might be imprecise.")
}
lambda[!Index] <- lambdaneu/sum(lambdaneu)
return(lambda)
}

## Evaluates the integral which gives the eigenvalues of the SSCM;
# input:
# evShape: eigenvalues of shape matrix
# wei: weights of Gauss-Jacobi quadrature
# kno: knots of Gauss-Jacobi quadrature
# output: defined integral


fuint <- function(lambda,wei=NULL,kno=NULL) {
pdata <- length(lambda)
pdata <- min(pdata,200)
p <- pdata/2+1
if(is.null(wei)|is.null(kno)) {
jacobiquad <- get(load(system.file("extdata", "jacobiquad", package = "sscor")))
kno <- jacobiquad[[1]][pdata,]
wei <- jacobiquad[[2]][pdata,]
}
wei <- wei*(1-kno)^(-p)
kno <- (1+kno)/(1-kno)
foo <- function(x){ 1/(sqrt(apply(1+outer(lambda,x),2,prod)))}
f1 <- foo(kno)
fd <- 1/(1+outer(lambda,kno))
w <- f1*t(fd)
bi <- 2*apply(wei*w,2,sum)
return(bi)
}


## Evaluates the integral which gives the eigenvalues of the SSCM;
# input:
# evShape: eigenvalues of shape matrix
# wei: weights of Gauss-Jacobi quadrature
# kno: knots of Gauss-Jacobi quadrature
# output: defined integral


fuint2 <- function(lambda,wei=NULL,kno=NULL) {
pdata <- length(lambda)
pdata <- min(pdata,200)
p <- pdata/2+1
if(is.null(wei)|is.null(kno)) {
jacobiquad <- get(load(system.file("extdata", "jacobiquad", package = "sscor")))
kno <- jacobiquad[[1]][pdata,]
wei <- jacobiquad[[2]][pdata,]
}
wei <- wei*(1-kno)^(-p)
kno <- (1+kno)/(1-kno)
foo <- function(x){ x/(sqrt(apply(1+outer(lambda,x),2,prod)))}
f1 <- foo(kno)
fd <- 1/(1+outer(lambda,kno))^2
w <- f1*t(fd)
bi <- 2*apply(wei*w,2,sum)
return(bi)
}


fuint3 <- function(lambda,wei=NULL,kno=NULL) {
pdata <- length(lambda)
pdata <- min(pdata,200)
p <- pdata/2+1
if(is.null(wei)|is.null(kno)) {
jacobiquad <- get(load(system.file("extdata", "jacobiquad", package = "sscor")))
kno <- jacobiquad[[1]][pdata,]
wei <- jacobiquad[[2]][pdata,]
}
wei <- wei*(1-kno)^(-p)
kno <- (1+kno)/(1-kno)
foo <- function(x){ x/(sqrt(apply(1+outer(lambda,x),2,prod)))}
f1 <- foo(kno)
fd <- 1/(1+outer(lambda,kno))
fd <- apply(fd,2,function(x) return(x%*%t(x)))
w <- f1*t(fd)
bi <- 2*apply(wei*w,2,sum)
return(matrix(bi,ncol=pdata,nrow=pdata))
}



fuint4 <- function(lambda,wei=NULL,kno=NULL) {
pdata <- length(lambda)
pdata <- min(pdata,200)
p <- pdata/2+1
if(is.null(wei)|is.null(kno)) {
jacobiquad <- get(load(system.file("extdata", "jacobiquad", package = "sscor")))
kno <- jacobiquad[[1]][pdata,]
wei <- jacobiquad[[2]][pdata,]
}
wei <- wei*(1-kno)^(-p)
kno <- (1+kno)/(1-kno)

foo <- function(x){ 1/(sqrt(apply(1+outer(lambda,x),2,prod)))}
f1 <- foo(kno)
fd <- 1/(1+outer(lambda,kno))
w <- f1*t(fd)
bi1 <- apply(wei*w,2,sum)

foo <- function(x){ x/(sqrt(apply(1+outer(lambda,x),2,prod)))}
f1 <- foo(kno)
fd <- 1/(1+outer(lambda,kno))
fd <- apply(fd,2,function(x) return(x%*%t(x)))
w <- f1*t(fd)
bi2 <- matrix(2*apply(wei*w,2,sum),ncol=pdata)
gew <- matrix(-1/2,ncol=pdata,nrow=pdata)
gew <- gew - diag(rep(1,pdata))
gew <- gew * lambda/2
bi2 <- bi2*gew

return(diag(bi1)+bi2)

}



