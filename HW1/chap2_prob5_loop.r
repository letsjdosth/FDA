# R file

# chap 2 problem 2.5.5

library(fields)
library(expm)
library(fda)

set.seed(2019311252)

num.sample = 50

#generate bump function
bump = function(x, c0, r0, a0){
    result = rep(0, length(x))
    for(i in 1:length(x)){
        if(abs(x[i] - c0) < r0){
            result[i] = a0*exp(-(1-((x[i]-c0)/r0)^2)^-1)
        }
    }
    return(result)
}

times = seq(0,0.99,0.01)
bumpval1 = bump(times, 3/8, 1/4, 5)
bumpval2 = bump(times, 5/8, 1/4, 5)

#generate Matern process
mt.sig2 = 1 #pointwise variance
mt.rho = 0.1 #range parameter (how quickly falls off)
mt.nu = 1/2 #smoothness parameter

loop.varpropvec = list()
loop.coeffvec = list()

for(i in 1:10){
    cat("iter: ",i,"\n")
    d_mat = abs(outer(times,times,"-"))
    C_1 = apply(d_mat, c(1,2), FUN=Matern, range=mt.rho, nu=mt.nu)
    C_1 = C_1*mt.sig2
    C_1_sq = sqrtm(C_1) #<-오래걸림

    MaternSamples = matrix(0, length(times), num.sample)
    for(i in 1:num.sample){
        Z<-rnorm(length(times))
        MaternSamples[,i] = C_1_sq%*%Z
    }

    MtSample.group1 = MaternSamples[, 1:(num.sample/2)]
    MtSample.group2 = MaternSamples[, (num.sample/2+1):num.sample]

    # mean adjust
    for(i in 1:(num.sample/2)){
        MtSample.group1[,i] = MtSample.group1[,i] + bumpval1
        MtSample.group2[,i] = MtSample.group2[,i] + bumpval2
    }
    adjMtSamples = cbind(MtSample.group1, MtSample.group2)

    # basic fitting to b-spline
    b_spline_basis = create.bspline.basis(c(0,1), nbasis=6)
    adjMtSample.f = smooth.basis(times, adjMtSamples, b_spline_basis)$fd #굳이 par로 smoothing?

    #apply continuous registration
    adjMtSample.f.reg = register.fd(adjMtSample.f) #<-오래걸림

    #pca
    adjMtSample.f.pca = pca.fd(adjMtSample.f, nharm=1)
    adjMtSample.f.reg.pca = pca.fd(adjMtSample.f.reg$regfd, nharm=1)

    
    loop.varpropvec = rbind(loop.varpropvec, c(adjMtSample.f.pca$varprop, adjMtSample.f.reg.pca$varprop))


    #regression (뭘 하라는거야??)
    dummy_var = c(rep(0,num.sample/2),rep(1,num.sample/2))
    un_align_score = adjMtSample.f.pca$score
    align_score = adjMtSample.f.reg.pca$score

    lm1= lm(un_align_score ~ dummy_var)
    lm2= lm(align_score ~ dummy_var)
    loop.coeffvec= rbind(loop.coeffvec, c(lm1$coeff[2], lm1$coeff[2]))
}