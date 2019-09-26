# R file

# chap 2 problem 2.5.5

library(fields)
library(expm)
library(fda)

set.seed(2019311252)

num.sample = 50

# (a)
# Simulate a functional sample over the unit interval each with a sample size of 50 from the Mat´ern process. 
# For the ﬁrst half of the sample, set the mean function equal to the the bump function with parameters (c0,r0,a0) = (3/8,1/4,5). 
# For the second half use (c0,r0,a0) = (5/8,1/4,5). 
# You may choose the values for the Mat´ern covariance function as well as the number of points sampled per curve. 
# Plot all of the curves and include a curve for the overall mean function.

#define bump function
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
# plot(times, bumpval2)

#generate Matern process
mt.sig2 = 1 #pointwise variance
mt.rho = 0.1 #range parameter (how quickly falls off)
mt.nu = 1/2 #smoothness parameter

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

#plot
#original matern samples
par(mfrow=c(1,1))
plot(times,MaternSamples[,1],type="l", ylim=c(-5,5)) 
for(i in 2:num.sample)
    points(times,MaternSamples[,i],type="l")
}

# mean adjust
for(i in 1:(num.sample/2)){
    MtSample.group1[,i] = MtSample.group1[,i] + bumpval1
    MtSample.group2[,i] = MtSample.group2[,i] + bumpval2
}
adjMtSamples = cbind(MtSample.group1, MtSample.group2)


#plot after adjusted by bump function
par(mfrow=c(1,3))
plot(times,MtSample.group1[,1],type="n", ylim=c(-3,4), ylab="value")
for(i in 1:(num.sample/2)){
    points(times, MtSample.group1[,i],type="l", col="red")
    points(times, MtSample.group2[,i], type="l", col="blue")
}
points(times, bumpval1, type="l", col="red", lwd=3)
points(times, bumpval2, type="l", col="blue", lwd=3)
points(times, rowMeans(adjMtSamples[,1:(num.sample/2)]), type="l", col="yellow", lwd=3)
points(times, rowMeans(adjMtSamples[,(num.sample/2+1):num.sample]), type="l", col="green", lwd=3)
#
plot(times, bumpval1, type="l", col="red", lwd=3, ylim=c(-3,4), ylab="bump function value & group sample mean value")
points(times, bumpval2, type="l", col="blue", lwd=3)
points(times, rowMeans(adjMtSamples[,1:(num.sample/2)]), type="l", col="yellow", lwd=3)
points(times, rowMeans(adjMtSamples[,(num.sample/2+1):num.sample]), type="l", col="green", lwd=3)
#
plot(times,MtSample.group1[,1],type="n", ylim=c(-3,4), ylab="mean of bump function value & sample mean value")
for(i in 1:(num.sample)){
    points(times, adjMtSamples[,i],type="l", col="grey")
}
points(times, 0.5*(bumpval1+bumpval2), type="l", col="black", lwd=3)
points(times, rowMeans(adjMtSamples), type="l", col="yellow", lwd=3)





# (b) Align the curves using continuous registration. 
# Plot the resulting curves and include a mean function. 
# Comment on any diﬀerences with (a) and if the registered curves exhibit any odd patterns. 

# basic fitting to b-spline
b_spline_basis = create.bspline.basis(c(0,1), nbasis=6)
adjMtSample.f = smooth.basis(times, adjMtSamples, b_spline_basis)$fd

#plot: curves before alignment
par(mfrow=c(1,3))
plot(adjMtSample.f, ylim=c(-3,4), col='grey'); plot(mean(adjMtSample.f), add=T, lwd=3)

#apply continuous registration
adjMtSample.f.reg = register.fd(adjMtSample.f) #<-take long~ time

#plot: aligned curves
plot(adjMtSample.f.reg$regfd, col='grey', ylim=c(-3,4))
plot(mean(adjMtSample.f.reg$regfd), add=T, lwd=3, col='red')

#plot: compare only mean
plot(mean(adjMtSample.f), lwd=3, col='black', ylim=c(-3,4))
plot(mean(adjMtSample.f.reg$regfd), add=T, lwd=3, col='red')

# (c) Carry out an FPCA with one PC on the unaligned and aligned curves separately. 
# For each, do a simple linear regression of the score onto a dummy variable (coded 0/1) 
# indicating which type of mean the function had 
# (i.e. is it from the ﬁrst or second half of the sample). 
# Calculate a p-value to determine if the estimated slope parameters you get are signiﬁcant.
# Compare with the aligned and unaligned curves. What did aligning do to the p-value? 
# You may want to rerun your simulations a few times to see how the p-values change. 
adjMtSample.f.pca = pca.fd(adjMtSample.f, nharm=1)
adjMtSample.f.reg.pca = pca.fd(adjMtSample.f.reg$regfd, nharm=1)


#pca component plot
par(mfrow=c(1,1))
plot(adjMtSample.f.pca$harmonics, ylim=c(-3,4))
plot(adjMtSample.f.reg.pca$harmonics, col='red', add=T)
adjMtSample.f.pca$varprop
adjMtSample.f.reg.pca$varprop #왜 줄어드냐ㅋㅋ


#regression
dummy_var = c(rep(0,num.sample/2),rep(1,num.sample/2))
un_align_score = adjMtSample.f.pca$score
align_score = adjMtSample.f.reg.pca$score

par(mfrow=c(1,2))
lm1_unaligned= lm(un_align_score ~ dummy_var)
plot(dummy_var, un_align_score, ylim=c(-1.5, 1))
abline(lm1_unaligned$coefficient, col='red')
summary(lm1_unaligned)

lm2_aligned= lm(align_score ~ dummy_var)
plot(dummy_var, align_score, ylim=c(-1.5, 1))
abline(lm2_aligned$coefficient, col='red')
summary(lm2_aligned)
