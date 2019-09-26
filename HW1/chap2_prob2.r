# R file

#chap2 problem 2.5.2.

library(fda)
library(fds)
# ?FedYieldcurve #<-monthly data

#(a)
# Smooth the interest rates (yields) in January 1982 using a B–spline basis with four basis functions. 
# Plot the raw and smoothed interest rates on one graph. 
yield.J1982 = FedYieldcurve$y[,1]
terms = FedYieldcurve$x

b_spline_basis4 = create.bspline.basis(c(1,120), nbasis=4)
yield.F4 = Data2fd(terms, yield.J1982, b_spline_basis4)
plot(yield.F4, col='black', xlim=c(0,120), ylim=c(12.5,15.5), xlab="maturity", ylab="yield rate")
points(terms, yield.J1982, type='l', col='gray')


#(b)
# Re–ﬁt the January 1982 yields using a penalized smoothing based on six basis functions (as many as data points)
# with the smoothing parameter λ = 1, and the second derivative as the penalty operator. 
# Add the smooth in red to the graph you obtained in part (a) and comment on the result.
b_spline_basis6 = create.bspline.basis(c(1,120), nbasis=6)
par.F6.S1 = fdPar(b_spline_basis6, Lfdobj=2, lambda=1)
yield.F6.S1 = smooth.basis(terms, yield.J1982, par.F6.S1)

par(new=T)
plot(yield.F6.S1, xlim=c(0,120), ylim=c(12.5,15.5), ylab="", xlab="", col='red')
yield.F6.S1$gcv
# plot(deriv(yield.F6.S1$fd))

#(c)
# Repeat part (b) with several other smoothing parameters λ. Which λ gives the most informative smooth curve?
candid_lambda = seq(1, 200, 1)
gcvs = rep(0, length(candid_lambda))
for(i in 1:length(candid_lambda)){
    par.F6.Si = fdPar(b_spline_basis6, Lfdobj=2, lambda=candid_lambda[i])
    yield.F6.Si = smooth.basis(terms, yield.J1982, par.F6.Si)
    gcvs[i]=mean(yield.F6.Si$gcv)
}
plot(candid_lambda, gcvs)
which.min(gcvs) #91


#take lambda=91!
par.F6.S91 = fdPar(b_spline_basis6, Lfdobj=2, lambda=91)
yield.F6.S91 = smooth.basis(terms, yield.J1982, par.F6.S91)

par(mfrow=c(1,2))
plot(yield.F4, col='black', xlim=c(0,120), ylim=c(12.5,15.5),, xlab="maturity", ylab="yield rate")
points(terms, yield.J1982, type='l', col='gray')
par(new=T)
plot(yield.F6.S1, xlim=c(0,120), ylim=c(12.5,15.5), ylab="", xlab="", col='red')
par(new=T)
plot(yield.F6.S91, xlim=c(0,120), ylim=c(12.5,15.5), ylab="", xlab="", col='blue')


#for 2nd derivative
plot(deriv.fd(yield.F4), col='black', ylim=c(-0.1,0.4), xlab="maturity", ylab="2nd derivative of yield")
par(new=T)
plot(deriv(yield.F6.S1$fd), col='red', ylim=c(-0.1,0.4), ylab="", xlab="")
par(new=T)
plot(deriv(yield.F6.S91$fd), col='blue', ylim=c(-0.1,0.4), ylab="", xlab="")