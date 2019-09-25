# R file

#chap2 problem 2.5.2.

library(fda)
library(fds)
# ?FedYieldcurve #<-monthly data

#(a)
# Smooth the interest rates (yields) in January 1982 using a B–spline basis with four basis functions. 
# Plot the raw and smoothed interest rates on one graph. 
yield1982 = FedYieldcurve$y[,1:12]
terms = FedYieldcurve$x

b_spline_basis = create.bspline.basis(c(1,120), nbasis=4)
yield.F = Data2fd(terms, yield1982, b_spline_basis)
plot(yield.F, col='black', xlim=c(0,120), ylim=c(7.5,15))
for(i in 1:12){
    points(terms, yield1982[,i], type='l', col='gray')
}


#(b)
# Re–ﬁt the January 1982 yields using a penalized smoothing based on six basis functions (as many as data points)
# with the smoothing parameter λ = 1, and the second derivative as the penalty operator. 
# Add the smooth in red to the graph you obtained in part (a) and comment on the result.
b_spline_basis6 = create.bspline.basis(c(1,120), nbasis=6)
par.F1 = fdPar(b_spline_basis6, Lfdobj=2, lambda=1)
yield.S1 = smooth.basis(terms, yield1982, par.F1)

par(new=T)
plot(yield.S1, xlim=c(0,120), ylim=c(7.5,15), ylab="", xlab="", col='red')
mean(yield.S1$gcv)
# plot(deriv(yield.S1$fd))

#(c)
# Repeat part (b) with several other smoothing parameters λ. Which λ gives the most informative smooth curve?
candid_lambda = seq(1, 100, 1)
gcvs = rep(0, length(candid_lambda))
for(i in 1:length(candid_lambda)){
    par.F2 = fdPar(b_spline_basis2, Lfdobj=2, lambda=candid_lambda[i])
    yield.S2 = smooth.basis(terms, yield1982, par.F2)
    gcvs[i]=mean(yield.S2$gcv)
}
plot(candid_lambda, gcvs)

#take lambda=80!
par.F80 = fdPar(b_spline_basis6, Lfdobj=2, lambda=80)
yield.S80 = smooth.basis(terms, yield1982, par.F80)
plot(yield.F, col='black', xlim=c(0,120), ylim=c(7.5,15))
for(i in 1:12){
    par(new=T)
    points(terms, yield1982[,i], type='l', col='gray')
}
par(new=T)
plot(yield.S1, xlim=c(0,120), ylim=c(7.5,15), ylab="", xlab="", col='red')
par(new=T)
plot(yield.S80, xlim=c(0,120), ylim=c(7.5,15), ylab="", xlab="", col='blue')
mean(yield.S80$gcv)



#for 2nd derivative
plot(deriv.fd(yield.F), col='black', ylim=c(-0.1,0.6))
par(new=T)
plot(deriv(yield.S1$fd), col='red', ylim=c(-0.1,0.6))
par(new=T)
plot(deriv(yield.S80$fd), col='blue', ylim=c(-0.1,0.6))