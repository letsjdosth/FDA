# R file

library(fds)
library(fda)

yield = FedYieldcurve
dim(yield$y) # 6 * 330. 한 curve마다 6개 point
terms = yield$x

# (a)
# On one graph, plot the interest rates x(tj) for January 1982 and for June 2009 against the maturity terms tj. 
# How do the interest rates in these two months compare?
plot(terms, yield$y[,1], pch=15, ylab="Yield", ylim=c(0,16)) 
points(terms, yield$y[,330], pch=16)


# (b)
# Convert the yield data to functional objects using bspline basis with four basis functions. 
# Calculate and plot the the mean yield function. 
# What is the average behavior of interest rates as a function of the maturity?
b_spline_basis = create.bspline.basis(c(1,120), nbasis=4)
yield.F = Data2fd(terms, yield$y, b_spline_basis)
plot(yield.F, col='gray')

yield.F.mean = mean(yield.F)
plot(yield.F.mean,lwd=4,add=TRUE)

# (c)
# Plot the ﬁrst principal component of the interest rate curves. 
# What percentage of variance does this component explain? 
# Interpret the plot and the percentage of variance.
yield.F.pca = pca.fd(yield.F, nharm=1)
plot(yield.F.pca$harmonics, lwd=3, ylim=c(-1,1))
yield.F.pca$varprop #0.999982

