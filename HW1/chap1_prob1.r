# R file

library(fda)
# ?pinch
length(pinch[,1]) #151
length(pinch[1,]) #20 : 한 열이 한 curve
ts.plot(pinch) # 개형
head(pinch)

# 1.6.1.a. fit the fd object
b_spline_basis<-create.bspline.basis(c(1,151), nbasis=15)
pinch.F<-Data2fd(1:151, pinch, b_spline_basis)
plot(pinch.F)

#for test
names(pinch.F)
pinch.F$coefs

# 1.6.1.b. mean&sd and add to the plot
pinch.F.mean <- mean(pinch.F)
pinch.F.mean$coefs #흠..정확한 의미? 노치마다?

pinch.F.std <- std.fd(pinch.F)
pinch.F.std$coefs

par(mar=c(4,4,1,1))
plot(pinch.F, col="grey")
plot(pinch.F.mean,lwd=4,add=TRUE)
plot(pinch.F.std, lwd=4, add=TRUE, col="red")
#pointwise라는게 spline마다를 의미하나??? 걍 pinch matrix 가져다놓고 깡으로 구해??

#1.6.1.c
pinch.F.cov <- var.fd(pinch.F)
dim(pinch.F.cov$coef) #15*15
grid <- 1:15
pinch.F.cov.mat = eval.bifd(grid, grid, pinch.F.cov) #<-뭐하는 코드인지?
persp(grid, grid, pinch.F.cov.mat, xlab="s", ylab="t", zlab="cov(s,t)")
contour(grid, grid, pinch.F.cov.mat)

#1.6.1.d
pinch.F.pca = pca.fd(pinch.F, nharm=4)
plot(pinch.F.pca$harmonics, lwd=3)
sum(pinch.F.pca$varprop) #0.98
# nharm : sum of varprop
# 4 : 0.987
# 3 : 0.967
# 2 : 0.921
# 1 : 0.672