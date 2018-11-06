
rm(list=ls())

n = 50

alpha = 0.05
alpha_tilde = floor((n+1)*alpha)/(n+1)




######################################
# generate mixture of normals data 

x = sort(runif(n, -1, 1))
a = 0.5
b = 1.2
y_true = exp(a + b * x)
shape = 10
y = rgamma(n, rate=shape/y_true, shape=shape)


######################################
# conformity score

# kde usse (ydata,ynew)
# evaluated at yat

resid = function(dat,xnew,ynew,yat) {
  obj = glm(c(dat$y,ynew) ~ c(dat$x,xnew), family=Gamma(link="log"))
  ypred = predict(obj, xnew=xnew)
  return(abs(ypred - yat))
}

######################################

pi_n = function(dat,xnew,ynew) {
  yscores = sapply(dat$y, function(yat)resid(dat,xnew,ynew,yat))
  yscore_here = resid(dat,xnew,ynew,ynew)
  return(mean(c(yscores,yscore_here) <= yscore_here))
}

######################
# plot setup

par(mfrow=c(3,1), mar=c(4,4,1,1), bty="n")
xlim=c(-1,1.3)

####################
# scatter plot 

plot(x,y, pch=16, xlim=xlim)

xn1 = -0.25
ystar = 2.5

abline(v=xn1, col="red")
points(xn1,ystar, pch=21, col="white", bg="red", cex=2, lwd=2)
text(xn1,ystar, pos=4, "y", col="red")
mtext(expression(x[n+1]), side=1, line=1, at=xn1, col="red")

xnew = c(x,xn1)
ynew = c(y,ystar)
obj = glm(ynew ~ xnew, family=Gamma(link="log"))
ypred = predict(obj, type="response")
ord = order(xnew)
lines(xnew[ord], ypred[ord], col="red")
text(xnew[n], ypred[n], expression(hat(mu)[y](x)), pos=4, col="red")

####################
# resid plot 

resids = abs(ypred - ynew)
hist(resids, main="", col="gray", freq=FALSE, xlab="Absolute Residuals", axes=FALSE)
axis(2)
axis(1, at=0:4, lab=0:4)
points(resids,rep(0,n+1), pch=4)
points(resids[n+1], 0, pch=16, col="red")
mtext(expression(abs(y-hat(mu)[y](x[n+1]))), 
     at=resids[n+1],
     side=1,
     line=1,
     col="red")


# pi's

dat = data.frame(x=x,y=y)

yseq = seq(1e-4, 7, by=0.01)

pi_ns = sapply(yseq, function(yi) pi_n(dat,xn1,yi))

pi_step = stepfun(yseq, c(0,pi_ns))

plot(yseq, pi_step(yseq), ylim=c(0,1), 
     xlim=range(yseq),
     type="l", 
     ylab=expression(pi_n (y)), 
     xlab="y")

abline(h=1-alpha_tilde, col="red", lty="dashed")

# intervals

switches = -1*c(0,diff(pi_ns > 1-alpha_tilde))

ups = which(switches==1)
downs = which(switches==-1)

if(length(downs)>length(ups)) {
  ups[1] = 1

}

for(i in 1:length(ups)) {
  rect(yseq[ups[i]], -1, yseq[downs[i]], 2, col=rgb(0,0,0,alpha=0.3))
  text((yseq[ups[i]] + yseq[downs[i]])/2, 0.6, expression(C^alpha(y)))
}



