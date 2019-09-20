

The following performs our Monte Carlo simulation of $B = 50$ iterations 
when $n = 150$.

<<regression.150.2, cache = TRUE>>=
set.seed(13)
beta <- c(2, 5)
n <- 150
bins <- 2
B <- 50
system.time(out.regression.150.2 <- do.call(cbind, 
  lapply(1:B, FUN = function(j){
    unlist(regression.simulator(beta = beta, n = n, 
      bins = bins))
})))
@

<<>>=
regression.150.2 <- cbind( 
  rowMeans(out.regression.150.2, na.rm = TRUE),  
  apply(out.regression.150.2, 1, 
  FUN = function(x){ 
    sds <- sd(x, na.rm = TRUE)
    lengths <- length(which(!is.na(x)))
    sds / sqrt(lengths)
  }))
@

The following performs our Monte Carlo simulation of $B = 50$ iterations 
when $n = 150$.

<<regression.250.3, cache = TRUE>>=
n <- 250
bins <- 3
system.time(out.regression.250.3 <- do.call(cbind, 
  lapply(1:B, FUN = function(j){
    unlist(regression.simulator(beta = beta, n = n, 
      bins = bins))
})))
@

<<>>=
regression.250.3 <- cbind( 
  rowMeans(out.regression.250.3, na.rm = TRUE),  
  apply(out.regression.250.3, 1, 
  FUN = function(x){ 
    sds <- sd(x, na.rm = TRUE)
    lengths <- length(which(!is.na(x)))
    sds / sqrt(lengths)
  }))
@

The following performs our Monte Carlo simulation of $B = 50$ iterations 
when $n = 500$.

<<regression.500.3, cache = TRUE>>=
n <- 500
system.time(out.regression.500.3 <- do.call(cbind, 
  lapply(1:B, FUN = function(j){
    unlist(regression.simulator(beta = beta, n = n, 
      bins = bins))
})))
@

<<>>=
regression.500.3 <- cbind( 
  rowMeans(out.regression.500.3, na.rm = TRUE),  
  apply(out.regression.500.3, 1, 
  FUN = function(x){ 
    sds <- sd(x, na.rm = TRUE)
    lengths <- length(which(!is.na(x)))
    sds / sqrt(lengths)
  }))
@



\subsection{Results}

Results form our simulations are depicted in 
Table~\ref{Tab:regression-results} and Figure~\ref{Fig:regresion.inx}.  
This table and figure depicts the estimatated area, prediction 
error, and local coverage probabilities for all five considered prediction 
regions. 

In these simulations, errors about the estimated mean function are symmetric 
and homogeneous across the support.  We therefore expect for the LS and LSLW 
conformal prediction regions to perform nearly as well as the oracle HD 
prediction region.  These prediction regions also exhibit finite-sample 
marginal validity, local validity with respect to binning, and near 
conditional validity across the support.  The parametric conformal prediction 
region is similar to the LS and LSLW conformal prediciton regions and the HD 
prediction region in area, prediction error, finite-sample coverage 
properties, and appearance.  However, the parametric conformal prediciton 
region is slightly larger and gives more conservative coverage that these 
other prediction regions.  The nonparametric conformal prediciton region 
is larger and gives larger prediciton errors than the other prediciton 
regions.  It also appears to not visually fit the data well while the others 
do as seen in Section~\ref{sec:regressionplots}.


<<echo = FALSE>>=
regression.150.2[, 1] <- round(regression.150.2[, 1], 3)
regression.250.3[, 1] <- round(regression.250.3[, 1], 3)
regression.500.3[, 1] <- round(regression.500.3[, 1], 3)
regression.150.2[, 2] <- round(regression.150.2[, 2], 4)
regression.250.3[, 2] <- round(regression.250.3[, 2], 4)
regression.500.3[, 2] <- round(regression.500.3[, 2], 4)
@


\begin{table}[h!]
\tiny
\begin{center}
\begin{tabular}{llccccc}
  & & parametric conformal & nonparametric conformal & LS conformal & 
    LSLW conformal & HD region \\ 
  $n = 150$
    & marginal coverage &
  $\Sexpr{regression.150.2[1, 1]} \; (\Sexpr{regression.150.2[1, 2]})$ & 
  $\Sexpr{regression.150.2[31, 1]} \; (\Sexpr{regression.150.2[31, 2]})$ & 
  $\Sexpr{regression.150.2[61, 1]} \; (\Sexpr{regression.150.2[61, 2]})$ & 
  $\Sexpr{regression.150.2[91, 1]}  (\Sexpr{regression.150.2[91, 2]})$ & 
  $\Sexpr{regression.150.2[121, 1]} \; (\Sexpr{regression.150.2[121, 2]})$ \\
    & local coverage when $0 < x < 1/2$ & 
  $\Sexpr{regression.150.2[2, 1]} \; (\Sexpr{regression.150.2[2, 2]})$ & 
  $\Sexpr{regression.150.2[32, 1]} \; (\Sexpr{regression.150.2[32, 2]})$ & 
  $\Sexpr{regression.150.2[62, 1]} \; (\Sexpr{regression.150.2[62, 2]})$ & 
  $\Sexpr{regression.150.2[92, 1]} \; (\Sexpr{regression.150.2[92, 2]})$ & 
  $\Sexpr{regression.150.2[122, 1]} \; (\Sexpr{regression.150.2[122, 2]})$ \\
    & local coverage when $1/2 \leq x < 1$ & 
  $\Sexpr{regression.150.2[3, 1]} \; (\Sexpr{regression.150.2[3, 2]})$ & 
  $\Sexpr{regression.150.2[33, 1]} \; (\Sexpr{regression.150.2[33, 2]})$ & 
  $\Sexpr{regression.150.2[63, 1]} \; (\Sexpr{regression.150.2[63, 2]})$ & 
  $\Sexpr{regression.150.2[93, 1]} \; (\Sexpr{regression.150.2[93, 2]})$ & 
  $\Sexpr{regression.150.2[123, 1]} \; (\Sexpr{regression.150.2[123, 2]})$ \\
    & area &
  $\Sexpr{regression.150.2[29, 1]} \; (\Sexpr{regression.150.2[29, 2]})$ & 
  $\Sexpr{regression.150.2[59, 1]} \; (\Sexpr{regression.150.2[59, 2]})$ & 
  $\Sexpr{regression.150.2[89, 1]} \; (\Sexpr{regression.150.2[89, 2]})$ & 
  $\Sexpr{regression.150.2[119, 1]} \; (\Sexpr{regression.150.2[119, 2]})$ & 
  $\Sexpr{regression.150.2[149, 1]} \; (\Sexpr{regression.150.2[149, 2]})$ \\
    & prediction error &
  $\Sexpr{regression.150.2[30, 1]} \; (\Sexpr{regression.150.2[30, 2]})$ & 
  $\Sexpr{regression.150.2[60, 1]} \; (\Sexpr{regression.150.2[60, 2]})$ & 
  $\Sexpr{regression.150.2[90, 1]} \; (\Sexpr{regression.150.2[90, 2]})$ & 
  $\Sexpr{regression.150.2[120, 1]} \; (\Sexpr{regression.150.2[120, 2]})$ & 
  $\Sexpr{regression.150.2[150, 1]} \; (\Sexpr{regression.150.2[150, 2]})$ \\
  \hline
  $n = 250$  
    & marginal coverage &
  $\Sexpr{regression.250.3[1, 1]} \; (\Sexpr{regression.250.3[1, 2]})$ & 
  $\Sexpr{regression.250.3[32, 1]} \; (\Sexpr{regression.250.3[32, 2]})$ & 
  $\Sexpr{regression.250.3[63, 1]} \; (\Sexpr{regression.250.3[63, 2]})$ & 
  $\Sexpr{regression.250.3[94, 1]} \; (\Sexpr{regression.250.3[94, 2]})$ & 
  $\Sexpr{regression.250.3[125, 1]} \; (\Sexpr{regression.250.3[125, 2]})$ \\
    & local coverage when $0 < x < 1/3$ & 
  $\Sexpr{regression.250.3[2, 1]} \; (\Sexpr{regression.250.3[2, 2]})$ & 
  $\Sexpr{regression.250.3[33, 1]} \; (\Sexpr{regression.250.3[33, 2]})$ & 
  $\Sexpr{regression.250.3[64, 1]} \; (\Sexpr{regression.250.3[64, 2]})$ & 
  $\Sexpr{regression.250.3[95, 1]} \; (\Sexpr{regression.250.3[95, 2]})$ & 
  $\Sexpr{regression.250.3[126, 1]} \; (\Sexpr{regression.250.3[126, 2]})$ \\
    & local coverage when $1/3 \leq x < 2/3$ & 
  $\Sexpr{regression.250.3[3, 1]} \; (\Sexpr{regression.250.3[3, 2]})$ & 
  $\Sexpr{regression.250.3[34, 1]} \; (\Sexpr{regression.250.3[34, 2]})$ & 
  $\Sexpr{regression.250.3[65, 1]} \; (\Sexpr{regression.250.3[65, 2]})$ & 
  $\Sexpr{regression.250.3[96, 1]} \; (\Sexpr{regression.250.3[96, 2]})$ & 
  $\Sexpr{regression.250.3[127, 1]} \; (\Sexpr{regression.250.3[127, 2]})$ \\
    & local coverage when $2/3 \leq x < 1$ &
  $\Sexpr{regression.250.3[4, 1]} \; (\Sexpr{regression.250.3[4, 2]})$ & 
  $\Sexpr{regression.250.3[35, 1]} \; (\Sexpr{regression.250.3[35, 2]})$ & 
  $\Sexpr{regression.250.3[66, 1]} \; (\Sexpr{regression.250.3[66, 2]})$ & 
  $\Sexpr{regression.250.3[97, 1]} \; (\Sexpr{regression.250.3[97, 2]})$ & 
  $\Sexpr{regression.250.3[128, 1]} \; (\Sexpr{regression.250.3[128, 2]})$ \\
    & area & 
  $\Sexpr{regression.250.3[30, 1]} \; (\Sexpr{regression.250.3[30, 2]})$ & 
  $\Sexpr{regression.250.3[61, 1]} \; (\Sexpr{regression.250.3[61, 2]})$ & 
  $\Sexpr{regression.250.3[92, 1]} \; (\Sexpr{regression.250.3[92, 2]})$ & 
  $\Sexpr{regression.250.3[123, 1]} \; (\Sexpr{regression.250.3[123, 2]})$ & 
  $\Sexpr{regression.250.3[154, 1]} \; (\Sexpr{regression.250.3[154, 2]})$ \\
    & prediction error & 
  $\Sexpr{regression.250.3[31, 1]} \; (\Sexpr{regression.250.3[31, 2]})$ & 
  $\Sexpr{regression.250.3[62, 1]} \; (\Sexpr{regression.250.3[62, 2]})$ & 
  $\Sexpr{regression.250.3[93, 1]} \; (\Sexpr{regression.250.3[93, 2]})$ & 
  $\Sexpr{regression.250.3[124, 1]} \; (\Sexpr{regression.250.3[124, 2]})$ & 
  $\Sexpr{regression.250.3[155, 1]} \; (\Sexpr{regression.250.3[155, 2]})$ \\
  \hline
  $n = 500$
    & marginal coverage & 
  $\Sexpr{regression.500.3[1, 1]} \; (\Sexpr{regression.500.3[1, 2]})$ & 
  $\Sexpr{regression.500.3[32, 1]} \; (\Sexpr{regression.500.3[32, 2]})$ & 
  $\Sexpr{regression.500.3[63, 1]} \; (\Sexpr{regression.500.3[63, 2]})$ & 
  $\Sexpr{regression.500.3[94, 1]} \; (\Sexpr{regression.500.3[94, 2]})$ & 
  $\Sexpr{regression.500.3[125, 1]} \; (\Sexpr{regression.500.3[125, 2]})$ \\
    & local coverage when $0 < x < 1/3$ & 
  $\Sexpr{regression.500.3[2, 1]} \; (\Sexpr{regression.500.3[2, 2]})$ & 
  $\Sexpr{regression.500.3[33, 1]} \; (\Sexpr{regression.500.3[33, 2]})$ & 
  $\Sexpr{regression.500.3[64, 1]} \; (\Sexpr{regression.500.3[64, 2]})$ & 
  $\Sexpr{regression.500.3[95, 1]} \; (\Sexpr{regression.500.3[95, 2]})$ & 
  $\Sexpr{regression.500.3[126, 1]} \; (\Sexpr{regression.500.3[126, 2]})$ \\
    & local coverage when $1/3 \leq x < 2/3$ & 
  $\Sexpr{regression.500.3[3, 1]} \; (\Sexpr{regression.500.3[3, 2]})$ & 
  $\Sexpr{regression.500.3[34, 1]} \; (\Sexpr{regression.500.3[34, 2]})$ & 
  $\Sexpr{regression.500.3[65, 1]} \; (\Sexpr{regression.500.3[65, 2]})$ & 
  $\Sexpr{regression.500.3[96, 1]} \; (\Sexpr{regression.500.3[96, 2]})$ & 
  $\Sexpr{regression.500.3[127, 1]} \; (\Sexpr{regression.500.3[127, 2]})$ \\
    & local coverage when $2/3 \leq x < 1$ & 
  $\Sexpr{regression.500.3[4, 1]} \; (\Sexpr{regression.500.3[4, 2]})$ & 
  $\Sexpr{regression.500.3[35, 1]} \; (\Sexpr{regression.500.3[35, 2]})$ & 
  $\Sexpr{regression.500.3[66, 1]} \; (\Sexpr{regression.500.3[66, 2]})$ & 
  $\Sexpr{regression.500.3[97, 1]} \; (\Sexpr{regression.500.3[97, 2]})$ & 
  $\Sexpr{regression.500.3[128, 1]} \; (\Sexpr{regression.500.3[128, 2]})$ \\
    & area & 
  $\Sexpr{regression.500.3[30, 1]} \; (\Sexpr{regression.500.3[30, 2]})$ & 
  $\Sexpr{regression.500.3[61, 1]} \; (\Sexpr{regression.500.3[61, 2]})$ & 
  $\Sexpr{regression.500.3[92, 1]} \; (\Sexpr{regression.500.3[92, 2]})$ & 
  $\Sexpr{regression.500.3[123, 1]} \; (\Sexpr{regression.500.3[123, 2]})$ & 
  $\Sexpr{regression.500.3[154, 1]} \; (\Sexpr{regression.500.3[154, 2]})$ \\
    & prediction error & 
  $\Sexpr{regression.500.3[31, 1]} \; (\Sexpr{regression.500.3[31, 2]})$ & 
  $\Sexpr{regression.500.3[62, 1]} \; (\Sexpr{regression.500.3[62, 2]})$ & 
  $\Sexpr{regression.500.3[93, 1]} \; (\Sexpr{regression.500.3[93, 2]})$ & 
  $\Sexpr{regression.500.3[124, 1]} \; (\Sexpr{regression.500.3[124, 2]})$ & 
  $\Sexpr{regression.500.3[155, 1]} \; (\Sexpr{regression.500.3[155, 2]})$ 
\end{tabular}
\end{center}
\caption{Diagnostics for conformal prediction regions for linear regression 
  models with normal errors and constant variance.  Local and marginal 
  coverage properties, areas, and prediction errors are presented for the 
    parametric conformal prediction region (third column),
    nonparametric conformal prediction region (fourth column),
    LS conformal prediction region (fifth column), 
    LSLW conformal prediction region (sixth column), and 
    HD prediction region (seventh column). Standard errors are in parentheses.}
\label{Tab:regression-results}
\end{table}





\newpage
\begin{figure}[h!]
\begin{center}
<<Fig-regression-inx-500, echo = FALSE, fig.height = 3>>=
par(mfrow = c(1,3), oma = c(4,4,1,0), mar = c(1,3,3,0))

plot.new()
plot.window(xlim = c(0, 1), ylim = c(0.60, 1.00))
lines(inx, regression.150.2[4:28, 1], col = "red", lty = 1, lwd = 1.5)
lines(inx, regression.150.2[34:58, 1], col = "green", lty = 2, lwd = 1.5)
lines(inx, regression.150.2[64:88, 1], col = "purple", lty = 1, lwd = 1.5)
lines(inx, regression.150.2[94:118, 1], col = "black", lty = 2, lwd = 1.5)
lines(inx, regression.150.2[124:148, 1], col = "blue", lty = 4, lwd = 1.5)
abline(h = 0.90, lty = 1, col = "grey", lwd = 1.5)
axis(1); axis(2)
mtext("n = 150", side = 3, cex = 1, line = 2)

plot.new()
plot.window(xlim = c(0, 1), ylim = c(0.60, 1.00))
lines(inx, regression.250.3[5:29, 1], col = "red", lty = 1, lwd = 1.5)
lines(inx, regression.250.3[36:60, 1], col = "green", lty = 2, lwd = 1.5)
lines(inx, regression.250.3[67:91, 1], col = "purple", lty = 1, lwd = 1.5)
lines(inx, regression.250.3[98:122, 1], col = "black", lty = 2, lwd = 1.5)
lines(inx, regression.250.3[129:153, 1], col = "blue", lty = 4, lwd = 1.5)
abline(h = 0.90, lty = 1, col = "grey", lwd = 1.5)
axis(1)
mtext("n = 250", side = 3, cex = 1, line = 2)

#plot.new()
#plot.window(xlim = c(0, 1), ylim = c(0.60, 1.00))
legend(0.15, 0.85, legend=c("parametric", 
  "nonparametric", "LSLW", "LS", "HD", "nominal coverage"), 
  col=c("red", "green", "purple", "black", "blue", "grey"), cex = 1,
  lty=c(1,2,1,2,4,1), bty = "n")

plot.new()
plot.window(xlim = c(0, 1), ylim = c(0.60, 1.00))
lines(inx, regression.500.3[5:29, 1], col = "red", lty = 1, lwd = 1.5)
lines(inx, regression.500.3[36:60, 1], col = "green", lty = 2, lwd = 1.5)
lines(inx, regression.500.3[67:91, 1], col = "purple", lty = 1, lwd = 1.5)
lines(inx, regression.500.3[98:122, 1], col = "black", lty = 2, lwd = 1.5)
lines(inx, regression.500.3[129:153, 1], col = "blue", lty = 4, lwd = 1.5)
abline(h = 0.90, lty = 1, col = "grey", lwd = 1.5)
axis(1)
mtext("n = 500", side = 3, cex = 1, line = 2)


mtext("x", side = 1, cex = 1, line = 2, outer = TRUE)
mtext("coverage", side = 2, cex = 1.25, line = 0, outer = TRUE)

@
\end{center}
\caption{Plot of the estimated coverage probabilities of prediction regions 
  across $x$ and sample sizes.}
\label{Fig:regresion.inx}
\end{figure}





\newpage
\section{Example plots of prediction regions}
\label{sec:plotsofregions}

In this section we construct prediction regions corresponding to the 
simulations and results in Sections~\ref{sec:Gamma} through 
\ref{sec:regressionplots}.  In Gamma analyses we generate a dataset for each 
shape parameter considered with $n = 150$ and in regression analyses we 
generate a dataset for all sample sizes considered.  For each of these datasets 
we depict the parametric, nonparametric, LS, LSLW conformal prediction regions 
over the observed data to visually assess the appropriateness of each 
prediction region.  The findings from these figures are consistent with the 
findings from our numerical diagnostics.  
The parametric conformal prediction region gives visually natural bounds for 
the observed data when the model is correctly specified in small to moderate 
sample sizes and is appropriate when the model is misspecified. 
The LSLW conformal prediction region gives visually natural bounds for the 
observed data when the model is correctly specified, is appropriate under mild 
model misspecification, and is ill-fitting in settings where deviations about 
an estimated mean function are clearly not symmetric.
The nonparametric conformal prediction region is coarse and is larger than 
necessary when the mean function is steep relative to its variability.  
The LS conformal prediction performs well when deviations about the estimated 
mean function are symmetirc and is sensitive to mild departures from that 
setting.  This prediciton region is seen to provide overcoverage 
(undercoverage) in regions of the predictor space where variability about the 
estimated mean function is relatively small (large).



\newpage
\begin{figure}[h!]
\begin{center}
<<conformal-plots, fig.height = 7, cache = TRUE, echo = FALSE>>=
par(mfrow = c(3, 4), oma = c(4,4,2,0), mar = c(1,1,1,1))

######### Sim setting A ####################
## generate random data
set.seed(13)
alpha <- 0.10
cores <- 6
n <- 150
bins <- 2
p <- k <- length(beta) - 1
x <- matrix(runif(n), ncol = p)

## generate gamma regression model 
beta <- c(1.25, -1); shape = 1
rate <- (cbind(1, x) %*% beta) * shape
y <- rgamma(n = n, shape = shape, rate = rate)
y <- y / sd(y)
data <- data.frame(y = y, x = x)
colnames(data)[2:(p+1)] <- paste("x", 1:p, sep = "")

## fit the Gamma regression model
fit <- glm(y ~ x1, family = "Gamma", data = data)
formula <- fit$formula
newdata <- data
respname <- all.vars(formula)[1]
newdata <- newdata[, !(colnames(data) %in% respname)]
newdata <- as.matrix(newdata)

## obtain parametric and nonparametric conformal 
## prediction regions
cpred <- conformal.glm(fit, parametric = TRUE, 
  nonparametric = FALSE, alpha = alpha,
  bins = bins, cores = cores)
paraCI <- cpred$paraconformal

cpred <- conformal.glm(fit, parametric = FALSE, 
  nonparametric = TRUE, alpha = alpha,
  bins = bins, cores = cores)
nonparaCI <- cpred$nonparaconformal

## obtain LS conformal prediction region
p1.tibs <- conformal.pred(x = cbind(x, x^2, x^3), y = y, 
  x0 = cbind(x, x^2, x^3), train.fun = train.fun, 
  predict.fun = predict.fun, alpha = alpha, 
  grid.method = "linear")
LSCI <- cbind(p1.tibs$lo, p1.tibs$up)

## obtain LSLW conformal prediction region
cubic.model <- lm(y ~ x + I(x^2) + I(x^3))
abs.resid <- abs(cubic.model$resid)
smooth.call <- smooth.spline(x, abs.resid, 
  nknots = 10)
lambda <- smooth.call$lambda
df <- smooth.call$df
mad.train.fun <- function(x, y, out = NULL){
  smooth.spline(x[, 1], y, lambda = lambda, 
  df = df, nknots = 10)
}
p2.tibs <- conformal.pred(x = cbind(x,x^2,x^3), y = y, 
  x0 = cbind(x,x^2,x^3), 
  train.fun = train.fun, predict.fun = predict.fun, 
  mad.train.fun = mad.train.fun,
  mad.predict.fun = mad.predict.fun,
  alpha = alpha)
LSLWCI <- cbind(p2.tibs$lo, p2.tibs$up)  

## create plots
ix <- sort(x, index.return = TRUE)$ix

## parametric conformal
plot.new()
plot.window(xlim = c(0,1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], paraCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], paraCI[ix, 2], col = "red", lwd = 1.5)
axis(2)
mtext("Parametric \n conformal", side = 3, cex = 1, line = 0)
mtext("Setting A", side = 2, cex = 1, line = 2)

## LSLW conformal
plot.new()
plot.window(xlim = c(0,1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSLWCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSLWCI[ix, 2], col = "red", lwd = 1.5)
mtext("LSLW \n conformal", side = 3, cex = 1, line = 0)

## nonparametric conformal
plot.nonparametric(region = nonparaCI, x = x, 
  y = y, bins = bins)
mtext("Nonparametric \n conformal", side = 3, cex = 1, line = 0)

## LS conformal
plot.new()
plot.window(xlim = c(0,1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSCI[ix, 2], col = "red", lwd = 1.5)
mtext("LS \n conformal", side = 3, cex = 1, line = 0)

######### Sim setting B ####################
## generate gamma regression model (misspecified parametric conformal)
beta <- c(0.5, 1); shape = 10
rate <- (cbind(1, x) %*% beta) * shape
y <- rgamma(n = n, shape = shape, rate = rate)
y <- y / sd(y)
data <- data.frame(y = y, x = x)
colnames(data)[2:(p+1)] <- paste("x", 1:p, sep = "")

## fit the misspecified regression model
fit <- glm(y ~ x1 + I(x1^2) + I(x1^3), 
  family = "gaussian", data = data)
formula <- fit$formula
newdata <- data
respname <- all.vars(formula)[1]
newdata <- newdata[, !(colnames(data) %in% respname)]
newdata <- as.matrix(newdata)

## obtain parametric and nonparametric conformal 
## prediction regions
cpred <- conformal.glm(fit, parametric = TRUE, 
  nonparametric = FALSE, alpha = alpha,
  bins = bins, cores = cores)
paraCI <- cpred$paraconformal

cpred <- conformal.glm(fit, parametric = FALSE, 
  nonparametric = TRUE, alpha = alpha,
  bins = bins, cores = cores)
nonparaCI <- cpred$nonparaconformal

## obtain LS conformal prediction region
p1.tibs <- conformal.pred(x = cbind(x, x^2, x^3), y = y, 
  x0 = cbind(x, x^2, x^3), train.fun = train.fun, 
  predict.fun = predict.fun, alpha = alpha, 
  grid.method = "linear")
LSCI <- cbind(p1.tibs$lo, p1.tibs$up)

## obtain LSLW conformal prediction region
cubic.model <- lm(y ~ x + I(x^2) + I(x^3))
abs.resid <- abs(cubic.model$resid)
smooth.call <- smooth.spline(x, abs.resid, 
  nknots = 10)
lambda <- smooth.call$lambda
df <- smooth.call$df
mad.train.fun <- function(x, y, out = NULL){
  smooth.spline(x[, 1], y, lambda = lambda, 
  df = df, nknots = 10)
}
p2.tibs <- conformal.pred(x = cbind(x,x^2,x^3), y = y, 
  x0 = cbind(x,x^2,x^3), 
  train.fun = train.fun, predict.fun = predict.fun, 
  mad.train.fun = mad.train.fun,
  mad.predict.fun = mad.predict.fun,
  alpha = alpha)
LSLWCI <- cbind(p2.tibs$lo, p2.tibs$up)  

## parametric conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], paraCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], paraCI[ix, 2], col = "red", lwd = 1.5)
axis(2)
mtext("Setting B", side = 2, cex = 1, line = 2)

## LSLW conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSLWCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSLWCI[ix, 2], col = "red", lwd = 1.5)

## nonparametric conformal
plot.nonparametric(region = nonparaCI, x = x, 
  y = y, bins = bins)

## LS conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSCI[ix, 2], col = "red", lwd = 1.5)

######### Sim setting C ####################
## generate linear regression model
beta <- c(2, 5)
mu <- (cbind(1, x) %*% beta) 
y <- rnorm(n = n, mean = mu, sd = 1)
y <- y / sd(y)
data <- data.frame(y = y, x = x)
colnames(data)[2:(p+1)] <- paste("x", 1:p, sep = "")

## fit the regression model
fit <- glm(y ~ x1, family = "gaussian", data = data)
formula <- fit$formula
newdata <- data
respname <- all.vars(formula)[1]
newdata <- newdata[, !(colnames(data) %in% respname)]
newdata <- as.matrix(newdata)

## obtain parametric and nonparametric conformal 
## prediction regions
cpred <- conformal.glm(fit, parametric = TRUE, 
  nonparametric = FALSE, alpha = alpha,
  bins = bins, cores = cores)
paraCI <- cpred$paraconformal

cpred <- conformal.glm(fit, parametric = FALSE, 
  nonparametric = TRUE, alpha = alpha,
  bins = bins, cores = cores)
nonparaCI <- cpred$nonparaconformal

## obtain LS conformal prediction region
p1.tibs <- conformal.pred(x = x, y = y, x0 = x, 
  train.fun = train.fun, 
  predict.fun = predict.fun, alpha = alpha, 
  grid.method = "linear")
LSCI <- cbind(p1.tibs$lo, p1.tibs$up)

## obtain LSLW conformal prediction region
fit <- lm(y ~ x1, data = data)
abs.resid <- abs(fit$resid)
smooth.call <- smooth.spline(x, abs.resid, 
  nknots = 10)
lambda <- smooth.call$lambda
df <- smooth.call$df
mad.train.fun <- function(x, y, out = NULL){
  smooth.spline(x[, 1], y, lambda = lambda, 
  df = df, nknots = 10)
}
p2.tibs <- conformal.pred(x = x, y = y, x0 = x, 
  train.fun = train.fun, predict.fun = predict.fun, 
  mad.train.fun = mad.train.fun,
  mad.predict.fun = mad.predict.fun,
  alpha = alpha)
LSLWCI <- cbind(p2.tibs$lo, p2.tibs$up)  

## parametric conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], paraCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], paraCI[ix, 2], col = "red", lwd = 1.5)
axis(1); axis(2)
mtext("Setting C", side = 2, cex = 1, line = 2)

## LSLW conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSLWCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSLWCI[ix, 2], col = "red", lwd = 1.5)
axis(1)

## nonparametric conformal
plot.nonparametric(region = nonparaCI, x = x, 
  y = y, bins = bins)
axis(1)

## LS conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSCI[ix, 2], col = "red", lwd = 1.5)
axis(1)

mtext("x", side = 1, cex = 1, line = 2, outer = TRUE)
mtext("y", side = 2, cex = 1, line = 3, outer = TRUE)
@
\end{center}
\caption{The depiction of conformal prediction regions when $n = 150$ that 
  appears in \citet{eck2019conformal}.  The rows display conformal prediction 
  regions across simulation settings.  The columns display the different 
  conformal prediction regions.  The top, middle, and bottom rows correspond 
  to simulation setting a with shape parameter equal to 1, simulation setting 
  b with shape parameter equal to 10, and simulation setting c respectively.  
  The first column displays the parametric conformal prediction region which 
  is misspecified in row 2, the second column displays the least squares 
  locally weighted conformal prediction region, the third column displays the 
  nonparametric conformal prediction region, and the fourth column displays 
  the least squares conformal prediction region.}
\label{conformal-plots}
\end{figure}











\newpage
\subsection{Plots corresponding to Section~\ref{sec:Gamma}}
\label{sec:gammaplots}

%Plot of conformal prediciton regions in sim setting A when n = 150
\begin{figure}[h!]
\begin{center}
<<conformal-plots-A-150, fig.height = 7, cache = TRUE, echo = FALSE>>=
par(mfrow = c(4, 4), oma = c(4,4,2,0), mar = c(1,2,1,1))

############ shape = 0.75 ################
## generate random data
set.seed(13)
alpha <- 0.10
cores <- 6
n <- 150
bins <- 2
p <- k <- length(beta) - 1
x <- matrix(runif(n), ncol = p)

## generate gamma regression model 
beta <- c(1.25, -1); shape = 0.75
rate <- (cbind(1, x) %*% beta) * shape
y <- rgamma(n = n, shape = shape, rate = rate)
y <- y / sd(y)
data <- data.frame(y = y, x = x)
colnames(data)[2:(p+1)] <- paste("x", 1:p, sep = "")

## fit the Gamma regression model
fit <- glm(y ~ x1, family = "Gamma", data = data)
formula <- fit$formula
newdata <- data
respname <- all.vars(formula)[1]
newdata <- newdata[, !(colnames(data) %in% respname)]
newdata <- as.matrix(newdata)

## obtain parametric and nonparametric conformal 
## prediction regions
cpred <- conformal.glm(fit, parametric = TRUE, 
  nonparametric = FALSE, alpha = alpha,
  bins = bins, cores = cores)
paraCI <- cpred$paraconformal

cpred <- conformal.glm(fit, parametric = FALSE, 
  nonparametric = TRUE, alpha = alpha,
  bins = bins, cores = cores)
nonparaCI <- cpred$nonparaconformal

## obtain LS conformal prediction region
p1.tibs <- conformal.pred(x = cbind(x, x^2, x^3), y = y, 
  x0 = cbind(x, x^2, x^3), train.fun = train.fun, 
  predict.fun = predict.fun, alpha = alpha, 
  grid.method = "linear")
LSCI <- cbind(p1.tibs$lo, p1.tibs$up)

## obtain LSLW conformal prediction region
cubic.model <- lm(y ~ x + I(x^2) + I(x^3))
abs.resid <- abs(cubic.model$resid)
smooth.call <- smooth.spline(x, abs.resid, 
  nknots = 10)
lambda <- smooth.call$lambda
df <- smooth.call$df
mad.train.fun <- function(x, y, out = NULL){
  smooth.spline(x[, 1], y, lambda = lambda, 
  df = df, nknots = 10)
}
p2.tibs <- conformal.pred(x = cbind(x,x^2,x^3), y = y, 
  x0 = cbind(x,x^2,x^3), 
  train.fun = train.fun, predict.fun = predict.fun, 
  mad.train.fun = mad.train.fun,
  mad.predict.fun = mad.predict.fun,
  alpha = alpha)
LSLWCI <- cbind(p2.tibs$lo, p2.tibs$up)  

## create plots
ix <- sort(x, index.return = TRUE)$ix

## parametric conformal
plot.new()
plot.window(xlim = c(0,1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], paraCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], paraCI[ix, 2], col = "red", lwd = 1.5)
axis(2)
mtext("Parametric \n conformal", side = 3, cex = 1, line = 0)
mtext("shape = 0.75", side = 2, cex = 1, line = 2)

## LSLW conformal
plot.new()
plot.window(xlim = c(0,1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSLWCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSLWCI[ix, 2], col = "red", lwd = 1.5)
mtext("LSLW\n conformal", side = 3, cex = 1, line = 0)

## nonparametric conformal
plot.nonparametric(region = nonparaCI, x = x, 
  y = y, bins = bins)
mtext("nonparametric\n conformal", side = 3, cex = 1, line = 0)

## LS conformal
plot.new()
plot.window(xlim = c(0,1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSCI[ix, 2], col = "red", lwd = 1.5)
mtext("LS\n conformal", side = 3, cex = 1, line = 0)


############ shape = 2 ################
## generate gamma regression model 
beta <- c(1.25, -1); shape = 2
rate <- (cbind(1, x) %*% beta) * shape
y <- rgamma(n = n, shape = shape, rate = rate)
y <- y / sd(y)
data <- data.frame(y = y, x = x)
colnames(data)[2:(p+1)] <- paste("x", 1:p, sep = "")

## fit the Gamma regression model
fit <- glm(y ~ x1, family = "Gamma", data = data)
formula <- fit$formula
newdata <- data
respname <- all.vars(formula)[1]
newdata <- newdata[, !(colnames(data) %in% respname)]
newdata <- as.matrix(newdata)

## obtain parametric and nonparametric conformal 
## prediction regions
cpred <- conformal.glm(fit, parametric = TRUE, 
  nonparametric = FALSE, alpha = alpha,
  bins = bins, cores = cores)
paraCI <- cpred$paraconformal

cpred <- conformal.glm(fit, parametric = FALSE, 
  nonparametric = TRUE, alpha = alpha,
  bins = bins, cores = cores)
nonparaCI <- cpred$nonparaconformal

## obtain LS conformal prediction region
p1.tibs <- conformal.pred(x = cbind(x, x^2, x^3), y = y, 
  x0 = cbind(x, x^2, x^3), train.fun = train.fun, 
  predict.fun = predict.fun, alpha = alpha, 
  grid.method = "linear")
LSCI <- cbind(p1.tibs$lo, p1.tibs$up)

## obtain LSLW conformal prediction region
cubic.model <- lm(y ~ x + I(x^2) + I(x^3))
abs.resid <- abs(cubic.model$resid)
smooth.call <- smooth.spline(x, abs.resid, 
  nknots = 10)
lambda <- smooth.call$lambda
df <- smooth.call$df
mad.train.fun <- function(x, y, out = NULL){
  smooth.spline(x[, 1], y, lambda = lambda, 
  df = df, nknots = 10)
}
p2.tibs <- conformal.pred(x = cbind(x,x^2,x^3), y = y, 
  x0 = cbind(x,x^2,x^3), 
  train.fun = train.fun, predict.fun = predict.fun, 
  mad.train.fun = mad.train.fun,
  mad.predict.fun = mad.predict.fun,
  alpha = alpha)
LSLWCI <- cbind(p2.tibs$lo, p2.tibs$up)  

## parametric conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], paraCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], paraCI[ix, 2], col = "red", lwd = 1.5)
axis(2)
mtext("shape = 2", side = 2, cex = 1, line = 2)

## LSLW conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSLWCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSLWCI[ix, 2], col = "red", lwd = 1.5)

## nonparametric conformal
plot.nonparametric(region = nonparaCI, x = x, 
  y = y, bins = bins)

## LS conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSCI[ix, 2], col = "red", lwd = 1.5)


############ shape = 10 ################
## generate gamma regression model 
beta <- c(1.25, -1); shape = 10
rate <- (cbind(1, x) %*% beta) * shape
y <- rgamma(n = n, shape = shape, rate = rate)
y <- y / sd(y)
data <- data.frame(y = y, x = x)
colnames(data)[2:(p+1)] <- paste("x", 1:p, sep = "")

## fit the Gamma regression model
fit <- glm(y ~ x1, family = "Gamma", data = data)
formula <- fit$formula
newdata <- data
respname <- all.vars(formula)[1]
newdata <- newdata[, !(colnames(data) %in% respname)]
newdata <- as.matrix(newdata)

## obtain parametric and nonparametric conformal 
## prediction regions
cpred <- conformal.glm(fit, parametric = TRUE, 
  nonparametric = FALSE, alpha = alpha,
  bins = bins, cores = cores)
paraCI <- cpred$paraconformal

cpred <- conformal.glm(fit, parametric = FALSE, 
  nonparametric = TRUE, alpha = alpha,
  bins = bins, cores = cores)
nonparaCI <- cpred$nonparaconformal

## obtain LS conformal prediction region
p1.tibs <- conformal.pred(x = cbind(x, x^2, x^3), y = y, 
  x0 = cbind(x, x^2, x^3), train.fun = train.fun, 
  predict.fun = predict.fun, alpha = alpha, 
  grid.method = "linear")
LSCI <- cbind(p1.tibs$lo, p1.tibs$up)

## obtain LSLW conformal prediction region
cubic.model <- lm(y ~ x + I(x^2) + I(x^3))
abs.resid <- abs(cubic.model$resid)
smooth.call <- smooth.spline(x, abs.resid, 
  nknots = 10)
lambda <- smooth.call$lambda
df <- smooth.call$df
mad.train.fun <- function(x, y, out = NULL){
  smooth.spline(x[, 1], y, lambda = lambda, 
  df = df, nknots = 10)
}
p2.tibs <- conformal.pred(x = cbind(x,x^2,x^3), y = y, 
  x0 = cbind(x,x^2,x^3), 
  train.fun = train.fun, predict.fun = predict.fun, 
  mad.train.fun = mad.train.fun,
  mad.predict.fun = mad.predict.fun,
  alpha = alpha)
LSLWCI <- cbind(p2.tibs$lo, p2.tibs$up)  

## parametric conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], paraCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], paraCI[ix, 2], col = "red", lwd = 1.5)
axis(2)
mtext("shape = 10", side = 2, cex = 1, line = 2)

## LSLW conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSLWCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSLWCI[ix, 2], col = "red", lwd = 1.5)

## nonparametric conformal
plot.nonparametric(region = nonparaCI, x = x, 
  y = y, bins = bins)

## LS conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSCI[ix, 2], col = "red", lwd = 1.5)


############ shape = 50 ################
## generate gamma regression model 
beta <- c(1.25, -1); shape = 50
rate <- (cbind(1, x) %*% beta) * shape
y <- rgamma(n = n, shape = shape, rate = rate)
y <- y / sd(y)
data <- data.frame(y = y, x = x)
colnames(data)[2:(p+1)] <- paste("x", 1:p, sep = "")

## fit the Gamma regression model
fit <- glm(y ~ x1, family = "Gamma", data = data)
formula <- fit$formula
newdata <- data
respname <- all.vars(formula)[1]
newdata <- newdata[, !(colnames(data) %in% respname)]
newdata <- as.matrix(newdata)

## obtain parametric and nonparametric conformal 
## prediction regions
cpred <- conformal.glm(fit, parametric = TRUE, 
  nonparametric = FALSE, alpha = alpha,
  bins = bins, cores = cores)
paraCI <- cpred$paraconformal

cpred <- conformal.glm(fit, parametric = FALSE, 
  nonparametric = TRUE, alpha = alpha,
  bins = bins, cores = cores)
nonparaCI <- cpred$nonparaconformal

## obtain LS conformal prediction region
p1.tibs <- conformal.pred(x = cbind(x, x^2, x^3), y = y, 
  x0 = cbind(x, x^2, x^3), train.fun = train.fun, 
  predict.fun = predict.fun, alpha = alpha, 
  grid.method = "linear")
LSCI <- cbind(p1.tibs$lo, p1.tibs$up)

## obtain LSLW conformal prediction region
cubic.model <- lm(y ~ x + I(x^2) + I(x^3))
abs.resid <- abs(cubic.model$resid)
smooth.call <- smooth.spline(x, abs.resid, 
  nknots = 10)
lambda <- smooth.call$lambda
df <- smooth.call$df
mad.train.fun <- function(x, y, out = NULL){
  smooth.spline(x[, 1], y, lambda = lambda, 
  df = df, nknots = 10)
}
p2.tibs <- conformal.pred(x = cbind(x,x^2,x^3), y = y, 
  x0 = cbind(x,x^2,x^3), 
  train.fun = train.fun, predict.fun = predict.fun, 
  mad.train.fun = mad.train.fun,
  mad.predict.fun = mad.predict.fun,
  alpha = alpha)
LSLWCI <- cbind(p2.tibs$lo, p2.tibs$up)  

## parametric conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], paraCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], paraCI[ix, 2], col = "red", lwd = 1.5)
axis(1);axis(2)
mtext("shape = 50", side = 2, cex = 1, line = 2)

## LSLW conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSLWCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSLWCI[ix, 2], col = "red", lwd = 1.5)
axis(1)

## nonparametric conformal
plot.nonparametric(region = nonparaCI, x = x, 
  y = y, bins = bins)
axis(1)

## LS conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSCI[ix, 2], col = "red", lwd = 1.5)
axis(1)

mtext("x", side = 1, cex = 1, line = 2, outer = TRUE)
mtext("y", side = 2, cex = 1, line = 2, outer = TRUE)
@
\end{center}
\caption{The depiction of conformal prediction regions under simulation 
  setting A when $n = 150$ and the number of bins equals 2.
}
\label{conformal-plots-A-150}
\end{figure}











\newpage
%Plot of conformal prediciton regions in sim setting A when n = 250
\begin{figure}[h!]
\begin{center}
<<conformal-plots-A-250, fig.height = 7, cache = TRUE, echo = FALSE>>=
par(mfrow = c(4, 4), oma = c(4,4,2,0), mar = c(1,2,1,1))

############ shape = 0.75 ################
## generate random data
set.seed(13)
alpha <- 0.10
cores <- 6
n <- 250
bins <- 3
p <- k <- length(beta) - 1
x <- matrix(runif(n), ncol = p)

## generate gamma regression model 
beta <- c(1.25, -1); shape = 0.75
rate <- (cbind(1, x) %*% beta) * shape
y <- rgamma(n = n, shape = shape, rate = rate)
y <- y / sd(y)
data <- data.frame(y = y, x = x)
colnames(data)[2:(p+1)] <- paste("x", 1:p, sep = "")

## fit the Gamma regression model
fit <- glm(y ~ x1, family = "Gamma", data = data)
formula <- fit$formula
newdata <- data
respname <- all.vars(formula)[1]
newdata <- newdata[, !(colnames(data) %in% respname)]
newdata <- as.matrix(newdata)

## obtain parametric and nonparametric conformal 
## prediction regions
cpred <- conformal.glm(fit, parametric = TRUE, 
  nonparametric = FALSE, alpha = alpha,
  bins = bins, cores = cores)
paraCI <- cpred$paraconformal

cpred <- conformal.glm(fit, parametric = FALSE, 
  nonparametric = TRUE, alpha = alpha,
  bins = bins, cores = cores)
nonparaCI <- cpred$nonparaconformal

## obtain LS conformal prediction region
p1.tibs <- conformal.pred(x = cbind(x, x^2, x^3), y = y, 
  x0 = cbind(x, x^2, x^3), train.fun = train.fun, 
  predict.fun = predict.fun, alpha = alpha, 
  grid.method = "linear")
LSCI <- cbind(p1.tibs$lo, p1.tibs$up)

## obtain LSLW conformal prediction region
cubic.model <- lm(y ~ x + I(x^2) + I(x^3))
abs.resid <- abs(cubic.model$resid)
smooth.call <- smooth.spline(x, abs.resid, 
  nknots = 10)
lambda <- smooth.call$lambda
df <- smooth.call$df
mad.train.fun <- function(x, y, out = NULL){
  smooth.spline(x[, 1], y, lambda = lambda, 
  df = df, nknots = 10)
}
p2.tibs <- conformal.pred(x = cbind(x,x^2,x^3), y = y, 
  x0 = cbind(x,x^2,x^3), 
  train.fun = train.fun, predict.fun = predict.fun, 
  mad.train.fun = mad.train.fun,
  mad.predict.fun = mad.predict.fun,
  alpha = alpha)
LSLWCI <- cbind(p2.tibs$lo, p2.tibs$up)  

## create plots
ix <- sort(x, index.return = TRUE)$ix

## parametric conformal
plot.new()
plot.window(xlim = c(0,1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], paraCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], paraCI[ix, 2], col = "red", lwd = 1.5)
axis(2)
mtext("Parametric\n conformal", side = 3, cex = 1, line = 0)
mtext("shape = 0.75", side = 2, cex = 1, line = 2)

## LSLW conformal
plot.new()
plot.window(xlim = c(0,1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSLWCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSLWCI[ix, 2], col = "red", lwd = 1.5)
mtext("LSLW\n conformal", side = 3, cex = 1, line = 0)

## nonparametric conformal
plot.nonparametric(region = nonparaCI, x = x, 
  y = y, bins = bins)
mtext("Nonparametric\n conformal", side = 3, cex = 1, line = 0)

## LS conformal
plot.new()
plot.window(xlim = c(0,1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSCI[ix, 2], col = "red", lwd = 1.5)
mtext("LS\n conformal", side = 3, cex = 1, line = 0)


############ shape = 2 ################
## generate gamma regression model 
beta <- c(1.25, -1); shape = 2
rate <- (cbind(1, x) %*% beta) * shape
y <- rgamma(n = n, shape = shape, rate = rate)
y <- y / sd(y)
data <- data.frame(y = y, x = x)
colnames(data)[2:(p+1)] <- paste("x", 1:p, sep = "")

## fit the Gamma regression model
fit <- glm(y ~ x1, family = "Gamma", data = data)
formula <- fit$formula
newdata <- data
respname <- all.vars(formula)[1]
newdata <- newdata[, !(colnames(data) %in% respname)]
newdata <- as.matrix(newdata)

## obtain parametric and nonparametric conformal 
## prediction regions
cpred <- conformal.glm(fit, parametric = TRUE, 
  nonparametric = FALSE, alpha = alpha,
  bins = bins, cores = cores)
paraCI <- cpred$paraconformal

cpred <- conformal.glm(fit, parametric = FALSE, 
  nonparametric = TRUE, alpha = alpha,
  bins = bins, cores = cores)
nonparaCI <- cpred$nonparaconformal

## obtain LS conformal prediction region
p1.tibs <- conformal.pred(x = cbind(x, x^2, x^3), y = y, 
  x0 = cbind(x, x^2, x^3), train.fun = train.fun, 
  predict.fun = predict.fun, alpha = alpha, 
  grid.method = "linear")
LSCI <- cbind(p1.tibs$lo, p1.tibs$up)

## obtain LSLW conformal prediction region
cubic.model <- lm(y ~ x + I(x^2) + I(x^3))
abs.resid <- abs(cubic.model$resid)
smooth.call <- smooth.spline(x, abs.resid, 
  nknots = 10)
lambda <- smooth.call$lambda
df <- smooth.call$df
mad.train.fun <- function(x, y, out = NULL){
  smooth.spline(x[, 1], y, lambda = lambda, 
  df = df, nknots = 10)
}
p2.tibs <- conformal.pred(x = cbind(x,x^2,x^3), y = y, 
  x0 = cbind(x,x^2,x^3), 
  train.fun = train.fun, predict.fun = predict.fun, 
  mad.train.fun = mad.train.fun,
  mad.predict.fun = mad.predict.fun,
  alpha = alpha)
LSLWCI <- cbind(p2.tibs$lo, p2.tibs$up)  

## parametric conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], paraCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], paraCI[ix, 2], col = "red", lwd = 1.5)
axis(2)
mtext("shape = 2", side = 2, cex = 1, line = 2)

## LSLW conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSLWCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSLWCI[ix, 2], col = "red", lwd = 1.5)

## nonparametric conformal
plot.nonparametric(region = nonparaCI, x = x, 
  y = y, bins = bins)

## LS conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSCI[ix, 2], col = "red", lwd = 1.5)


############ shape = 10 ################
## generate gamma regression model 
beta <- c(1.25, -1); shape = 10
rate <- (cbind(1, x) %*% beta) * shape
y <- rgamma(n = n, shape = shape, rate = rate)
y <- y / sd(y)
data <- data.frame(y = y, x = x)
colnames(data)[2:(p+1)] <- paste("x", 1:p, sep = "")

## fit the Gamma regression model
fit <- glm(y ~ x1, family = "Gamma", data = data)
formula <- fit$formula
newdata <- data
respname <- all.vars(formula)[1]
newdata <- newdata[, !(colnames(data) %in% respname)]
newdata <- as.matrix(newdata)

## obtain parametric and nonparametric conformal 
## prediction regions
cpred <- conformal.glm(fit, parametric = TRUE, 
  nonparametric = FALSE, alpha = alpha,
  bins = bins, cores = cores)
paraCI <- cpred$paraconformal

cpred <- conformal.glm(fit, parametric = FALSE, 
  nonparametric = TRUE, alpha = alpha,
  bins = bins, cores = cores)
nonparaCI <- cpred$nonparaconformal

## obtain LS conformal prediction region
p1.tibs <- conformal.pred(x = cbind(x, x^2, x^3), y = y, 
  x0 = cbind(x, x^2, x^3), train.fun = train.fun, 
  predict.fun = predict.fun, alpha = alpha, 
  grid.method = "linear")
LSCI <- cbind(p1.tibs$lo, p1.tibs$up)

## obtain LSLW conformal prediction region
cubic.model <- lm(y ~ x + I(x^2) + I(x^3))
abs.resid <- abs(cubic.model$resid)
smooth.call <- smooth.spline(x, abs.resid, 
  nknots = 10)
lambda <- smooth.call$lambda
df <- smooth.call$df
mad.train.fun <- function(x, y, out = NULL){
  smooth.spline(x[, 1], y, lambda = lambda, 
  df = df, nknots = 10)
}
p2.tibs <- conformal.pred(x = cbind(x,x^2,x^3), y = y, 
  x0 = cbind(x,x^2,x^3), 
  train.fun = train.fun, predict.fun = predict.fun, 
  mad.train.fun = mad.train.fun,
  mad.predict.fun = mad.predict.fun,
  alpha = alpha)
LSLWCI <- cbind(p2.tibs$lo, p2.tibs$up)  

## parametric conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], paraCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], paraCI[ix, 2], col = "red", lwd = 1.5)
axis(2)
mtext("shape = 10", side = 2, cex = 1, line = 2)

## LSLW conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSLWCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSLWCI[ix, 2], col = "red", lwd = 1.5)

## nonparametric conformal
plot.nonparametric(region = nonparaCI, x = x, 
  y = y, bins = bins)

## LS conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSCI[ix, 2], col = "red", lwd = 1.5)


############ shape = 50 ################
## generate gamma regression model 
beta <- c(1.25, -1); shape = 50
rate <- (cbind(1, x) %*% beta) * shape
y <- rgamma(n = n, shape = shape, rate = rate)
y <- y / sd(y)
data <- data.frame(y = y, x = x)
colnames(data)[2:(p+1)] <- paste("x", 1:p, sep = "")

## fit the Gamma regression model
fit <- glm(y ~ x1, family = "Gamma", data = data)
formula <- fit$formula
newdata <- data
respname <- all.vars(formula)[1]
newdata <- newdata[, !(colnames(data) %in% respname)]
newdata <- as.matrix(newdata)

## obtain parametric and nonparametric conformal 
## prediction regions
cpred <- conformal.glm(fit, parametric = TRUE, 
  nonparametric = FALSE, alpha = alpha,
  bins = bins, cores = cores)
paraCI <- cpred$paraconformal

cpred <- conformal.glm(fit, parametric = FALSE, 
  nonparametric = TRUE, alpha = alpha,
  bins = bins, cores = cores)
nonparaCI <- cpred$nonparaconformal

## obtain LS conformal prediction region
p1.tibs <- conformal.pred(x = cbind(x, x^2, x^3), y = y, 
  x0 = cbind(x, x^2, x^3), train.fun = train.fun, 
  predict.fun = predict.fun, alpha = alpha, 
  grid.method = "linear")
LSCI <- cbind(p1.tibs$lo, p1.tibs$up)

## obtain LSLW conformal prediction region
cubic.model <- lm(y ~ x + I(x^2) + I(x^3))
abs.resid <- abs(cubic.model$resid)
smooth.call <- smooth.spline(x, abs.resid, 
  nknots = 10)
lambda <- smooth.call$lambda
df <- smooth.call$df
mad.train.fun <- function(x, y, out = NULL){
  smooth.spline(x[, 1], y, lambda = lambda, 
  df = df, nknots = 10)
}
p2.tibs <- conformal.pred(x = cbind(x,x^2,x^3), y = y, 
  x0 = cbind(x,x^2,x^3), 
  train.fun = train.fun, predict.fun = predict.fun, 
  mad.train.fun = mad.train.fun,
  mad.predict.fun = mad.predict.fun,
  alpha = alpha)
LSLWCI <- cbind(p2.tibs$lo, p2.tibs$up)  

## parametric conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], paraCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], paraCI[ix, 2], col = "red", lwd = 1.5)
axis(1);axis(2)
mtext("shape = 50", side = 2, cex = 1, line = 2)

## LSLW conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSLWCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSLWCI[ix, 2], col = "red", lwd = 1.5)
axis(1)

## nonparametric conformal
plot.nonparametric(region = nonparaCI, x = x, 
  y = y, bins = bins)
axis(1)

## LS conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSCI[ix, 2], col = "red", lwd = 1.5)
axis(1)

mtext("x", side = 1, cex = 1, line = 2, outer = TRUE)
mtext("y", side = 2, cex = 1, line = 2, outer = TRUE)
@
\end{center}
\caption{The depiction of conformal prediction regions under simulation 
  setting A when $n = 250$ and the number of bins equals 3.
}
\label{conformal-plots-A-250}
\end{figure}








\newpage
%Plot of conformal prediciton regions in sim setting A when n = 500
\begin{figure}[h!]
\begin{center}
<<conformal-plots-A-500, fig.height = 7, cache = TRUE, echo = FALSE>>=
par(mfrow = c(4, 4), oma = c(4,4,2,0), mar = c(1,2,1,1))

############ shape = 0.75 ################
## generate random data
set.seed(13)
alpha <- 0.10
cores <- 6
n <- 500
bins <- 3
p <- k <- length(beta) - 1
x <- matrix(runif(n), ncol = p)

## generate gamma regression model 
beta <- c(1.25, -1); shape = 0.75
rate <- (cbind(1, x) %*% beta) * shape
y <- rgamma(n = n, shape = shape, rate = rate)
y <- y / sd(y)
data <- data.frame(y = y, x = x)
colnames(data)[2:(p+1)] <- paste("x", 1:p, sep = "")

## fit the Gamma regression model
fit <- glm(y ~ x1, family = "Gamma", data = data)
formula <- fit$formula
newdata <- data
respname <- all.vars(formula)[1]
newdata <- newdata[, !(colnames(data) %in% respname)]
newdata <- as.matrix(newdata)

## obtain parametric and nonparametric conformal 
## prediction regions
cpred <- conformal.glm(fit, parametric = TRUE, 
  nonparametric = FALSE, alpha = alpha,
  bins = bins, cores = cores)
paraCI <- cpred$paraconformal

cpred <- conformal.glm(fit, parametric = FALSE, 
  nonparametric = TRUE, alpha = alpha,
  bins = bins, cores = cores)
nonparaCI <- cpred$nonparaconformal

## obtain LS conformal prediction region
p1.tibs <- conformal.pred(x = cbind(x, x^2, x^3), y = y, 
  x0 = cbind(x, x^2, x^3), train.fun = train.fun, 
  predict.fun = predict.fun, alpha = alpha, 
  grid.method = "linear")
LSCI <- cbind(p1.tibs$lo, p1.tibs$up)

## obtain LSLW conformal prediction region
cubic.model <- lm(y ~ x + I(x^2) + I(x^3))
abs.resid <- abs(cubic.model$resid)
smooth.call <- smooth.spline(x, abs.resid, 
  nknots = 10)
lambda <- smooth.call$lambda
df <- smooth.call$df
mad.train.fun <- function(x, y, out = NULL){
  smooth.spline(x[, 1], y, lambda = lambda, 
  df = df, nknots = 10)
}
p2.tibs <- conformal.pred(x = cbind(x,x^2,x^3), y = y, 
  x0 = cbind(x,x^2,x^3), 
  train.fun = train.fun, predict.fun = predict.fun, 
  mad.train.fun = mad.train.fun,
  mad.predict.fun = mad.predict.fun,
  alpha = alpha)
LSLWCI <- cbind(p2.tibs$lo, p2.tibs$up)  

## create plots
ix <- sort(x, index.return = TRUE)$ix

## parametric conformal
plot.new()
plot.window(xlim = c(0,1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], paraCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], paraCI[ix, 2], col = "red", lwd = 1.5)
axis(2)
mtext("Parametric\n conformal", side = 3, cex = 1, line = 0)
mtext("shape = 0.75", side = 2, cex = 1, line = 2)

## LSLW conformal
plot.new()
plot.window(xlim = c(0,1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSLWCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSLWCI[ix, 2], col = "red", lwd = 1.5)
mtext("LSLW\n conformal", side = 3, cex = 1, line = 0)

## nonparametric conformal
plot.nonparametric(region = nonparaCI, x = x, 
  y = y, bins = bins)
mtext("Nonparametric\n conformal", side = 3, cex = 1, line = 0)

## LS conformal
plot.new()
plot.window(xlim = c(0,1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSCI[ix, 2], col = "red", lwd = 1.5)
mtext("LS\n conformal", side = 3, cex = 1, line = 0)


############ shape = 2 ################
## generate gamma regression model 
beta <- c(1.25, -1); shape = 2
rate <- (cbind(1, x) %*% beta) * shape
y <- rgamma(n = n, shape = shape, rate = rate)
y <- y / sd(y)
data <- data.frame(y = y, x = x)
colnames(data)[2:(p+1)] <- paste("x", 1:p, sep = "")

## fit the Gamma regression model
fit <- glm(y ~ x1, family = "Gamma", data = data)
formula <- fit$formula
newdata <- data
respname <- all.vars(formula)[1]
newdata <- newdata[, !(colnames(data) %in% respname)]
newdata <- as.matrix(newdata)

## obtain parametric and nonparametric conformal 
## prediction regions
cpred <- conformal.glm(fit, parametric = TRUE, 
  nonparametric = FALSE, alpha = alpha,
  bins = bins, cores = cores)
paraCI <- cpred$paraconformal

cpred <- conformal.glm(fit, parametric = FALSE, 
  nonparametric = TRUE, alpha = alpha,
  bins = bins, cores = cores)
nonparaCI <- cpred$nonparaconformal

## obtain LS conformal prediction region
p1.tibs <- conformal.pred(x = cbind(x, x^2, x^3), y = y, 
  x0 = cbind(x, x^2, x^3), train.fun = train.fun, 
  predict.fun = predict.fun, alpha = alpha, 
  grid.method = "linear")
LSCI <- cbind(p1.tibs$lo, p1.tibs$up)

## obtain LSLW conformal prediction region
cubic.model <- lm(y ~ x + I(x^2) + I(x^3))
abs.resid <- abs(cubic.model$resid)
smooth.call <- smooth.spline(x, abs.resid, 
  nknots = 10)
lambda <- smooth.call$lambda
df <- smooth.call$df
mad.train.fun <- function(x, y, out = NULL){
  smooth.spline(x[, 1], y, lambda = lambda, 
  df = df, nknots = 10)
}
p2.tibs <- conformal.pred(x = cbind(x,x^2,x^3), y = y, 
  x0 = cbind(x,x^2,x^3), 
  train.fun = train.fun, predict.fun = predict.fun, 
  mad.train.fun = mad.train.fun,
  mad.predict.fun = mad.predict.fun,
  alpha = alpha)
LSLWCI <- cbind(p2.tibs$lo, p2.tibs$up)  

## parametric conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], paraCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], paraCI[ix, 2], col = "red", lwd = 1.5)
axis(2)
mtext("shape = 2", side = 2, cex = 1, line = 2)

## LSLW conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSLWCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSLWCI[ix, 2], col = "red", lwd = 1.5)

## nonparametric conformal
plot.nonparametric(region = nonparaCI, x = x, 
  y = y, bins = bins)

## LS conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSCI[ix, 2], col = "red", lwd = 1.5)


############ shape = 10 ################
## generate gamma regression model 
beta <- c(1.25, -1); shape = 10
rate <- (cbind(1, x) %*% beta) * shape
y <- rgamma(n = n, shape = shape, rate = rate)
y <- y / sd(y)
data <- data.frame(y = y, x = x)
colnames(data)[2:(p+1)] <- paste("x", 1:p, sep = "")

## fit the Gamma regression model
fit <- glm(y ~ x1, family = "Gamma", data = data)
formula <- fit$formula
newdata <- data
respname <- all.vars(formula)[1]
newdata <- newdata[, !(colnames(data) %in% respname)]
newdata <- as.matrix(newdata)

## obtain parametric and nonparametric conformal 
## prediction regions
cpred <- conformal.glm(fit, parametric = TRUE, 
  nonparametric = FALSE, alpha = alpha,
  bins = bins, cores = cores)
paraCI <- cpred$paraconformal

cpred <- conformal.glm(fit, parametric = FALSE, 
  nonparametric = TRUE, alpha = alpha,
  bins = bins, cores = cores)
nonparaCI <- cpred$nonparaconformal

## obtain LS conformal prediction region
p1.tibs <- conformal.pred(x = cbind(x, x^2, x^3), y = y, 
  x0 = cbind(x, x^2, x^3), train.fun = train.fun, 
  predict.fun = predict.fun, alpha = alpha, 
  grid.method = "linear")
LSCI <- cbind(p1.tibs$lo, p1.tibs$up)

## obtain LSLW conformal prediction region
cubic.model <- lm(y ~ x + I(x^2) + I(x^3))
abs.resid <- abs(cubic.model$resid)
smooth.call <- smooth.spline(x, abs.resid, 
  nknots = 10)
lambda <- smooth.call$lambda
df <- smooth.call$df
mad.train.fun <- function(x, y, out = NULL){
  smooth.spline(x[, 1], y, lambda = lambda, 
  df = df, nknots = 10)
}
p2.tibs <- conformal.pred(x = cbind(x,x^2,x^3), y = y, 
  x0 = cbind(x,x^2,x^3), 
  train.fun = train.fun, predict.fun = predict.fun, 
  mad.train.fun = mad.train.fun,
  mad.predict.fun = mad.predict.fun,
  alpha = alpha)
LSLWCI <- cbind(p2.tibs$lo, p2.tibs$up)  

## parametric conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], paraCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], paraCI[ix, 2], col = "red", lwd = 1.5)
axis(2)
mtext("shape = 10", side = 2, cex = 1, line = 2)

## LSLW conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSLWCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSLWCI[ix, 2], col = "red", lwd = 1.5)

## nonparametric conformal
plot.nonparametric(region = nonparaCI, x = x, 
  y = y, bins = bins)

## LS conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSCI[ix, 2], col = "red", lwd = 1.5)


############ shape = 50 ################
## generate gamma regression model 
beta <- c(1.25, -1); shape = 50
rate <- (cbind(1, x) %*% beta) * shape
y <- rgamma(n = n, shape = shape, rate = rate)
y <- y / sd(y)
data <- data.frame(y = y, x = x)
colnames(data)[2:(p+1)] <- paste("x", 1:p, sep = "")

## fit the Gamma regression model
fit <- glm(y ~ x1, family = "Gamma", data = data)
formula <- fit$formula
newdata <- data
respname <- all.vars(formula)[1]
newdata <- newdata[, !(colnames(data) %in% respname)]
newdata <- as.matrix(newdata)

## obtain parametric and nonparametric conformal 
## prediction regions
cpred <- conformal.glm(fit, parametric = TRUE, 
  nonparametric = FALSE, alpha = alpha,
  bins = bins, cores = cores)
paraCI <- cpred$paraconformal

cpred <- conformal.glm(fit, parametric = FALSE, 
  nonparametric = TRUE, alpha = alpha,
  bins = bins, cores = cores)
nonparaCI <- cpred$nonparaconformal

## obtain LS conformal prediction region
p1.tibs <- conformal.pred(x = cbind(x, x^2, x^3), y = y, 
  x0 = cbind(x, x^2, x^3), train.fun = train.fun, 
  predict.fun = predict.fun, alpha = alpha, 
  grid.method = "linear")
LSCI <- cbind(p1.tibs$lo, p1.tibs$up)

## obtain LSLW conformal prediction region
cubic.model <- lm(y ~ x + I(x^2) + I(x^3))
abs.resid <- abs(cubic.model$resid)
smooth.call <- smooth.spline(x, abs.resid, 
  nknots = 10)
lambda <- smooth.call$lambda
df <- smooth.call$df
mad.train.fun <- function(x, y, out = NULL){
  smooth.spline(x[, 1], y, lambda = lambda, 
  df = df, nknots = 10)
}
p2.tibs <- conformal.pred(x = cbind(x,x^2,x^3), y = y, 
  x0 = cbind(x,x^2,x^3), 
  train.fun = train.fun, predict.fun = predict.fun, 
  mad.train.fun = mad.train.fun,
  mad.predict.fun = mad.predict.fun,
  alpha = alpha)
LSLWCI <- cbind(p2.tibs$lo, p2.tibs$up)  

## parametric conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], paraCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], paraCI[ix, 2], col = "red", lwd = 1.5)
axis(1);axis(2)
mtext("shape = 50", side = 2, cex = 1, line = 2)

## LSLW conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSLWCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSLWCI[ix, 2], col = "red", lwd = 1.5)
axis(1)

## nonparametric conformal
plot.nonparametric(region = nonparaCI, x = x, 
  y = y, bins = bins)
axis(1)

## LS conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSCI[ix, 2], col = "red", lwd = 1.5)
axis(1)

mtext("x", side = 1, cex = 1, line = 2, outer = TRUE)
mtext("y", side = 2, cex = 1, line = 2, outer = TRUE)
@
\end{center}
\caption{The depiction of conformal prediction regions under simulation 
  setting A when $n = 500$ and the number of bins equals 3.
}
\label{conformal-plots-A-500}
\end{figure}







\newpage
\subsection{Plots corresponding to Section~\ref{sec:misspec}}
\label{sec:misspecplots}


%Plot of conformal prediciton regions in sim setting B when n = 150
\begin{figure}[h!]
\begin{center}
<<conformal-plots-B-150, fig.height = 7, cache = TRUE, echo = FALSE>>=
par(mfrow = c(4, 4), oma = c(4,4,2,0), mar = c(1,2,1,1))

############ shape = 0.75 ################
## generate random data
set.seed(13)
alpha <- 0.10
cores <- 6
n <- 150
bins <- 2
p <- k <- length(beta) - 1
x <- matrix(runif(n), ncol = p)

## generate gamma regression model 
beta <- c(0.5, 1); shape = 0.75
rate <- (cbind(1, x) %*% beta) * shape
y <- rgamma(n = n, shape = shape, rate = rate)
y <- y / sd(y)
data <- data.frame(y = y, x = x)
colnames(data)[2:(p+1)] <- paste("x", 1:p, sep = "")

## fit the misspecified regression model
fit <- glm(y ~ x1 + I(x1^2) + I(x1^3), 
  family = "gaussian", data = data)
formula <- fit$formula
newdata <- data
respname <- all.vars(formula)[1]
newdata <- newdata[, !(colnames(data) %in% respname)]
newdata <- as.matrix(newdata)

## obtain parametric and nonparametric conformal 
## prediction regions
cpred <- conformal.glm(fit, parametric = TRUE, 
  nonparametric = FALSE, alpha = alpha,
  bins = bins, cores = cores)
paraCI <- cpred$paraconformal

cpred <- conformal.glm(fit, parametric = FALSE, 
  nonparametric = TRUE, alpha = alpha,
  bins = bins, cores = cores)
nonparaCI <- cpred$nonparaconformal

## obtain LS conformal prediction region
p1.tibs <- conformal.pred(x = cbind(x, x^2, x^3), y = y, 
  x0 = cbind(x, x^2, x^3), train.fun = train.fun, 
  predict.fun = predict.fun, alpha = alpha, 
  grid.method = "linear")
LSCI <- cbind(p1.tibs$lo, p1.tibs$up)

## obtain LSLW conformal prediction region
cubic.model <- lm(y ~ x + I(x^2) + I(x^3))
abs.resid <- abs(cubic.model$resid)
smooth.call <- smooth.spline(x, abs.resid, 
  nknots = 10)
lambda <- smooth.call$lambda
df <- smooth.call$df
mad.train.fun <- function(x, y, out = NULL){
  smooth.spline(x[, 1], y, lambda = lambda, 
  df = df, nknots = 10)
}
p2.tibs <- conformal.pred(x = cbind(x,x^2,x^3), y = y, 
  x0 = cbind(x,x^2,x^3), 
  train.fun = train.fun, predict.fun = predict.fun, 
  mad.train.fun = mad.train.fun,
  mad.predict.fun = mad.predict.fun,
  alpha = alpha)
LSLWCI <- cbind(p2.tibs$lo, p2.tibs$up)  

## create plots
ix <- sort(x, index.return = TRUE)$ix

## parametric conformal
plot.new()
plot.window(xlim = c(0,1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], paraCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], paraCI[ix, 2], col = "red", lwd = 1.5)
axis(2)
mtext("Parametric\n conformal", side = 3, cex = 1, line = 0)
mtext("shape = 0.75", side = 2, cex = 1, line = 2)

## LSLW conformal
plot.new()
plot.window(xlim = c(0,1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSLWCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSLWCI[ix, 2], col = "red", lwd = 1.5)
mtext("LSLW\n conformal", side = 3, cex = 1, line = 0)

## nonparametric conformal
plot.nonparametric(region = nonparaCI, x = x, 
  y = y, bins = bins)
mtext("Nonparametric\n conformal", side = 3, cex = 1, line = 0)

## LS conformal
plot.new()
plot.window(xlim = c(0,1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSCI[ix, 2], col = "red", lwd = 1.5)
mtext("LS\n conformal", side = 3, cex = 1, line = 0)


############ shape = 2 ################
## generate gamma regression model 
beta <- c(0.5, 1); shape = 2
rate <- (cbind(1, x) %*% beta) * shape
y <- rgamma(n = n, shape = shape, rate = rate)
y <- y / sd(y)
data <- data.frame(y = y, x = x)
colnames(data)[2:(p+1)] <- paste("x", 1:p, sep = "")

## fit the misspecified regression model
fit <- glm(y ~ x1 + I(x1^2) + I(x1^3), 
  family = "gaussian", data = data)
formula <- fit$formula
newdata <- data
respname <- all.vars(formula)[1]
newdata <- newdata[, !(colnames(data) %in% respname)]
newdata <- as.matrix(newdata)

## obtain parametric and nonparametric conformal 
## prediction regions
cpred <- conformal.glm(fit, parametric = TRUE, 
  nonparametric = FALSE, alpha = alpha,
  bins = bins, cores = cores)
paraCI <- cpred$paraconformal

cpred <- conformal.glm(fit, parametric = FALSE, 
  nonparametric = TRUE, alpha = alpha,
  bins = bins, cores = cores)
nonparaCI <- cpred$nonparaconformal

## obtain LS conformal prediction region
p1.tibs <- conformal.pred(x = cbind(x, x^2, x^3), y = y, 
  x0 = cbind(x, x^2, x^3), train.fun = train.fun, 
  predict.fun = predict.fun, alpha = alpha, 
  grid.method = "linear")
LSCI <- cbind(p1.tibs$lo, p1.tibs$up)

## obtain LSLW conformal prediction region
cubic.model <- lm(y ~ x + I(x^2) + I(x^3))
abs.resid <- abs(cubic.model$resid)
smooth.call <- smooth.spline(x, abs.resid, 
  nknots = 10)
lambda <- smooth.call$lambda
df <- smooth.call$df
mad.train.fun <- function(x, y, out = NULL){
  smooth.spline(x[, 1], y, lambda = lambda, 
  df = df, nknots = 10)
}
p2.tibs <- conformal.pred(x = cbind(x,x^2,x^3), y = y, 
  x0 = cbind(x,x^2,x^3), 
  train.fun = train.fun, predict.fun = predict.fun, 
  mad.train.fun = mad.train.fun,
  mad.predict.fun = mad.predict.fun,
  alpha = alpha)
LSLWCI <- cbind(p2.tibs$lo, p2.tibs$up)  

## create plots
ix <- sort(x, index.return = TRUE)$ix

## parametric conformal
plot.new()
plot.window(xlim = c(0,1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], paraCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], paraCI[ix, 2], col = "red", lwd = 1.5)
axis(2)
mtext("shape = 2", side = 2, cex = 1, line = 2)

## LSLW conformal
plot.new()
plot.window(xlim = c(0,1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSLWCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSLWCI[ix, 2], col = "red", lwd = 1.5)

## nonparametric conformal
plot.nonparametric(region = nonparaCI, x = x, 
  y = y, bins = bins)

## LS conformal
plot.new()
plot.window(xlim = c(0,1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSCI[ix, 2], col = "red", lwd = 1.5)


############ shape = 10 ################
## generate gamma regression model 
beta <- c(0.5, 1); shape = 10
rate <- (cbind(1, x) %*% beta) * shape
y <- rgamma(n = n, shape = shape, rate = rate)
y <- y / sd(y)
data <- data.frame(y = y, x = x)
colnames(data)[2:(p+1)] <- paste("x", 1:p, sep = "")

## fit the misspecified regression model
fit <- glm(y ~ x1 + I(x1^2) + I(x1^3), 
  family = "gaussian", data = data)
formula <- fit$formula
newdata <- data
respname <- all.vars(formula)[1]
newdata <- newdata[, !(colnames(data) %in% respname)]
newdata <- as.matrix(newdata)

## obtain parametric and nonparametric conformal 
## prediction regions
cpred <- conformal.glm(fit, parametric = TRUE, 
  nonparametric = FALSE, alpha = alpha,
  bins = bins, cores = cores)
paraCI <- cpred$paraconformal

cpred <- conformal.glm(fit, parametric = FALSE, 
  nonparametric = TRUE, alpha = alpha,
  bins = bins, cores = cores)
nonparaCI <- cpred$nonparaconformal

## obtain LS conformal prediction region
p1.tibs <- conformal.pred(x = cbind(x, x^2, x^3), y = y, 
  x0 = cbind(x, x^2, x^3), train.fun = train.fun, 
  predict.fun = predict.fun, alpha = alpha, 
  grid.method = "linear")
LSCI <- cbind(p1.tibs$lo, p1.tibs$up)

## obtain LSLW conformal prediction region
cubic.model <- lm(y ~ x + I(x^2) + I(x^3))
abs.resid <- abs(cubic.model$resid)
smooth.call <- smooth.spline(x, abs.resid, 
  nknots = 10)
lambda <- smooth.call$lambda
df <- smooth.call$df
mad.train.fun <- function(x, y, out = NULL){
  smooth.spline(x[, 1], y, lambda = lambda, 
  df = df, nknots = 10)
}
p2.tibs <- conformal.pred(x = cbind(x,x^2,x^3), y = y, 
  x0 = cbind(x,x^2,x^3), 
  train.fun = train.fun, predict.fun = predict.fun, 
  mad.train.fun = mad.train.fun,
  mad.predict.fun = mad.predict.fun,
  alpha = alpha)
LSLWCI <- cbind(p2.tibs$lo, p2.tibs$up)  

## create plots
ix <- sort(x, index.return = TRUE)$ix

## parametric conformal
plot.new()
plot.window(xlim = c(0,1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], paraCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], paraCI[ix, 2], col = "red", lwd = 1.5)
axis(2)
mtext("shape = 10", side = 2, cex = 1, line = 2)

## LSLW conformal
plot.new()
plot.window(xlim = c(0,1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSLWCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSLWCI[ix, 2], col = "red", lwd = 1.5)

## nonparametric conformal
plot.nonparametric(region = nonparaCI, x = x, 
  y = y, bins = bins)

## LS conformal
plot.new()
plot.window(xlim = c(0,1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSCI[ix, 2], col = "red", lwd = 1.5)


############ shape = 50 ################
## generate gamma regression model 
beta <- c(0.5, 1); shape = 50
rate <- (cbind(1, x) %*% beta) * shape
y <- rgamma(n = n, shape = shape, rate = rate)
y <- y / sd(y)
data <- data.frame(y = y, x = x)
colnames(data)[2:(p+1)] <- paste("x", 1:p, sep = "")

## fit the misspecified regression model
fit <- glm(y ~ x1 + I(x1^2) + I(x1^3), 
  family = "gaussian", data = data)
formula <- fit$formula
newdata <- data
respname <- all.vars(formula)[1]
newdata <- newdata[, !(colnames(data) %in% respname)]
newdata <- as.matrix(newdata)

## obtain parametric and nonparametric conformal 
## prediction regions
cpred <- conformal.glm(fit, parametric = TRUE, 
  nonparametric = FALSE, alpha = alpha,
  bins = bins, cores = cores)
paraCI <- cpred$paraconformal

cpred <- conformal.glm(fit, parametric = FALSE, 
  nonparametric = TRUE, alpha = alpha,
  bins = bins, cores = cores)
nonparaCI <- cpred$nonparaconformal

## obtain LS conformal prediction region
p1.tibs <- conformal.pred(x = cbind(x, x^2, x^3), y = y, 
  x0 = cbind(x, x^2, x^3), train.fun = train.fun, 
  predict.fun = predict.fun, alpha = alpha, 
  grid.method = "linear")
LSCI <- cbind(p1.tibs$lo, p1.tibs$up)

## obtain LSLW conformal prediction region
cubic.model <- lm(y ~ x + I(x^2) + I(x^3))
abs.resid <- abs(cubic.model$resid)
smooth.call <- smooth.spline(x, abs.resid, 
  nknots = 10)
lambda <- smooth.call$lambda
df <- smooth.call$df
mad.train.fun <- function(x, y, out = NULL){
  smooth.spline(x[, 1], y, lambda = lambda, 
  df = df, nknots = 10)
}
p2.tibs <- conformal.pred(x = cbind(x,x^2,x^3), y = y, 
  x0 = cbind(x,x^2,x^3), 
  train.fun = train.fun, predict.fun = predict.fun, 
  mad.train.fun = mad.train.fun,
  mad.predict.fun = mad.predict.fun,
  alpha = alpha)
LSLWCI <- cbind(p2.tibs$lo, p2.tibs$up)  

## parametric conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], paraCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], paraCI[ix, 2], col = "red", lwd = 1.5)
axis(1);axis(2)
mtext("shape = 50", side = 2, cex = 1, line = 2)

## LSLW conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSLWCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSLWCI[ix, 2], col = "red", lwd = 1.5)
axis(1)

## nonparametric conformal
plot.nonparametric(region = nonparaCI, x = x, 
  y = y, bins = bins)
axis(1)

## LS conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSCI[ix, 2], col = "red", lwd = 1.5)
axis(1)

mtext("x", side = 1, cex = 1, line = 2, outer = TRUE)
mtext("y", side = 2, cex = 1, line = 2, outer = TRUE)
@
\end{center}
\caption{The depiction of conformal prediction regions under simulation 
  setting B when $n = 150$ and the number of bins equals 2.
}
\label{conformal-plots-B-150}
\end{figure}










\newpage
%Plot of conformal prediciton regions in sim setting B when n = 250
\begin{figure}[h!]
\begin{center}
<<conformal-plots-B-250, fig.height = 7, cache = TRUE, echo = FALSE>>=
par(mfrow = c(4, 4), oma = c(4,4,2,0), mar = c(1,2,1,1))

############ shape = 0.75 ################
## generate random data
set.seed(13)
alpha <- 0.10
cores <- 6
n <- 250
bins <- 3
p <- k <- length(beta) - 1
x <- matrix(runif(n), ncol = p)

## generate gamma regression model 
beta <- c(0.5, 1); shape = 0.75
rate <- (cbind(1, x) %*% beta) * shape
y <- rgamma(n = n, shape = shape, rate = rate)
y <- y / sd(y)
data <- data.frame(y = y, x = x)
colnames(data)[2:(p+1)] <- paste("x", 1:p, sep = "")

## fit the misspecified regression model
fit <- glm(y ~ x1 + I(x1^2) + I(x1^3), 
  family = "gaussian", data = data)
formula <- fit$formula
newdata <- data
respname <- all.vars(formula)[1]
newdata <- newdata[, !(colnames(data) %in% respname)]
newdata <- as.matrix(newdata)

## obtain parametric and nonparametric conformal 
## prediction regions
cpred <- conformal.glm(fit, parametric = TRUE, 
  nonparametric = FALSE, alpha = alpha,
  bins = bins, cores = cores)
paraCI <- cpred$paraconformal

cpred <- conformal.glm(fit, parametric = FALSE, 
  nonparametric = TRUE, alpha = alpha,
  bins = bins, cores = cores)
nonparaCI <- cpred$nonparaconformal

## obtain LS conformal prediction region
p1.tibs <- conformal.pred(x = cbind(x, x^2, x^3), y = y, 
  x0 = cbind(x, x^2, x^3), train.fun = train.fun, 
  predict.fun = predict.fun, alpha = alpha, 
  grid.method = "linear")
LSCI <- cbind(p1.tibs$lo, p1.tibs$up)

## obtain LSLW conformal prediction region
cubic.model <- lm(y ~ x + I(x^2) + I(x^3))
abs.resid <- abs(cubic.model$resid)
smooth.call <- smooth.spline(x, abs.resid, 
  nknots = 10)
lambda <- smooth.call$lambda
df <- smooth.call$df
mad.train.fun <- function(x, y, out = NULL){
  smooth.spline(x[, 1], y, lambda = lambda, 
  df = df, nknots = 10)
}
p2.tibs <- conformal.pred(x = cbind(x,x^2,x^3), y = y, 
  x0 = cbind(x,x^2,x^3), 
  train.fun = train.fun, predict.fun = predict.fun, 
  mad.train.fun = mad.train.fun,
  mad.predict.fun = mad.predict.fun,
  alpha = alpha)
LSLWCI <- cbind(p2.tibs$lo, p2.tibs$up)  

## create plots
ix <- sort(x, index.return = TRUE)$ix

## parametric conformal
plot.new()
plot.window(xlim = c(0,1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], paraCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], paraCI[ix, 2], col = "red", lwd = 1.5)
axis(2)
mtext("Parametric\n conformal", side = 3, cex = 1, line = 0)
mtext("shape = 0.75", side = 2, cex = 1, line = 2)

## LSLW conformal
plot.new()
plot.window(xlim = c(0,1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSLWCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSLWCI[ix, 2], col = "red", lwd = 1.5)
mtext("LSLW\n conformal", side = 3, cex = 1, line = 0)

## nonparametric conformal
plot.nonparametric(region = nonparaCI, x = x, 
  y = y, bins = bins)
mtext("Nonparametric\n conformal", side = 3, cex = 1, line = 0)

## LS conformal
plot.new()
plot.window(xlim = c(0,1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSCI[ix, 2], col = "red", lwd = 1.5)
mtext("LS\n conformal", side = 3, cex = 1, line = 0)


############ shape = 2 ################
## generate gamma regression model 
beta <- c(0.5, 1); shape = 2
rate <- (cbind(1, x) %*% beta) * shape
y <- rgamma(n = n, shape = shape, rate = rate)
y <- y / sd(y)
data <- data.frame(y = y, x = x)
colnames(data)[2:(p+1)] <- paste("x", 1:p, sep = "")

## fit the misspecified regression model
fit <- glm(y ~ x1 + I(x1^2) + I(x1^3), 
  family = "gaussian", data = data)
formula <- fit$formula
newdata <- data
respname <- all.vars(formula)[1]
newdata <- newdata[, !(colnames(data) %in% respname)]
newdata <- as.matrix(newdata)

## obtain parametric and nonparametric conformal 
## prediction regions
cpred <- conformal.glm(fit, parametric = TRUE, 
  nonparametric = FALSE, alpha = alpha,
  bins = bins, cores = cores)
paraCI <- cpred$paraconformal

cpred <- conformal.glm(fit, parametric = FALSE, 
  nonparametric = TRUE, alpha = alpha,
  bins = bins, cores = cores)
nonparaCI <- cpred$nonparaconformal

## obtain LS conformal prediction region
p1.tibs <- conformal.pred(x = cbind(x, x^2, x^3), y = y, 
  x0 = cbind(x, x^2, x^3), train.fun = train.fun, 
  predict.fun = predict.fun, alpha = alpha, 
  grid.method = "linear")
LSCI <- cbind(p1.tibs$lo, p1.tibs$up)

## obtain LSLW conformal prediction region
cubic.model <- lm(y ~ x + I(x^2) + I(x^3))
abs.resid <- abs(cubic.model$resid)
smooth.call <- smooth.spline(x, abs.resid, 
  nknots = 10)
lambda <- smooth.call$lambda
df <- smooth.call$df
mad.train.fun <- function(x, y, out = NULL){
  smooth.spline(x[, 1], y, lambda = lambda, 
  df = df, nknots = 10)
}
p2.tibs <- conformal.pred(x = cbind(x,x^2,x^3), y = y, 
  x0 = cbind(x,x^2,x^3), 
  train.fun = train.fun, predict.fun = predict.fun, 
  mad.train.fun = mad.train.fun,
  mad.predict.fun = mad.predict.fun,
  alpha = alpha)
LSLWCI <- cbind(p2.tibs$lo, p2.tibs$up)  

## create plots
ix <- sort(x, index.return = TRUE)$ix

## parametric conformal
plot.new()
plot.window(xlim = c(0,1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], paraCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], paraCI[ix, 2], col = "red", lwd = 1.5)
axis(2)
mtext("shape = 2", side = 2, cex = 1, line = 2)

## LSLW conformal
plot.new()
plot.window(xlim = c(0,1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSLWCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSLWCI[ix, 2], col = "red", lwd = 1.5)

## nonparametric conformal
plot.nonparametric(region = nonparaCI, x = x, 
  y = y, bins = bins)

## LS conformal
plot.new()
plot.window(xlim = c(0,1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSCI[ix, 2], col = "red", lwd = 1.5)


############ shape = 10 ################
## generate gamma regression model 
beta <- c(0.5, 1); shape = 10
rate <- (cbind(1, x) %*% beta) * shape
y <- rgamma(n = n, shape = shape, rate = rate)
y <- y / sd(y)
data <- data.frame(y = y, x = x)
colnames(data)[2:(p+1)] <- paste("x", 1:p, sep = "")

## fit the misspecified regression model
fit <- glm(y ~ x1 + I(x1^2) + I(x1^3), 
  family = "gaussian", data = data)
formula <- fit$formula
newdata <- data
respname <- all.vars(formula)[1]
newdata <- newdata[, !(colnames(data) %in% respname)]
newdata <- as.matrix(newdata)

## obtain parametric and nonparametric conformal 
## prediction regions
cpred <- conformal.glm(fit, parametric = TRUE, 
  nonparametric = FALSE, alpha = alpha,
  bins = bins, cores = cores)
paraCI <- cpred$paraconformal

cpred <- conformal.glm(fit, parametric = FALSE, 
  nonparametric = TRUE, alpha = alpha,
  bins = bins, cores = cores)
nonparaCI <- cpred$nonparaconformal

## obtain LS conformal prediction region
p1.tibs <- conformal.pred(x = cbind(x, x^2, x^3), y = y, 
  x0 = cbind(x, x^2, x^3), train.fun = train.fun, 
  predict.fun = predict.fun, alpha = alpha, 
  grid.method = "linear")
LSCI <- cbind(p1.tibs$lo, p1.tibs$up)

## obtain LSLW conformal prediction region
cubic.model <- lm(y ~ x + I(x^2) + I(x^3))
abs.resid <- abs(cubic.model$resid)
smooth.call <- smooth.spline(x, abs.resid, 
  nknots = 10)
lambda <- smooth.call$lambda
df <- smooth.call$df
mad.train.fun <- function(x, y, out = NULL){
  smooth.spline(x[, 1], y, lambda = lambda, 
  df = df, nknots = 10)
}
p2.tibs <- conformal.pred(x = cbind(x,x^2,x^3), y = y, 
  x0 = cbind(x,x^2,x^3), 
  train.fun = train.fun, predict.fun = predict.fun, 
  mad.train.fun = mad.train.fun,
  mad.predict.fun = mad.predict.fun,
  alpha = alpha)
LSLWCI <- cbind(p2.tibs$lo, p2.tibs$up)  

## create plots
ix <- sort(x, index.return = TRUE)$ix

## parametric conformal
plot.new()
plot.window(xlim = c(0,1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], paraCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], paraCI[ix, 2], col = "red", lwd = 1.5)
axis(2)
mtext("shape = 10", side = 2, cex = 1, line = 2)

## LSLW conformal
plot.new()
plot.window(xlim = c(0,1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSLWCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSLWCI[ix, 2], col = "red", lwd = 1.5)

## nonparametric conformal
plot.nonparametric(region = nonparaCI, x = x, 
  y = y, bins = bins)

## LS conformal
plot.new()
plot.window(xlim = c(0,1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSCI[ix, 2], col = "red", lwd = 1.5)


############ shape = 50 ################
## generate gamma regression model 
beta <- c(0.5, 1); shape = 50
rate <- (cbind(1, x) %*% beta) * shape
y <- rgamma(n = n, shape = shape, rate = rate)
y <- y / sd(y)
data <- data.frame(y = y, x = x)
colnames(data)[2:(p+1)] <- paste("x", 1:p, sep = "")

## fit the misspecified regression model
fit <- glm(y ~ x1 + I(x1^2) + I(x1^3), 
  family = "gaussian", data = data)
formula <- fit$formula
newdata <- data
respname <- all.vars(formula)[1]
newdata <- newdata[, !(colnames(data) %in% respname)]
newdata <- as.matrix(newdata)

## obtain parametric and nonparametric conformal 
## prediction regions
cpred <- conformal.glm(fit, parametric = TRUE, 
  nonparametric = FALSE, alpha = alpha,
  bins = bins, cores = cores)
paraCI <- cpred$paraconformal

cpred <- conformal.glm(fit, parametric = FALSE, 
  nonparametric = TRUE, alpha = alpha,
  bins = bins, cores = cores)
nonparaCI <- cpred$nonparaconformal

## obtain LS conformal prediction region
p1.tibs <- conformal.pred(x = cbind(x, x^2, x^3), y = y, 
  x0 = cbind(x, x^2, x^3), train.fun = train.fun, 
  predict.fun = predict.fun, alpha = alpha, 
  grid.method = "linear")
LSCI <- cbind(p1.tibs$lo, p1.tibs$up)

## obtain LSLW conformal prediction region
cubic.model <- lm(y ~ x + I(x^2) + I(x^3))
abs.resid <- abs(cubic.model$resid)
smooth.call <- smooth.spline(x, abs.resid, 
  nknots = 10)
lambda <- smooth.call$lambda
df <- smooth.call$df
mad.train.fun <- function(x, y, out = NULL){
  smooth.spline(x[, 1], y, lambda = lambda, 
  df = df, nknots = 10)
}
p2.tibs <- conformal.pred(x = cbind(x,x^2,x^3), y = y, 
  x0 = cbind(x,x^2,x^3), 
  train.fun = train.fun, predict.fun = predict.fun, 
  mad.train.fun = mad.train.fun,
  mad.predict.fun = mad.predict.fun,
  alpha = alpha)
LSLWCI <- cbind(p2.tibs$lo, p2.tibs$up)  

## parametric conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], paraCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], paraCI[ix, 2], col = "red", lwd = 1.5)
axis(1);axis(2)
mtext("shape = 50", side = 2, cex = 1, line = 2)

## LSLW conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSLWCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSLWCI[ix, 2], col = "red", lwd = 1.5)
axis(1)

## nonparametric conformal
plot.nonparametric(region = nonparaCI, x = x, 
  y = y, bins = bins)
axis(1)

## LS conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSCI[ix, 2], col = "red", lwd = 1.5)
axis(1)

mtext("x", side = 1, cex = 1, line = 2, outer = TRUE)
mtext("y", side = 2, cex = 1, line = 2, outer = TRUE)
@
\end{center}
\caption{The depiction of conformal prediction regions under simulation 
  setting B when $n = 250$ and the number of bins equals 3.
}
\label{conformal-plots-B-250}
\end{figure}









\newpage
%Plot of conformal prediciton regions in sim setting B when n = 500
\begin{figure}[h!]
\begin{center}
<<conformal-plots-B-500, fig.height = 7, cache = TRUE, echo = FALSE>>=
par(mfrow = c(4, 4), oma = c(4,4,2,0), mar = c(1,2,1,1))

############ shape = 0.75 ################
## generate random data
set.seed(13)
alpha <- 0.10
cores <- 6
n <- 500
bins <- 3
p <- k <- length(beta) - 1
x <- matrix(runif(n), ncol = p)

## generate gamma regression model 
beta <- c(0.5, 1); shape = 0.75
rate <- (cbind(1, x) %*% beta) * shape
y <- rgamma(n = n, shape = shape, rate = rate)
y <- y / sd(y)
data <- data.frame(y = y, x = x)
colnames(data)[2:(p+1)] <- paste("x", 1:p, sep = "")

## fit the misspecified regression model
fit <- glm(y ~ x1 + I(x1^2) + I(x1^3), 
  family = "gaussian", data = data)
formula <- fit$formula
newdata <- data
respname <- all.vars(formula)[1]
newdata <- newdata[, !(colnames(data) %in% respname)]
newdata <- as.matrix(newdata)

## obtain parametric and nonparametric conformal 
## prediction regions
cpred <- conformal.glm(fit, parametric = TRUE, 
  nonparametric = FALSE, alpha = alpha,
  bins = bins, cores = cores)
paraCI <- cpred$paraconformal

cpred <- conformal.glm(fit, parametric = FALSE, 
  nonparametric = TRUE, alpha = alpha,
  bins = bins, cores = cores)
nonparaCI <- cpred$nonparaconformal

## obtain LS conformal prediction region
p1.tibs <- conformal.pred(x = cbind(x, x^2, x^3), y = y, 
  x0 = cbind(x, x^2, x^3), train.fun = train.fun, 
  predict.fun = predict.fun, alpha = alpha, 
  grid.method = "linear")
LSCI <- cbind(p1.tibs$lo, p1.tibs$up)

## obtain LSLW conformal prediction region
cubic.model <- lm(y ~ x + I(x^2) + I(x^3))
abs.resid <- abs(cubic.model$resid)
smooth.call <- smooth.spline(x, abs.resid, 
  nknots = 10)
lambda <- smooth.call$lambda
df <- smooth.call$df
mad.train.fun <- function(x, y, out = NULL){
  smooth.spline(x[, 1], y, lambda = lambda, 
  df = df, nknots = 10)
}
p2.tibs <- conformal.pred(x = cbind(x,x^2,x^3), y = y, 
  x0 = cbind(x,x^2,x^3), 
  train.fun = train.fun, predict.fun = predict.fun, 
  mad.train.fun = mad.train.fun,
  mad.predict.fun = mad.predict.fun,
  alpha = alpha)
LSLWCI <- cbind(p2.tibs$lo, p2.tibs$up)  

## create plots
ix <- sort(x, index.return = TRUE)$ix

## parametric conformal
plot.new()
plot.window(xlim = c(0,1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], paraCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], paraCI[ix, 2], col = "red", lwd = 1.5)
axis(2)
mtext("Parametric\n conformal", side = 3, cex = 1, line = 0)
mtext("shape = 0.75", side = 2, cex = 1, line = 2)

## LSLW conformal
plot.new()
plot.window(xlim = c(0,1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSLWCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSLWCI[ix, 2], col = "red", lwd = 1.5)
mtext("LSLW\n conformal", side = 3, cex = 1, line = 0)

## nonparametric conformal
plot.nonparametric(region = nonparaCI, x = x, 
  y = y, bins = bins)
mtext("Nonparametric\n conformal", side = 3, cex = 1, line = 0)

## LS conformal
plot.new()
plot.window(xlim = c(0,1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSCI[ix, 2], col = "red", lwd = 1.5)
mtext("LS\n conformal", side = 3, cex = 1, line = 0)


############ shape = 2 ################
## generate gamma regression model 
beta <- c(0.5, 1); shape = 2
rate <- (cbind(1, x) %*% beta) * shape
y <- rgamma(n = n, shape = shape, rate = rate)
y <- y / sd(y)
data <- data.frame(y = y, x = x)
colnames(data)[2:(p+1)] <- paste("x", 1:p, sep = "")

## fit the misspecified regression model
fit <- glm(y ~ x1 + I(x1^2) + I(x1^3), 
  family = "gaussian", data = data)
formula <- fit$formula
newdata <- data
respname <- all.vars(formula)[1]
newdata <- newdata[, !(colnames(data) %in% respname)]
newdata <- as.matrix(newdata)

## obtain parametric and nonparametric conformal 
## prediction regions
cpred <- conformal.glm(fit, parametric = TRUE, 
  nonparametric = FALSE, alpha = alpha,
  bins = bins, cores = cores)
paraCI <- cpred$paraconformal

cpred <- conformal.glm(fit, parametric = FALSE, 
  nonparametric = TRUE, alpha = alpha,
  bins = bins, cores = cores)
nonparaCI <- cpred$nonparaconformal

## obtain LS conformal prediction region
p1.tibs <- conformal.pred(x = cbind(x, x^2, x^3), y = y, 
  x0 = cbind(x, x^2, x^3), train.fun = train.fun, 
  predict.fun = predict.fun, alpha = alpha, 
  grid.method = "linear")
LSCI <- cbind(p1.tibs$lo, p1.tibs$up)

## obtain LSLW conformal prediction region
cubic.model <- lm(y ~ x + I(x^2) + I(x^3))
abs.resid <- abs(cubic.model$resid)
smooth.call <- smooth.spline(x, abs.resid, 
  nknots = 10)
lambda <- smooth.call$lambda
df <- smooth.call$df
mad.train.fun <- function(x, y, out = NULL){
  smooth.spline(x[, 1], y, lambda = lambda, 
  df = df, nknots = 10)
}
p2.tibs <- conformal.pred(x = cbind(x,x^2,x^3), y = y, 
  x0 = cbind(x,x^2,x^3), 
  train.fun = train.fun, predict.fun = predict.fun, 
  mad.train.fun = mad.train.fun,
  mad.predict.fun = mad.predict.fun,
  alpha = alpha)
LSLWCI <- cbind(p2.tibs$lo, p2.tibs$up)  

## create plots
ix <- sort(x, index.return = TRUE)$ix

## parametric conformal
plot.new()
plot.window(xlim = c(0,1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], paraCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], paraCI[ix, 2], col = "red", lwd = 1.5)
axis(2)
mtext("shape = 2", side = 2, cex = 1, line = 2)

## LSLW conformal
plot.new()
plot.window(xlim = c(0,1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSLWCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSLWCI[ix, 2], col = "red", lwd = 1.5)

## nonparametric conformal
plot.nonparametric(region = nonparaCI, x = x, 
  y = y, bins = bins)

## LS conformal
plot.new()
plot.window(xlim = c(0,1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSCI[ix, 2], col = "red", lwd = 1.5)


############ shape = 10 ################
## generate gamma regression model 
beta <- c(0.5, 1); shape = 10
rate <- (cbind(1, x) %*% beta) * shape
y <- rgamma(n = n, shape = shape, rate = rate)
y <- y / sd(y)
data <- data.frame(y = y, x = x)
colnames(data)[2:(p+1)] <- paste("x", 1:p, sep = "")

## fit the misspecified regression model
fit <- glm(y ~ x1 + I(x1^2) + I(x1^3), 
  family = "gaussian", data = data)
formula <- fit$formula
newdata <- data
respname <- all.vars(formula)[1]
newdata <- newdata[, !(colnames(data) %in% respname)]
newdata <- as.matrix(newdata)

## obtain parametric and nonparametric conformal 
## prediction regions
cpred <- conformal.glm(fit, parametric = TRUE, 
  nonparametric = FALSE, alpha = alpha,
  bins = bins, cores = cores)
paraCI <- cpred$paraconformal

cpred <- conformal.glm(fit, parametric = FALSE, 
  nonparametric = TRUE, alpha = alpha,
  bins = bins, cores = cores)
nonparaCI <- cpred$nonparaconformal

## obtain LS conformal prediction region
p1.tibs <- conformal.pred(x = cbind(x, x^2, x^3), y = y, 
  x0 = cbind(x, x^2, x^3), train.fun = train.fun, 
  predict.fun = predict.fun, alpha = alpha, 
  grid.method = "linear")
LSCI <- cbind(p1.tibs$lo, p1.tibs$up)

## obtain LSLW conformal prediction region
cubic.model <- lm(y ~ x + I(x^2) + I(x^3))
abs.resid <- abs(cubic.model$resid)
smooth.call <- smooth.spline(x, abs.resid, 
  nknots = 10)
lambda <- smooth.call$lambda
df <- smooth.call$df
mad.train.fun <- function(x, y, out = NULL){
  smooth.spline(x[, 1], y, lambda = lambda, 
  df = df, nknots = 10)
}
p2.tibs <- conformal.pred(x = cbind(x,x^2,x^3), y = y, 
  x0 = cbind(x,x^2,x^3), 
  train.fun = train.fun, predict.fun = predict.fun, 
  mad.train.fun = mad.train.fun,
  mad.predict.fun = mad.predict.fun,
  alpha = alpha)
LSLWCI <- cbind(p2.tibs$lo, p2.tibs$up)  

## create plots
ix <- sort(x, index.return = TRUE)$ix

## parametric conformal
plot.new()
plot.window(xlim = c(0,1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], paraCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], paraCI[ix, 2], col = "red", lwd = 1.5)
axis(2)
mtext("shape = 10", side = 2, cex = 1, line = 2)

## LSLW conformal
plot.new()
plot.window(xlim = c(0,1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSLWCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSLWCI[ix, 2], col = "red", lwd = 1.5)

## nonparametric conformal
plot.nonparametric(region = nonparaCI, x = x, 
  y = y, bins = bins)

## LS conformal
plot.new()
plot.window(xlim = c(0,1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSCI[ix, 2], col = "red", lwd = 1.5)


############ shape = 50 ################
## generate gamma regression model 
beta <- c(0.5, 1); shape = 50
rate <- (cbind(1, x) %*% beta) * shape
y <- rgamma(n = n, shape = shape, rate = rate)
y <- y / sd(y)
data <- data.frame(y = y, x = x)
colnames(data)[2:(p+1)] <- paste("x", 1:p, sep = "")

## fit the misspecified regression model
fit <- glm(y ~ x1 + I(x1^2) + I(x1^3), 
  family = "gaussian", data = data)
formula <- fit$formula
newdata <- data
respname <- all.vars(formula)[1]
newdata <- newdata[, !(colnames(data) %in% respname)]
newdata <- as.matrix(newdata)

## obtain parametric and nonparametric conformal 
## prediction regions
cpred <- conformal.glm(fit, parametric = TRUE, 
  nonparametric = FALSE, alpha = alpha,
  bins = bins, cores = cores)
paraCI <- cpred$paraconformal

cpred <- conformal.glm(fit, parametric = FALSE, 
  nonparametric = TRUE, alpha = alpha,
  bins = bins, cores = cores)
nonparaCI <- cpred$nonparaconformal

## obtain LS conformal prediction region
p1.tibs <- conformal.pred(x = cbind(x, x^2, x^3), y = y, 
  x0 = cbind(x, x^2, x^3), train.fun = train.fun, 
  predict.fun = predict.fun, alpha = alpha, 
  grid.method = "linear")
LSCI <- cbind(p1.tibs$lo, p1.tibs$up)

## obtain LSLW conformal prediction region
cubic.model <- lm(y ~ x + I(x^2) + I(x^3))
abs.resid <- abs(cubic.model$resid)
smooth.call <- smooth.spline(x, abs.resid, 
  nknots = 10)
lambda <- smooth.call$lambda
df <- smooth.call$df
mad.train.fun <- function(x, y, out = NULL){
  smooth.spline(x[, 1], y, lambda = lambda, 
  df = df, nknots = 10)
}
p2.tibs <- conformal.pred(x = cbind(x,x^2,x^3), y = y, 
  x0 = cbind(x,x^2,x^3), 
  train.fun = train.fun, predict.fun = predict.fun, 
  mad.train.fun = mad.train.fun,
  mad.predict.fun = mad.predict.fun,
  alpha = alpha)
LSLWCI <- cbind(p2.tibs$lo, p2.tibs$up)  

## parametric conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], paraCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], paraCI[ix, 2], col = "red", lwd = 1.5)
axis(1);axis(2)
mtext("shape = 50", side = 2, cex = 1, line = 2)

## LSLW conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSLWCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSLWCI[ix, 2], col = "red", lwd = 1.5)
axis(1)

## nonparametric conformal
plot.nonparametric(region = nonparaCI, x = x, 
  y = y, bins = bins)
axis(1)

## LS conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSCI[ix, 2], col = "red", lwd = 1.5)
axis(1)

mtext("x", side = 1, cex = 1, line = 2, outer = TRUE)
mtext("y", side = 2, cex = 1, line = 2, outer = TRUE)
@
\end{center}
\caption{The depiction of conformal prediction regions under simulation 
  setting B when $n = 500$ and the number of bins equals 3.
}
\label{conformal-plots-B-500}
\end{figure}









\newpage
\subsection{Plots corresponding to Section~\ref{sec:regression}}
\label{sec:regressionplots}

%Plot of conformal prediciton regions in sim setting C
\begin{figure}[h!]
\begin{center}
<<conformal-plots-C, fig.height = 7, cache = TRUE, echo = FALSE>>=
par(mfrow = c(3, 4), oma = c(4,4,2,0), mar = c(1,2,1,1))

######### n = 150 ####################
## generate random data
set.seed(13)
alpha <- 0.10
cores <- 6
n <- 150
bins <- 2
p <- k <- length(beta) - 1
x <- matrix(runif(n), ncol = p)

## generate linear regression model
beta <- c(2, 5)
mu <- (cbind(1, x) %*% beta) 
y <- rnorm(n = n, mean = mu, sd = 1)
y <- y / sd(y)
data <- data.frame(y = y, x = x)
colnames(data)[2:(p+1)] <- paste("x", 1:p, sep = "")

## fit the regression model
fit <- glm(y ~ x1, family = "gaussian", data = data)
formula <- fit$formula
newdata <- data
respname <- all.vars(formula)[1]
newdata <- newdata[, !(colnames(data) %in% respname)]
newdata <- as.matrix(newdata)

## obtain parametric and nonparametric conformal 
## prediction regions
cpred <- conformal.glm(fit, parametric = TRUE, 
  nonparametric = FALSE, alpha = alpha,
  bins = bins, cores = cores)
paraCI <- cpred$paraconformal

cpred <- conformal.glm(fit, parametric = FALSE, 
  nonparametric = TRUE, alpha = alpha,
  bins = bins, cores = cores)
nonparaCI <- cpred$nonparaconformal

## obtain LS conformal prediction region
p1.tibs <- conformal.pred(x = x, y = y, x0 = x, 
  train.fun = train.fun, 
  predict.fun = predict.fun, alpha = alpha, 
  grid.method = "linear")
LSCI <- cbind(p1.tibs$lo, p1.tibs$up)

## obtain LSLW conformal prediction region
fit <- lm(y ~ x1, data = data)
abs.resid <- abs(fit$resid)
smooth.call <- smooth.spline(x, abs.resid, 
  nknots = 10)
lambda <- smooth.call$lambda
df <- smooth.call$df
mad.train.fun <- function(x, y, out = NULL){
  smooth.spline(x[, 1], y, lambda = lambda, 
  df = df, nknots = 10)
}
p2.tibs <- conformal.pred(x = x, y = y, x0 = x, 
  train.fun = train.fun, predict.fun = predict.fun, 
  mad.train.fun = mad.train.fun,
  mad.predict.fun = mad.predict.fun,
  alpha = alpha)
LSLWCI <- cbind(p2.tibs$lo, p2.tibs$up)  

## create plots
ix <- sort(x, index.return = TRUE)$ix

## parametric conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], paraCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], paraCI[ix, 2], col = "red", lwd = 1.5)
axis(2)
mtext("n = 150", side = 2, cex = 1, line = 2)
mtext("Parametric\n conformal", side = 3, cex = 1, line = 0)

## LSLW conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSLWCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSLWCI[ix, 2], col = "red", lwd = 1.5)
mtext("LSLW\n conformal", side = 3, cex = 1, line = 0)

## nonparametric conformal
plot.nonparametric(region = nonparaCI, x = x, 
  y = y, bins = bins)
mtext("Nonparametric\n conformal", side = 3, cex = 1, line = 0)

## LS conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSCI[ix, 2], col = "red", lwd = 1.5)
mtext("LS\n conformal", side = 3, cex = 1, line = 0)


######### n = 250 ####################
## generate random data
n <- 250
bins <- 3
p <- k <- length(beta) - 1
x <- matrix(runif(n), ncol = p)

## generate linear regression model
beta <- c(2, 5)
mu <- (cbind(1, x) %*% beta) 
y <- rnorm(n = n, mean = mu, sd = 1)
y <- y / sd(y)
data <- data.frame(y = y, x = x)
colnames(data)[2:(p+1)] <- paste("x", 1:p, sep = "")

## fit the regression model
fit <- glm(y ~ x1, family = "gaussian", data = data)
formula <- fit$formula
newdata <- data
respname <- all.vars(formula)[1]
newdata <- newdata[, !(colnames(data) %in% respname)]
newdata <- as.matrix(newdata)

## obtain parametric and nonparametric conformal 
## prediction regions
cpred <- conformal.glm(fit, parametric = TRUE, 
  nonparametric = FALSE, alpha = alpha,
  bins = bins, cores = cores)
paraCI <- cpred$paraconformal

cpred <- conformal.glm(fit, parametric = FALSE, 
  nonparametric = TRUE, alpha = alpha,
  bins = bins, cores = cores)
nonparaCI <- cpred$nonparaconformal

## obtain LS conformal prediction region
p1.tibs <- conformal.pred(x = x, y = y, x0 = x, 
  train.fun = train.fun, 
  predict.fun = predict.fun, alpha = alpha, 
  grid.method = "linear")
LSCI <- cbind(p1.tibs$lo, p1.tibs$up)

## obtain LSLW conformal prediction region
fit <- lm(y ~ x1, data = data)
abs.resid <- abs(fit$resid)
smooth.call <- smooth.spline(x, abs.resid, 
  nknots = 10)
lambda <- smooth.call$lambda
df <- smooth.call$df
mad.train.fun <- function(x, y, out = NULL){
  smooth.spline(x[, 1], y, lambda = lambda, 
  df = df, nknots = 10)
}
p2.tibs <- conformal.pred(x = x, y = y, x0 = x, 
  train.fun = train.fun, predict.fun = predict.fun, 
  mad.train.fun = mad.train.fun,
  mad.predict.fun = mad.predict.fun,
  alpha = alpha)
LSLWCI <- cbind(p2.tibs$lo, p2.tibs$up)  

## create plots
ix <- sort(x, index.return = TRUE)$ix

## parametric conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], paraCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], paraCI[ix, 2], col = "red", lwd = 1.5)
axis(2)
mtext("n = 250", side = 2, cex = 1, line = 2)

## LSLW conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSLWCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSLWCI[ix, 2], col = "red", lwd = 1.5)

## nonparametric conformal
plot.nonparametric(region = nonparaCI, x = x, 
  y = y, bins = bins)

## LS conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSCI[ix, 2], col = "red", lwd = 1.5)


######### n = 500 ####################
## generate random data
n <- 500
bins <- 3
p <- k <- length(beta) - 1
x <- matrix(runif(n), ncol = p)

## generate linear regression model
beta <- c(2, 5)
mu <- (cbind(1, x) %*% beta) 
y <- rnorm(n = n, mean = mu, sd = 1)
y <- y / sd(y)
data <- data.frame(y = y, x = x)
colnames(data)[2:(p+1)] <- paste("x", 1:p, sep = "")

## fit the regression model
fit <- glm(y ~ x1, family = "gaussian", data = data)
formula <- fit$formula
newdata <- data
respname <- all.vars(formula)[1]
newdata <- newdata[, !(colnames(data) %in% respname)]
newdata <- as.matrix(newdata)

## obtain parametric and nonparametric conformal 
## prediction regions
cpred <- conformal.glm(fit, parametric = TRUE, 
  nonparametric = FALSE, alpha = alpha,
  bins = bins, cores = cores)
paraCI <- cpred$paraconformal

cpred <- conformal.glm(fit, parametric = FALSE, 
  nonparametric = TRUE, alpha = alpha,
  bins = bins, cores = cores)
nonparaCI <- cpred$nonparaconformal

## obtain LS conformal prediction region
p1.tibs <- conformal.pred(x = x, y = y, x0 = x, 
  train.fun = train.fun, 
  predict.fun = predict.fun, alpha = alpha, 
  grid.method = "linear")
LSCI <- cbind(p1.tibs$lo, p1.tibs$up)

## obtain LSLW conformal prediction region
fit <- lm(y ~ x1, data = data)
abs.resid <- abs(fit$resid)
smooth.call <- smooth.spline(x, abs.resid, 
  nknots = 10)
lambda <- smooth.call$lambda
df <- smooth.call$df
mad.train.fun <- function(x, y, out = NULL){
  smooth.spline(x[, 1], y, lambda = lambda, 
  df = df, nknots = 10)
}
p2.tibs <- conformal.pred(x = x, y = y, x0 = x, 
  train.fun = train.fun, predict.fun = predict.fun, 
  mad.train.fun = mad.train.fun,
  mad.predict.fun = mad.predict.fun,
  alpha = alpha)
LSLWCI <- cbind(p2.tibs$lo, p2.tibs$up)  

## create plots
ix <- sort(x, index.return = TRUE)$ix

## parametric conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], paraCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], paraCI[ix, 2], col = "red", lwd = 1.5)
axis(1); axis(2)
mtext("n = 500", side = 2, cex = 1, line = 2)

## LSLW conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSLWCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSLWCI[ix, 2], col = "red", lwd = 1.5)
axis(1)

## nonparametric conformal
plot.nonparametric(region = nonparaCI, x = x, 
  y = y, bins = bins)
axis(1)

## LS conformal
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(y)-1, max(y)+1))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSCI[ix, 1], col = "red", lwd = 1.5)
lines(x[ix], LSCI[ix, 2], col = "red", lwd = 1.5)
axis(1)

mtext("x", side = 1, cex = 1, line = 2, outer = TRUE)
mtext("y", side = 2, cex = 1, line = 2, outer = TRUE)
@
\end{center}
\caption{The depiction of conformal prediction regions under simulation 
  setting C.
}
\label{conformal-plots-C}
\end{figure}

