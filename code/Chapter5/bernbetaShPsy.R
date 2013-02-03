# The functions here are adapted from code by John Kruschke:
# BernBeta.R
# available at http://www.indiana.edu/~kruschke/DoingBayesianDataAnalysis/Programs/
# This code authored by Scott Jackson
# available on GitHub at https://github.com/shoestringpsycholing/kruschke_remix
bern.beta <- function(prior.shape, data, cred.mass = 0.95, tol = 1e-8) {
# Bayesian updating for Bernoulli likelihood and beta prior.
# Input arguments:
#   prior.shape = shape parameters for prior beta distribution, e.g., c(1, 1)
#   data        = vector of 1's and 0's
# Output:
#   S3 object of type "bernbetapost" with:
#   $prior.shape   = the shape parameters of the prior
#   $data.summary  = hits and number of observations
#   $post.shape    = shape parameters for the posterior
#   $evidence   = p(data)

# needs hdi.icdf function
  source("book_code/HDIofICDF.R")
  
# Check for errors in input arguments:
  if(any(data != 1 & data != 0)) { stop("dataVec must be a vector of 1s and 0s.") }

  a <- prior.shape[1]
  b <- prior.shape[2]
  z <- sum(data == 1) # number of 1's in data
  N <- length(data)   # number of observations
# Compute the posterior shape parameters:
  post.shape <- c( a+z , b+N-z )
# Compute the evidence, p(D):
  p.data <- beta(z+a, N-z+b) / beta(a, b)

  bern.beta.post <- list()
  bern.beta.post[["prior.shape"]] <- list(a = a, b = b)
  bern.beta.post[["data.summary"]] <- list(hits = z, n.obs = N)
  bern.beta.post[["post.shape"]] <- list(a = post.shape[1], b = post.shape[2])
  bern.beta.post[["evidence"]] <- p.data
  bern.beta.post[["hdi"]] <- list()
  bern.beta.post[["hdi"]][["cred.mass"]] <- cred.mass
  bern.beta.post[["hdi"]][["interval"]] <- HDIofICDF(ICDFname = qbeta, credMass = cred.mass, tol = tol, shape1 = post.shape[1], shape2 = post.shape[2])
  names(bern.beta.post$hdi$interval) <- c("low", "hi")
  class(bern.beta.post) <- "bern.beta.post"
  return(bern.beta.post)
}

# Plot function for use with bern.beta():
#   Creates a three-panel graph of prior, likelihood, and posterior
#   with highest posterior density interval.

plot.bern.beta.post <- function(bern.beta.post, cred.mass=0.95, binwidth = 0.005) {
  source("code/Chapter5/HDIofICDFShPsy.R")
  a <- bern.beta.post$prior.shape$a
  b <- bern.beta.post$prior.shape$b
  z <- bern.beta.post$data.summary$hits
  N <- bern.beta.post$data.summary$n.obs
  p.data <- bern.beta.post$evidence
  hdi <- bern.beta.post$hdi$interval
  cred.mass <- bern.beta.post$hdi$cred.mass
  # Construct grid of theta values, used for graphing.
  Theta <- seq(from = binwidth/2, to = 1-(binwidth/2), by = binwidth)
  # Compute the prior at each value of theta.
  pTheta <- dbeta(Theta, a, b)
  # Compute the likelihood of the data at each value of theta.
  pDataGivenTheta <- Theta^z * (1-Theta)^(N-z)
  # Compute the posterior at each value of theta.
  pThetaGivenData <- dbeta(Theta, a+z, b+N-z)
  # Open a window with three panels.
  layout(matrix( c( 1,2,3 ) ,nrow=3 ,ncol=1 ,byrow=FALSE)) # 3x1 panels
  par(mar=c(3,3,1,0) , mgp=c(2,1,0) , mai=c(0.5,0.5,0.3,0.1)) # margin specs
  maxY <- max(c(pTheta,pThetaGivenData)) # max y for plotting
  # Plot the prior.
  plot(Theta, pTheta ,type="l", lwd=3,
      xlim=c(0,1), ylim=c(0,maxY), cex.axis=1.2,
      xlab=bquote(theta), ylab=bquote(p(theta)), cex.lab=1.5,
      main="Prior", cex.main=1.5)
if (a > b) { textx <- 0 ; textadj <- c(0,1) } 
else { textx <- 1 ; textadj <- c(1,1) }
text( textx, 1.0*max(pThetaGivenData) ,
      bquote("beta(" * theta * "|" * .(a) * "," * .(b) * ")" ),
      cex=2.0, adj=textadj)
# Plot the likelihood: p(data|theta)
plot(Theta, pDataGivenTheta ,type="l" ,lwd=3 ,
      xlim=c(0,1), cex.axis=1.2, xlab=bquote(theta),
      ylim=c(0,1.1*max(pDataGivenTheta)),
      ylab=bquote( "p(D|" * theta * ")" ),
      cex.lab=1.5, main="Likelihood", cex.main=1.5)
if(z > .5*N) { textx <- 0 ; textadj <- c(0,1) }
else { textx <- 1 ; textadj <- c(1,1) }
text(textx, 1.0*max(pDataGivenTheta), cex=2.0,
      bquote("Data: z=" * .(z) * ",N=" * .(N)), adj=textadj)
# Plot the posterior.
plot(Theta, pThetaGivenData, type="l", lwd=3,
      xlim=c(0,1), ylim=c(0,maxY), cex.axis=1.2,
      xlab=bquote(theta), ylab=bquote("p(" * theta * "|D)"),
      cex.lab=1.5, main="Posterior", cex.main=1.5)
if(a+z > b+N-z) { textx <- 0 ; textadj <- c(0,1) }
else { textx <- 1 ; textadj <- c(1,1) }
text(textx, 1.00*max(pThetaGivenData), cex=2.0,
      bquote("beta(" * theta * "|" * .(a+z) * "," * .(b+N-z) * ")"),
      adj=textadj)
text(textx, 0.75*max(pThetaGivenData), cex=2.0,
      bquote("p(D)=" * .(signif(p.data,3))), adj=textadj)
# Mark the HDI in the posterior.
hpdHt <- mean(c(dbeta(hdi[["low"]],a+z,b+N-z), dbeta(hdi[["hi"]],a+z,b+N-z)))
lines(c(hdi[["low"]],hdi[["low"]]), c(-0.5,hpdHt), type="l", lty=2, lwd=1.5)
lines(c(hdi[["hi"]],hdi[["hi"]]), c(-0.5,hpdHt), type="l", lty=2, lwd=1.5)
lines(c(hdi[["low"]], hdi[["hi"]]), c(hpdHt,hpdHt), type="l", lwd=2)
text(mean(c(hdi[["low"]], hdi[["hi"]])), hpdHt, bquote(.(100*cred.mass) * "% HDI"),
      adj=c(0.5,-1.0), cex=2.0)
text(hdi[["low"]], hpdHt, bquote(.(round(hdi[["low"]],3))),
      adj=c(1.1,-0.1), cex=1.2)
text(hdi[["hi"]], hpdHt, bquote(.(round(hdi[["hi"]],3))),
      adj=c(-0.1,-0.1), cex=1.2)
}



