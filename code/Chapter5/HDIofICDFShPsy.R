hdi.icdf.bernbeta <- function(bernbetapost, cred.mass=0.95, tol=1e-8, ... ) {
   # Arguments:
   #   ICDFname is R's name for the inverse cumulative density function
   #     of the distribution.
   #   cred.mass is the desired mass of the HDI region.
   #   tol is passed to R's optimize function.
   # Return value:
   #   Highest density iterval (HDI) limits in a vector.
   # Example of use: For determining HDI of a beta(30,12) distribution, type
   #   HDIofICDF(qbeta, shape1 = 30 , shape2 = 12 )
   #   Notice that the parameters of the ICDFname must be explicitly named;
   #   e.g., HDIofICDF( qbeta , 30 , 12 ) does not work.
   # Adapted and corrected from Greg Snow's TeachingDemos package.
   if(class(bernbetapost) == "bernbetapost") {
     ICDFname <- qbeta
     shape1 <- bernbetapost$post.shape$a
     shape2 <- bernbetapost$post.shape$b
   } else { stop(paste("function 'hdi.icdf.bernbeta' needs an object of class 'bernbetapost', not", as.character(class(bernbetapost)))) }
   incred.mass <- 1.0 - cred.mass
   intervalWidth <- function(lowTailPr, ICDFname, cred.mass, ... ) {
      ICDFname(cred.mass + lowTailPr, shape1 = shape1, shape2 = shape2, ... ) - ICDFname(lowTailPr, shape1 = shape1, shape2 = shape2, ... )
   }
   optInfo <- optimize(intervalWidth, c(0, incred.mass), ICDFname=ICDFname,
                        cred.mass=cred.mass, tol=tol, ... )
   HDIlowTailPr = optInfo$minimum
   hdi <- list()
   hdi["low"] <- ICDFname(HDIlowTailPr, shape1 = shape1, shape2 = shape2, ... )
   hdi["hi"] <- ICDFname(cred.mass + HDIlowTailPr, shape1 = shape1, shape2 = shape2, ... )
   hdi["cred.mass"] <- cred.mass
   class(hdi) <- "hdi"
   return(hdi)
}
