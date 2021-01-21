#------------------------------------------------------------------------------#
#                                                                              #
#     Functions generating examples of Random Field/Functional data            #
#                                                                              #
#------------------------------------------------------------------------------#
# Required packages:
#   - abind
#   - spatstat
#   - MASS
#
# Contained functions:
#   - ArbCovNoise (rdy)
#   - SinCosSumNoise
#   - GaussDensitySumNoise
#   - OUNoise
#   - DegrasNonGaussNoise
#   - BernsteinSumNoise
#   - HermiteSumNoise
#   - SquaredExp2DNoise
#   - GaussDensitySum2DNoise
#
#------------------------------------------------------------------------------#
# Developer notes:
# - Style guideline NOT included
#------------------------------------------------------------------------------#

#' Samples from a process having a given covariance structure
#'
#' Creates sample paths from a 1D random field having an arbitrary covariance
#' function. The default is gthe square exponential covariance function:
#'  C(h)=exp( -h^2/(4*nu^2) )
#'
#' @param N Integer amount of realisations of the random field to be generated.
#' @param x Vector locations at which the random field is evaluated.
#' @param cov.fun vectorized function which computes the covariance between two
#' points in the vector x.
#' @param ... additional arguments to cov.fun
#' @return Matrix containing the realisations of the random field as columns. The rows are the locations.
#'
#' @export
ArbCovProcess <- function( N,
                           x = seq( 0, 1, length.out = 100 ),
                           covf = covf.square.exp, ...
){
  # Get the covariance matrix
  Sigma = outer( x, x, FUN = Vectorize( function( y, z ) covf( y, z, ... ) ) )

  # Simulate realizations from the Gaussian process with the given covariance
  # function
  Y = t( MASS::mvrnorm( N, mu = rep( 0, length( x ) ), Sigma = Sigma ) )
}


#' Samples from Random Linear Combination of Given Basis Functions
#'
#' Creates sample paths from a random field generated as a random sum of given
#' basis functions.
#'
#' @param N Integer amount of realisations of the random field to be generated.
#' @param x Vector locations at which the random field is evaluated.
#' @param randf Function generating a vector of random numbers with mean zero and variance 1. Default is rnorm().
#' @return Matrix containing the realisations of the random field as columns. The rows are the locations.
#' @export
RandomLinearSum <- function( N,
                             x = seq( 0, 1, length.out = 100 ),
                             basisf,
                             randf = rnorm,
                             normalize = TRUE ){
  # Get the basis functions evaluated on x.
  f <- basisf( x )
  # Get the number of basis functions
  nBasis <- dim( f )[2]

  # Normalize the functions such that the resulting processes will have variance
  # one if the random multipliers have.
  if( normalize ){
    fSqSum <- apply( f^2, 1, sum )
    f <- f / sqrt( fSqSum )
  }

  return( f %*% matrix( randNumber( nBasis * N ), nBasis, N ) )
}


#' Samples from the Cosine Field
#'
#' Creates sample paths from a 1D Gaussian fields generated as a random sum of sine and cosine functions with frequency pi/2 and standard Gaussians as coefficients.
#' @param N Integer amount of realisations of the random field to be generated.
#' @param x Vector locations at which the random field is evaluated.
#' @return Matrix containing the realisations of the random field as columns. The rows are the locations.
#' @export
SinCosSumNoise <- function( N,
                            x     = seq( 0, 1, length.out = 100 ),
                           ){

  return( ( sin( pi/2 * x ) %*% t( rnorm( N, 0, 1 ) ) +
              cos( pi/2 * x ) %*% t( rnorm( N, 0, 1 ) ) ) )
}

#' Samples from Random Linear Combination of Gaussian Densities
#'
#' Creates sample paths from a 1D field generated as a random sum of Gaussian densities with different means and variances and random coefficients.
#' @param N Integer amount of realisations of the random field to be generated.
#' @param x Vector locations at which the random field is evaluated.
#' @param randNumber Function generating a vector of random numbers with mean zero and variance 1. Default is rnorm().
#' @param anchorPoints Vector containing the locations of the mean of the Gaussian densities in the random sum. Default value are 15 equidistant knots of the interval specifying the domain.
#' @param anchorSd Vector containing the standard deviations of the stand
#' @return Matrix containing the realisations of the random field as columns. The rows are the locations.
#' @export
GaussDensitySumNoise <- function( N,
                                  x          = seq( 0, 1, length.out = 100 ),
                                  randNumber = rnorm,
                                  anchorPoints = NULL,
                                  anchorSd   = NULL){
  if(is.null(anchorPoints)){
    anchorPoints <- seq( from = min(x), to = max(x), length.out = 15 )
  }
  nAnchorPoints = length(anchorPoints)
  if(is.null(anchorSd)){
    anchorSd <- rep(diff(range(x))/nAnchorPoints, nAnchorPoints)
  }
  f <- sapply(1:length(anchorPoints),
              function(pt) dnorm(x, mean=anchorPoints[pt],
                                 sd=anchorSd[pt]))
  fSqSum <- apply(f^2, 1, sum)
  fNorm <- f / sqrt(fSqSum)
  return( (fNorm %*% matrix( randNumber( nAnchorPoints * N ),
                             nAnchorPoints,
                             N ) ) )
}

#' Samples from Ornstein-Uhlenbeck process
#'
#' Creates sample paths from a 1D Ornstein Uhlenbeck process.
#' @param N Integer amount of realisations of the random field to be generated.
#' @param x Vector locations at which the random field is evaluated.
#' @param alphaOU Numeric value of the retraction parameter in an OU process. Default is 5.
#' @param sigmaOU Numeric value of the sigma parameter in an OU process. Default is sqrt(10).
#' @param ganna Vector of the same length as x giving the mean of the OU process.
#' @return Matrix containing the realisations of the random field as columns. The rows are the locations.
#' @export
OUNoise <- function( N,
                     x       = seq( 0, 1, length.out = 100 ),
                     alphaOU = 5,
                     sigmaOU = sqrt( 10 ),
                     gamma   = NULL ){
  ## Initialize matrix containing as rows a path of the process
  Y  = matrix( NA, length(x), N )

  ## Mean Curve vector
  if( is.null( gamma ) ){
    gamma <- rep( 0, length( x ) )
  }

  ## Get starting values
  Y[1,] <- rnorm(N, mean = 0, sd = sigmaOU / sqrt( 2*alphaOU )) + gamma[1]

  ## Compute the Ornstein-Uhlenbeck trajectory
  for (n in 2:length(x) ){
    dt <- x[n]-x[n-1]
    Y[n,] <- vapply( Y[n-1,], function(y) rnorm(1, mean = (y - gamma[n-1])
                                                * exp(-alphaOU * dt) + gamma[n], sd = sigmaOU
                                                * sqrt((1 - exp(-2 * alphaOU * dt)) / 2 / alphaOU) )
                     , FUN.VALUE=1)
  }
  return( Y )
}


#' Samples from a Non-Gaussian process used in Degras(2011)
#'
#' Creates sample paths from a 1D non Gaussian field as used in Degras (2011, Simultaneous confidence bands for nonparametric regression with functional data. Statistica Sinica, 1735-1765.).
#' @param N Integer amount of realisations of the random field to be generated.
#' @param x Vector locations at which the random field is evaluated.
#' @param Function computing the pointwise variance of the field. Default value is unit variance everywhere.
#' @return Matrix containing the realisations of the random field as columns. The rows are the locations.
#' @export
DegrasNonGaussNoise <- function( N,
                                 x     = seq( 0, 1, length.out = 100)
                                 ){

  rcoeff <- cbind( rchisq(N, 1, ncp = 0), rexp(N) )

  vapply( 1:N,
          function(l) sqrt( 2 ) / 6 * ( rcoeff[ l, 1 ] - 1 ) * sin( pi * x ) +
            2 / 3 * ( rcoeff[ l, 2 ] - 1 ) * ( x - 0.5 ),
          FUN.VALUE = rep( 0, length( x ) ) )
}


#' Samples from Random Linear Combination of Bernstein Polynomials
#'
#' Creates sample paths from a 1D random field generated as a random sum ofthe first seven Bernstein polynomials.
#' @param N Integer amount of realisations of the random field to be generated.
#' @param x Vector locations at which the random field is evaluated.
#' @param randNumber Function generating a vector of random numbers with mean zero and variance 1. Default is rnorm().
#' @return Matrix containing the realisations of the random field as columns. The rows are the locations.
#' @export
BernsteinSumNoise <- function( N,
                               x          = seq( 0, 1, length.out = 100 ),
                               randNumber = rnorm ){

  f <- cbind( ( 1 - x )^6,
              6 * x * ( 1 - x )^5,
              15 * x^2 * ( 1 - x )^4,
              20 * x^3 * ( 1 - x )^3,
              15 * x^4 * ( 1 - x )^2,
              6 * x^5 * ( 1 - x ),
              x^6 )
  fSqSum <- apply( f^2, 1, sum )
  fNorm <- f / sqrt( fSqSum )
  nBasis <- dim(f)[2]
  return( fNorm %*% matrix( randNumber( nBasis * N), nBasis, N ) )
}

#' Samples from Random Linear Combination of Hermite Polynomials
#'
#' Creates sample paths from a 1D random field generated as a random sum of the first five Hermite polynomials.
#' @param N Integer amount of realisations of the random field to be generated.
#' @param x Vector locations at which the random field is evaluated.
#' @param randNumber Function generating a vector of random numbers with mean zero and variance 1. Default is rnorm().
#' @return Matrix containing the realisations of the random field as columns. The rows are the locations.
#' @export
HermiteSumNoise <- function( N,
                             x          = seq( 0, 1, length.out = 100 ),
                             randNumber = rnorm
){

  f <- cbind( 1,
              6 * x,
              4 * ( 3 * x )^2 - 2,
              8 * ( 3 * x )^3 - 12 * ( 3 * x ),
              16 * ( 3 * x )^4 - 48 * ( 3 * x)^2 + 12 )
  fSqSum <- apply( f^2, 1, sum )
  fNorm <- f / sqrt( fSqSum )
  nBasis <- dim(f)[2]
  return( fNorm %*% matrix( rnorm( nBasis * N), nBasis, N ) )
}

#' Samples from Random Linear Combination of Gaussian Densities
#'
#'Creates a sample from an Gaussian isotropic 2D random field having a square exponential function: C(h)=exp( -h^2/(4*nu^2) ) as covariance structure by smoothing independent noise with a Gaussian kernel.
#' @param N Integer amount of realisations of the random field to be generated.
#' @param x Vector locations at which the random field is evaluated.
#' @param nu Numeric positiv real number. Smoothing bandwidth of the Gaussian kernel
#' @return Array containing the realisations of the random field. Last dimension indexes the realisations.
#'
#' @export
SquaredExp2DNoise <- function( N,
                               x     = seq( 1, 50, by = 1),
                               nu    = 5
){
  a = 4 * nu
  K <- L <- length( x )

  ### Generate fields of white noise
  Y = array( rnorm( N*(K+2*a)*(L+2*a) ) , c((K+2*a), (L+2*a), N) )

  ### Smooth the white noise field with a Gaussian kernel
  Y = array( apply(Y, 3, function(I) as.vector(as.matrix(spatstat::Smooth(spatstat::as.im(I),normalise=FALSE,sigma=nu)) )), dim=c((K+2*a), (L+2*a), N) )

  ### Scale the variance of the random fields and return the realisations
  return( Y[ (a+1):(K+a), (a+1):(L+a),]  * 2* nu *sqrt(pi) )
}

#'Creates a sample from an 2D random field, which is the sum of 2-D Gaussian densities..
#' @param N Integer amount of realisations of the random field to be generated.
#' @param x Vector locations at which the random field is evaluated.
#' @param randNumber Function generating a vector of random numbers with mean zero and variance 1. Default is rnorm().
#' @param M Integer positiv real number. Square root of Gaussian kernels used for generation of the noise.
#' @param bdx Numeric postiv real number. Possibility to center the Gaussian kernels outside the range of x. Bdx simply extends the range by +-bdx and then a equidistant grid having M points from the interval range(x)-+bdx is formed giving the means for the Gaussian kernel bumps.
#' @param Fnorm Array containing the scaling to pointwise unit variance, if it is NULL this is done internally. This option is only included to speed up simulations.
#' @return Array containing the realisations of the random field. Last dimension indexes the realisations.
#'
#' @export
GaussDensitySum2DNoise <- function(N,
                                   x = seq( 0, 1, length.out = 50),
                                   randNumber=rnorm, M=6, bdx=0.02, Fnorm=NULL){
  if(is.null(Fnorm)){
    # Compute vector with gridpoints for the mean of the random Gaussians
    nT <- nT1 <- nT2 <- length(x)
    rx     = range(x)
    grd    = seq(rx[1]+bdx,rx[2]-bdx,length.out=M)
    ptGrid = as.matrix(expand.grid(grd,grd))

    FF=NULL
    for( k in 1:dim(ptGrid)[1] ){
      F1 = outer(x,x, FUN=Vectorize(function(s,t){ dnorm( sqrt(sum((c(s,t)-ptGrid[k,])^2)), mean=0, sd=0.1) }, vectorize.args = c("s","t")) )
      FF=abind::abind(FF,F1, along=3)
    }
    FFsqSum = sqrt(apply(FF^2, MARGIN=c(1, 2), sum))
    Fnorm   =  matrix( FF / array(rep(FFsqSum, M^2), dim=c(nT,nT, M^2) ), nT^2,M^2)
    rm(FF, FFsqSum)
  }else{
    dimF = dim(Fnorm)
    M    = sqrt(dimF[3])
    nT1  = dimF[1]
    nT2  = dimF[2]
    Fnorm = matrix(Fnorm, nT1*nT2,M^2)
  }

  rcoef <- matrix( randNumber(M^2*N), M^2, N )

  ### Scale the variance of the random fields and return the realisations
  return( array( Fnorm%*%rcoef, dim=c(nT1,nT2, N) ) )
}
