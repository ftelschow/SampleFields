#------------------------------------------------------------------------------#
#                                                                              #
#     Functions generating examples of Random Field/Functional data            #
#                                                                              #
#------------------------------------------------------------------------------#
# Required packages:
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
# - a function only allowing vector input for x is called process
# - a fucntion allowing an arbitrary matrix input for x is called field
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# Definitions of 1D processes ( x only vector valued )
#------------------------------------------------------------------------------#
#' Samples from a process having a pre-defined covariance structure
#'
#' Creates sample paths from a 1D random field having an arbitrary covariance
#' function. The default is the square exponential covariance function:
#'  C(h)=exp( -h^2/(4*nu^2) )
#'
#' @param N Integer amount of realisations of the random field to be generated.
#' @param x Vector locations at which the random process is evaluated.
#' @param cov.fun vectorized function which computes the covariance between two
#' points in the vector x.
#' @param ... additional arguments to cov.fun
#' @return RandomField class object containing the realisations of the random
#' field and the locations.
#'
#' @export
ArbCovProcess <- function( N,
                           x = seq( 0, 1, length.out = 100 ),
                           covf = covf.square.exp, ...
){
  # Get the covariance matrix
  Sigma = outer( x, x, FUN = Vectorize( function( y, z ) covf( y, z, ... ) ) )

  # Simulate realizations from the Gaussian process with the given covariance
  # function and output as RandomField object
  samp.field = list( values = t( MASS::mvrnorm( N,
                                                mu = rep( 0, length( x ) ),
                                                Sigma = Sigma ) ),
                     locations = x )
  class( samp.field ) = "RandomField"

  return( samp.field )
}


#' Samples from a Non-Gaussian process used in Degras(2011)
#'
#' Creates sample paths from a 1D non Gaussian field as used in Degras
#' (2011, Simultaneous confidence bands for nonparametric regression with
#' functional data. Statistica Sinica, 1735-1765.).
#'
#' @inheritParams ArbCovProcess
#' @inherit ArbCovProcess return
#'
#' @export
DegrasNonGaussProcess <- function( N, x = seq( 0, 1, length.out = 100) ){

  rcoeff <- cbind( rchisq(N, 1, ncp = 0), rexp(N) )

  samp.field = vapply( 1:N,
                       function( l ){ sqrt( 2 ) / 6 * ( rcoeff[ l, 1 ] - 1 ) *
                           sin( pi * x ) + 2 / 3 *
                           ( rcoeff[ l, 2 ] - 1 ) * ( x - 0.5 )
                       },
                       FUN.VALUE = rep( 0, length( x ) ) )

  # Create the output class object
  samp.field = list( values = samp.field,
                     locations = x )
  class( samp.field ) = "RandomField"

  return( samp.field )
}


#' Samples from Ornstein-Uhlenbeck process
#'
#' Creates sample paths from a 1D Ornstein Uhlenbeck process.
#'
#' @inheritParams ArbCovProcess
#' @inherit ArbCovProcess return
#'
#' @param alphaOU Numeric value of the retraction parameter in an OU process. Default is 5.
#' @param sigmaOU Numeric value of the sigma parameter in an OU process. Default is sqrt(10).
#' @param gamma Vector of the same length as x giving the mean of the OU process.
#' @return RandomField class object containing the realisations of the random
#' field and the locations.
#' @export
OUProcess <- function( N,
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


#------------------------------------------------------------------------------#
# Definitions of random fields ( x might be matrix valued )
#------------------------------------------------------------------------------#

#' Samples from Random Linear Combination of Given Basis Functions
#'
#' Creates sample paths from a random field generated as a random sum of given
#' basis functions.
#'
#' @inheritParams ArbCovProcess
#' @inherit ArbCovProcess return
#' @param x Vector or matrix of dimension T x D, where T is the amount of
#' locations and D the dimension of the domain, defining the coordinates of
#' the locations at which the random field is evaluated.
#' @param basisf Function evaluating the basis functions on x and returning a
#' matrix of dimension T x K, where K is the number of basis functions.
#' @param randf Function generating a vector of random numbers as coefficients
#' in the linear combination. Default is rnorm().
#' @param normalize boolean If true the resulting random field is scaled to have
#' variance one.
#'
#' @export
RandomBasisSum <- function( N,
                            x = seq( 0, 1, length.out = 100 ),
                            basisf,
                            randf = rnorm.mod,
                            normalize = TRUE ){
  # Get the basis functions evaluated on x.
  f <- basisf( x )
  # Get the number of basis functions
  nBasis <- dim( f )[2]

  # Normalize the functions such that the resulting processes will have variance
  # one if the random multipliers have.
  if( normalize ){
    f <- f / sqrt( rowSums( f^2 ) )
  }

  # Create the output class object
  samp.field = list( values =  f %*% randf( nBasis, N ),
                     locations = x )
  class( samp.field ) = "RandomField"

  return( samp.field )
}


#' Samples from Random Linear Combination of First Six Bernstein Polynomials
#'
#' Creates sample paths from a 1D random field generated as a random sum of the
#' first seven Bernstein polynomials.
#'
#' @inheritParams RandomBasisSum
#' @inherit RandomBasisSum return
#'
#' @export
BernsteinSumNoise <- function( N,
                               x = seq( 0, 1, length.out = 100 ),
                               randf = rnorm.mod,
                               normalize = TRUE ){

  # define the basis functions as Bernstein polynomials
  basisf <- function(x){ cbind( ( 1 - x )^6,
                                6 * x * ( 1 - x )^5,
                                15 * x^2 * ( 1 - x )^4,
                                20 * x^3 * ( 1 - x )^3,
                                15 * x^4 * ( 1 - x )^2,
                                6 * x^5 * ( 1 - x ),
                                x^6 )
  }

  # Output a random linear combination of the first 6 Bernstein polynomials
  RandomBasisSum( N = N,
                  x = x,
                  basisf,
                  randf = randf,
                  normalize = normalize )
}


#' Samples from Random Linear Combination of First Five Hermite Polynomials
#'
#' Creates sample paths from a 1D random field generated as a random sum of the
#' first five Hermite polynomials.
#'
#' @inheritParams RandomBasisSum
#' @inherit RandomBasisSum return
#'
#' @export
HermiteSumNoise <- function( N,
                             x = seq( 0, 1, length.out = 100 ),
                             randf = rnorm.mod,
                             normalize = TRUE ){

  # define the basis functions as Hermite polynomials
  basisf <- function( x ){ cbind( 1,
                                  6 * x,
                                  4 * ( 3 * x )^2 - 2,
                                  8 * ( 3 * x )^3 - 12 * ( 3 * x ),
                                  16 * ( 3 * x )^4 - 48 * ( 3 * x)^2 + 12 )
                        }

  # Output a random linear combination of the first 6 Bernstein polynomials
  RandomBasisSum( N = N,
                  x = x,
                  basisf,
                  randf = randf,
                  normalize = normalize )
}


#' Samples from Random Linear Combination of Gaussian Densities
#'
#' Creates sample paths from a 1D field generated as a random sum of Gaussian
#' densities with different means and variances and random coefficients.
#'
#' @inheritParams RandomBasisSum
#' @inherit RandomBasisSum return
#' @param means Matrix of dimension K x D indicating the mean locations of the
#' K multivariate Gaussian densities used as the basis functions.
#' @param sigmas Matrix of dimension K x D*D containing the covariance matrix
#' as vectors, i.e. \code{matrix(sigmas[k,], D, D)} must be the covariance of
#' the k-th kernel.
#'
#' @export
RandomNormalSumField <- function( N,
                                  x = seq( 0, 1, length.out = 100 ),
                                  randf = rnorm.mod,
                                  means = NULL,
                                  sigmas = NULL,
                                  normalize = TRUE ){
  # Get the dimension of the domain
  dimx = dim( x )

  if( is.null( dimx ) ){
    # Get the default values of means
    if( is.null( means ) ){
      means <- seq( from = min( x ), to = max( x ), length.out = 15 )
    }
    # Amount of means
    Nmeans = length( means )

    # Get the default values for the sds
    if( is.null( sigmas ) ){
      sigmas <- rep( diff( range( x ) ) / Nmeans, Nmeans )
    }

    f <- sapply( 1:Nmeans,
                function( pt ) dnorm( x,
                                      mean = means[pt],
                                      sd = sigmas[pt] ) )
  }else{
    # Get the dimension of the domain
    dimB   = dim( means )
    Nmeans = dimB[1]
    D      = dimB[2]

    # Get the basis function
    f = vapply( 1:Nmeans,
                function( l ){
                   mvtnorm::dmvnorm( x, means[l,], matrix( sigmas[l,], D, D ) )
                },
                FUN.VALUE = seq( 0, pi, length.out = dim( x )[1] )
            )
  }

  # Normalize the functions such that the resulting processes will have variance
  # one if the random multipliers have.
  if( normalize ){
    f <- f / sqrt( rowSums( f^2 ) )
  }

  # Create the output class object
  samp.field = list( values = f %*% randf( Nmeans, N ),
                     locations = x )
  class( samp.field ) = "RandomField"

  return( samp.field )

}


#' Samples from the Cosine Field
#'
#' Creates sample paths from a 1D Gaussian fields generated as a random sum of
#' sine and cosine functions with frequency pi/2 and standard Gaussians as
#' coefficients.
#'
#' @inheritParams RandomBasisSum
#' @inherit RandomBasisSum return
#'
#' @export
CosineField <- function( N, x = seq( 0, 1, length.out = 100 ) ){

  samp.field = sin( pi/2 * x ) %*% t( rnorm( N, 0, 1 ) ) +
               cos( pi/2 * x ) %*% t( rnorm( N, 0, 1 ) )

  # Create the output class object
  samp.field = list( values = samp.field,
                     locations = x )
  class( samp.field ) = "RandomField"

  return( samp.field )

}


#' Creates a sample from an Gaussian isotropic 2D random field having a square
#' exponential function: C(h)=exp( -h^2/(4*nu^2) ) as covariance structure by
#' smoothing independent noise with a Gaussian kernel.
#'
#' @inheritParams aws::kernsm
#' @inheritParams RandomBasisSum
#' @inherit RandomBasisSum return
#'
#' @export
KernelSmoothedField <- function( N,
                                 x = seq( 1, 50, by = 1 ),
                                 randf = rnorm.mod,
                                 h = 5,
                                 kern  = "Gaussian",
                                 unit  = c( "SD", "FWHM" )
                                  ){
  # Get the dimension of the input coordinates
  dimx = dim( x )
  if( is.null( dimx ) ){
    Nx = length( x )
  }else{
    D  = dimx[2]
    Nx = dimx[1]
  }

  # Create the array of noise to be smoothed
  y = array( randf( dim( x ) * N ), dim = c( dimx, N ) )

  samp.field = aws::kernsm( y, h = h, kern = kern, unit = unit )

  # Create the output class object
  samp.field = list( values = samp.field@yhat,
                     locations = x )
  class( samp.field ) = "RandomField"

  return( samp.field )

}

