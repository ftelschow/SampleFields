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
#  - Only 1D
#   - ArbCovProcess
#   - DegrasNonGaussProcess
#   - OUProcess
#
#  - arbitrary dimension
#   - RandomBasisSum
#   - RandomBernsteinSum
#   - RandomHermiteSum
#   - RandomNormalSum
#   - SquaredExp2DNoise
#   - GaussDensitySum2DNoise
#   - CosineField
#   - KernelSmoothedField
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
#' @param covf Function computing the covariance between two points in the
#' vector x.
#' @param ... additional arguments to covf
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

  # Add additional important variables to RandomField class
  samp.field$dim    = dim( samp.field$values )
  gDim              = getDim( samp.field$locations )
  samp.field$D      = gDim$D
  samp.field$nloc   = gDim$nloc
  samp.field$N      = N

  # Make samp.field a S3 class
  class( samp.field ) = "RandomField"

  return( samp.field )
}

#' Creates sample paths from a 1D random field having at each location an
#' iid sample from a pre-defined probability distribution.
#' The default is a standard zero-mean, unit-variance Gaussian.
#'
#' @param N Integer amount of realisations of the random field to be generated.
#' @param x Vector locations at which the random process is evaluated.
#' @param rdistr Function computing random samples from a distribution. Note
#' that the function must have an paramter "n" for the sample size drawn from
#' this distribution.
#' @param ... additional arguments to rdistr
#' @return RandomField class object containing the realisations of the random
#' field and the locations.
#'
#' @export
IIDProcess <- function( N,
                        x = seq( 0, 1, length.out = 100 ),
                        rdistr = rnorm, ...
){
  # Get the covariance matrix
  lx = length(x)

  # Obtain the values of the field
  fieldValues = matrix(rdistr(n = N * lx, ...), lx, N)

  # Simulate realizations from the Gaussian process with the given covariance
  # function and output as RandomField object
  samp.field = list( values    = fieldValues,
                     locations = x )

  # Add additional important variables to RandomField class
  samp.field$dim    = dim( samp.field$values )
  gDim              = getDim( samp.field$locations )
  samp.field$D      = gDim$D
  samp.field$nloc   = gDim$nloc
  samp.field$N      = N

  # Make samp.field a S3 class
  class( samp.field ) = "RandomField"

  return( samp.field )
}


#' Samples from a Non-Gaussian process used in Degras(2011). The only difference
#' is that we scaled the process to have variance 1 at each location.
#'
#' Creates sample paths from a 1D non Gaussian field as used in Degras
#' (2011, Simultaneous confidence bands for nonparametric regression with
#' functional data. Statistica Sinica, 1735-1765.).
#'
#' @inheritParams ArbCovProcess
#' @inherit ArbCovProcess return
#'
#' @export
DegrasNonGaussProcess <- function( N, x = seq( 0, 1, length.out = 100 ) ){

  rcoeff <- cbind( rchisq( N, 1, ncp = 0 ), rexp( N ) )

  samp.field = vapply( 1:N,
                       function( l ){ sqrt( 2 ) / 6 * ( rcoeff[ l, 1 ] - 1 ) *
                           sin( pi * x ) + 2 / 3 *
                           ( rcoeff[ l, 2 ] - 1 ) * ( x - 0.5 )
                       },
                       FUN.VALUE = rep( 0, length( x ) ) )

  # Return RandomField object
  return( RandomField( field = samp.field / sqrt(sin(pi * x)^2 / 9 + 4 / 9 * (x - 0.5)^2),
                       locations = x ) )
}


#' Samples from Ornstein-Uhlenbeck process
#'
#' Creates sample paths from a 1D Ornstein Uhlenbeck process.
#'
#' @inheritParams ArbCovProcess
#' @inherit ArbCovProcess return
#'
#' @param alpha Numeric value of the retraction parameter in an OU process. Default is 5.
#' @param sigma Numeric value of the sigma parameter in an OU process. Default is sqrt(10).
#' @param gamma Vector of the same length as x giving the mean of the OU process.
#'
#' @export
OUProcess <- function( N,
                       x       = seq( 0, 1, length.out = 100 ),
                       alpha = 5,
                       sigma = sqrt( 10 ) ){

  # Initialize matrix containing as columns a path of the process
  Y = matrix( NA, length( x ), N )

  # Get starting values
  Y[ 1, ] <- rnorm( N, mean = 0, sd = sigma / sqrt( 2*alpha ) )

  # Compute the Ornstein-Uhlenbeck sample path
  for( n in 2:length( x ) ){
    dt <- x[ n ] - x[ n - 1 ]
    Y[ n, ] <- vapply( Y[ n - 1, ],
                       function( y ){
                         rnorm( 1,
                                mean = ( y - gamma[ n - 1 ] ) * exp( -alpha * dt),
                                sd   = sigma * sqrt( (1 - exp( -2 * alpha * dt ) ) / 2 / alpha )
                                )
                       }
                     , FUN.VALUE=1)
  }


  # Return RandomField object
  return( RandomField( field = Y, locations = x ) )
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
#'
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
  if(is.function(basisf)){
    # Get the basis functions evaluated on x.
    f <- basisf(x)
  }else if(is.array(basisf)){
    if(dim(basisf)[1] == length(x)){
      f <- basisf
    }else{
      stop("If basisf is an array it needs to have length(x) rows!")
    }
  }else{
    stop("basisf needs to be either a function or an array with length(x) rows!")
  }

  # Get the number of basis functions
  nBasis <- dim( f )[2]

  # Normalize the functions such that the resulting processes will have variance
  # one if the random multipliers have.
  if( normalize ){
    f <- f / sqrt( rowSums( f^2 ) )
  }

  # Return RandomField object
  return( RandomField( field = f %*% randf( nBasis, N ), locations = x ) )
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
RandomBernsteinSum <- function( N,
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
RandomHermiteSum <- function( N,
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
#'
#' @param means Matrix of dimension K x D indicating the mean locations of the
#' K multivariate Gaussian densities used as the basis functions.
#' @param sigmas Matrix of dimension K x D*D containing the covariance matrix
#' as vectors, i.e. \code{matrix(sigmas[k,], D, D)} must be the covariance of
#' the k-th kernel.
#'
#' @export
RandomNormalSum <- function( N,
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
    # Default value for the mean locations
    if( is.null( means ) | length( means ) == 1 ){
      if( length( means ) == 1 ){
        Nmeans = means
      }else{
        Nmeans = ceiling( dimx[1] * 0.3 )
      }

      means  = x[ seq( 1, dimx[1], length.out = Nmeans ), ]
      rownames(means) <- NULL
    }

    # Get the dimension of the domain
    dimB   = dim( means )
    if( is.null( dimB ) ){
      means = matrix( means, 1, length(means) )
    }else{
      Nmeans = dimB[1]
      D      = dimB[2]
    }

    if( is.null(Nmeans) ){
      Nmeans = 1
      D = length(means)
      means = matrix( means, 1, D )
    }

    # Default value for the standard deviations
    if( is.null( sigmas ) | length(sigmas) == 1 ){
      if( is.null( sigmas ) ){
        lambda = 0.001
      }else{
        lambda = sigmas
      }

      # Get the sigmas matrix
      sigmas = matrix( 0, Nmeans, D^2 )

      # Fill as diagonals depending on dimension
      sigmas[ , 1 ]   <- rep( 1, Nmeans )
      sigmas[ , D^2 ] <- rep( 1, Nmeans )

      if( D == 3 ){
        sigmas[ , 5 ] <- rep( 1, Nmeans )
      }
      # Scale the sigma matrix appropriately
      sigmas = sigmas * diff( range( x ) ) * lambda
    }

    # Get the basis function
    f = vapply( 1:Nmeans,
                function( l ){
                   mvtnorm::dmvnorm( x,
                                     mean = t( means[l,] ),
                                     sigma = matrix( sigmas[l,], D, D ) )
                },
                FUN.VALUE = seq( 0, pi, length.out = dim( x )[1] )
            )
  }

  # Normalize the functions such that the resulting processes will have variance
  # one if the random multipliers have.
  if( normalize ){
    f <- f / sqrt( rowSums( f^2 ) )
  }

  # Return RandomField object
  return( RandomField( field = f %*% randf( Nmeans, N ), locations = x ) )

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

  samp.field = sin( pi / 2 * x ) %*% t( rnorm( N, 0, 1 ) ) +
               cos( pi / 2 * x ) %*% t( rnorm( N, 0, 1 ) )

  # Return RandomField object
  return( RandomField( field = samp.field, locations = x ) )

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
    Nx <- dimx <- length( x )
  }else{
    # get the grid from the coordinates
    grid.vals = coords2grid( x )
    D  = length( grid.vals )
    Nx = dimx[1]

    # Get the dimension of the grid
    dimx = unlist( lapply( grid.vals, length ) )
  }

  # Create the array of noise to be smoothed
  y = array( randf( Nx, N ), dim = c( dimx, N ) )

  # Create the sample fields by smoothing the noise field
  samp.field = aws::kernsm( y, h = h, kern = kern, unit = unit )

  # Return RandomField object
  return( RandomField( field = samp.field@yhat, locations = x ) )

}

