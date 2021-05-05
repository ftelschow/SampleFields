#------------------------------------------------------------------------------#
#                                                                              #
#     Functions generating samples of different Models for random fields       #
#                                                                              #
#------------------------------------------------------------------------------#
# Required packages:
#
# Contained functions:
#   - SignalPlusNoise
#
#------------------------------------------------------------------------------#
# Developer notes:
# - write documentation, especially figure out how to get formulas
#------------------------------------------------------------------------------#

#' Generate Samples from Signal-plus-Noise Model
#'
#' SignalPlusNoise() returns samples from a signal-plus-noise model
#' $$ Y(x) = \mu(x) + \sigma(x) Z(x) $$
#' are computed.
#' The population mean, variance, error field and observation noise can be
#' seperately defined as functions or predefined functions from this package
#' can be used.
#'
#' @inheritParams RandomBasisSum
#' @inherit RandomBasisSum return
#'
#' @param mu Function computing the population mean of the signal plus noise
#' model from x.
#' @param noise Function computing a sample of the noise field Z.
#' @param sigma Function computing the standard deviation from x.
#' @param sd_ObsNoise Numeric positiv real number giving the standard deviation
#' for independent additiv Gaussian observation noise.
#' @param ... additional parameters for the 'noise' function
#'
#' @export
SignalPlusNoise <- function( N,
                             x     = seq( 0, 1, length.out = 100 ),
                             mu    = NULL,
                             noise = RandomNormalSum,
                             sigma = NULL,
                             obs.noise = 0,
                             ... ){
  dimx = dim( x )
  if( is.null( dimx ) ){
    lx = length( x )
  }else{
    lx = dimx[1]
  }

  #---- If mu is NULL use the function returning 0
  if( is.null( mu ) ){
    mu <- function( x ){ rep( 0, lx )  }
  }else if( !is.function( mu ) ){
    stop( "sigma must be a function. See help for more information." )
  }

  #---- If sigma is NULL use the function returning 1
  if( is.null( sigma ) ){
    sigma <- function( x ){ rep( 1, lx ) }
  }else if( !is.function( sigma ) ){
    stop( "sigma must be a function. See help for more information." )
  }

  #---- If obs.noise is a numeric define as default independent normal
  #     observation noise
  if( is.numeric( obs.noise ) ){
    if( obs.noise > 0 ){
      sd.obs <- obs.noise
      obs.noise <- function( N, x ){
        array( rnorm( lx * N, mean = 0, sd = sd.obs ),
               dim = c( lx, N ) )
      }
    }else if( obs.noise == 0 ){
      obs.noise <- function( N, x ){ 0 }
    }else{
        stop( "If obs.noise is a numeric, it must be non-negative." )
    }

  }else if( !is.function( obs.noise ) ){
      stop( "obs.noise must be a non-negative numeric or a function." )
  }

  # Output the samples
  samp.f = noise( N = N, x = x, ... )

  samp.f$values =  mu( x ) + sigma( x ) * samp.f$values + obs.noise( N = N, x = x )

  return( samp.f )
}


#' Generate Samples from Gaussian related fields
#'
#' GRF() returns samples from a Gaussian related field
#' $$ Y(x) = F( X_1,...,X_N ) $$
#' are computed.
#' The population mean, variance, error field and observation noise can be
#' seperately defined as functions or predefined functions from this package
#' can be used.
#'
#' @inheritParams RandomBasisSum
#' @inherit RandomBasisSum return
#'
#' @param mu Function computing the population mean of the signal plus noise
#' model from x.
#' @param noise Function computing a sample of the noise field Z.
#' @param sigma Function computing the standard deviation from x.
#' @param sd_ObsNoise Numeric positiv real number giving the standard deviation
#' for independent additiv Gaussian observation noise.
#' @param ... additional parameters for the 'noise' function
#'
#' @export
chi2_Field <- function( N,
                  x     = seq( 0, 1, length.out = 100 ),
                  df,
                  mu    = NULL,
                  noise = RandomNormalSum,
                  sigma = NULL,
                  obs.noise = 0,
                  ... ){

  values = matrix( NA, getDim( x )$nloc, N )

  # Loop over number of samples
  for( n in 1:N ){
    # get df independent samples
    Y = SignalPlusNoise( N = df, x = x, mu = mu, noise = noise,
                         sigma = sigma, obs.noise = obs.noise, ... )
    # Compute the chi2 field from the sample
    values[,n] = rowSums( Y$values^2 )
  }

  return( RandomField( field = values, locations = x ) )
}


#' Generate Samples from Gaussian related fields
#'
#' GRF() returns samples from a Gaussian related field
#' $$ Y(x) = F( X_1,...,X_N ) $$
#' are computed.
#' The population mean, variance, error field and observation noise can be
#' seperately defined as functions or predefined functions from this package
#' can be used.
#'
#' @inheritParams RandomBasisSum
#' @inherit RandomBasisSum return
#'
#' @param mu Function computing the population mean of the signal plus noise
#' model from x.
#' @param noise Function computing a sample of the noise field Z.
#' @param sigma Function computing the standard deviation from x.
#' @param sd_ObsNoise Numeric positiv real number giving the standard deviation
#' for independent additiv Gaussian observation noise.
#' @param ... additional parameters for the 'noise' function
#'
#' @export
t_Field <- function( N,
                    x     = seq( 0, 1, length.out = 100 ),
                    df,
                    mu    = NULL,
                    noise = RandomNormalSum,
                    sigma = NULL,
                    obs.noise = 0,
                    ... ){

  values = matrix( NA, getDim( x )$nloc, N )

  # Loop over number of samples
  for( n in 1:N ){
    # get df independent samples
    Y = SignalPlusNoise( N = df+1, x = x, mu = mu, noise = noise,
                         sigma = sigma, obs.noise = obs.noise, ... )
    # Compute the t field from the sample
    values[,n] = sqrt( df ) * Y$values[,1] / sqrt( rowSums( Y$values[,-1]^2 ) )
  }

  return( RandomField( field = values, locations = x ) )
}

#' Generate Samples from Gaussian related fields
#'
#' GRF() returns samples from a Gaussian related field
#' $$ Y(x) = F( X_1,...,X_N ) $$
#' are computed.
#' The population mean, variance, error field and observation noise can be
#' seperately defined as functions or predefined functions from this package
#' can be used.
#'
#' @inheritParams RandomBasisSum
#' @inherit RandomBasisSum return
#'
#' @param mu Function computing the population mean of the signal plus noise
#' model from x.
#' @param noise Function computing a sample of the noise field Z.
#' @param sigma Function computing the standard deviation from x.
#' @param sd_ObsNoise Numeric positiv real number giving the standard deviation
#' for independent additiv Gaussian observation noise.
#' @param ... additional parameters for the 'noise' function
#'
#' @export
F_Field <- function( N,
                    x     = seq( 0, 1, length.out = 100 ),
                    df,
                    mu    = NULL,
                    noise = RandomNormalSum,
                    sigma = NULL,
                    obs.noise = 0,
                    ... ){

  values = matrix( NA, getDim( x )$nloc, N )

  # Loop over number of samples
  for( n in 1:N ){
    # get df independent samples
    Y1 = SignalPlusNoise( N = df[1], x = x, mu = mu, noise = noise,
                          sigma = sigma, obs.noise = obs.noise, ... )
    Y2 = SignalPlusNoise( N = df[2], x = x, mu = mu, noise = noise,
                          sigma = sigma, obs.noise = obs.noise, ... )
    # Compute the t field from the sample
    values[,n] = df[2] / df[1] * rowSums( Y1$values^2 ) / rowSums( Y2$values^2 )
  }

  return( RandomField( field = values, locations = x ) )
}
