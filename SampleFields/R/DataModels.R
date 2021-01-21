#------------------------------------------------------------------------------#
#                                                                              #
#     Functions generating samples of different random fields                  #
#                                                                              #
#------------------------------------------------------------------------------#
# Required packages:
#   - abind
#   - spatstat
#   - MASS
#
# Contained functions:
#   - SignalPlusNoise
#
#------------------------------------------------------------------------------#
# Developer notes:
# - Style guideline NOT included
#------------------------------------------------------------------------------#

#' Generate samples from signal-plus-noise model
#'
#' This function generates samples from a signal-plus-noise model
#' *insert formula*
#' The population mean, variance, error field and observation noise can be
#' seperately defined as functions or predefined functions from this package
#' can be used.
#' @param N Integer amount of realisations of the random field to be generated.
#' @param x Vector/Array locations at which the random field is evaluated. Must be a valid input for the functions 'mu', 'noise' and 'sigma'.
#' @param mu Function computing the population mean of the signal plus noise model from a vector of locations.
#' @param noise Function computing a sample of random fields.
#' @param sigma Function computing the mean from a vector.
#' @param sd_ObsNoise Numeric positiv real number giving the standard deviation for independent additiv Gaussian observation noise.
#' @param ... additional parameters for the 'noise' function
#'
#' @export
SignalPlusNoise <- function( N,
                             x     = seq( 0, 1, length.out = 100 ),
                             mu    = NULL,
                             noise = GaussDensitySumNoise,
                             sigma = NULL,
                             obs.noise = 0,
                             ... ){

  #---- If mu is NULL use the function returning 0
  if( is.null( mu ) ){
    mu <- function( N, x ){ 0 }
  }else if( !is.function( mu ) ){
    stop( "sigma must be a function. See help for more information." )
  }

  #---- If sigma is NULL use the function returning 1
  if( is.null( sigma ) ){
    sigma <- function( N, x ){ 1 }
  }else if( !is.function( sigma ) ){
    stop( "sigma must be a function. See help for more information." )
  }

  #---- If obs.noise is a numeric define as default independent normal
  #     observation noise
  if( is.numeric( obs.noise ) ){
    if( obs.noise > 0 ){
      sd.obs <- obs.noise
      obs.noise <- function( N, x ){
        array( rnorm( N * prod( mdim ), mean = 0, sd = sd.obs ),
               dim = c( mdim, N ) )
      }
    }else if( obs.noise == 0 ){
      obs.noise <- function( N, x ){ 0 }
    }else{
        stop( "If obs.noise is a numeric, it must be non-negative." )
    }

  }else if( !is.function( obs.noise ) ){
      stop( "obs.noise must be a non-negative numeric or a function." )
  }

  # Get the population mean
  m = mu( x )
  if( is.vector( m ) ){
    mdim = length( x )
    m    = matrix( m, mdim, N )
  }else{
    mdim = dim( m )
    m    = array( rep( m, N ), c( mdim, N ) )
  }

  # Output the samples
  return( m +
          sigma( N = N, x = x ) * noise( N = N, x = x, ... ) +
          obs.noise( N = N, x = x ) )
}
