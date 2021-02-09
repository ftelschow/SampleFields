#------------------------------------------------------------------------------#
#                                                                              #
#     Different covariance functions of Random Processes                       #
#                                                                              #
#------------------------------------------------------------------------------#
# Contained functions:
#  - covf.square.exp
#  - covf.nonst.matern
#
#------------------------------------------------------------------------------#
# Developer notes:
# - add the descriptions
#------------------------------------------------------------------------------#

#' Covariance function: Square exponential
#'
#' This implements the square exponential covariance function, i.e.,
#' *c( x, y ) = exp( -( x - y )^2 / 4 / h^2 )*
#'
#' @param x Numeric. First argument of cov( x, y ).
#' @param y Numeric. First argument of cov( x, y ).
#' @param h Numeric. Scaling parameter.
#'
#' @return Value of the covariance function evaluated at (x, y).
#'
#' @export
covf.square.exp <- function( x, y, h = 0.01, FWHM = FALSE ){
  if( FWHM ){
    h = h / sqrt( 8*log(2) )
  }
  exp( -( x - y )^2 / 4 / h^2 )
}


#' Modified Matern Covariance Function (varying roughness parameter)
#'
#' This implements the square exponential covariance function, i.e.,
#' *c( x, y ) = exp( -( x - y )^2 / 4 / h^2 )*
#'
#' @inheritParams covf.square.exp
#' @inherit covf.square.exp return
#'
#' @param params Vector with 3 elements: c(nu1, nu2, sigma).
#'
#' @export
covf.nonst.matern <- function( x,
                               y,
                               params = c( 3 / 2, 1 / 2, 1 )
                              ){

  nu    <- params[1] + sqrt( max( x, y ) ) * ( params[2] - params[1] )
  sigma <- params[3]
  l     <- 1
  d     <- sqrt( 2 * nu ) * abs( x - y ) / l
  if( d > 0 ){
    sigma^2 * 2^( 1 - nu ) / gamma( nu ) * d^nu * besselK( d, nu )
  }else{
    sigma^2
  }
}

