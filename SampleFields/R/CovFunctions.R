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

#' Covariance function: Square exponential variable smoothness
#'
#' This implements the square exponential covariance function, i.e.,
#' *c( x, y ) = exp( -( x - y )^2 / 4 / h(t) / h(s) )*
#'
#' @param x Numeric. First argument of cov( x, y ).
#' @param y Numeric. First argument of cov( x, y ).
#' @param h Numeric. Scaling parameter.
#'
#' @return Value of the covariance function evaluated at (x, y).
#'
#' @export
covf.square.exp.nonstat <- function( x, y, h = c(0.1, 0.01),
                                     FWHM = FALSE, f1=F ){
  if( FWHM ){
    h = h / sqrt( 8*log(2) )
  }

  # exp( -( x - y )^2 / 4 / (h[1] + x*diff(h)) / (h[1] + y*diff(h)) )
  if(f1){
    return(exp( -( x - y )^2 / 4 / exp(log(h[2])*x+(1-x)*log(h[1]) + log(h[2])*y+(1-y)*log(h[1])) ))
  }else{
    hx = h[2] + (x-0.5)^2 * (h[1]-h[2]) * 4
    hy = h[2] + (y-0.5)^2 * (h[1]-h[2]) * 4
#    hy = (h[1]-h[2])*(cos(pi*y+pi)+1.)/2 + h[2]
    return(exp( -( x - y )^2 / 4 / hx / hy ))
  }
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

