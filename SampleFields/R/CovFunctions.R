#------------------------------------------------------------------------------#
#                                                                              #
#     Different covariance functions to be included in generation of           #
#     functional data samples                                                  #
#                                                                              #
#------------------------------------------------------------------------------#
# Contained functions:
#  - covf.nonst.matern
#
#------------------------------------------------------------------------------#
# Developer notes:
# - Style guideline included
#------------------------------------------------------------------------------#

#' Modified Matern Covariance Function (varying roughness parameter)
#'
#' @param x1 First argument of cov(x1,x2). Caution: It is assumed that 0<=x1<=1.
#' @param x2 Second argument of cov(x1,x2). Caution: It is assumed that 0<=x2<=1.
#' @param params Covariance function parameters: params=c(nu1, nu2, sigma).
#' @export
covf.SquareExp <- function( x, y, h = 0.01 ){
  exp( -( x - y )^2 / 4 / h^2 )
}


#' Modified Matern Covariance Function (varying roughness parameter)
#'
#' @param x1 First argument of cov(x1,x2). Caution: It is assumed that 0<=x1<=1.
#' @param x2 Second argument of cov(x1,x2). Caution: It is assumed that 0<=x2<=1.
#' @param params Covariance function parameters: params=c(nu1, nu2, sigma).
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

