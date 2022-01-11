#------------------------------------------------------------------------------#
#                                                                              #
#     Functions generating examples of mean functions                          #
#                                                                              #
#------------------------------------------------------------------------------#
# Required packages:
#
# Contained functions:
# - mu.example1.1d
# - mu.prod
#------------------------------------------------------------------------------#
# Developer notes:
# - need update of documentation
#------------------------------------------------------------------------------#

#' Example mean function in 1D
#'
#' Function computing the population mean curve. The function is
#' *insert formula*
#' and is used in Telschow & Schwartzman(2021).
#'
#' @inheritParams ArbCovProcess
#' @inherit ArbCovProcess return
#'
#' @export
mu.example1.1d <- function( x = seq( 0, 1, length.out = 100 ) ){
  if( !is.null( dim( x ) ) ){
    break( "For this function x must be a vector." )
  }
  sin( 8 * pi * x ) * exp( -3 * x )
}

#' Example mean function in 2D
#'
#' Function computing the population mean curve. The function is
#' *insert formula*
#' and is used in Telschow & Schwartzman(2021).
#'
#' @param x Integer amount of realisations of the random field to be generated.
#'
#' @return Array containing the realisations of the random field.
#'
#' @export
mu.prod <- function( x = expand.grid( seq( 0, 1, length.out = 50 ),
                                      seq( 0, 1, length.out = 50 ) ) ){

  # Compute the mean function
  apply( x, 1, function( y ) prod( y ) )

}


#' Example mean function in 2D given by Gaussian bumbs
#'
#' Function computing the population mean curve. The function is
#' *insert formula*
#' and is used in Telschow & Schwartzman(2021).
#'
#' @param x Integer amount of realisations of the random field to be generated.
#'
#' @return Array containing the realisations of the random field.
#'
#' @export
mu.prod <- function( x = expand.grid(seq(0, 1, length.out = 50),
                                     seq(0, 1, length.out = 50)),
                     locs   = rbind(expand.grid(c(1/6, 3/6, 5/6), c(1/6, 3/6, 5/6))),
                     spread = t(apply(locs, 1, function(x) c(1, 0, 1)))
                     ){

  function(x, x0, Sigma){
    exp( t(x - x0) %*% Sigma %*% (x - x0) )
  }
  mvtnorm::dmvnorm()
  # Compute the mean function
  apply( x, 1, function( y ) prod( y ) )

}
