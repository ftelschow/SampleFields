#------------------------------------------------------------------------------#
#                                                                              #
#     Functions generating examples of variance functions                      #
#                                                                              #
#------------------------------------------------------------------------------#
# Required packages:
#
# Contained functions:
# - var.example1.1d
# - var.example1.2d
#
#------------------------------------------------------------------------------#
# Developer notes:
# - Style guideline NOT included
#------------------------------------------------------------------------------#

#' Example variance function in 1D
#'
#' Function computing the variance. The function is
#' *insert formula*
#' and is used in Telschow & Schwartzman(2021).
#'
#' @inheritParams ArbCovProcess
#'
#' @param h Numeric a scaling constant for the variance.
#' @return Vector containing the the variance values for x.
#'
#' @export
var.example1.1d <- function( x = seq( 0, 1, length.out = 100 ), h = 1 / 6 ){
  abs( h ) * ( ( 0.6 - x )^2 + 1 )
}

#' Example variance function in 2D
#'
#' Function computing the variance. The function is
#' *insert formula*
#' and is used in Telschow & Schwartzman(2021).
#'
#' @param x Integer amount of realisations of the random field to be generated.
#' @param h Vector locations at which the random field is evaluated.
#' @return Array containing the realisations of the random field.
#'
#' @export
var.example1.2d <- function( x = expand.grid( seq( 0, 1, length.out = 50 ),
                                              seq( 0, 1, length.out = 50 ) ),
                             h = 1 ){

  if( dim( x )[2] != 2 ){
    break( "Error, this function requires x to be a Tx2 matrix." )
  }

  # compute the variance function
  h * apply( x, 1,  FUN = function( r ){ ( r[1] + 1 ) / ( r[2]^2 + 1 ) / 3 } )

}
