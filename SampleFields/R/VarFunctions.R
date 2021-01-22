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
#' @param x Integer amount of realisations of the random field to be generated.
#' @param h Vector locations at which the random field is evaluated.
#' @return Array containing the realisations of the random field.
#'
#' @export
var.example1.1d <- function( x = seq( 0, 1, length.out = 100 ), h = 1/6 ){
  h * ( ( 0.6 - x )^2 + 1 )
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
var.example1.2d <- function( x = seq( 0, 1, length.out = 100 ), y = NULL, h = 1 ){

  # use the x-values as y-values if nothing else specified
  if( is.null(y) ){
    y = x
  }

  # compute the variance function
  h * array( outer( x,
                    x,
                    FUN = function( s, t ) ( s + 1 ) / ( t^2 + 1 ) ) / 3,
             c( rep( length(x), 2 ), 1 ) )
}
