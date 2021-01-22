#------------------------------------------------------------------------------#
#                                                                              #
#     Functions generating examples of mean functions                          #
#                                                                              #
#------------------------------------------------------------------------------#
# Required packages:
#
# Contained functions:
#
#------------------------------------------------------------------------------#
# Developer notes:
# - Style guideline NOT included
#------------------------------------------------------------------------------#

#' Example mean function in 1D
#'
#' Function computing the population mean curve. The function is
#' *insert formula*
#' and is used in Telschow & Schwartzman(2021).
#'
#' @param x Integer amount of realisations of the random field to be generated.
#' @return Array containing the realisations of the random field.
#'
#' @export
mu.example1.1d <- function( x = seq( 0, 1, length.out = 100 ) ){
  sin( 8 * pi * x ) * exp( -3 * x )
}

#' Example mean function in 2D
#'
#' Function computing the population mean curve. The function is
#' *insert formula*
#' and is used in Telschow & Schwartzman(2021).
#'
#' @param x Integer amount of realisations of the random field to be generated.
#' @param h Vector locations at which the random field is evaluated.
#' @return Array containing the realisations of the random field.
#'
#' @export
mu.prod.2d <- function( x = seq( 0, 1, length.out = 100 ), y = NULL ){

  # use the x-values as y-values if nothing else specified
  if( is.null(y) ){
    y = x
  }

  # Compute the mean function
  outer(x,x, FUN = "*" )

}
