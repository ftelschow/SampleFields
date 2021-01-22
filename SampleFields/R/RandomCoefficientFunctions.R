#------------------------------------------------------------------------------#
#                                                                              #
#     Functions generating examples of random coefficient functions            #
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
rnorm.mod <- function( k, N ){
  return( rnorm( k * N, k, N ) )
}
