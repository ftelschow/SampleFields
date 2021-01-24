#------------------------------------------------------------------------------#
#                                                                              #
#     Functions generating examples of random coefficient functions            #
#                                                                              #
#------------------------------------------------------------------------------#
# Required packages:
#
# Contained functions:
# - rnorm.mod
#
#------------------------------------------------------------------------------#
# Developer notes:
# - need documention
#------------------------------------------------------------------------------#

#' Example mean function in 1D
#'
#' Function computing the population mean curve. The function is
#' *insert formula*
#' and is used in Telschow & Schwartzman(2021).
#'
#' @param K Integer: output number of rows.
#' @param N Integer: output number of columns.
#'
#' @return Matrix of dimension K x N containing standard normal random
#'  variables.
#'
#' @export
rnorm.mod <- function( K, N ){
  return( matrix( rnorm( K * N ), K, N ) )
}
