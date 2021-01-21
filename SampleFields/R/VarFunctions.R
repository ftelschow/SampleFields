#------------------------------------------------------------------------------#
#                                                                              #
#     Functions generating examples of variance functions                      #
#                                                                              #
#------------------------------------------------------------------------------#
# Required packages:
#
# Contained functions:
# - var.example1.2d
#
#------------------------------------------------------------------------------#
# Developer notes:
# - Style guideline NOT included
#------------------------------------------------------------------------------#

#'Creates a sample from an 2D random field, which is the sum of 2-D Gaussian densities..
#' @param x Integer amount of realisations of the random field to be generated.
#' @param h Vector locations at which the random field is evaluated.
#' @return Array containing the realisations of the random field.
#'
#' @export
var.example1.2d <- function( x, h = 1 ){
  h * array( outer( x,
                    x,
                    FUN = function( s, t ) ( s + 1 ) / ( t^2 + 1 ) ) / 3,
             c( rep( length(x), 2 ), 1 ) )
}
