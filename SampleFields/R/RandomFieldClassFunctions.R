#------------------------------------------------------------------------------#
#                                                                              #
#     Functions for RandomField class                                          #
#                                                                              #
#------------------------------------------------------------------------------#
# Required packages:
#
# Contained functions:
#  - Create random field
#
#------------------------------------------------------------------------------#
# Developer notes:
# -
# -
#------------------------------------------------------------------------------#
#' Regular grid from coordinates
#'
#' This function takes a T x D matrix of coordinates and checks whether the
#' points lie on a regular grid. If so it outputs the coordinate vectors for the
#' grid.
#'
#' @param x Matrix of dimension T x D containing coordinates. D > 1 required.
#'
#' @return list containing the coordinates in each direction.
#'
#'
#' @export
getDim <- function( x ){
  # Get the dimension of x
  dimx <- dim( x )

  #
  if( is.null( dimx ) ){
    D    = 1
    nloc = length( x )
  }else{
    D    = dimx[2]
    nloc = dimx[1]
  }

  # return a list containing the dimension and number of locations
  return( list( D = D, nloc = nloc ) )

}
#' RandomField constructs an object of class RandomField from input data
#'
#' Function computing the population mean curve. The function is
#' *insert formula*
#' and is used in Telschow & Schwartzman(2021).
#'
#' @param field Integer: output number of rows.
#' @param location Integer: output number of columns.
#'
#' @return RandomField class object
#'
#' @export
RandomField <- function( field, locations ){

  # Get the dimension of the field and locations matrix
  dimRF = dim( field )
  gDim  = getDim( locations )

  # Check whether the input fits together
  if( dimRF[1] != gDim$nloc ){
    break( "Error: First dimension of field must agree with number of locations in the domain." )
  }

  # Create the output class object
  RF = list( values = field, locations = locations )

  # Add additional important variables to RandomField class
  RF$dim    = dimRF
  RF$D      = gDim$D
  RF$nloc   = gDim$nloc
  RF$N      = dimRF[2]

  # Make RF a S3 class
  class( RF ) = "RandomField"

  return( RF )
}
