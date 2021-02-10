#------------------------------------------------------------------------------#
#                                                                              #
#     Miscaleaneous functions                                                  #
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
#' @export
coords2grid <- function( x ){
  # Check the input of x is a matrix
  if( !is.matrix( x ) && !is.list( x ) ){
    stop( "x must be a matrix" )
  }

  # Get the dimension of x
  dimx <- dim( x )
  D = dimx[ 2 ]

  # Pre-allocate the the ve
  X <- dX <- list()

  # get the coordinates for each dimension
  for( i in 1:D ){
    X[[ i ]]  <- unique( x[ ,i ] )
    dX[[ i ]] <- diff( X[[ i ]] )

    # Don't proceed if the coordinates are not regularly spaced
    if( any( abs( dX[[ i ]] - dX[[ i ]][1] ) > 1e-12 ) ){
      stop( "These coordinates do not define a regular grid." )
    }
  }

  # expand the unique vectors to a grid
  if( D == 2 ){
    xx = expand.grid( X[[ 1 ]], X[[ 2 ]] )
  }else{
    xx = expand.grid( X[[ 1 ]], X[[ 2 ]], X[[ 3 ]] )
  }

  if( all( unlist( abs( xx - x ) ) < 1e-12 ) ){
    return( X )
  }else{
    stop( "These coordinates do not define a regular grid." )
  }

}


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
