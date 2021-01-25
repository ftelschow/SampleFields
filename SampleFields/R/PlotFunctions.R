#------------------------------------------------------------------------------#
#                                                                              #
#     Functions creating plots of RandomFields                                 #
#                                                                              #
#------------------------------------------------------------------------------#
# Required packages:
# - plotly
# - ggplot
#
# Contained functions:
# -
#------------------------------------------------------------------------------#
# Developer notes:
# - need update of documentation
#------------------------------------------------------------------------------#

#' Computes simultaneous confidence bands for the mean of a sample from a one dimensional
#' functional signal plus noise model. It is possible to choose between different estimators
#' for the quantile.
#'
#' @param Y array of dimension K_1 x ... x K_d x N containing N-realizations of
#' a random field over a d-dimensional domain.
#' @param x Numeric the targeted covering probability. Must be strictly between 0 and 1.
#'
#' @return list with elements
#'  \itemize{
#'   \item hatmean pointwise sample mean
#'   \item scb list containing the upper and lower bounds of the simultaneous confidence band
#'   \item level targeted covering probability
#'   \item q quantile of the maximum of the residual field
#' }
#' @export
sample2tibble <- function( Y, x ){
  N = dim( Y )[2]
  Y.tibble <- as_tibble( Y ) %>% add_column( location = x )
  pivot_longer( Y.tibble,
                cols      = names( Y.tibble )[1:N],
                names_to  = "Sample",
                values_to = "value" )
}


#' Multiple plot function
#'
#' ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
#' - cols:   Number of columns in layout
#' - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#'
#' If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
#' then plot 1 will go in the upper left, 2 will go in the upper right, and
#' 3 will go all the way across the bottom.
#' @export
multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c( list(...), plotlist )

  numPlots = length( plots )

  # If layout is NULL, then use 'cols' to determine layout
  if( is.null( layout ) ){
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix( seq(1, cols * ceiling( numPlots / cols ) ),
                      ncol = cols, nrow = ceiling( numPlots / cols ) )
  }

  if( numPlots == 1 ){
    print( plots[[1]] )

  } else {
    # Set up the page
    grid.newpage()
    pushViewport( viewport( layout = grid.layout( nrow( layout ),
                                                  ncol( layout ) )
    ) )

    # Make each plot, in the correct location
    for( i in 1:numPlots ){
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame( which( layout == i, arr.ind = TRUE ) )

      print( plots[[i]], vp = viewport( layout.pos.row = matchidx$row,
                                        layout.pos.col = matchidx$col ) )
    }
  }
}

#' Covariance function: Square exponential
#'
#' This implements the square exponential covariance function, i.e.,
#' *c( x, y ) = exp( -( x - y )^2 / 4 / h^2 )*
#'
#' @param x Numeric. First argument of cov( x, y ).
#' @param y Numeric. First argument of cov( x, y ).
#' @param h Numeric. Scaling parameter.
#'
#' @return Value of the covariance function evaluated at (x, y).
#'
plot <- function( x, ... ){
  UseMethod( "plot" )
}

#' Covariance function: Square exponential
#'
#' This implements the square exponential covariance function, i.e.,
#' *c( x, y ) = exp( -( x - y )^2 / 4 / h^2 )*
#'
#' @param x Numeric. First argument of cov( x, y ).
#' @param y Numeric. First argument of cov( x, y ).
#' @param h Numeric. Scaling parameter.
#'
#' @return Value of the covariance function evaluated at (x, y).
#'
plot.RandomField <- function( rf, surface = FALSE, ... ){
  if( rf$D == 1 ){
    rf.tib = sample2tibble( rf$values, rf$locations ) %>% group_by( Sample )

    ggplot( rf.tib ) + geom_line( aes( x = location, y = value, color = Sample ) )

  }else if( rf$D == 2 ){

    coords = coords2grid( rf$locations )

    z.vals = matrix( rf$values[,1],
                     length( coords[[1]] ),
                     length( coords[[2]] ) )

    fig <- plotly::plot_ly( z = ~ z.vals )

    if( surface ){
      fig <- fig %>% plotly::add_surface()
    }

    fig

  }
}
