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

#' redefine * for input of RandomFields
#'
#' GRF() returns samples from a Gaussian related field
#' $$ Y(x) = F( X_1,...,X_N ) $$
#' are computed.
#' The population mean, variance, error field and observation noise can be
#' seperately defined as functions or predefined functions from this package
#' can be used.
#'
#' @inheritParams RandomBasisSum
#' @inherit RandomBasisSum return
#'
#' @param mu Function computing the population mean of the signal plus noise
#' model from x.
#' @param noise Function computing a sample of the noise field Z.
#' @param sigma Function computing the standard deviation from x.
#' @param sd_ObsNoise Numeric positiv real number giving the standard deviation
#' for independent additiv Gaussian observation noise.
#' @param ... additional parameters for the 'noise' function
#'
#' @export
"%+%" = function( f, g ){
  if( is( f, "RandomField" ) & is( g, "RandomField" ) ){
    out = f
    out$values = f$values + g$values
  }else if( is( f, "RandomField" ) & is.numeric(g) ){
    out = f
    out$values = f$values + g
  }else if( is.numeric(f) & is( g, "RandomField" ) ){
    out = g
    out$values =  f + g$values
  }

  return( out )
}

#' redefine * for input of RandomFields
#'
#' GRF() returns samples from a Gaussian related field
#' $$ Y(x) = F( X_1,...,X_N ) $$
#' are computed.
#' The population mean, variance, error field and observation noise can be
#' seperately defined as functions or predefined functions from this package
#' can be used.
#'
#' @inheritParams RandomBasisSum
#' @inherit RandomBasisSum return
#'
#' @param mu Function computing the population mean of the signal plus noise
#' model from x.
#' @param noise Function computing a sample of the noise field Z.
#' @param sigma Function computing the standard deviation from x.
#' @param sd_ObsNoise Numeric positiv real number giving the standard deviation
#' for independent additiv Gaussian observation noise.
#' @param ... additional parameters for the 'noise' function
#'
#' @export
"%**%" = function( f, g ){
  if( is( f, "RandomField" ) & is( g, "RandomField" ) ){
    out = f
    out$values = f$values * g$values
  }else if( is( f, "RandomField" ) & is.numeric(g) ){
    out = f
    out$values = f$values * g
  }else if( is.numeric(f) & is( g, "RandomField" ) ){
    out = g
    out$values =  f * g$values
  }

  return( out )
}
