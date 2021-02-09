#------------------------------------------------------------------------------#
#                                                                              #
#     Testing the SampleFields package                                         #
#                                                                              #
#------------------------------------------------------------------------------#
# load package
require( SampleFields )
require( plotly )
require( tidyverse )

#----- Test signal plus noise model
N = 10
# Default 1D
x = seq( 0, 2, length.out = 100 )
Y <- SignalPlusNoise( N = N, x = x )
plot( Y )

# Default 2D
x = seq( 0, 2, length.out = 50 )
x = expand.grid( x, x )
Y <- SignalPlusNoise( N = N, x = x )
plot( Y )
plot( Y, surface = TRUE )

# Default 1D with predefined mean and sigma
x = seq( 0, 2, length.out = 100 )
Y <- SignalPlusNoise( N = N,
                      x = x,
                      mu = mu.example1.1d,
                      sigma = var.example1.1d,
                      obs.noise = 2 )
plot( Y )

# Default 2D with predefined mean and sigma
x = seq( 0, 2, length.out = 50 )
x = expand.grid( x, x )
Y <- SignalPlusNoise( N = N,
                      x = x,
                      mu = mu.prod,
                      sigma = var.example1.2d,
                      obs.noise = 2 )
fig = plot( Y, ncols = 2 )

plot( Y, surface = TRUE )


#----- Test some of the error process models
# square exp
x = seq( 0, 1, length.out = 200 )
Y = ArbCovProcess( N = N,
                   x = x,
                   covf = covf.square.exp,
                   h = 0.01 )
plot( Y )


# matern covariance
x = seq( 0, 1, length.out = 200 )
Y = ArbCovProcess( N = N,
                   x = x,
                   covf = covf.nonst.matern,
                   params = c( 3 / 2, 1 / 2, 1 ) )
plot( Y )


# 2D RandomNormalSum
x  = seq( 0, 1, length.out = 50 )
grid = expand.grid(x,x)
lx = dim(grid)[1]
Y = RandomNormalSum( N = N,
                     expand.grid(x,x),
                     normalize = TRUE,
                     means = lx*0.3,
                     sigmas = 0.002 )
plot( Y, ncols = 2 )
plot( Y, surface = TRUE )
