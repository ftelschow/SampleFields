#------------------------------------------------------------------------------#
#                                                                              #
#     Testing the SampleFields package                                         #
#                                                                              #
#------------------------------------------------------------------------------#
# load package
require( SampleFields )

#----- Test signal plus noise model
Y <- SignalPlusNoise( N = 20, x = seq( 0, 2, length.out = 200 ) )
dim( Y )

#----- Test signal plus noise model
