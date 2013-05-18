import numpy
from mayavi.mlab import *
import pickle
import prettytable
global redBasis





def plane_z( x, y ):
	plane_normal = numpy.cross( redBasis[ :, 0 ] , redBasis[ : , 1 ] )	
	return (- plane_normal[ 0 ] * x - plane_normal[ 1 ] * y )/( plane_normal[ 2 ])

def split_cords( given ):
	x = [ given[ i ][ 0 ] for i in range( len( given ) ) ]
	y = [ given[ i ][ 1 ] for i in range( len( given ) ) ]
	z = [ given[ i ][ 2 ] for i in range( len( given ) ) ]
	return x,y,z


file_read = open( './solnPWL_AOI.txt','r')
solnPWL, redBasis, linP, linPRed = pickle.load( file_read )
file_read.close()
x, y,z = split_cords( solnPWL )
TPWL = plot3d(x, y, z, tube_radius = 0.005,color = (1,0,0) )
