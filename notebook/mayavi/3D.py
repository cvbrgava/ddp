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

#file_read = open('./dump.txt','r')
#solncheck = pickle.load( file_read )
#file_read.close()

#x_check ,y_check , z_check  = split_cords ( solncheck )

#x_MOR ,y_MOR , z_MOR  = split_cords ( yact )

x_linP ,y_linP, z_linP = split_cords( linP )

x_MOR = [ linPRed[ i ][ 0 ] for i in range( len( linPRed ) ) ]
y_MOR = [ linPRed[ i ][ 1 ] for i in range( len( linPRed ) ) ]
MOR_points = [ X_ * redBasis[ :, 0 ] + Y_ * redBasis[:,1 ] for X_,Y_ in zip( x_MOR,y_MOR ) ]

linPRed_x, linPRed_y, linPRed_z = split_cords( MOR_points )

error = [ numpy.linalg.norm(numpy.array( (x_linP[ i ]  - linPRed_x[ i ] , y_linP[ i ]  - linPRed_y[ i ] ,z_linP[ i ]  - linPRed_z[ i ]  )) ) for i in range( len( x_linP )  )]
normLinp = [ numpy.linalg.norm( numpy.array( (x_linP[ i ] , y_linP[ i ] , z_linP[ i ] ) ) ) for i in range( len( x_linP ) ) ] 
percError = [ error[ i ] * 100 / normLinp[ i ]  for i in range( len ( error ) ) ]



#table = prettytable.PrettyTable( ["Point", "Percentage Error"] )
#for linpt in range( len( x_linP) -1) :
#    table.add_row( [linpt+1 , percError[ linpt ] ] )
#print table

x_plane, y_plane = numpy.mgrid[ 1:5,3:6 ]
planeMOR = surf( x_plane, y_plane, plane_z , color = ( 0.3 ,0.3 , 1 ) ,opacity = 0.3)
axes()
#ax.label_text_property.font_family = 'times'
#ax.label_text_property.font_size = 3
#text('Ideal Differential Amplifier 3D state space' )
outline( planeMOR )

TPWL = plot3d(x, y, z, tube_radius = 0.005,color = (1,0,0) )
#plot3d(x_check ,y_check , z_check, tube_radius = 0.00015 )

points3d( linPRed_x[ : 50], linPRed_y[ : 50], linPRed_z[ : 50],color = ( 0,1,0), scale_factor = 0.025 , opacity = 0.5)
points3d( x_linP[ : 50] ,y_linP[ : 50], z_linP[ : 50] , color = (1, 1 , 0 ), scale_factor = 0.025)


show()

