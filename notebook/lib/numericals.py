import numpy
import sympy

from sympy.parsing.sympy_parser import parse_expr
from parsers import *
from config import *
from getparameters import Sym2NumArray
from progressbar import *

class Point(ProgressBarWidget):
	def update(self,pbar):
		return '%2d' % pbar.currval
		
class matrices( object ) :
	def __init__( self , name ):
		self.name = name
		self.linPJac = sympy.zeros( (order,order) )
		self.linPVal = sympy.zeros( (order, 1) )
	
	def call( self ):
		print self.name



def init_statematrices(count,order,input_list,output_list, stateorder,state):
	''' Initialze all the state space matrices needed.
	Input : 
	
	count = number of linearization points
	order = order of the system
	input_list = list of the inputs to the system
	ouput_list = list of the outputs of the system
	stateorder = symbolic matrix with the state space ordering
	Returns : 
	
	linPJac = Jacobian of the system at some linearization point
	B = input matrix
	C = Output matrix
	linPVal = Value of the system at some linearization point 
	inputm = value of the systems's input at some linearizaton point'''
	C = numpy.zeros( ( len( output_list ), order ) )
	C_ones = [i for i in range( len( stateorder ) ) if str( stateorder[ i ] ) in output_list ]
	for i , j in zip( range( len( output_list ) ) , C_ones) :
		C[ i , j ] = 1
	D = [ 0 ]

	inputm = sympy.Matrix( numpy.zeros( (3*len(input_list),1 ) ))
	for k in range( len(input_list) ) :
		for i in range(  3 ) :
        		inputm[ i + ( 3 * k ) ] = state[ str( input_list[ k ] ) ] ** (i+1)
	#offset=numpy.ones((count,len(input_list)*3,1))
	#diff=numpy.ones((order,1))
	#proj=numpy.zeros((order,order))
	#linPJac = numpy.zeros( ( count, order, order) )
	#linPVal = numpy.zeros( ( count, order, 1) )
	B = numpy.zeros( (count, order, 3*len( input_list ) ) )
	
	return  C , B , D , inputm


def get_statematrices(  B, C, D, inputm, count, order, regions, state, regexp, datapoints, stateorder, inputorder, Vth): 
	''' Initialze all the state space matrices needed.
	Input : 
		count = number of linearization points
		order = order of the system
		regions = list of all the regions of operation of the transistor
		datapoints = snapshots of the sytem at various linearization points
		stateorder = symbolic matrix with the state space order to be followed
	Returns : 
		linPJac = Jacobian of the system at some linearization point
		B = input matrix
		C = Output matrix
		linPVal = Value of the system at some linearization point '''
	widgets = ['Building State matrices @ Point',Point(),Bar('.')]
	pbar = ProgressBar(widgets = widgets,maxval = count,term_width = 60).start()	
	
	Matrix = []
	for linpt in range(count):
		pbar.update( linpt + 1 )
		Matrix.append( matrices( linpt ) )
		eqs = get_nonlinear_matrix(state,regions[linpt],Vth)
		#Matrix[ linpt ].linPVal = eqs
		for i in range(order):
			B[linpt][i], string = parse_nonlinear( str( eqs[ i ] ), regexp) 
			eqs[i]=sympy.parsing.sympy_parser.parse_expr(string, state)
		
		#sub1={state[i]:abs(float(datapoints[linpt][i])) for i in state}
		
	    
		# dx/dt = Jacobian * ( state ) + B * ( input ) + ( Val @ state - Jacobian * linPoint )
	
	   
	
		#linPVal[ linpt ] = Sym2NumArray( eqs.evalf(subs = sub1) ) 
        	Matrix[ linpt ].linPVal = eqs

		jack = eqs.jacobian( stateorder )
		#linPJac[ linpt ] =  jack.evalf( subs = sub1 )  )        
		Matrix[ linpt ].linPJac =  jack
		#inputm[ linpt ] = Sym2NumArray( inputorder.evalf( subs = sub1 ) )
		#inputm[ linpt ] = inputorder
		#temp=numpy.matrix(B[linpt])
		#proj=numpy.linalg.pinv(temp)
	
		#diff=(linPVal[linpt]-numpy.dot(linPJac[linpt],linP[linpt]))
		#offset[linpt]=proj*diff
	pbar.finish()
	return Matrix, B, C, D

# The following functions are used for integration of the PWL model

def get_parameters_integration( initialcond, intg_end, stateorder ):
	print "PWL Integration initialized...."
	y0 = numpy.array( [ (float( initialcond[ str( i ) ]) )for i in stateorder ])
	time  = numpy.linspace(0, intg_end, 100000)
	return y0, time

