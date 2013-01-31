import numpy
import sympy

from sympy.parsing.sympy_parser import parse_expr
from parsers import *
from config import *
from getparameters import Sym2NumArray

def init_statematrices(count,order,input_list,output_list, stateorder):
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
	#offset=numpy.ones((count,len(input_list)*3,1))
	#diff=numpy.ones((order,1))
	#proj=numpy.zeros((order,order))
	linPJac = numpy.zeros( ( count, order, order) )
	linPVal = numpy.zeros( ( count, order, 1) )
	B = numpy.zeros( (count, order, 3*len( input_list ) ) )
	inputm = numpy.zeros( (count, 3*len(input_list),1 ) )
	return linPJac, B, C, D, linPVal, inputm


def get_statematrices( linPJac, B, C, D, linPVal, inputm, count, order, regions, state, regexp, datapoints, stateorder, inputorder, Vth): 
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
	for linpt in range(count):
		print "Building matrices for point",linpt+1,"...."
		eqs = get_nonlinear_matrix(state,regions[linpt],Vth)
		
		for i in range(order):
			B[linpt][i], string = parse_nonlinear( str( eqs[ i ] ), regexp) 
			eqs[i]=sympy.parsing.sympy_parser.parse_expr(string, state)
		
		sub1={state[i]:abs(float(datapoints[linpt][i])) for i in state}
		
	    
		# dx/dt = Jacobian * ( state ) + B * ( input ) + ( Val @ state - Jacobian * linPoint )
	
	   
	
		linPVal[ linpt ] = Sym2NumArray( eqs.evalf(subs = sub1) ) 
        
		jack = eqs.jacobian( stateorder )
		linPJac[ linpt ] = Sym2NumArray( jack.evalf( subs = sub1 )  )        

		inputm[ linpt ] = Sym2NumArray( inputorder.evalf( subs = sub1 ) )
		#temp=numpy.matrix(B[linpt])
		#proj=numpy.linalg.pinv(temp)
	
		#diff=(linPVal[linpt]-numpy.dot(linPJac[linpt],linP[linpt]))
		#offset[linpt]=proj*diff
	return linPJac, B, C, D, linPVal, inputm 

# The following functions are used for integration of the PWL model

def get_parameters_integration( initialcond, intg_end, stateorder ):
	print "PWL Integration initialized...."
	y0 = numpy.array( [ (float( initialcond[ str( i ) ]) )for i in stateorder ])
	time  = numpy.linspace(0, intg_end, 10)
	return y0, time

