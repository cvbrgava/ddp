import prettytable
import numpy
from getparameters import *
from config import *

global ratio ,moments,linP,state
 
ratio = 100
moments = 2 
redOrder = 2




def y0MOR(y,t):
	'''Generates the lower dimension initial conditions'''
	global state

	input_signals = get_input_signals( t ) 

	sub1 = { state[ str( input_list[ k ] ) ] : input_signals[ k ] for k in range( len( input_list ) ) }

	redBasis , linPValRed, linPJacRed = get_redBasis(  t , moments , order,redOrder )
        
	
	redBasisPro = numpy.dot( numpy.linalg.inv(numpy.dot( redBasis.T, redBasis ) ) , redBasis.T ) 

	return numpy.dot( redBasisPro , y )

def LD2HD ( t, y ):
	'''Used to convert from Lower dimension to Higher dimension'''

	redBasis, linPValRed, linPJacRed = get_redBasis(  t , moments , order,redOrder)

	return numpy.dot( redBasis, y )

def get_redBasis(  t , moment, order,redOrder):
	'''Returns basis after Moment Matching '''

	spanAgg = numpy.concatenate( linP , 1 )

	input_signals = get_input_signals( t ) 

	sub1 = { state[ str( input_list[ k ] ) ] : input_signals[ k ] for k in range( len( input_list ) ) }

	linPVal=[ Sym2NumArray( ( Matrix[ i ].linPVal ).evalf( subs = sub1 ) )  for i in range( count )]

	linPJac=[ Sym2NumArray( ( Matrix[ i ].linPJac ).evalf( subs = sub1 ) ) for i in range( count ) ]

	for linpt in range( count ):
		spanV1 = [ numpy.dot( numpy.linalg.inv( linPJac[ linpt ] ),B[ linpt ]) ]
        
		for k in range( moment ):
			spanV1.append( numpy.dot( numpy.linalg.inv( linPJac[ linpt ] ),spanV1[ -1 ]) )

		spanV1 = numpy.concatenate( spanV1, 1 )
        

		spanV2 = [ numpy.dot( numpy.linalg.inv( linPJac[ linpt ] ),linPVal[ linpt ] - (numpy.dot( linPJac[ linpt ], linP[ linpt ] )).reshape( order ,1 ) )]
		for k in range( moment ):
			spanV2.append( numpy.dot( numpy.linalg.inv( linPJac[ linpt ] ),spanV2[ -1 ]) )
		spanV2 = numpy.concatenate( spanV2, 1 )	
		
		unitSpanV1 , R = numpy.linalg.qr( spanV1 )

		unitSpanV2 , R = numpy.linalg.qr( spanV2 )

	spanAgg = numpy.concatenate( ( spanAgg, unitSpanV1 , unitSpanV2) , 1 )

	eigenC , eigenV , eigenR = numpy.linalg.svd(spanAgg)
        
	redBasis = eigenC[ : , :redOrder ]
	
	linPValRed = [numpy.dot( redBasis.T,linPVal[ k ]) for k in range( count ) ] 
	
	linPJacRed = [numpy.dot( redBasis.T, numpy.dot( linPJac[ k ] , redBasis ) ) for k in range( count )]
	
	return redBasis, linPValRed, linPJacRed

def dervPWLMOR(yNew,t):
	''' Produces reduced order derivative at any point in time'''

	temp = 0

	spanAgg = numpy.concatenate( linP , 1 )

	input_signals = get_input_signals( t ) 

	sub1 = { state[ str( input_list[ k ] ) ] : input_signals[ k ] for k in range( len( input_list ) ) }

	redBasis, redOrder , linPValRed, linPJacRed = get_redBasis(  t , moments , order )
	
	redBasisPro = numpy.dot( numpy.linalg.inv(numpy.dot( redBasis.T, redBasis ) ) , redBasis.T ) 
	
	#yNew = numpy.dot( redBasisPro, y )
        linPRed = [ numpy.dot( redBasisPro, linP[ i ] ) for i in range( count ) ]

        BRed = [numpy.dot( redBasis.T, B[ k ] ) for k in range( count )]

	inputM = Sym2NumArray(inputm.evalf( subs = sub1 ))
        
	norm = numpy.array([ 10.0**(-20*numpy.linalg.norm( linPRed[ i ] - yNew.reshape(redOrder,1) )) for i in range(count) ])

        weights = norm/norm.sum()
        #print weights
        
           
    	for i in range(count):
		temp = temp +  weights[ i ] * ( linPValRed[ i ] + (numpy.dot( linPJacRed[ i ], yNew )).reshape( redOrder ,1 ) - (numpy.dot( linPJacRed[ i ], linPRed[ i ] )).reshape( redOrder ,1 ) + (numpy.dot( BRed[ i ], inputM )) )
    
	return  (temp.reshape(1, redOrder))[0] 




def check_stability_MOR(moments, count, state):
	''' Checks the stability of the Jacobians created after reduction with a certain number of moments matched '''
	global linP

	x = prettytable.PrettyTable( )

	x.add_column("Linearization Point",[ k for k in range( count ) ] )
	#redOrder = 10
	#moments = 2

	for red_order in range( 1 , order +1 ):

		Jaccheck = []

		for linpt in range( count  ):

			redOrder = red_order

			redBasis, linPValRed, linPJacRed = get_redBasis( time[ linpt ] - sim_begin , moments ,order ,redOrder)
	
			redBasisPro = numpy.dot( numpy.linalg.inv(numpy.dot( redBasis.T, redBasis ) ) , redBasis.T ) 
	
			linPRed = [ numpy.dot( redBasisPro, linP[ i ] ) for i in range( count ) ]
	    
			norm = numpy.array([ 10.0**(-20*numpy.linalg.norm( linPRed[ i ] - linPRed[ linpt] )) for i in range(count) ])

			weights = norm/norm.sum()

			temp = 0 
	   
			for j in range(count):

				temp = temp +  weights[ j ] * (  linPJacRed[ j ] )

			Jaccheck.append( temp )

		x.add_column("Order:"+str(redOrder),[ numpy.real(numpy.linalg.eigvals(Jaccheck[ k ] )).max() for k in range( count ) ])

	print moments," moments matched:"

	print x

