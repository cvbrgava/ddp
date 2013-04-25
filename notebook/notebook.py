import sympy
import os
from getparameters import *
from currenteqs import *
from parsers import *
from config import *
from numericals import *
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# This file is the topmost layer of the project. Only function level calls and data transfer is shown here. 
global linPJac, B, C, D, linPVal, inputm, count, linP, state
#
# These functions are used in the integration process
# Derivative of the function at any point is given by
# dx/dt = Summation( weight*(Jacobian * ( state ) + B * ( input ) + ( Val @ state - Jacobian * linPoint )))
# Calculates the weighing function depending on the proximity to linearizatoin points
def normcalc( y ) :
	norm = numpy.array([ 10.0**(-20*numpy.linalg.norm( linP[ i ] - y.reshape(order,1) )) for i in range(count) ])
	return norm/norm.sum()

def dervPWL(y,t):
	#order = len(y)
	weights = normcalc( y )
	#count = len(weights)
	temp = 0
    	
	linPVal = []
	linPJac = []
	
	for i in range( count ):
		input_signals = get_input_signals( t ) 
		sub1 = {state[ k ]: (float(datapoints[ i ][ k ])) if str( state[ k ] ) not in constant_dict.keys() else constant_dict[ str( k ) ] for k in state}
		
		for k in range( len( input_list ) ) :
			sub1[ state[ str( input_list[ k ] ) ] ] = input_signals[ k ]
		
		linPVal.append( Sym2NumArray( ( Matrix[ i ].linPVal ).evalf( subs = sub1 ) ) )
		linPJac.append( Sym2NumArray( ( Matrix[ i ].linPJac ).evalf( subs = sub1 ) ) )
	inputM = Sym2NumArray(inputm.evalf( subs = sub1 ))
        
    
	for i in range(count):
		temp = temp +  weights[ i ] * ( linPVal[ i ] + (numpy.dot( linPJac[ i ], y )).reshape( order ,1 ) - (numpy.dot( linPJac[ i ], linP[ i ] )).reshape( order ,1 ) + (numpy.dot( B[ i ], inputM )) )
	return (temp.reshape(1, order))[0]
#------------------------------------------------------------------------ARNOLDI--------------------------------------------------------
def y0MOR(y,t):
	global state

	input_signals = get_input_signals( t ) 

	sub1 = { state[ str( input_list[ k ] ) ] : input_signals[ k ] for k in range( len( input_list ) ) }

	redBasis , linPValRed, linPJacRed = get_redBasis(  t , moments , order,redOrder )
        
	
	redBasisPro = numpy.dot( numpy.linalg.inv(numpy.dot( redBasis.T, redBasis ) ) , redBasis.T ) 

	return numpy.dot( redBasisPro , y )

def LD2HD ( t, y ):
	

	redBasis, linPValRed, linPJacRed = get_redBasis(  t , moments , order,redOrder)

	return numpy.dot( redBasis, y )

def get_redBasis(  t , moment, order,redOrder):
	

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
	

	temp = 0

	spanAgg = numpy.concatenate( linP , 1 )

	input_signals = get_input_signals( t ) 

	sub1 = { state[ str( input_list[ k ] ) ] : input_signals[ k ] for k in range( len( input_list ) ) }

	redBasis, linPValRed, linPJacRed = get_redBasis(  t , moments , order ,redOrder)
	
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




def check_stability_MOR(count, state):
	global linP

	x = prettytable.PrettyTable( )

	x.add_column("Linearization Point",[ k for k in range( count ) ])

	for moments in range( 1,order ):
		for red_order in range( 2 , order +1 ):

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

			eigenV = [ numpy.real(numpy.linalg.eigvals(Jaccheck[ k ] )).max() for k in range( count ) ]
			checkStab = numpy.array( [ 1 if egV < 0 else 0 for egV in eigenV ] )
			if checkStab.sum() == count:
				stab = 'Stable'
				redOrder = red_order
				moment = moments
			else:
				stab = 'Unstable'
		    

			x.add_column("Order:"+str(redOrder)+" -- "+str(stab),eigenV)
			if stab == 'Stable' and red_order < order:
			
				break
		if stab == 'Stable' and red_order < order:
			break
		

	if stab == 'Stable' and red_order < order:
		print x
		print "Reducing order to", redOrder ,"with ",moments," moments matched"
		
	else :
		print "No STABLE reduced order model exists"
		redOrder = order
		moments = order
	return redOrder, moments

#------------------------------------------------------------------------ARNOLDI--------------------------------------------------------

def initialize():
    get_newdata_heading( file_voltage )
    if os.path.isfile('./config.py'):
        print "Config file in ", os.path.abspath('./config.py')
        print "-----------------------------------------------"
    else:
        print "Config file in PYTHONPATH"
        print "-----------------------------------------------"
    

initialize()

initialcond = get_initialcond(file_voltage, sim_begin)

state = get_states(file_voltage,sim_begin)

stateorder = get_stateorder(state)

steadystate = get_steadystate(file_voltage,initialcond,stateorder, sim_end,order)

delta = get_linpdiff(initialcond,steadystate,stateorder,order)

print "Simulation begins at ", sim_begin

count,time = get_linPoints(file_voltage,initialcond,delta/denominator,order,stateorder,sim_begin)

datapoints = get_datapoints(file_voltage, count , time)

currents = get_currents(file_current,count,time)

regions = get_region_operation(datapoints,file_netlist,count)

Vth = get_vth()

regexp = get_regexp_eqs(input_list)

inputorder = get_inputorder(state, input_list )

linP = get_num_linP(count,order,datapoints,stateorder)

C, B , D , inputm = init_statematrices(count,order,input_list,output_list,stateorder,state)

Matrix, B, C, D = get_statematrices( B, C, D, inputm, count, order, regions, state, regexp, datapoints, stateorder, inputorder, Vth)

y0 , time_intg = get_parameters_integration( initialcond, intg_end, stateorder )

time_check = numpy.linspace( 0, 1e-7 , 100)

#----------------------------------------------------------------------------------
print "Checking for stable reduced order models"

redOrder ,moments = check_stability_MOR(count, state) # Add error checking to make it complete. 

if redOrder < order :
	print "TPWL + MOR integration intiated....."
	solnPWLMOR = odeint(dervPWLMOR, y0MOR( y0 , 0  ), time_check)
	solnPWL =  [ LD2HD( time_check[ i ] , solnPWLMOR[ i ] ) for i in range( len( time_check) ) ]
else :
	print "TPWL integration intiated....."
	solnPWL = odeint(dervPWL, y0, time_check)	
#---------------------------------------------------------------------------------

print "Plotting results....."	
time_linp = [ float( i ) - sim_begin for i in time] 
C_ones = [i for i in range( len( stateorder ) ) if str( stateorder[ i ] ) in output_list ]
intg_end = time_check[ -1 ]
for k in range(len(C_ones)):
    plt.figure(k)
    calc, = plt.plot(time_check, numpy.array([ abs(float(solnPWL[i][C_ones[k]])) for i in range(len(time_check) )] ), '--')
    plt.ylabel(r'$V_{'+str(stateorder[C_ones[k]])+'}$' )
    plt.xlabel(r'$time$')
    x , y  = [], []
    for i in import_text(file_voltage,"\t"):
        if (float( i[ 'time' ] ) < intg_end+sim_begin) and (float( i[ 'time' ]  )> sim_begin) :
            y.append( i[ str( stateorder[ C_ones[ k ] ] ) ] ), x.append( float(i[ 'time' ])-sim_begin )
    plt.ylabel(r'$V_{'+str(stateorder[C_ones[k]])+'}$' )
    plt.xlabel(r'$time$')
    plt.title(str(cirname)+' output\n MOR from '+str(order)+' to '+str(redOrder)+', Moments matched '+str(moments))
    spice, = plt.plot(x,y)
    #plt.figsize(10,10)
    plt.figlegend([calc,spice],('Calculated','Spice'),'lower right')
    x_time , y_val = [], []
    for i in range( count   ):
        x_time.append(float( time_linp[ i ] ) )
        y_val.append( float(datapoints[ i ] [str(stateorder[ C_ones[ k ] ] )]) )
    plt.plot( x_time,y_val,'ro' )
plt.show()


#----------------------------------------------------------------------------------






