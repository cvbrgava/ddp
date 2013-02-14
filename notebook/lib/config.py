import sympy
from currenteqs import current

# Add the name of the circuit being simulated here
cirname = 'Differential Amplifier'

# The functions in this library change with the circuit being used. 
# All the state space equations go into the get_nonlinear_matrix() function
# The state space ordering is done by the get_stateorder() function
print "Working on",cirname

# The following parameters are to be customized for the circuit being simulated
# These parameters are the bridge between the project and SPICE.
#-------------------------------------------------------------------------------------------- 
order = 13
sim_begin = 2e-6
sim_end = 4.796677620601137e-006
file_voltage = './Data/voltage_diffamp.txt'
file_current = './Data/currents.txt'
file_netlist = './Data/diffamp.net'
input_list = ['innegative','inpositive'] 
output_list = ['out1pos','out1neg']
denominator = 4
intg_end = 2e-7
#--------------------------------------------------------------------------------------------


def get_nonlinear_matrix(state,regions,Vth):
	''' 
	1.This function is instance specific. 
	2.All the state equations are defined in here
	Input:
		
	state = symbolic vector with states defined it
	regions = contains information regarding each transistor with the region of operation specified
	Vth = dictionary with threshold voltages of each transistor
	
	Returns :
	The symbolic expression with non-linear state space equations in it.'''
	
	# Define all the parasitic capacitors in here. Generally we assume that they are constant. 
	# Name the parasitics accordingly
	#-----------------------------------------------------------------------------------------------------
	Cpara_vdo=1e-15
	Cpara_vd12=1e-15
	Cpara_vd11=1e-15
	Cpara_vcmfb2=1e-15
	Cpara_vom=1e-15
	Cpara_vop=1e-15
	Cpara_vdc3=1e-15
	Cpara_vdc5=1e-15
	Cpara_vcmfb1=1e-15
	Cpara_vd1=1e-15
	Cpara_vd2=1e-15
	#-----------------------------------------------------------------------------------------------------	
	# This provides a dictionary with
	# KEY	:	Name of the transistor
	# VALUE	:	Polynomial expansion of current equation at a particular linearization point 
	# To access any current equation use I[ 'name of the transistor' ] 
	
	I={i: current(regions[i],state,Vth[i]) for i in regions.keys()}

	
	# The equations from here are all circuit dependant. Add your equations from here
	# NOTE:
	
	# To access any current through a transistor use I[ 'name of transistor' ], we assume this current is the DRAIN current
	#-----------------------------------------------------------------------------------------------------
	
#    IIII = ( ( I[ 'm15' ] - I[ 'm13' ] ) / Cpara_vcmfb2 ) - ( ( I[ 'm10' ] - I[ 'm9' ] ) / Cpara_vom )
#    II = ( ( I [ 'm15' ] - I[ 'm13' ] ) / Cpara_vcmfb2 ) - ( ( I[ 'm8' ] - I[ 'm7' ] ) / Cpara_vop ) 
#    cur=sympy.Matrix( ( ( IIII ) , ( II ) ) )
#    cap=sympy.Matrix( ( ( ( 1/1e-12 + 1/Cpara_vcmfb2 + 1/Cpara_vom ) , ( 1/Cpara_vcmfb2 ) ) , ( 1/Cpara_vcmfb2 , 1/1e-12 + 1/Cpara_vcmfb2 + 1/Cpara_vop ) ) )
	[IC10,IC9]=[0,0]#cap.inv()*cur
	
#    III = ( ( I['mc2'] - I[ 'mc4' ] ) / Cpara_vcmfb1 ) - ( ( state['ig9'] + I[ 'm4' ] - I[ 'm2' ] ) / Cpara_vd2 ) 
#    IV = ( ( I[ 'mc2' ] - I[ 'mc4' ] ) / Cpara_vcmfb1 ) - ( ( state[ 'ig7' ] + I[ 'm3' ] - I[ 'm1' ] ) / Cpara_vd1 )
#    cap=sympy.Matrix( ( ( 1/7e-12 + 1/Cpara_vd2 + 1/Cpara_vcmfb1 , 1/Cpara_vcmfb1 ) , (Cpara_vcmfb1 ,1/7e-12 + 1/Cpara_vd1 + 1/Cpara_vcmfb1 ) ) )
#    cur=sympy.Matrix( ( ( III ) , ( IV ) ) )
	[IC4,IC2]=[0,0]#cap.inv()*cur
	
	vd0dot=( - I[ 'm1' ] - I[ 'm2' ] + I[ 'm0' ] ) / Cpara_vdo   #V Im0 V Im1 VIm2 chkd
	vd12dot=( -I[ 'm12' ] + I[ 'm14' ] ) / Cpara_vd12  #V Im14 V im12 chkd
	vd11dot=( I[ 'm12' ] + I[ 'm13' ] - I[ 'm11' ] ) / Cpara_vd11    #V Im12 V Im13 V Im11 chkd
	vcmfb2dot=( I[ 'm15' ] - I[ 'm13' ] - IC10 - IC9 ) / Cpara_vcmfb2  # V Im15 <IC9 >IC10 V Im13
	vomdot=( I[ 'm10' ] + IC10 - I[ 'm9' ] - state[ 'ig9' ] ) / Cpara_vom   # > IC10 V Im10 V Im9  < ig9 chkd
	vopdot=( I[ 'm8' ] + IC9 - I[ 'm7' ] - state[ 'ig7' ] ) / Cpara_vop          # < IC9 V Im8 V Im7  > ig7 chkd
	vdc3dot=( I[ 'mc1' ] + I[ 'm6' ] - I[ 'mc3' ] ) / Cpara_vdc3                 # V Imc1  V Im6 V Imc3 
	vdc5dot = ( - I[ 'mc1' ] - I[ 'm6' ] - I[ 'mc2' ] + I[ 'mc5' ] ) / Cpara_vdc5   # V Imc5 V Imc1 V Im6 V Imc2
	vcmfb1dot = ( I[ 'mc2' ] - I[ 'mc4' ] - IC2 - IC4 ) / Cpara_vcmfb1 #V Imc2 < IC2 >IC4 V Imc4
	vd1dot = ( state[ 'ig7' ] - I[ 'm3' ] + I[ 'm1' ] + IC2 ) / Cpara_vd1   # >ig7 V Im1 V Im3 < IC2
	vd2dot = ( state[ 'ig9' ] - I[ 'm4' ] + I[ 'm2' ] + IC4 ) / Cpara_vd2  # < ig9 V Im2 V Im4 > IC4
	ig9dot = 0.001 * ( vomdot - ( state[ 'ig9' ] / 7e-12 ) - vd2dot ) # <ig9
	ig7dot = 0.001 * ( vopdot - ( state[ 'ig7' ] / 7e-12 ) - vd1dot ) # > ig7
	
	#-----------------------------------------------------------------------------------------------------
	
	# The next line packs all the equations in a certain order which is binding. This order decides the state space co-ordinates
	# Ensure that the same ordering is maintained with the stateorder function.
	
	eqs=sympy.Matrix( [(vd0dot),(vd12dot),(vd11dot),(vcmfb2dot),(vomdot),(vopdot),(vdc3dot),(vdc5dot),(vcmfb1dot),(vd1dot),(vd2dot),(ig9dot),(ig7dot)] )
	    
	return eqs


# Define the time dependant input signal. 
def get_input_signals( t ):
	''' Returns:
	The values of all the inputs at any given time instant'''

	return [1.3,0.5]
		

def get_stateorder(state):
	''' Returns stateorder, which is necessary for the construction of the state 
	Check with the state space equations and ensure that the same ordering is being followed'''
	
	stateorder=sympy.Matrix([(state['n001']),(state['n002']),(state['n005']),(state['cfmb2']),(state['vom']),(state['vop']),(state['n009']),(state['n007']),(state['cmfb1']),(state['out1neg']),(state['out1pos']),(state['ig9']),(state['ig7'])])
	return stateorder

