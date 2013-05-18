import sympy
from currenteqs import current

# Add the name of the circuit being simulated here
cirname = '4 input NAND-gate '

# The functions in this library change with the circuit being used. 
# All the state space equations go into the get_nonlinear_matrix() function
# The state space ordering is done by the get_stateorder() function
print "Working on",cirname

# The following parameters are to be customized for the circuit being simulated
# These parameters are the bridge between the project and SPICE.
#-------------------------------------------------------------------------------------------- 
order = 4
sim_begin = 0.0
sim_end = 3.500000000000000e-006
file_voltage = './Data/4inpNAND_voltage.txt'
file_current = './Data/4inpNAND_current.txt'
file_netlist = './Data/4inpNAND.net'
input_list = ['ainp']
output_list = ['n001']
denominator = 40
constant_dict = {'vdd':5 , '0' :0}
intg_end = 2e-6

#--------------------------------------------------------------------------------------------
# Check-list before running the program : 
#	1. All the details pertaining to the SPICE file
#	2. All the constant voltages in the SPICE simulation
#	3. State equations for the current circuit
#	4. State order for the current state space
#	5. Inputs specific to the current scenario ( in sync with the input list mentioned above )	


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
	Cpara_out = 1e-15
	Cpara_M1 = 1e-15
	Cpara_M2 = 1e-15
	Cpara_M3 = 1e-15

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
	   
	voutdot = (I[ 'm6' ] + I[ 'm5' ] + I[ 'm7' ] + I[ 'm8' ] - I[ 'm1' ] )	/ Cpara_out
	v2dot = ( I[ 'm1' ] - I[ 'm2' ] ) / Cpara_M1
	v3dot = ( I[ 'm2' ] - I[ 'm3' ] ) / Cpara_M2
	v4dot = ( I[ 'm3' ] - I[ 'm4' ] ) / Cpara_M3
	#-----------------------------------------------------------------------------------------------------
	
	# The next line packs all the equations in a certain order which is binding. This order decides the state space co-ordinates
	# Ensure that the same ordering is maintained with the stateorder function.
	
#	eqs=sympy.Matrix( [(vd0dot),(vd12dot),(vd11dot),(vcmfb2dot),(vomdot),(vopdot),(vdc3dot),(vdc5dot),(vcmfb1dot),(vd1dot),(vd2dot),(ig9dot),(ig7dot)] )
	eqs = sympy.Matrix( [(voutdot),(v2dot),(v3dot),(v4dot)] )	    
	return eqs


# Define the time dependant input signal. 
def get_input_signals( t ):
	''' Returns:
	The values of all the inputs at any given time instant'''
	if t <= 1e-6:
		return [ 5 - 4e6*t ]
	if t > 1e-6:
		return [ 1 ]

		

def get_stateorder(state):
	''' Returns stateorder, which is necessary for the construction of the state 
	Check with the state space equations and ensure that the same ordering is being followed'''
	stateorder=sympy.Matrix([(state['n001']),(state['n002']),(state['n003']),(state['n004'])])
	#stateorder=sympy.Matrix([(state['n001']),(state['n002']),(state['n005']),(state['cfmb2']),(state['vom']),(state['vop']),(state['n009']),(state['n007']),(state['cmfb1']),(state['out1neg']),(state['out1pos']),(state['ig9']),(state['ig7'])])
	return stateorder

