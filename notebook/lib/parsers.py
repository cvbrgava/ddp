import re
import sympy
from currenteqs import region
from config import input_list, constant_dict 
# There are two parsers used in the project. These parsers use regular expressions for sequence matching followed by manipulation.
# 1. Net-list parser :
# 	Used to find where transistors are connected in the circuit. This information can be used to find the regions of operations. Which wil be used while substituting 
#	appropriate current equations in the state space matrix
# 2. Equation parser :
# 	Used to parse through the state space equations and separate the non-linear input terms in it. These terms are then stored in a different matrix ( input matrix B )




def get_region_operation( datapoints, file_netlist, count):
	''' Parses through the netlist and finds transistors in the circuit

	Calculates the regions of operations of these transistors at linearization points

	Returns:

	list of dictionary with each transistor as key and it's region of operation as value'''
	print "Parsing netlist...."
	regions = [ {} for i in range(count) ]	
	nl = open(file_netlist)
	x = re.compile('M(?#comment)')
	for point in range(count):
		a = 0
		while True:
			i = nl.readline()    
			if re.match(x,i):
				a = a+1
				inf = i.split()
				sma = [low.lower() for low in i.split()]
				r = region(sma[5].upper(),float(datapoints[point][sma[1]]),float(datapoints[point][sma[2]]),float(datapoints[point][sma[3]]))
				regions[point][sma[0]] = [r,sma[1],sma[2],sma[3],sma[5].upper(),sma[7],sma[6]]
			if not i:           
				nl.seek(0)
				break
	print "This circuit has",a,"transistors"
	return regions

def get_regexp_eqs( inp_list ):
	''' For a given set of inputs of the system, returns a list of regular expression rules which can be used to identify and isolate non linear input terms 
	returns:	
	
	[input[0]**1,input[0]**2,input[0]**3,input[1].....,input[2]...] '''
	regexp = []
	for i in inp_list :
		regexp.append( re.compile(r"([-+]?[ ]?\d+\.\d+)."+i+" "))
		regexp.append( re.compile(r"([-+]?[ ]?\d+\.\d+)."+i+"\*\*2 "))
		regexp.append( re.compile(r"([-+]?[ ]?\d+\.\d+)."+i+"\*\*3 "))
	return regexp		

def parse_nonlinear(string, regexp):
	''' 1. Parses the non-linear equation 	
		2. Searches for non-linear input terms
		3. Isolates their coefficients
	Returns :
		The row of B matrix and the remainder of the equation. '''
	B_row = [ 0 ]*len(regexp)
	for i in range( len( regexp ) ):
		if re.search( regexp[ i ], string) :
			temp = re.search( regexp[ i ], string)
# To check while debugging uncomment the next line to see which coefficient is being removed						
			#print temp.group(1)            
			B_row[ i ] = float((temp.group(1)).replace(" ",""))
			string = string.replace((temp.group()),"")
	return B_row, string


def parse_within( linpt, eqs, order, state, datapoints ) :
	''' This function extracts the input terms from the Matrices. This makes sure there are no redundant calculations during the Integration '''
	for Order in range( order ):
		polyDict = sympy.collect( eqs[ Order ], [ state[ z ] for z in input_list] ,evaluate = False)
		
		sub1 = {state[ k ]: (float(datapoints[ linpt ][ k ])) if str( state[ k ] ) not in constant_dict.keys() else constant_dict[ str( k ) ] for k in state}
	       	
		polyDict = {k:polyDict[ k ].evalf( subs = sub1 ) for k in polyDict.keys()} 

    		temp = 0
    		for k in polyDict.keys():
        		temp +=  k*polyDict[ k ] 
		eqs[ Order ] = temp
    	
	return eqs


