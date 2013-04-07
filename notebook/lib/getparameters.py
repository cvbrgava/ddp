import csv
import sympy
import numpy
import prettytable
import re
import fileinput
from currenteqs import region

def get_newdata_heading( filename ) :
	for line in fileinput.input(filename,inplace =1 ):
		if fileinput.isfirstline():
			print ( get_newhead(line))[:-1]
		else :
			print ( line )[:-1]
	fileinput.close()

def get_newhead( line ):
	check = re.compile(r'(V\().*(\))')
	
	if re.search( check , line ):
	        temp = re.search( check, line )
		line = line.replace(temp.group(1),"")
        	line = line.replace(temp.group(2),"")
	return line

def import_text(filename, separator):
    ''' Import data from wave-form files exported from LTSPICE.
	Generates a python dictionary for each row it reads with the column-name as the key.
	Usage:
		import_text(<absolute_path>,'\t') use this for tab separated columns '''
	
    reader= csv.DictReader(open(filename), delimiter=separator, skipinitialspace=True)
    for line in reader:
        if line:
            yield line

def get_steadystate(filename,initialcond,stateorder,final_time,order):
    ''' Provides the steady-state value for the given waveforms.
	
	current version however just takes the final-time ( obtained by visual inspection ) and provides the values at that time instant. 
	
	Intended use of the function is this

	if |x0 - xnew |/|x0| < error :
		count += 1
	if count >= some number :
		return x0 
	'''
	
    print "Getting steady state ...."
    pre=[float(initialcond[str(stateorder[k])]) for k in range(order)]
    count=0
    for data in import_text(filename,'\t'):
        curt=[float(data[str(stateorder[k])]) for k in range(order)]
        diff=numpy.linalg.norm(numpy.array(pre)-numpy.array(curt))
        pre=curt
        if(float(data['time'])== final_time):
            return curt
       
     
def get_states(filename,sim_begin):
	''' Given absolute path to the waveform files provides states needed for the system '''
	print "Getting states...."
	initialcond=get_initialcond(filename,sim_begin)
	state={}
	state['0']=sympy.Symbol('0')
	for ind in initialcond.keys():
		state[ind]=sympy.Symbol(ind)
	return state

def get_initialcond(filename, sim_begin):
	''' Returns:
	The initial condition given the absolute path to wave form files '''
	print "Getting initial conditions..."
	for data in import_text(filename,'\t'):
		if(float(data['time']) == sim_begin):
			initial = data
			initial['0']= 0
			return initial	

def get_linpdiff(initial,steady,stateorder,order):
    ''' Provides the norm of the difference between the steady state and initial condition in state space 
	
		Returns |x(initial) - x(steadystate)| = delta '''
    init=[float(initial[str(stateorder[k])]) for k in range(order) ]

    return numpy.linalg.norm(numpy.array(init)-steady)/numpy.linalg.norm(init)

def get_linPoints(filename,initialcond,delta,order,stateorder,sim_begin):
	'''Choice of the points of linearization is done as follows
	if |xnew - xi|/|xi| < delta :
		continue
	else :
	    x(i+1) = xnew
	Returns the number of linearization points and the time's at which they were chosen
	'''
	pre = numpy.array( [ float( initialcond[ str(stateorder[k]) ] ) for k in range(order) ] )
	count = 1
	time = []
	time.append( float( sim_begin ) )
	table = prettytable.PrettyTable( [ "Point", "|xnew - xi|/|xi|" ,"Time" ] )
	for i in import_text(filename,'\t'):
		cur = numpy.array( [ float(i[ str( stateorder[k] ) ] ) for k in range(order) ] )
		diff = abs( numpy.linalg.norm(pre-cur) / numpy.linalg.norm(pre) )
		if(diff >= delta):
			pre = cur
			count = count+1
			table.add_row( [count, diff, i['time'] ] )
			time.append( float(i['time']) )
	print "The Linearization points chosen if |xnew - xi|/|xi| > ", delta	
	print "DC operating point is chosen as the first linearization point"
	print table
	return count,time
        
def Sym2NumArray(F):
	''' After evaluation of an expression in sympy, convert that array into a numpy readable array '''
	shapeF = F.shape
	B=numpy.zeros( shapeF )
	for i in range( 0 , shapeF[0] ):
		for j in range(0 , shapeF[1] ):
			B[i,j] = sympy.N( F[ i, j] )
	return B

def get_datapoints( file_voltage, count, time):
    ''' Get snapshot of the system at required time instants 
    Returns :
        Dictionary with state variable as key and value at the time instant as value'''
    datapoints=[{} for i in range(count)]
    for data in import_text(file_voltage,'\t'):
        for i in range(count):
            if(float(data['time'])==time[i]):
                datapoints[i]=data
	for i in range(count):
		datapoints[ i ]['0'] = 0
    return datapoints

def get_currents( file_current, count ,time):
	''' Get snapshots of the currents through vairous transistors in the system.
	Returns:
		Dictionary containing transistor label as key and current at the time instant as the value'''
	currents = [ {} for i in range(count) ]
	for data in import_text(file_current,'\t'):
		for i in range(count):
			if(float(data['time'])==time[i]):
				currents[i]=data
	return currents

def get_inputorder( state,inp_list):
	'''Returns:
	symbolic non-linear input vector 
	inp[0]
	inp[0]**2
	inp[0]**3
	inp[1]
	inp[1]**2
	.
	.	
	.	
	'''
	check = sympy.Matrix( numpy.zeros( (3*len(inp_list),1) ) )
	for i in range(len(inp_list) ):
		check[3*i] = state[inp_list[i]]**1
		check[3*i + 1] = state[inp_list[i]]**2
		check[3*i + 2] = state[inp_list[i]]**3
	return check

def get_num_linP(count,order,datapoints,stateorder):
	'''Returns:
	numpy compatible array for the relavent state space values at linearization points'''
	linP=numpy.zeros((count,order,1))
	for i in range(count):
		for k in range(order):
			linP[i][k]=(datapoints[i][str(stateorder[k])])          
	return linP

