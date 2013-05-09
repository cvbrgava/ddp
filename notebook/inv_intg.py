from scipy.integrate import odeint
from currenteqs import *
import numpy
import matplotlib.pyplot as plt
def input_signal( t ):
	if t < 1e-6:
		return 5e6*t
	if t > 1e-6:
		return 5

def current1(information,vth):
    """Current equations verified for level 1 PMOS and NMOS
	The equations in this library have been incorporated after comparing results with LTSPICE models
	Any change in these equations has to be verified thoroughly with appropriate LTSPICE simulations.
	
	Returns:
	Symbolic equations with polynomial expansion """
    d=information[1]
    g=information[2]
    s=information[3]
    w=information[5]
    l=information[6]
    vth=float(vth)


    
    if(information[4]=='CMOSP'):
        if(information[0]=='c'):
            return 0
        elif(information[0]=='s'):
            cc=(43e-6) * (w / l) * ((g**2) + (1 - 2 * vth * 0.09) * (s**2) + (vth**2) + ((0.09 * vth - 1) * 2 * g * s) + (2 * g * vth) + (-2 * vth + 0.09 * vth**2) * s - (0.09 * d * vth**2) - (0.09 * d * g**2) -  (0.09 * d * s**2) + ( 2 * g * s * 0.09 * d) - (2 * g * vth * 0.09 * d) + ( 2* s * vth * 0.09 * d ) + (0.09 * s * g**2) + (0.09 * s**3) - (0.09 * 2 * g * s**2))
            return cc
        elif(information[0]=='l'):
            cc=(86e-6) * (w / l) * ((0.5 * s**2) - (0.5 * d**2) - (g * s) - (s * vth) + (g * d) + (d * vth) + (0.09 * 0.5 * s**3) - (0.09 * s * 0.5 * d**2) - (0.09 * g * s**2) - (s**2 * 0.09 * vth) + (2 * 0.09 * s * g * d) + (2 * 0.09 * s * d * vth) - (0.5 * d * 0.09 * s**2) + (0.5 * 0.09 * d**3) - (0.09 * g * d**2) - (0.09 * d**2 *vth))
            return  cc
        
    if(information[4]=='CMOSN'):
        if(information[0]=='c'):
            return 0
        elif(information[0]=='s'):
            cc = (0.5 * 285e-6) * (w / l) * ((g**2) + (1-2 * vth * 0.14) * (s**2) + (vth**2) + ((0.14 * vth - 1) * 2 * g * s) - (2 * g * vth) + (2 * vth - 0.14 * vth**2) * s + (0.14 * d * vth**2) + (0.14 * d * g**2) + (0.14 * d * s**2) - (2 * g * s * 0.14 * d) - (2 * g * vth * 0.14 * d) + (2 * s * vth * 0.14 * d) - (0.14 * s * g**2) - (0.14 * s**3) + (0.14 * 2 * g * s**2))
            return cc
        elif(information[0]=='l'):
            cc=(285e-6) * (w / l) * ((g * d) - (s * g) + (0.5 * s**2) - (vth * d) + (vth * s) - (0.5 * d**2) + (0.14 * g * d**2) - (2 * 0.14 * d * s * g) + (0.5 * 0.14 * d * s**2) - (0.14 * vth * d**2) + (2 * 0.14 * d * vth * s) - (0.5 * 0.14 * d**3) + (s**2 *g * 0.14) - (0.5 * 0.14 * s**3) - (0.14 * s**2 * vth) + (0.5 * 0.14 * s * d**2))
            return cc


def derv( y ,t ):
	I = {}

	regm1 = region( 'CMOSN',y,input_signal( t ),0 )
	regm2 = region( 'CMOSP',y,input_signal( t ),5 )

	I[ 'm1' ] = current1( [ regm1 ,y,input_signal( t ), 0 , 'CMOSN',0.72,0.36] ,  0.477)
	I[ 'm2' ] = current1( [ regm2 ,y,input_signal( t ), 5 , 'CMOSP',1.44,0.36]  , 0.477)
	
	return ( I[ 'm2' ] - I[ 'm1' ] - ( y / 100e3 ) ) * 1e15



# initial conditions

y0 = 5

t = numpy.linspace(0,10e-6,10)
soln = odeint( derv , y0 ,t )


plt.plot( t, soln)
plt.savefig('output.png',format='png')
