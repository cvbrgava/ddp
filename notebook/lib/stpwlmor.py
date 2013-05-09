from numpy import dot as nD
import numpy
import scipy

global redOrder 
redOrder = 3

def get_redBasis( moment, order, redOrder):
	'''Returns basis after Moment Matching '''

	spanAgg = numpy.concatenate( linP , 1 )

	#input_signals = get_input_signals( t ) 

	#sub1 = { state[ str( input_list[ k ] ) ] : input_signals[ k ] for k in range( len( input_list ) ) }

	linPVal=[  Matrix[ i ].linPVal for i in range( count )]

	linPJac=[ Matrix[ i ].linPJac for i in range( count ) ]

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
	
	return redBasis#, linPValRed, linPJacRed

def get_U_linp( Matrix ):
	''' Assumes stable system and returns a stable rowspace for MOR using lyapunov stability criterion. Builds lyapunov function using -I '''
	V = get_redBasis( 2 , order , redOrder )
	P = [ scipy.linalg.solve_sylvester( ( Matrix[ linpt ].linPJac ).T , ( Matrix[ linpt ].linPJac ) , -numpy.eye( order ) ) if max( numpy.array ( numpy.linalg.eigvals( Matrix[ linpt ].linPJac) ) ) < 0 else numpy.eye( order ,order ) for linpt in range( count ) ]
	U = [ ( nD( numpy.linalg.inv(nD( nD( V.T, P[ linpt ] ) , V ) ) , nD( V.T , P[ linpt ] ) ) ).T for linpt in range( count ) ]
	return U

def obliqueprojection( rowspace, columnspace, points, fullD = False):
	''' Returns the result of the oblique projection with a given rowspace and columnpsace of the projection. The input can be a single vector or a list of vectors. 
	The result is returned in lower dimension unless full dimension is asked for'''
	proj = nD( numpy.linalg.inv( nD( rowspace.T , columnspace ) ) , rowspace.T )
	if ( type( points ) == list ) :
		if not fullD:
			return [ nD( proj , point ) for point in points ]
		else :
			return [ nD( columnspace, nD( proj , point ) ) for point in points ]
	else:
		if not fullD:
			return nD( proj , points )
		else :
			return nD( columnspace, nD(proj , points ) )
        
def get_mu( y , U , V , linP, redOrder , count ):
	''' Obliquely projects all the linearization points on the same subspace and calculated mu '''
	linPRed = obliqueprojection( U , V , list( linP ) )
	norm = numpy.array([ 10.0**(-20*numpy.linalg.norm( linPRed[ i ] - y.reshape(redOrder,1) )) for i in range(count) ])
	mu = norm/norm.sum()
	return mu


def get_U( mu , U ) :
	''' Returns cummulative U from Piecewise U's obtained '''
	temp = 0
	for linpt in range( count ):
		temp = temp + mu[ linpt ] * U[ linpt ]
	return temp

def check_stab( Matrix ):
	''' Checks the stability of the Piecewise jacobians '''

	if 1.0 in [ sign ( max( numpy.array (numpy.linalg.eigvals( Matrix[ i ].linPJac ) ) )  ) for i in range( count ) ] :
		print " Unstable"
		return False
	else:
		return True


def get_y0_MOR( y ):
	''' Function to get the initial condition in the required basis'''
	UlinP = ( get_U_linp( Matrix ) )[ 0 ]
	VlinP = get_redBasis( 2 , order , redOrder )    
	return obliqueprojection( UlinP , VlinP , y )

def get_MOR_Jac( Matrix , B ,linP,  count):
	''' Returns the order reduced Jacobian and new B '''
	Ulinp = get_U_linp( Matrix )
	V = get_redBasis( 2 , order , redOrder )

	linPJac_Red = [ [ nD( nD((Ulinp[k]).T, Matrix[ i ].linPJac ) , V )for i in range( count  ) ] for k in range( count ) ]

	Bnew = [ numpy.concatenate(( B[ linpt ] , Matrix[ linpt ].linPVal - (numpy.dot( Matrix[ i ].linPJac, linP[ i ] )).reshape( order ,1 ) ),1 ) for linpt in range( count ) ]

	BRed = [ [  nD( ( Ulinp[ k ] ).T, Bnew[ i ] ) for i in range( count  ) ] for k in range( count ) ]
	return linPJac_Red , BRed




class STPWLMOR( object ):
	'''Contains all the information needed to carry out integration using STPWL-MOR algorithm on all Hurwitz system. '''
	def __init__( self ):        
		self.newBasis = numpy.matrix( numpy.eye(order,redOrder) )
		self.UlinP = get_U_linp( Matrix )
		self.oldBasis = self.UlinP[ 0 ] 
		self.V = get_redBasis( 2 , order , redOrder )
		self.JacRed , self.BRed = get_MOR_Jac(Matrix, B , linP, count )
		self.linPRed = [ obliqueprojection( self.UlinP[ linpt ] , self.V , linP[ linpt ] ) for linpt in range( count ) ]
    
	def derv(self, y , t ):
		''' Returns derivative of the state space at a given time-instant '''
        
		mu = get_mu(y ,self.oldBasis, self.V , linP ,redOrder , count )
		U = get_U( mu , self.UlinP )
		#print mu
		self.oldBasis = U
		#linPRed = obliqueprojection( U , V , list( linP ), fullD = True )
		#yNew = obliqueprojection( U , V , y , fullD = True)
		temp = 0 
		sub1 = { state[ str( input_list[ k ] ) ]  : (get_input_signals( t ))[ k ] for k in range( len( input_list ) ) }
		inputM = numpy.concatenate( (Sym2NumArray( inputm.evalf( subs = sub1) ) , numpy.matrix( ([ 1.0 ]) ) ) , 0) 
		#print inputM
		weight = self.weights( y, mu )
		mutemp = 0 
		for Mu in range( count ):
			temp = 0
			for wei in range(count):
				temp = temp +  weight[ wei ] * ( (numpy.dot( self.JacRed[ Mu ] [ wei ] , y.reshape( redOrder , 1 ) )).reshape( redOrder ,1 )  + (numpy.dot( self.BRed[ Mu ][ wei ], inputM ) ) )
                
			mutemp = mutemp + mu[ Mu ] * numpy.array( temp )
        
		return  ( ( nD( numpy.linalg.inv( nD( U.T , self.V ) ) , mutemp ) ).reshape( 1, redOrder )  ) [ 0 ]
    
	def check_stab( self ):
		''' Checks if the MOR of the system provides us with stable Jacobians'''
		if 1 in [  sign( max( numpy.array( numpy.linalg.eigvals(self.JacRed[ k ][ i ] ) ) ) )for i in range( count ) for k in range( count ) ]:
			#print "There is an unstable MOR Jacobian "
			return False
		else:
			#print " System is stable "
			return True  
    
	def get_STPWL_weights_Bound( self, mu ):
		'''Calculates the numerical bound on the Weighing functions'''
		P = [ scipy.linalg.solve_sylvester( ( Matrix[ linpt ].linPJac ).T , ( Matrix[ linpt ].linPJac ) , -numpy.eye( order ) ) for linpt in range( count ) ]
		Phat = [ nD( nD( (self.V).T, P[ linpt ] ) , self.V ) for linpt in range( count ) ]
		MaxTerm = []
		MinTerm = []
		wStar = []
		for m in range( count ):
			rowSpan = range( count ) 
			rowSpan.remove( m )
			singluarValues = numpy.array( [ ( numpy.array( [ max( numpy.linalg.svd( nD( Phat[ m ] , self.JacRed[ k ] [ i ] ) , compute_uv= False ) ) for i in rowSpan ]) ).sum() for k in range( count ) ] )
			MaxTerm.append( numpy.dot ( numpy.array( mu ) , singluarValues ) )
			MinTerm.append( mu[ m ] * min ( numpy.linalg.svd( nD( Phat[ m ] , self.JacRed[ m ][ m ] ) , compute_uv=False ) ) )
			wStar.append( MaxTerm[ m ]  / ( MaxTerm[ m ]  + MinTerm[ m ]  ) )
			#print MinTerm [ m ] 
		return wStar
        
	def get_Beta( self, y , mu ):
		''' Returns the Beta function required to construct the weighting function'''
        
		wStar = self.get_STPWL_weights_Bound( mu )
		norm = numpy.array([ numpy.linalg.norm( self.linPRed[ i ] - y.reshape(redOrder,1) ) for i in range(count) ])
		for Norm in range( len( norm ) ):
			if min( norm ) == norm[ Norm ]:
				m = Norm
		deltaM = 0.5 * max( numpy.array( [numpy.linalg.norm(   self.linPRed[ m ] - self.linPRed[ linpt ] ) for linpt in range( count ) ] ) )
		return log( wStar[ m ] )/ deltaM**2
    
	def weights( self, y, mu  ):
		''' Returns stabilizing weights for the STPWL '''
		beta = self.get_Beta( y, mu  )
		norm = numpy.array([ numpy.e**(beta*numpy.linalg.norm( self.linPRed[ i ] - y.reshape(redOrder,1) )) for i in range(count) ])
		return norm / norm.sum()
                
