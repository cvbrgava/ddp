import numpy
from numpy import concatenate as nC
from numpy import dot as nD
from numpy.linalg import inv as nI

# This file contains all the functions generally used in control system analysis. 
# The Kalman Decomposition of the system to obtain minimal realization is also provided.
# Certain linear algebraic functions such as finding Null-space , Columnspace have been added for readability of code
# Special functions such as finding intersection of basis and removal of basis have been introduced.




def observability( C , A  ):
	''' Checks for the observability of the system. 
	Returns:
	Observability matrix 
	True : if full rank
	False : if rank deficient '''

	start = [ C ]
	order = (shape( A ))[ 0 ]
	for i in range( order - 1 ) :
		start.append( nD( start[ -1 ], A ) )
	obs = nC( start , 0 )
	if( matrix_rank( obs) < order ) :
		return obs, False
	else:
		return obs, True

def controlability( B, A ):
	''' Returns the controllability matrix and also checks its rank
	returns 
	False : if rank deficient
	True : if full rank '''
	start = [ B ] 
	order = (shape( A ))[ 0 ]
	for i in range( order -1 ):
		start.append( nD( A, start[ -1 ] ) )
	cont = nC( start,1 )
    	if ( matrix_rank( cont ) < order ):
		return cont, False
	else :
		return cont, True

def columnProj( column , vect ,ret ):
	''' Generates the column-space projection of the vector onto the space spanned by columns of column.
	 ret = 1 : returns the coefficients of the projection
	 ret = 0 : returns the projection in standard basis'''
	if ret == 0:
		return numpy.dot( column, numpy.dot(numpy.linalg.inv(numpy.dot(column.T, column )), numpy.dot( column.T,vect ) ) )
	else:
		return numpy.dot(numpy.linalg.inv(numpy.dot(column.T, column )), numpy.dot( column.T,vect ) ) 

def columspace( mat ):
	''' Generates the Unitary basis vectors which span the column space of the given matrix'''
	col,eigen,row = numpy.linalg.svd(mat)
	return col[:,:matrix_rank(mat)]

def nullspace( mat ):
	'''Generates the Unitary basis for the null-space of the given matrix'''
	col, eigen, row = numpy.linalg.svd( mat )    
	return (row.T)[:,(matrix_rank(mat)):]

def KalmanBasis( check_A, check_B, check_C ):
	''' This function also produces Kalman basis but the generation of basis is slightly abmiguous, use KalmanBasisNew for the latest methods for basis generation'''
	check_cont , bool_cont = controlability( check_B , check_A )
	check_obs, bool_obs = observability( check_C, check_A )

	range_cont = columspace( check_cont ) 
	null_obs = nullspace( check_obs )
	inter = numpy.concatenate(( numpy.eye(matrix_rank(range_cont) ), numpy.zeros( (matrix_rank( range_cont) , matrix_rank( null_obs ) ) ) ) , 1)
	interBasis = numpy.dot( numpy.dot( range_cont, inter ), nullspace( numpy.concatenate( (range_cont,null_obs) , 1 ) ) )

	dim_int = min(shape(interBasis))
	dim_C = min( shape( range_cont ) )
	dim_N = min( shape( null_obs ) )

	T1 = range_cont[:,:(dim_C - dim_int)]
	T2 = interBasis[: , :dim_int]
	T4 = null_obs[: , dim_int: dim_N]
	T3 = nullspace( ( numpy.concatenate( (T1, T2,T4 ), 1 ) ).T)

	KBasis = numpy.concatenate( (T1,T2,T3,T4) , 1 )
	
	return KBasis, min( shape( T1 ) ) , min( shape( T2 ) ), min( shape( T3 ) ), min( shape( T4 ) ) 

def removeBasis( space, basis ):
	''' Removes the mentioned basis and provides the remaining basis from the subspace.
	returns basis which is unitary in the remaining space'''    
	to_remove = columnProj( space, basis, 1 ) 
	null = nullspace( to_remove.T )
	return numpy.dot( space, null )


def intersection( space1 , space2 ):	
	''' Finds the basis for the intersection of two spaces'''
	inter = numpy.concatenate(( numpy.eye(matrix_rank(space1) ), numpy.zeros( (matrix_rank( space1) , matrix_rank( space2 ) ) ) ) , 1)
	interBasis = numpy.dot( numpy.dot( space1, inter ), nullspace( numpy.concatenate( (space1,space2) , 1 ) ) )
	return interBasis

def KalmanBasisNew( check_A, check_B, check_C ):
	''' This produces the Kalman Basis needed for a kalman decomposition of the system, which is used for a minimal realization of the state space
	+---------------------+---------------+-------------------+-----------+--------------+
	|      Conditions     |       T1      |         T2        |     T3    |      T4      |
	+---------------------+---------------+-------------------+-----------+--------------+
	|        R, !N        |     R( C )    |         *         |  ~R( C )  |      *       |
	|  ( R N ) full rank  |     R( C )    |         *         |     *     |    N( O )    |
	|   ( R N ) low rank  |     R( C )    |         *         |  ~( R N ) |    N( O )    |
	|       Kalman        |  R( C ) - T2  | R( C ) int N( O ) |  ~( R N ) | N( O ) - T2  |
	+---------------------+---------------+-------------------+-----------+--------------+	'''

	check_cont , bool_cont = controlability( check_B , check_A )
	check_obs, bool_obs = observability( check_C, check_A )

	range_cont = columspace( check_cont ) 
	null_obs = nullspace( check_obs )


	# The system has R( Cont ) but has not N( Obs ) => Fully observable but not controllable. Similarity Transform into Cont basis can finish the job.
	
	if not bool_obs :
		interBasis = intersection( range_cont , null_obs )
	else :
		print "Completely Observable system + Uncontrollable system"
		new_A , new_B, new_C, KBasis = ContSpace( check_A, check_B, check_C , False )
		return KBasis, matrix_rank( check_cont ) , 0,0, matrix_rank( check_obs ) - matrix_rank( check_cont )
	
	
	dim_int = min(shape(interBasis))
	dim_C = min( shape( range_cont ) )
	dim_N = min( shape( null_obs ) )
		
	T2 = interBasis[: , 0:dim_int]

	if min( shape( T2 ) ) == 0:
		T1 = range_cont
		T4 = null_obs
		if matrix_rank( nC( (T1,T4), 1 ) ) == max( shape( check_A ) ):
			return nC( ( T1, T4 ) ,1 ), min( shape( T1 ) ), 0, 0, min(shape(T4) ) 
		else :
			T3 = nullspace ( ( nC( ( T1,T4) , 1 ) ).T )
			return nC( (T1,T3,T4) , 1 ) , min( shape( T1 ) ), 0 ,min( shape( T3 ) ),min( shape( T4 ) )
	
	else :	

	# The system has R( Cont ) and N( Obs ) => Usual Kalman Decomposition can take place.
		T1 = removeBasis( range_cont, T2 )
		T4 = removeBasis( null_obs , T2 )    
		T3 = nullspace( ( numpy.concatenate( (T1, T2,T4 ), 1 ) ).T)



		KBasis = numpy.concatenate( (T1,T2,T3,T4) , 1 )

		return KBasis, min( shape( T1 ) ) , min( shape( T2 ) ), min( shape( T3 ) ), min( shape( T4 ) ) 


def KBasisContObs( check_A, check_B, check_C ):
	''' This produces the Kalman Basis needed for a kalman decomposition of the system, which is used for a minimal realization of the state space Using a different algorithm than the one previously mentioned. This algorithm checks for controllability and then moves over to observability '''
	


	check_cont , bool_cont = controlability( check_B , check_A )
	check_obs, bool_obs = observability( check_C, check_A )

	range_cont = columspace( check_cont ) 
	null_obs = nullspace( check_obs )

	if bool_obs and bool_cont:
		print " Minimal realization already"
		return numpy.eye( max( shape( check_A ) ) ), max( shape( check_A ) ) , max( shape( check_A ) ),max( shape( check_A ) ),max( shape( check_A ) )

	else:
		Cont_A, Cont_B , Cont_C , Cont_col = ContSpace( check_A, check_B, check_C,False)
	    
	#print Cont_A, Cont_B , Cont_C

	# Stage I : Checking for Controllablity
	#-------------------------------------

	C = matrix_rank( range_cont )

	ACont,BCont, CCont = Cont_A[ :C , :C] , Cont_B[ :C , : ],Cont_C[ : , :C ]

	AUncont, BUncont, CUncont = Cont_A[ C: , C: ] ,Cont_B[ C: , : ], Cont_C[ :, C: ]

#	print AUncont, BUncont, CUncont

	# Stage II : checking for Observability
	#-------------------------------------

		# Stage II ( a ) : Checking for Observablity of Controllable states
		#--------------------------------------------------------------------

	check_obs , bool_obs = observability ( CCont, ACont )

	ObsACont, ObsBCont, ObsCCont , ObsCont_row = ObsSpace( ACont, BCont , CCont ,False ) # Observable from Controllable

	#print ObsACont, ObsBCont, ObsCCont
	
	O = matrix_rank( check_obs )

	AObsCont,BObsCont, CObsCont = ObsACont[ :O , :O],ObsBCont[ :O , : ] , ObsCCont[ : , :O ]

		# Stage II ( b ) : Checking for Observability of Uncontrollable states 
		#---------------------------------------------------------------------

	check_obs , bool_obs = observability ( CUncont, AUncont ) 

	ObsAUncont, ObsBUncont, ObsCUncont , ObsUncont_row = ObsSpace( AUncont, BUncont , CUncont ,False ) # Observable from Uncontrollable

#	print ObsAUncont, ObsBUncont, ObsCUncont

	O = matrix_rank( check_obs )

	AObsUncont, BOBsUncont, CObsUncont = ObsAUncont[ :O , :O] ,ObsBUncont[ :O , : ] , ObsCUncont[ : , :O ]

	# Stage III : Putting it all together, finding the basis of transformation
	#---------------------------------------------------------------------------

	Z1 = numpy.matrix( numpy.zeros( ( shape( nI(ObsCont_row) ) [ 0 ] , shape( nI(ObsUncont_row) ) [ 1 ]) ) )
	Z2 = numpy.matrix( numpy.zeros( ( shape( nI(ObsUncont_row) ) [ 0 ] , shape( nI(ObsCont_row) ) [ 1 ]) ) )

	PO = nC( ( nC( (nI(ObsCont_row),Z1) ,1 ) , nC( (Z2,nI(ObsUncont_row)) ,1 ) ),0 )
	Kbasis = nD( PO , nI(Cont_col) )
	# Truncate output : AObsCont, BObsCont, CObsCont
	return Kbasis.I , max( shape( AObsCont ) ) , max( shape(BObsCont) ), max( shape( CObsCont ) )


def BasisTrans( A,B, C, Basis ):
	''' Produces the state space after a similarity transformation'''
	new_A = numpy.dot( numpy.dot(numpy.linalg.inv( Basis ), A ), Basis )
	new_B = numpy.dot( numpy.linalg.inv(Basis) , B )
	new_C = numpy.dot( C , Basis )
	return new_A , new_B , new_C 

def minReal( check_A , check_B, check_C , method = 'k'):
	''' Produces the minimal realization of the state space using Kalman decomposition procedure
	Returns only the space which is both controllable and observable'''

	if method == 'k':
		KBasis , CO , C_O , _CO , _C_O = KalmanBasisNew( check_A , check_B, check_C )
	else :
		KBasis , CO , C_O , _CO  = KBasisContObs( check_A , check_B, check_C )
#        print matrix_rank( KBasis ) 
	new_A, new_B , new_C = BasisTrans( check_A , check_B, check_C , KBasis )
	minR_A = new_A[:CO, :CO]
	minR_B = new_B[ :CO, : ]
	minR_C = new_C[ :, :CO ]
	return minR_A, minR_B, minR_C

def ObsSpace( A, B, C,conv=True ):
	'''Does a basis transformation into the space where Observable and Unobservable states are orthogonal '''
	obs, obs_bool = observability( C ,A )
	col, eigen, row = numpy.linalg.svd( obs )
	new_A , new_B, new_C = BasisTrans( A, B,C, row.T )
	if conv :
		return new_A , new_B, new_C
	else :
		return new_A , new_B, new_C,row.T

def ContSpace( A , B ,C ,conv=True):
	'''Does a basis transformation into the space where Controllable and Uncontrollable states are orthogonal '''
	cont, cont_bool = controlability( B , A )
	col , eigen , row = numpy.linalg.svd( cont )
	new_A , new_B , new_C = BasisTrans(A , B , C , col )
	if conv :
		return new_A, new_B, new_C
	else :
		return new_A, new_B, new_C, col
