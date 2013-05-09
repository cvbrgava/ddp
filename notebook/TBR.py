
from scipy.linalg import *

for linpt in [1]:
        t = time[ linpt ] - sim_begin
	input_signals = get_input_signals( t ) 
	sub1 = { state[ str( input_list[ k ] ) ] : input_signals[ k ] for k in range( len( input_list ) ) }
# Develop all the state matrices

	val = Sym2NumArray( (Matrix[ linpt ].linPVal).evalf( subs = sub1 ))
	jac = Sym2NumArray( ( Matrix[ linpt ].linPJac).evalf( subs = sub1 ) )	
        

	offset = val  - (numpy.dot( jac, linP[ linpt ] )).reshape( order ,1 )

	Bnew = numpy.concatenate( (offset , B[ linpt ]) , 1 )

# Minimal realization does a kalman decomposition of the state matrices
	minR_A, minR_B, minR_C = minReal( jac , Bnew, C )

#Calculate the controllability gramian
	contGram = solve_lyapunov( minR_A,-numpy.dot( minR_B, transpose( minR_B ) ) )

# Calculate the Observability gramian

	obsGram = solve_lyapunov( transpose( minR_A ) , -numpy.dot( transpose( minR_C ), minR_C  ) )

#        print "Positive def of OBS GRAM",(numpy.array( [ i for i in numpy.linalg.eigvals(obsGram ) ])).min()

# Cholesky Decomposition on ContGram and Obs Gram

	LCGram = numpy.linalg.cholesky( contGram )
	LOGram = numpy.linalg.cholesky ( obsGram )
	U, eigenV , V = svd( numpy.dot( transpose( LOGram ) , LCGram ) )

        CO_order = max( shape( minR_A ) )
        #if CO_order < order :
           # x.add_row( [linpt, True ] )
        #else:
            #x.add_row( [linpt, False ] )
    
        truncOrder = 1;
	for truncOrder in range( len( eigenV) ):
		if eigenV[ 0 ] / eigenV[ truncOrder ] > 10 :
			truncOrder = truncOrder+1
			break 

# Truncation depending upon the Hankel Singular Values

	Sigma = numpy.diag( numpy.array ( [ eigenV[ i ]**-0.5 for i in range( len( eigenV ) )  if i < truncOrder  ] ) )
	print  truncOrder
	I1 = numpy.eye(truncOrder)
	I2 = numpy.zeros( (truncOrder, CO_order - truncOrder) ) 
	Trunc = transpose( nC( (I1,I2 ) , 1) )

# Basis of projection 

	Vbasis = nD( LCGram, nD( V , nD( Trunc,  Sigma )) )
	Wbasis = nD( LOGram , nD( U ,nD(Trunc, Sigma ) )) 
	Agg  = nC( (Vbasis, Wbasis) , 1 )
		
	Col,eigenV, R = svd( Agg )

