# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import sympy
import numpy
from scipy.integrate import odeint
from scipy import exp
import matplotlib.pyplot as plt

#Function to convert symbolic expression with numerical data to numpy array 
def Sym2NumArray(F):
    shapeF=F.shape
    B=numpy.zeros(shapeF)
    for i in range(0,shapeF[0]):
        for j in range(0,shapeF[1]):
            B[i,j]=sympy.N(F[i,j])
    return B

#Function for evaluating input signal
def signal_in(t):
    if t>0.2:
        out=1.0
    else:
        out=0.0
    return out

def signal_new(t):
    if t>0.2:
        out=numpy.sin(2*t)
    else:
        out=0.0
    return out


# <headingcell level=4>

# Calculates all the non linear Matrices needed to define the system

# <codecell>

ord=10
A=numpy.diag(numpy.ones(ord)*-2)+numpy.diag(numpy.ones(ord-1),1)+numpy.diag(numpy.ones(ord-1),-1)
B=numpy.zeros((ord,1));C=numpy.zeros((ord));
B[0]=1;C[ord-1]=1;D=numpy.zeros((1));
A[ord-1][ord-1]=-1

x=sympy.symbols("x1:11");
nonl=sympy.Matrix(numpy.zeros((ord,1)))

for i in range(1,ord-1):
    nonl[i]=sympy.exp(10*x[i-1]-10*x[i])-sympy.exp(10*x[i]-10*x[i+1])
nonl[0]=2-sympy.exp(10*x[0])-sympy.exp(10*x[0]-10*x[1])
nonl[ord-1]=sympy.exp(10*x[ord-2]-10*x[ord-1])-1

xMat=sympy.Matrix(x)
Anew=A*(xMat)+nonl
AJac=Anew.jacobian(xMat)


# <headingcell level=4>

# Solves the system for a 1st Order approximation of the Actual non Linear sytem at one point, here the Origin

# <codecell>

AJacOrig=AJac.evalf(subs={x[0]:0.,x[1]:0.,x[2]:0.,x[3]:0.,x[4]:0.,x[5]:0.,x[6]:0.,x[7]:0.,x[8]:0.,x[9]:0.})
AJacON=Sym2NumArray(AJacOrig)

def dervL1(y,t):
    temp=numpy.dot(AJacON,y)
    sigIn=numpy.dot(B,(signal_new(t)))
    sigIn_R=sigIn.reshape(1,10)
    tempL1=temp+sigIn_R[0]
    return tempL1

y0=numpy.zeros(ord)
time  = numpy.linspace(0, 30., 100)
solnL1 = odeint(dervL1, y0, time)

# <headingcell level=4>

# Plot of the system linearized at one point

# <codecell>

plt.figure(0)
plt.plot(time,solnL1[:,9])
plt.xlabel('Time')
plt.ylabel('Voltage')
plt.title('Solution obtained from model Linearized at Origin')


# <headingcell level=4>

# Solving the Actual Non Linear System for Training input

# <codecell>

def dervNL(y,t):
    tempObj=Anew.evalf(subs={x[0]:y[0],x[1]:y[1],x[2]:y[2],x[3]:y[3],x[4]:y[4],x[5]:y[5],x[6]:y[6],x[7]:y[7],x[8]:y[8],x[9]:y[9]})+B*signal_in(t)
    temp=Sym2NumArray(tempObj).reshape(1,10)
    return temp[0]
    
y0=numpy.zeros(ord)
timeNL = numpy.linspace(0, 30., 1000)
solnNL = odeint(dervNL, y0, timeNL)   

# <headingcell level=4>

# Plot of the solution for Non Linear System for the Training input

# <codecell>

plt.figure(1)
plt.plot(timeNL,solnNL[:,9])
plt.xlabel('Time')
plt.ylabel('Voltage')
plt.title('Solution from Training input Non Linear Model')



# <headingcell level=4>

# Solving for the non linear system for the actual input

# <codecell>

def dervNLActual(y,t):
    tempObj=Anew.evalf(subs={x[0]:y[0],x[1]:y[1],x[2]:y[2],x[3]:y[3],x[4]:y[4],x[5]:y[5],x[6]:y[6],x[7]:y[7],x[8]:y[8],x[9]:y[9]})+B*signal_new(t)
    temp=Sym2NumArray(tempObj).reshape(1,10)
    return temp[0]
    
y0=numpy.zeros(ord)
timeNL = numpy.linspace(0, 30., 1000)
solnNLActual = odeint(dervNLActual, y0, timeNL)   

plt.figure(2)
plt.plot(timeNL,solnNLActual[:,9])
plt.xlabel('Time')
plt.ylabel('Voltage')
plt.title('Solution from actual Non-Linear Model')


# <codecell>

plt.figure(3)
Train=plt.plot(timeNL,solnNL[:,9])
Actual=plt.plot(timeNL,solnNLActual[:,9])
plt.figlegend([Train,Actual],('Training input Trajectory','Actual input response'),'lower right')
plt.xlabel('Time')
plt.ylabel('Voltage')
plt.title('Comparision between the Training and Actual input')


# <headingcell level=4>

# Selection of Linearization Points

# <codecell>

points=[10,50,130,170,700]
linP=numpy.zeros((ord,len(points)))
linPJac=numpy.zeros((5,10,10))
offset=numpy.zeros((5))
for i in range(0,5):
    linP[:,i]=solnNL[points[i],:]

# <headingcell level=4>

# Evaluating Jacobian and other parameters needed

# <codecell>

for i in range(0,5):
    linPV=Anew.evalf(subs={x[0]:linP[0,i],x[1]:linP[1,i],x[2]:linP[2,i],x[3]:linP[3,i],x[4]:linP[4,i],x[5]:linP[5,i],x[6]:linP[6,i],x[7]:linP[7,i],x[8]:linP[8,i],x[9]:linP[9,i]})
    linPV_N=Sym2NumArray(linPV)
    linPJ=AJac.evalf(subs={x[0]:linP[0,i],x[1]:linP[1,i],x[2]:linP[2,i],x[3]:linP[3,i],x[4]:linP[4,i],x[5]:linP[5,i],x[6]:linP[6,i],x[7]:linP[7,i],x[8]:linP[8,i],x[9]:linP[9,i]})
    linPJ_N=Sym2NumArray(linPJ)
    linPJac[i]=linPJ_N
    temp=linPV_N.reshape((1,10))-numpy.dot(linPJ_N,linP[:,i])
    offset[i]=temp[0][0]

# <codecell>

def normcalc(y):
    normV=numpy.zeros(5)
    for i in range(0,5):
        temp=y-linP[:,i]
        normV[i]=numpy.linalg.norm(temp)
    temp=10**(-20*normV*normV/normV.min())       
    return temp/temp.sum()

# <headingcell level=4>

# Solving the PWL Model

# <codecell>

def dervPWL(y,t):
    weight=normcalc(y)
    temp=numpy.dot(linPJac[0],y)*weight[0]+numpy.dot(linPJac[1],y)*weight[1]+numpy.dot(linPJac[2],y)*weight[2]+numpy.dot(linPJac[3],y)*weight[3]+numpy.dot(linPJac[4],y)*weight[4]
    sigIn=numpy.dot(B,(signal_new(t)+numpy.dot(offset,weight)))
    sigIn_R=sigIn.reshape(1,10)
    tempPWL=temp+sigIn_R[0]
    return tempPWL

y0=numpy.zeros(ord)
time  = numpy.linspace(0, 30., 100)
solnPWL = odeint(dervPWL, y0, time)

# <headingcell level=4>

# Plot for the solution of the PWL Model

# <codecell>

plt.figure(4)
plt.plot(time,solnPWL[:,9])
plt.xlabel('Time')
plt.ylabel('Voltage')
plt.title('Solution obtained from full order PWL Model')


# <headingcell level=4>

# Comparing the solution obtained from Non Linear, Linear at one point and PWL Models 

# <codecell>

plt.figure(5)
PWL=plt.plot(time,solnPWL[:,9])
NL=plt.plot(timeNL,solnNL[:,9])
L1=plt.plot(time,solnL1[:,9])
NLActual=plt.plot(timeNL,solnNLActual[:,9])
plt.title('Comparision b/w solutions')
plt.figlegend([PWL,NL,L1,NLActual],('Piece Wise Linear','Training input Trajectory','Linear at Origin','Actual solution'),'lower right')
plt.xlabel('Time')
plt.ylabel('Voltage')


# <headingcell level=4>

# Arnoldi Basis evaluation

# <codecell>

def ArnoldiBasis_QR(Jac,B,ord,red_ord):
    spanB=numpy.zeros((ord,red_ord))
    A=numpy.linalg.inv(Jac)
    temp=numpy.dot(A,B)
    spanB[:,0]=temp.reshape(1,ord)
    for i in range(1,red_ord):
        spanB[:,i]=numpy.dot(A,spanB[:,i-1])
    unitB,R=numpy.linalg.qr(spanB)
    return unitB
    

# <headingcell level=4>

# Creating Orthonormal basis at each linear point and creating the Aggregate space

# <codecell>

red_ord=3
AggB=numpy.copy(linP)
for i in range(0,5):
    unitB=ArnoldiBasis_QR(linPJac[i],B,ord,red_ord)
    unitBoff=ArnoldiBasis_QR(linPJac[i],numpy.dot(B,offset[i]),ord,red_ord)
    AggB=numpy.concatenate((unitB,unitBoff,AggB),1)

# <headingcell level=4>

# SVD to find the Orthonormal vectors and significant eigen values

# <codecell>

eigenC,eigenV,eigenR=numpy.linalg.svd(AggB)
for i in range(1,ord):
    if  (eigenV[0]/eigenV[i])>100:
        red_ordNew=i-1
        break

# <codecell>

redBasis=eigenC[:,:red_ordNew]

# <codecell>

newReptemp=numpy.dot(numpy.linalg.inv(eigenC),linP)
newRep=newReptemp[:red_ordNew,:]

# <codecell>

def normcalcMOR(y):
    normV=numpy.zeros(5)
    for i in range(0,5):
        temp=y-newRep[:,i]
        normV[i]=numpy.linalg.norm(temp)
    temp=10**(-20*normV*normV/normV.min())
    return temp/(temp.sum())
    

# <codecell>


linPJacMOR=numpy.zeros((5,red_ordNew,red_ordNew))
for i in range(5):
    linPJacMOR[i]=numpy.dot(redBasis.T,numpy.dot(linPJac[i],redBasis))
    

# <headingcell level=4>

# Solving the PWL MOR system

# <codecell>

def dervPWLMOR(y,t):
    weight=normcalcMOR(y)
    temp=numpy.dot(linPJacMOR[0],y)*weight[0]+numpy.dot(linPJacMOR[1],y)*weight[1]+numpy.dot(linPJacMOR[2],y)*weight[2]+numpy.dot(linPJacMOR[3],y)*weight[3]+numpy.dot(linPJacMOR[4],y)*weight[4]
    sigIn=numpy.dot(redBasis.T,numpy.dot(B,(signal_new(t)+numpy.dot(offset,weight))))
    sigIn_R=sigIn.reshape(1,red_ordNew)
    tempPWLMOR=temp+sigIn_R[0]
    return tempPWLMOR

y0=numpy.zeros(red_ordNew)
time  = numpy.linspace(0, 30., 100)
solnMOR = odeint(dervPWLMOR, y0, time)
solnMOR_out=numpy.zeros(len(time))
for i in range(len(time)):
    solnMOR_out[i]=numpy.dot(redBasis[9,:],solnMOR[i,:])

# <headingcell level=4>

# Plotting the solution for PWL+MOR system 

# <codecell>

plt.figure(6)
plt.plot(time,solnMOR_out)
plt.title('Solution from PWL+MOR')
plt.xlabel('Time')
plt.ylabel('Voltage')


# <headingcell level=4>

# Comparing the solutions from various methods used

# <codecell>

plt.figure(7)
PWL=plt.plot(time,solnPWL[:,9])
NL=plt.plot(timeNL,solnNL[:,9])
PWLMOR=plt.plot(time,solnMOR_out)
L1=plt.plot(time,solnL1[:,9])
Actual=plt.plot(timeNL,solnNLActual[:,9])
plt.title('Comparision b/w solutions')
plt.figlegend([PWL,NL,PWLMOR,L1,Actual],('Piece Wise Linear','Training Input','PWLMOR','Linear at Origin','Actual NonLinear'),'lower right')
plt.xlabel('Time')
plt.ylabel('Voltage')
plt.show()

# <codecell>


# <codecell>


