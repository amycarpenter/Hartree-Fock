########################################################################
# The Helium Atom
########################################################################

from numpy import empty,ones,pi,exp,sqrt,dot,linspace
from numpy.linalg import inv,eig
import matplotlib.pyplot as plt

########################################################################
#  Constants
########################################################################
a = (0.298073, 1.242567, 5.782948, 38.474970) # alpha's for Gaussian basis functions
M = len(a)

########################################################################
#  Create the matricies S, h, and Q
########################################################################

### Sij
S = empty([M,M],float)
for i in range(M):
	for j in range(i,M):
		S[i,j] = S[j,i] = (pi/(a[i]+a[j]))**(1.5)
Sin = inv(S)

### hij
h = empty([M,M],float)

def funh(i,j):
	return 3*a[i]*a[j]*(pi**1.5)/(a[i]+a[j])**2.5 - (4*pi/(a[i]+a[j]))

for i in range(M):
	for j in range(i,M):
		h[i,j]=h[j,i]=funh(i,j)
	
### Qijkl
Q = empty([M,M,M,M],float)

def funQ(i,j,k,l):
	return 2*(pi**2.5)/((a[i]+a[k])*(a[j]+a[l])*(a[i]+a[j]+a[k]+a[l])**0.5)	
			
for i in range(M):
	for j in range(M):
		for k in range(M):
			for l in range(M):
				Q[i,j,k,l]=funQ(i,j,k,l)

########################################################################
#  Functions
########################################################################

# Normalize the wavefunction
def normlz(C):
	b=0.0
	for i in range(len(C)):
		for j in range(len(C)):
			b+=C[i]*S[i,j]*C[j]
	C*=b**(-0.5)
	return C
	
# Finding E from C values
def Energy(C):
	E = 0.0
	for i in range(M):
		for j in range(M):
			for k in range(M):
				for l in range(M):
					E += C[i]*C[j]*C[k]*C[l]*Q[i,k,j,l]
			E += 2*C[i]*C[j]*h[i,j]
	return E
    
########################################################################
#  Functions for graphs
########################################################################
def fun(r,i):
	return exp(-a[i]*r**2)
    
# Calculate radial part squared
def radialprob(r,C):
	value=0.0
	for i in range(M):
		for j in range(M):
			value+=C[i]*C[j]*fun(r,i)*fun(r,j)
	value*=4*pi*r**2
	return value
	
def EnergyDistr(r,C):
	E = 0.0
	for i in range(M):
		for j in range(M):
			for k in range(M):
				for l in range(M):
					E += C[i]*C[j]*C[k]*C[l]*fun(r,i)*fun(r,k)*fun(r,j)*fun(r,l)
			E += 2*C[i]*C[j]*h[i,j]
	return E

########################################################################
#  Finding C using the self-consistent field method/gradient iteration method
########################################################################

# Initialize C vector (and normalize)
C = ones([M,1],float)
C = normlz(C)

epsilon = 1E-8
Epold = 0.0 # just making up starting values
Epnew = 1.0

# For graphs
N=100 # use e.g. 1000 for higher resolution
r1points = linspace(0,3,N)
iterations = 0

while abs(Epnew - Epold) > epsilon:
	iterations+=1
	
	if iterations < 4:
		plt.plot(r1points,radialprob(r1points,C),label=('Step'+str(iterations)))
		
	F = empty([M,M],float)
	value=0.0
	for i in range(M):
		for j in range(M):
			for k in range(M):
				for l in range(M):
					value += Q[i,k,j,l]*C[k]*C[l]
			F[i,j] = value + h[i,j]
			value = 0.0
	
	# Solve the eigenvalue problem S-1FC=E+C
	A = dot(Sin,F)
	eigvals,eigvects = eig(A)
	
	# find the smallest eigenvalue and corresponding eigenvector
	sol=eigvals[0]
	sol_index=0
	for i in range(M):
		if eigvals[i]<sol:
			sol=eigvals[i]
			sol_index=i
	Epold = Epnew
	Epnew = eigvals[sol_index]
	C = eigvects[:,sol_index]
	
	# Normalize new C's
	C = normlz(C)
    
	# Calculate new Energy
	print("E:",Energy(C))
    
# Convergence Graph
plt.plot(r1points,radialprob(r1points,C),label=('Final'))
plt.xlabel("Distance from Nucleus")
plt.ylabel("Radial Probability Distribution")
plt.title("Convergence of Helium Wavefunction with Gradient Iteration Method")
plt.legend(loc="best")
plt.show()

# Radial Probability Density Graph
R=2
x1points = linspace(-R,R,N)
y1points = linspace(-R,R,N)
mesh = empty([N,N],float)
heatpoints = empty([N,N],float)
	
for x in range(0,N):
	for y in range(0,N):
		r = sqrt(x1points[x]**2+y1points[y]**2)
		mesh[x,y] = r
		heatpoints[x,y] = radialprob(r,C)
        
plt.title("Radial Probability Distribution")
plt.xlabel("Distance from Nucleus")
plt.ylabel("Distance from Nucleus")
plt.pcolor(x1points,y1points,heatpoints,cmap="Purples")
#plt.colorbar()
plt.show()
