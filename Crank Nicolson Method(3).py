import numpy as np
import matplotlib.pyplot as plt
import math
import numpy.linalg


Nt=20; # Try also Nt=80
Nx=39;
T0=0; Tf=1; # T0<=t<=Tf
a=-1; b=1; # a<=x<=b
tau=Tf/Nt; # Step for t
h=(b-a)/(Nx); # Step for x
m=tau/(h**2); # Courant

t=np.linspace(T0,Tf,Nt); # t span
x=np.linspace(a,b,Nx); # x span
U=np.zeros((Nt,Nx)); # Iinitialize U to save approximations

# Initial conditions
U[0,:]=np.sin((math.pi)*x);

# Boundary conditions
U[:,0]=0;
U[:,Nx-1]=0;

# Initialize matrices for a linear system D*X=B and the coefficients alpha
D=np.zeros((Nx-2,Nx-2));
B=np.zeros((Nx-2));
alpha=np.zeros((Nx-2));

# Define all coefficients alpha
for i in range(0,Nx-2):
    alpha[i]=1+(x[i+1]**2);

gama=-1 # Gama is a constant in the method

# Insert values in the tridiagonal matrix D, according to Crank-Nickolson Method
for j in range(0,Nx-3):
    D[j+1,j]=-m*alpha[j+1];
    D[j,j+1]=-m*alpha[j];
for i in range(0,Nx-2):    
    D[i,i]=2+2*m*alpha[i]+tau*gama;
    
# Save in matrix C all the elements of matrix U, except from first and last columns
C=np.zeros((Nt,Nx-2));
for n in range(0,Nt):
    for i in range(0,Nx-2):
        C[n,i]=U[n,i+1];
        
TC=np.transpose(C); # Transpose matrix C

# Solve the linear system D*TC=B for every step n
for n in range(0,Nt-1):
    for i in range(1,Nx-3):
        B[i]=(2-2*m*alpha[i]-tau*gama)*TC[i,n]+(m*alpha[i])*TC[i+1,n]+m*alpha[i]*TC[i-1,n];
    B[0]=(2-2*m*alpha[0]-tau*gama)*TC[0,n]+(m*alpha[0])*TC[1,n];
    B[Nx-3]=(2-2*m*alpha[Nx-3]-tau*gama)*TC[Nx-3,n]+m*alpha[Nx-3]*TC[Nx-4,n];
    TC[:,n+1]=numpy.linalg.solve(D,B)

# Save solutions in matrix U    
for n in range(0,Nt):
    for i in range(0,Nx-2):
        U[n,i+1]=C[n,i];

# Generate Mesh    
T,X=np.meshgrid(t,x);

# for i in range(0,19):
#     plt.plot(X,U[i,:]);
# plt.show()

# Approximation at n=5
plt.plot(X,U[5,:]);
plt.show()
# Approximation at n=10
plt.plot(X,U[10,:]);
plt.show()
# Approximation at n=19
plt.plot(X,U[19,:]);
plt.show()
