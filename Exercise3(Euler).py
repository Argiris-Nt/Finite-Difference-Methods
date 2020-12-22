import numpy as np
import matplotlib.pyplot as plt
import math
import numpy.linalg


Nx=39; Nt=400;
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

# Implicit Euler Method
# U[n+1,i]=(1/(1+2*m*(1+x[i]**2)))*(m*(1+x[i]**2)*(U[n+1,i+1]+U[n+1,i-1])+(1+tau)*U[n,i]);

# Initialize matrices for a linear system D*X=B
D=np.zeros((Nx-2,Nx-2));
B=np.zeros((Nx-2));

# Insert values in the tridiagonal matrix D, according to Implicit Euler Method
for j in range(0,Nx-3):
    D[j+1,j]=-m*(1+x[j+1]**2);
    D[j,j+1]=-m*(1+x[j]**2);
for i in range(0,Nx-2):    
    D[i,i]=1+2*m*(1+x[i]**2);

# Save in matrix C all the elements of matrix U, except from first and last columns
C=np.zeros((Nt,Nx-2));
for n in range(0,Nt):
    for i in range(0,Nx-2):
        C[n,i]=U[n,i+1];
        
TC=np.transpose(C); # Transpose matrix C

# Solve the linear system D*TC=B for every step n
for n in range(0,Nt-1): 
    for i in range(0,Nx-2):
        B[i]=(1+tau)*TC[i,n];  
    TC[:,n+1]=numpy.linalg.solve(D,B)
 
# Save solutions in matrix U 
for n in range(0,Nt):
    for i in range(0,Nx-2):
        U[n,i+1]=C[n,i];
        
# Generate Mesh        
T,X=np.meshgrid(t,x);

# for i in range(0,50):
#     plt.plot(X,U[i,:]);
# plt.show()

# Approximation at n=100
plt.plot(X,U[100,:]);
plt.show()
# Approximation at n=200
plt.plot(X,U[200,:]);
plt.show()
# Approximation at n=399
plt.plot(X,U[399,:]);
plt.show()
