import numpy as np
import matplotlib.pyplot as plt


Nt=250; # Try also Nt=300
Nx=100; # Try also Nx=150
T0=0; Tf=1; # T0<=t<=Tf
a=0; b=1; # a<=x<=b
tau=Tf/Nt; # Step for t
h=(b-a)/(Nx); # Step for x
ni=tau/h; # Courant
alpha=2; # Wind

t=np.linspace(T0,Tf,Nt); # t span
x=np.linspace(a,b,Nx); # x span
U=np.zeros((Nt,Nx)); # Iinitialize U to save approximations

# Initial conditions
U[0,:]=np.where((x>0.1)&(x<0.2),(10**4)*((0.1-x)**2)*(0.2-x)**2,0);

# Boundary conditions
U[:,0]=0;

# L-W Method
for n in range(0,Nt-1):
    for i in range(1,Nx-1):
        U[n+1,i]=((alpha*ni)/2)*(1+alpha*ni)*U[n,i-1]+(1-(alpha**2)*(ni**2))*U[n,i]-((alpha*ni)/2)*(1-alpha*ni)*U[n,i+1];

# Generate Mesh         
T,X=np.meshgrid(t,x);

# Approximation at n=0
plt.plot(X,U[0,:])
# Exact solution
plt.plot(x+2*T[:,0],U[0,:]);
plt.show()

# Approximation at n=35
plt.plot(X,U[35,:]);
# Exact solution
plt.plot(x+2*T[:,35],U[0,:]);
plt.show()

# Approximation at n=50
plt.plot(X,U[50,:]);
# Exact solution
plt.plot(x+2*T[:,50],U[0,:]);
plt.show()

# Approximation at n=80
plt.plot(X,U[80,:]);
# Exact solution
plt.plot(x+2*T[:,80],U[0,:]);
plt.show()
