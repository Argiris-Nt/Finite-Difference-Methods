import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg


Nt=100; # Try also Nt=150
Nx=100; # Try also Nx=115
T0=0; Tf=1; # T0<=t<=Tf
a=0; b=1; # a<=x<=b
tau=Tf/Nt; # Step for t
h=(b-a)/(Nx); # Step for x
ni=tau/h; # Courant

t=np.linspace(T0,Tf,Nt); # t span
x=np.linspace(a,b,Nx); # x span
U=np.zeros((Nt,Nx)); # Iinitialize U to save approximations

# Initial conditions
U[0,:]=np.exp(-10*(4*x-1)**2);

# L-W Method
for n in range(0,Nt-1):
    for i in range(1,Nx-1):
        U[n+1,i]=U[n,i]-(ni/4)*(U[n,i+1]**2-U[n,i-1]**2)+((ni**2)/8)*(((U[n,i+1]+U[n,i])**2)*(U[n,i+1]-U[n,i])-((U[n,i]+U[n,i-1])**2)*(U[n,i]-U[n,i-1]));

# Generate Mesh          
T,X=np.meshgrid(t,x);

# for i in range(0,70):
#     plt.plot(X,U[i,:]);
# plt.show()

# Approximation at n=0
plt.plot(X,U[0,:]);
plt.show()
# Approximation at n=10
plt.plot(X,U[10,:]);
plt.show()
# Approximation at n=20
plt.plot(X,U[20,:]);
plt.show()
# Approximation at n=30
plt.plot(X,U[30,:]);
plt.show()
# Approximation at n=40
plt.plot(X,U[40,:]);
plt.show()
# Approximation at n=60
plt.plot(X,U[60,:]);
plt.show()
