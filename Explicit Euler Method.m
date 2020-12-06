Nt=800; % Try also Nt=770
Nx=19;
T0=0; Tf=1; % T0<=t<=Tf
a=0; b=1; % a<=x<=b
tau=Tf/Nt; % Step for t
h=(b-a)/(Nx+1); % Step for x
m=tau/(h^2); % Courant

x=linspace(a,b,Nx); % t span
t=linspace(T0,Tf,Nt); % x span

U=zeros(Nt,Nx); % Iinitialize U to save approximations

% Initial conditions
B=intersect(find(x>=0),find(x<=1/2));
for i=1:length(B)
    U(1,i)=2*x(B(i));
end
C=intersect(find(x>1/2),find(x<=1));
for i=1:length(C)
    U(1,length(B)+i)=2-2*x(C(i));
end 

% Boundary conditions
for k=1:Nt
    U(k,1)=0;
    U(k,Nx)=0;
end

% Explicit Euler Method
for n=1:Nt-1
    for i=2:Nx-1
        U(n+1,i)=m*U(n,i+1)+(1-2*m)*U(n,i)+m*U(n,i-1);
    end
end

% Generate mesh
[T,X] = meshgrid(t, x);

% Show some approximations
plot(X,U(1,:));
title('Nx=19 and Nt=800')
hold all
for i=2:25:Nt
    plot(X,U(i,:));
    hold all
end
