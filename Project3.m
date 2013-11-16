%Fluids Project 3

%% Problem Setup

%solve for u(y,t)

%     du       d^2 u
%     --  = nu -----
%     dt        dy^2

nu=1; %viscosity
U=10; %Plate velocity


%% Explicit Solution Params

%time params
t0=0;
dt=1;
tf=100;

t=t0:dt:tf;
Nt=numel(t);

%space params
y0=0;
dy=1;
yf=100;

y=y0:dy:yf;
Ny=numel(y);

%% Explicit Solution

%Assume Solution
u=ones(Ny,Nt);


% Enforce Boundary Conditions
u(1:Ny,1)=0; % initial condition t=0
u(1,1:Nt)=U; % Boundary condition
u(Ny,1:Nt)=0; % Far Field

Gamma=nu*dt/dy^2;
for n=1:Nt-1
    for j=2:Ny-1
        u(j,n+1)= u(j,n) + Gamma*(u(j+1,n)-2*u(j,n)+u(j-1,n));
    end
end



%% Implicit Solution Params

%% Implicit Solution


%% Error Analysis