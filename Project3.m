%Fluids Project 3

close all 
clearvars
clc

%% Settings
doPlot=false;

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
dt=.25;
tf=1;

t=t0:dt:tf;
Nt=numel(t);

%space params
y0=0;
dy=1;
yf=10;

y=y0:dy:yf;
Ny=numel(y);

%% Explicit Solution

%Assume Solution
ue=ones(Ny,Nt);


% Enforce Boundary Conditions
ue(1:Ny,1)=0; % initial condition t=0
ue(1,1:Nt)=U; % Boundary condition
ue(Ny,1:Nt)=0; % Far Field

Gamma=nu*dt/dy^2;
disp(Gamma);
for n=1:Nt-1
    for j=2:Ny-1
        ue(j,n+1)= ue(j,n) + Gamma*(ue(j+1,n)-2*ue(j,n)+ue(j-1,n));
    end
end

if doPlot
    figure(1)
    plot(ue);
end

%% Implicit Solution Params
% if nothing is here they are the same as the Explicit parameters

%% Implicit Solution

%Assume Solution
ui=ones(Ny,Nt);

% Enforce Boundary Conditions
ui(1:Ny,1)=0; % initial condition t=0
ui(1,1:Nt)=U; % Boundary condition
ui(Ny,1:Nt)=0; % Far Field


A_mat=(1+2*Gamma)*eye(Ny-2);

for i=1:Ny-3
    A_mat(i+1,i)=-Gamma;
    A_mat(i,i+1)=-Gamma;
end

for n=1:Nt-1
   b=ui(2:end-1,n);
   b(1)=b(1)+Gamma*ui(1,n);
   b(end)=b(end)+Gamma*ui(end,n);
   x=tri_diag_solver(A_mat,b);
   ui(2:Ny-1,n+1)= A_mat\b; % TODO tri-diagonal solver
end

if doPlot
    figure(2)
    plot(ui);
end

%% Error Analysis