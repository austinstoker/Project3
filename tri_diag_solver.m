function [x] = tri_diag_solver(A,b)
% Works by turning Ax=b into LUx=b Ly=b (y=Ux)

N=size(A,1);
x=zeros(N,1);
temp=zeros(N,1);

[L1,U1] = lu(A);



%get L and U
b=A(1,1);
a=A(2,1);
c=A(1,2);

d(1) = b;
u(1) = c;

for i=1:N-2
    l(i)=a/d(i);
    u(i)=c;
    d(i+1)=b-u(i)*l(i);
end
l(N-1)=a/d(N-2);
u(N-1)=c;
d(N-1)=b-u(N-1)*l(N-1);


beta = A(1);
x(1) = b(1)/beta;
% forward subs
for i = 2:N
    temp(i) = U(i-1)/beta;
    beta = A(i)-L(i)*temp(i);
 x(i) = (A(i)-L(i)*x(i-1))/beta;
end
% back subs
for j = (N-1):-1:1
    %k = N-j;
    x(j) = x(j)-gamma(j+1)*x(j+1);
end
end