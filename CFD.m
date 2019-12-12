function [x,theta]=CFD(func,n)
f=@(eta) ppval(func,eta);
x=linspace(0,10,n);
h=10/(n-1);
A=zeros(n-2,n-2);
for i=2:n-2
    c=0.7*f(x(i))*h/2;
    A(i-1,i-1)=-2;
    A(i-1,i)=1+c;
    A(i,i-1)=1-c;
    A(i,i)=-2;
end
c=0.7*f(x(1))*h/2;
z=zeros(n-2,1);
z(1)=(-1+c)*1;
theta=A\z;
theta=[1;theta;0];
