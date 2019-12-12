function [x,theta]=CFD(func,n)
f=@(eta) ppval(func,eta);
x=linspace(0,10,n);
h=10/(n-1);
A=zeros(n-2,n-2);
Prfx=0.7*f(x)*h/2;
for i=1:n-2-1
    A(i,i)=-2;
    A(i,i+1)=1+Prfx(i+1);
    A(i+1,i)=1-Prfx(i+2);%for next row, use next f
end
A(end,end)=-2;
Prfx=0.7*f(x(2))*h/2;
z=zeros(n-2,1);
z(1)=-(1-Prfx)*1;
%z(end)=-(1+c)*0;
theta=A\z;
theta=[1;theta;0];