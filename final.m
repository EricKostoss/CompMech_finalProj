%% ME 3255 - Final
% * *Name:* Alden Lamp
% * *Date:* 12/__/2019

%% I
nu=1.5e-5; rho=1.2;%[m^2/s]; [kg/m^3] with laminar valued Reynolds Number
%nondimensional y: eta=y/x*(Re^.25); velocity: df/d'eta = u/Uo
%Uo/Ui = 4/(Re^.5) to ODE  d3f/d'eta^3 + f*d2f/d'eta^2 + 2*(df/d'eta)^2 = 0
%for eta=0, f=df/d'eta=0;  for eta=H=10, df/d'eta=d2f/d'eta^2 = 0

%% I.a
clc;clear;
% implicit solution given as
eta=@(f) log(sqrt(1+sqrt(f)+f)./(1-sqrt(f)))+sqrt(3)*atan(sqrt(3.*f)./(2+sqrt(f)));
f=0:0.0005:1;
etaf=eta(f);
figure(1);
plot(etaf,f);
xlim([0,10]);
xlabel('eta');
ylabel('script f')

%% I.b
% derivative
detaf=diffc2(etaf,0.0005);
df=1./detaf;
figure(2);
plot(etaf,df);
xlim([0,10]);
xlabel('eta');
ylabel('derivative script f')

%% I.c
% derivative and proving the thing
dtf=diffc2(df,0.0005);
figure(3);
plot(etaf,dtf);
xlim([0,10]);
xlabel('eta');
ylabel('second derivative script f')
disp(dtf(1));disp(1.778/8)

%% I.d
% maximum velocity
[fx,index]=goldmin_array(df);
ply=polyfit(etaf((index-5):(index+5)),df((index-5):(index+5)),2);
figure(2);
hold on;
plot(etaf((index-5):(index+5)),polyval(ply,etaf((index-5):(index+5))))

%% 
