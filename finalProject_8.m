%% Final Project
% Alden "Mac" Lamp, Eric Kostoss, Nathan Orsini; 11-25-2019
%% init
clear all
%% I
nu=1.5e-5; rho=1.2;%[m^2/s]; [kg/m^3] with laminar valued Reynolds Number
Ui=1; Re=@(x) Ui*x/nu;% [m/s]; []
%nondimensional y: eta=y/x*(Re^.25); velocity: df/d'eta = u/Uo
%Uo/Ui = 4/(Re^.5) to ODE  d3f/d'eta^3 + f*d2f/d'eta^2 + 2*(df/d'eta)^2 = 0
%for eta=0, f=df/d'eta=0;  for eta=H=10, df/d'eta=d2f/d'eta^2 = 0
%% I.a
%   implicit solution given as
etaf=@(f,n) log(sqrt(1+sqrt(f)+f)./(1-sqrt(f)))...
    +sqrt(3)*atan(sqrt(3.*f)./(2+sqrt(f)))-n;
eta=[0:0.05:10]'; f=zeros((10/.05)+1,1);
for i=1:length(eta)
    f(i)=bisectE(@(f) etaf(f,eta(i)),0,1,1e-8);
end
figure(1); hold off; plot(eta,f); xlabel('eta'); ylabel('script f')

%% I.b
% derivative
df=diffc2(f,0.05);
figure(2); hold off; plot(eta,df);
xlabel('eta'); ylabel('derivative script f')

%% I.c
% derivative and proving the thing
d2f=diffc2(df,0.05);
figure(3); plot(eta,d2f);
xlabel('eta'); ylabel('second derivative script f')
fprintf('Cf*(Re^(5/4))/8: %5.4f, compare %5.4f\n',d2f(1),1.778/8)

%% I.d
% maximum velocity
dfetaSpln=spline(eta,df);
fetaReg=@(eta,neg,dwn) ((-2*neg)+1)*(ppval(dfetaSpln,eta)-(dwn*.01));
etaMax=goldmin(fetaReg,0,10,1e-6,9999,true,false);
figure(2); hold on;
plot((etaMax-.1):.01:(etaMax+.1),...
    ppval(dfetaSpln,(etaMax-.1):.01:(etaMax+.1)),...
    etaMax,ppval(dfetaSpln,etaMax),'*')
fprintf('n [eta]: %5.4f, compare 2.029\n',etaMax)
fprintf('f''(n) =  %5.4f, compare %5.4f\n',...
    ppval(dfetaSpln,etaMax),2^(-5/3))

%% I.e
% wall jet momentum flux
Z=trapz(eta,df.^2)*trapz(eta,df);
fprintf('[flux]/(Ui*nu^2*4^3): %5.4f, compare %5.4f\n',Z,128/(9*(4^3)))

%% I.f
% Shear thickness location
delta=bisectE(fetaReg,etaMax,10,1e-8,9999,false,true);
fprintf(strcat('shear layer thickness delta1\n',...
    'delta = %5.4f, compare 6.72\n'),delta)

%% I.g
% Velocity at eta=H
v_H=(3*eta(end)*df(end)-f(end))*Ui;%eta(end)==10
fprintf('Analytical: %5.4f\n',v_H)

%% I.h
% Solve Eq.(1) and compare to f(n)
zah=fzero(@res,0,[],1);
[etaODE,fODE]=ode45(@dydx,[0 10],[0 0 zah],[],1);
figure(1); hold on; plot(etaODE,fODE(:,1),'y--');

%% II
Pr=.7;
%% II.i
% Solve ODE
zai=fzero(@res,0,[],2,fODE(:,1),etaODE);
[n,theta]=ode45(@dydx,[0 10],[1 zai],[],2,fODE(:,1),etaODE);
figure(4); hold off; plot(theta(:,1),n); xlim([0,1]);
xlabel('theta'); ylabel('eta')
%there's a theta value of -1.687e-16, so we squish the plot display limit.

%% II.j
% Thermal boundary location
thetaEtaSpln=spline(n,theta(:,1));
thEta=@(eta) ppval(thetaEtaSpln,eta)-.01;
deltaT=bisectE(thEta,0,10,1e-8,9999);
fprintf(strcat('thermal boundary layer thickness deltaT\n',...
    'deltaT = %5.4f, compare %5.4f\n'),deltaT,(delta*(Pr^(-1/3))))

%% II.k
% Tempature gradient
dtheta=diffc2(theta(:,1),n(2)-n(1));
fprintf('theta''(0): %5.4f, compare %5.4f\n',-1*dtheta(1),.235*(Pr^(1/3)))

%% II.l
% Solve Eq.(4) and compare to theta(n)
fetaSpln=spline(eta,f);
[x,theta]=CFD(fetaSpln,5);
%f=@(eta) ppval(fetaSpln,eta);
figure(4); hold on; plot(theta,x)
[xx,thetaa]=CFD(fetaSpln,11); [dad,thetaaa]=CFD(fetaSpln,33);
plot(thetaa,xx,thetaaa,dad,'--'); legend('I.i','5node','11node','33node')