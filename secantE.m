function [root,errr,trys]=secantE(fun,d,guess,varargin)
% 'secant' root finding method; modified newtraph.m
%   supply function, perturbation, and guess. optional Approximate Relative
%   Error 'ARE' name-value, good digits 'digits' name-value, and maximum
%   iterations 'tries' name-value.
ip2=inputParser;%get the optional parameters for tolerance, iterations
addRequired(ip2,'fun'); addRequired(ip2,'d'); addRequired(ip2,'guess');
addParameter(ip2,'ARE',1e-4,@isscalar);
addParameter(ip2,'digits',0,@(x) x==abs(round(x)));
addParameter(ip2,'tries',9999,@(x) x==abs(round(x)));
parse(ip2,fun,d,guess,varargin{:})%ip.Results.[var] to call arguments
if ip2.Results.digits==0
    error=ip2.Results.ARE;
else
    error=.5*10^(2-ip2.Results.digits);
end
root=ip2.Results.guess; trys = 0;
while true
    if fun(root)==0
        errr=0; break
    end
    rootNew=root-((fun(root)*d*root)/(fun(root+(d*root))-fun(root)));
    errr=abs((rootNew-root)/rootNew);
    root=rootNew;
    trys=trys+1;
    if (errr<error)
        break
    elseif (trys>ip2.Results.tries)
        disp('failed; too many tries'); break
    end
end