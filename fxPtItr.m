function [root,errr,trys] = fxPtItr(fun,varargin)
%Fixed-point iteration root finding.
ip4=inputParser;%get the optional parameters for guess, tolerance
addRequired(ip4,'fun'); addOptional(ip4,'guess',-1,@isscalar)
addParameter(ip4,'ARE',1e-4,@isscalar);
addParameter(ip4,'digits',0,@(x) x==abs(round(x)));
parse(ip4,fun,varargin{:})%ip.Results.[var] to call arguments
if ip4.Results.digits==0
    error=ip4.Results.ARE;
else
    error=.5*10^(2-ip4.Results.digits);
end
root=ip4.Results.guess; trys=0;
while true
    if ip4.Results.guess<(((4^(-4/3))-10)/5)
        disp('this won''t necessarily converge');
    end
    rootNew=fun(root);
    errr=abs((rootNew-root)/rootNew);
    root=rootNew;
    trys=trys+1;
    if (errr<error)
        break
    elseif (trys>99999)
        disp('failed; too many tries'); break
    end
end