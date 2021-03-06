function [root,fx,ARE,trys]=bisectE(fun,xl,xu,ea,tries,varargin)
% bisect: root location zeroes
% [root,fx,ea,iter]=bisect(func,xl,xu,es,maxit,p1,p2,...):
% uses bisection method to find the root of func
% input:
% fun= name of function
% xl, xu = lower and upper guesses
% ea = desired relative error (default = 0.0001%)
% tries = maximum allowable iterations (default = 50)
% p1,p2,... = additional parameters used by func
% output:
% root = real root
% fx = function value at root
% ARE = approximate relative error
% trys = number of iterations
if nargin<3, error('bisect(): more arguments please'), end
test=fun(xl,varargin{:})*fun(xu,varargin{:});
if test>0, error('no sign change'), end
if fun(xl,varargin{:})==0, root=xl; ARE=0; go=false;
elseif fun(xu,varargin{:})==0, root=xu; ARE=0; go=false;
else, root=xl; ARE=1; go=true;
end%check the lower and upper points for zero
if (nargin<4)||(isempty(ea)), ea=1e-6; end
if (nargin<5)||(isempty(tries)), tries=50; end
trys=0;
while go
  rtold=root;
  root=(xl+xu)/2;
  trys=trys+1;
  if root~=0, ARE=abs((root-rtold)/root); end
  test=fun(xl,varargin{:})*fun(root,varargin{:});
  if test<0
    xu=root;
  elseif test>0
    xl=root;
  else
    if fun(root,varargin{:})==0, ARE=0;
    else, root=xl; ARE=0;
    end%redundantly check lower point for case of upper bound 'Inf'
  end
  if (ARE<=ea)||(trys>=tries), break; end
end
fx = fun(root, varargin{:});