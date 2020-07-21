function xprime = mod_cyclic(t,X,U) 
%MOD_CYCLIC Summary of this function goes here
% 
% [OUTPUTARGS] = MOD_CYCLIC(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2018/01/17 21:58:16 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2018

%% Cyclic model
% parameter values
global delta;

% state equations
x1     = X(1);
x2     = X(2);
x3     = X(3);
u1     = U;

x1dot  = 1 - 1e4*x1^2*exp(-1/x3) - 400*x1*exp(-delta/x3) - x1;
x2dot  = 1e4*x1^2*exp(-1/x3) - x2;
x3dot  = u1 - x3;

% Output
xprime=[x1dot;x2dot;x3dot];

end
