function xprime = cyclic_ode(t,X)
%CYCLIC_ODE Summary of this function goes here
%
% ODE driver for Cyclic (nondissipative) process 
%
% [OUTPUTARGS] = CYCLIC_ODE(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2018/01/17 21:57:58 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2018

% global u NT %Make perturbed inputs/disturbances available to model
global uc;

% Store all inputs and disturbances
u_all = uc;

xprime=mod_cyclic(t,X,u_all);

end
