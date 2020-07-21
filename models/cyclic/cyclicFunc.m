function [t,states,xdot,inputs] = cyclicFunc()
%CYCLICFUNC Summary of this function goes here
% 
% [OUTPUTARGS] = CYCLICFUNC(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2018/01/17 22:23:45 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2018

import casadi.* 
global delta;

% symbolic primitives
t  = SX.sym('t');
x1 = SX.sym('x1');  
x2 = SX.sym('x2');
x3 = SX.sym('x3');
states = [x1;x2;x3];
u      = SX.sym('u');
inputs = u;

% Cyclic model
x1dot  = 1 - 1e4*x1^2*exp(-1/x3) - 400*x1*exp(-delta/x3) - x1;
x2dot  = 1e4*x1^2*exp(-1/x3) - x2;
x3dot  = u - x3;
xdot   = [x1dot;x2dot;x3dot];

end
