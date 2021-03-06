function [Tall, xmeasureAll, uAll, ObjVal, runtime] = sNmpcDistColA(optProblem, system, mpciterations, N, T, tmeasure, xmeasure, u0, scr, varargin)
%SNMPCDISTCOLA Summary of this function goes here
% 
% [OUTPUTARGS] = SNMPCDISTCOLA(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2018/12/27 23:48:32 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2018

import casadi.*
%% ITERATE ON MPC RUNNING
% outer loop of the scenario
global nx mpciter xf uf check scenario;
xf      = xmeasure;
uf      = u0(1);
mpciter = 1;
z1      = xmeasure;
% data initialization
Tall    = [];
Xall    = zeros(mpciterations,size(xmeasure,1));
Uall    = zeros(mpciterations,size(u0,1));
xmeasureAll = [];
uAll        = [];
runtime     = [];
u_pf_opt    = [];
x_pf_opt    = [];

while (mpciter <= mpciterations)
    %% printing iteration
    fprintf('-----------------------------\n');
    fprintf('MPC iteration: %d\n', mpciter);
     
    % obtain new initial value
    t0 = tmeasure;
    x0 = xmeasure;
    
    % multistage MPC
    [u_opt, x_opt, t0, xmeasure, elapsedqp] = msNmpcDistColA(optProblem, N, u0, tmeasure, xmeasure, scr);

    
    %% INJECT TO PLANT
    %x0 = x0_measure; 
    %t0 = t0_measure;
    x0 = xmeasure;
    %[tmeasure, xmeasure] = applyControl(system, T, t0, x0, u_opt(1,:));
    [tmeasure, xmeasure] = applyControl(system, T, t0, x0, u_opt(:,1));
    
    %% GET CLOSED-LOOP RESPONSE
    %ObjVal(mpciter)      = computeObjFuncCstr(u_pf_opt(1,1),xmeasure); % CSTR only
    ObjVal(mpciter) = 0; % for time being (temporary)
    
    if (mpciter == 50)
        u_optIdeal50 = u_opt;
        x_optIdeal50 = x_opt;
        save i-msNMPCIteration50.mat u_optIdeal50 x_optIdeal50;
    end
    
    if (mpciter == 100)
        u_optIdeal100 = u_opt;
        x_optIdeal100 = x_opt;
        save i-msNMPCIteration100.mat u_optIdeal100 x_optIdeal100;
    end
    
    %% COLLECT CLOSED-LOOP DATA
    Tall                = [Tall t0];
    Xall(mpciter+1,:)   = x0';
    Uall(mpciter+1,:)   = u0(:,1);
    xmeasureAll = [xmeasureAll;xmeasure];
    %uAll        = [uAll;u_opt(1,1)];
    uAll        = [uAll;u_opt(:,1)];
    runtime     = [runtime;elapsedqp];
    
    %% SHIFT CONTROL INPUT
    %u0 = shiftHorizon(u_opt(1,:));
    u0 = shiftHorizon(u_opt);
    
    mpciter = mpciter+1;
end
xmeasureAll = reshape(xmeasureAll,nx,mpciterations);

end
