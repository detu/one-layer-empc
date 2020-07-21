function [Tall, xmeasureAll, uAll, ObjVal, runtime] = spfNmpcCyclic(optProblem, system, mpciterations, N, T, tmeasure, xmeasure, u0, scr, varargin)
%SPFNMPCCYCLIC Summary of this function goes here
% 
% [OUTPUTARGS] = SPFNMPCCYCLIC(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2018/11/06 23:16:16 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2018
import casadi.*
%% ITERATE ON MPC RUNNING
% outer loop of the scenario
global nx delta xf uf check scenario;
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
    
    % change parameter value
    if mpciter == 1
        check = 0;
        delta = 0.55;
        keyboard;
    end
    
    if mpciter == 25
        check = 0;
        delta = 1.1;
        keyboard;
    end
    
    if mpciter == 80
        check = 0;
        delta = 2.2;
        keyboard;
    end
    
    % obtain new initial value
    t0 = tmeasure;
    x0 = xmeasure;
    
    %% FIRST ITERATION CALL ideal-msNMPC
   if mpciter == 1 || mpciter == 25 || mpciter == 80        
        [u_opt,x_opt, t0,xmeasure,elapsedqp] = msNmpcCyclic(optProblem, N, u0, tmeasure, xmeasure, scr);
   else
        % call advanced-step-multistage-NMPC
        z1 = xmeasure;
        %% BACKGROUND STEP
        for i=1:scr.numScr
            % call MPC with v0=uk and z0=xk
            scenario = i;
            [primalNLP(:,i), dualNLP(:,i), lb, ub, ~, params] = solveMsOCP(optProblem, system, N, t0, x0, u0, T, mpciter, u_pf_opt, x_pf_opt, z1);
        end
        
        % BEFORE SOLVING NLP-SENSITIVITY CHOOSE THE CLOSEST SCENARIO FROM MEASUREMENT
        distance = zeros(scr.numScr,1);
        for i=1:scr.numScr
            distance(i) = norm(primalNLP(9:11,i) - x0,2); % hardcode position 9 until 11
        end
        [~,index] = min(distance);
        scenario  = index;
        fprintf('index number = %d \n', index);
        pfNLP     = primalNLP(:,index);
        dfNLP     = dualNLP(:,index);
        %% ONLINE STEP
        
        % APPLY PATH-FOLLOWING ALGORITHM HERE!
        % re-arrange NLP solutions
        [~, x_nlp_opt] = arrangeOptResults(pfNLP, lb, ub, N);
        
        p_init  = pfNLP(9:11);
        p_final = x0;
        xstart  = pfNLP;
        ystart  = dfNLP;
        
        % choose number of path-following step
        %delta_t = 0.5;
        delta_t = 1;
        
        lb_init = lb;
        ub_init = ub;
        
        % NLP sensitivity (predictor-corrector)
        [primalPF, ~, elapsedqp] = jpredictor_licq_pure_3(@(p)cyclicMultistageDerivatives(p), p_init, p_final, xstart, ystart, delta_t, lb_init, ub_init, 0, N);
        
        [u_pf_opt, x_pf_opt] = arrangeOptResults(primalPF, lb, ub, N);
        z1 = x_pf_opt(1:nx,5); % 4 = (d+1) + 1 (d=number of collocation point)
    end
    
    %% INJECT TO PLANT
    %x0 = x0_measure; 
    %t0 = t0_measure;
    x0 = xmeasure;
    %[tmeasure, xmeasure] = applyControl(system, T, t0, x0, u_pf_opt(1,:));  
    [tmeasure, xmeasure] = applyControl(system, T, t0, x0, u_opt(1,:));
    
    %% GET CLOSED-LOOP RESPONSE
    %ObjVal(mpciter)      = computeObjFuncCstr(u_pf_opt(1,1),xmeasure); % CSTR only
    ObjVal(mpciter) = 0; % for time being (temporary)
    
    %% COLLECT CLOSED-LOOP DATA
    Tall                = [Tall t0];
    Xall(mpciter+1,:)   = x0';
    Uall(mpciter+1,:)   = u0(:,1);
    xmeasureAll = [xmeasureAll;xmeasure];
    %uAll        = [uAll;u_pf_opt(1,1)];
    uAll        = [uAll;u_opt(1,1)];
    runtime     = [runtime;elapsedqp];
    
    %% SHIFT CONTROL INPUT
    %u0 = shiftHorizon(u_pf_opt(1,:));
    u0 = shiftHorizon(u_opt(1,:));
    
    mpciter = mpciter+1;
end
xmeasureAll = reshape(xmeasureAll,nx,mpciterations);
end
