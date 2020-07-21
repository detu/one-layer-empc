function [Tall, xmeasureAll, uAll, ObjVal, runtime] = spfNmpcCyclicWithoutINmpc(optProblem, system, mpciterations, N, T, tmeasure, xmeasure, u0, scr, varargin)
%SPFNMPCCYCLICWITHOUTINMPC Summary of this function goes here
% 
% [OUTPUTARGS] = SPFNMPCCYCLICWITHOUTINMPC(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2018/11/28 15:38:28 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2018


% to do:
% 1. remove ideal multistage NMPC
% 2. in the offline step: solve ideal multistage.. 

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

    
    % include parameter to measurement data
    xmeasure = [xmeasure;delta];
    
    % obtain new initial value
    t0 = tmeasure;
    x0 = xmeasure;
    
    %% FIRST ITERATION CALL ideal-msNMPC
    %if mpciter == 1 || mpciter == 25 || mpciter == 80
    if mpciter == 1 % run ideal ms-NMPC only at 1st iteration
        [u_pf_opt,primal,t0,xmeasure,elapsedqp] = msNmpcCyclicNew(optProblem, N, u0, tmeasure, xmeasure, scr);
        %z1 = [primal(9:11,1);delta];
    else
        %call advanced-step-multistage-NMPC
        % CHANGE THE OFFLINE STEP with idea-msNMPC 
        z1       = xmeasure;
        %z1 = [z1;delta];
        %% BACKGROUND STEP
%         for i=1:scr.numScr
%             % call MPC with v0=uk and z0=xk
%             scenario = i;
%             [primalNLP(:,i), dualNLP(:,i), lb, ub, ~, params] = solveMsOCP(optProblem, system, N, t0, x0, u0, T, mpciter, u_pf_opt, x_pf_opt, z1);
%         end
        %[primalNLP, dualNLP, lb, ub] = OfflineMsNmpcCyclic(optProblem, N, u0, tmeasure, xmeasure, scr);
        [primalNLP, dualNLP, lb, ub] = OfflineMsNmpcCyclic(optProblem, N, u0, tmeasure, z1, scr);
        
        % BEFORE SOLVING NLP-SENSITIVITY CHOOSE THE CLOSEST SCENARIO FROM MEASUREMENT
        distance = zeros(scr.numScr,1);
        for i=1:scr.numScr
            %distance(i) = norm(primalNLP(9:11,i) - x0(1:nx),2); % hardcode position 9 until 11
            distance(i) = norm(primalNLP(10:12,i) - x0(1:nx),2);
        end
        [~,index] = min(distance);
        scenario  = index;
        fprintf('index number = %d \n', index);
        pfNLP     = primalNLP(:,index);
        %dfNLP     = dualNLP(:,index);
        dfNLP.lam_x  = dualNLP.lam_x(:,index);
        dfNLP.lam_g  = dualNLP.lam_g(:,index);
        %% ONLINE STEP
        
        % APPLY PATH-FOLLOWING ALGORITHM HERE!
        % to-do: ADD THE TIME-VARYING PARAMETER AS PARAMETER !
        
        % re-arrange NLP solutions
        [~, x_nlp_opt] = arrangeOptResultsNew(pfNLP, lb, ub, N);
        
        p_init  = [pfNLP(10:12);delta];
        p_final = x0;
        xstart  = pfNLP;
        ystart  = dfNLP;
        
        % choose number of path-following step
        %delta_t = 0.25;
        delta_t = 1;
        
        lb_init = lb;
        ub_init = ub;
        
        % NLP sensitivity (predictor-corrector)
        %[primalPF, ~, elapsedqp] = jpredictor_licq_pure_3(@(p)cyclicMultistageDerivatives(p), p_init, p_final, xstart, ystart, delta_t, lb_init, ub_init, 0, N);
        [primalPF, ~, elapsedqp] = jpredictor_licq_pure_3(@(p)NewCyclicMultistageDerivatives(p), p_init, p_final, xstart, ystart, delta_t, lb_init, ub_init, 0, N);
        
        [u_pf_opt, x_pf_opt] = arrangeOptResultsNew(primalPF, lb, ub, N);
        z1 = x_pf_opt(1:nx,5); % 4 = (d+1) + 1 (d=number of collocation point)
    end
    
    %% INJECT TO PLANT
    %x0 = x0_measure; 
    %t0 = t0_measure;
    x0 = xmeasure;
    [tmeasure, xmeasure] = applyControl(system, T, t0, x0(1:nx), u_pf_opt(1,:));     
    
    %% GET CLOSED-LOOP RESPONSE
    %ObjVal(mpciter)      = computeObjFuncCstr(u_pf_opt(1,1),xmeasure); % CSTR only
    ObjVal(mpciter) = 0; % for time being (temporary)
    
    %% COLLECT CLOSED-LOOP DATA
    Tall                = [Tall t0];
    Xall(mpciter+1,:)   = x0(1:nx)';
    Uall(mpciter+1,:)   = u0(:,1);
    xmeasureAll = [xmeasureAll;xmeasure];
    uAll        = [uAll;u_pf_opt(1,1)];
    runtime     = [runtime;elapsedqp];
    
    %% SHIFT CONTROL INPUT
    u0 = shiftHorizon(u_pf_opt(1,:));
    
    mpciter = mpciter+1;
end
xmeasureAll = reshape(xmeasureAll,nx,mpciterations);
end

