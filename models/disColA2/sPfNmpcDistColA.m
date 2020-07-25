function [Tall, xmeasureAll, uAll, ObjVal, runtime] = sPfNmpcDistColA(optProblem, system, mpciterations, N, T, tmeasure, xmeasure, u0, scr, varargin)
%SPFNMPCDISTCOLA Summary of this function goes here
% 
% [OUTPUTARGS] = SPFNMPCDISTCOLA(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2019/01/15 23:27:02 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2019

import casadi.*
%% ITERATE ON MPC RUNNING
% outer loop of the scenario
global nx xf uf mpciter scenario;
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
    
%     % include parameter to measurement data
%     xmeasure = [xmeasure;delta];
    
    % obtain new initial value
    t0 = tmeasure;
    x0 = xmeasure;
    
    %% FIRST ITERATION CALL ideal-msNMPC
    if mpciter == 1 || mpciter == 50 || mpciter == 100 % run ideal ms-NMPC only at 1st iteration
        
        [u_pf_opt, primal, t0, xmeasure, elapsedqp] = msNmpcDistColA(optProblem, N, u0, tmeasure, xmeasure, scr);
        
        %z1 = [primal(9:11,1);delta];
%         if mpciter == 50
%             save pf-msNMPCIteration50.mat u_pf_opt primal;
%         end
%         if mpciter == 100
%             save pf-msNMPCIteration100.mat u_pf_opt primal;
%         end
    else
        %call advanced-step-multistage-NMPC
        % CHANGE THE OFFLINE STEP with idea-msNMPC 
        z1       = xmeasure;

        %% BACKGROUND STEP
        [primalNLP, dualNLP, lb, ub] = OfflineMsNmpcDistColA(optProblem, N, u0, tmeasure, z1, scr);
        
        % BEFORE SOLVING NLP-SENSITIVITY CHOOSE THE CLOSEST SCENARIO FROM MEASUREMENT
        distance = zeros(scr.numScr,1);
        for i=1:scr.numScr
            distance(i) = norm(primalNLP(520:603,i) - x0,2);
        end
        [~,index] = min(distance);
        scenario  = index;
        fprintf('index number = %d \n', index);
        pfNLP        = primalNLP(:,index);
        dfNLP.lam_x  = dualNLP.lam_x(:,index);
        dfNLP.lam_g  = dualNLP.lam_g(:,index);
        %% ONLINE STEP
        
        % APPLY PATH-FOLLOWING ALGORITHM HERE!
        % to-do: ADD THE TIME-VARYING PARAMETER AS PARAMETER !
        
        % re-arrange NLP solutions
        [~, x_nlp_opt] = arrangeOptResultsNew(pfNLP, lb, ub, N);
        
        p_init  = pfNLP(520:603);
        p_final = x0;
        xstart  = pfNLP;
        ystart  = dfNLP;
        
        % choose number of path-following step
        %delta_t = 0.25;
        delta_t = 1;
        
        lb_init = lb;
        ub_init = ub;
        
        % NLP sensitivity (predictor-corrector)
        [primalPF, ~, elapsedqp] = jpredictor_licq_pure_3(@(x,y,p,N)DistColAMultistageDerivatives(x,y,p,N), p_init, p_final, xstart, ystart, delta_t, lb_init, ub_init, 0, N);
        
        [u_pf_opt, x_pf_opt] = arrangeOptResultsNew(primalPF, lb, ub, N);
        z1 = x_pf_opt(1:nx,5); % 4 = (d+1) + 1 (d=number of collocation point)
    end
    
    %% INJECT TO PLANT
    x0 = xmeasure;    
    [tmeasure, xmeasure] = applyControl(system, T, t0, x0, u_pf_opt(:,1));
    
    %% GET CLOSED-LOOP RESPONSE
    %ObjVal(mpciter)      = computeObjFuncCstr(u_pf_opt(1,1),xmeasure); % CSTR only
    ObjVal(mpciter) = 0; % for time being (temporary)
    
    %% COLLECT CLOSED-LOOP DATA
    Tall                = [Tall t0];
    Xall(mpciter+1,:)   = x0';
    Uall(mpciter+1,:)   = u0(:,1);
    xmeasureAll = [xmeasureAll;xmeasure];
    uAll        = [uAll;u_pf_opt(:,1)];
    runtime     = [runtime;elapsedqp];
    
    %% SHIFT CONTROL INPUT
    u0 = shiftHorizon(u_pf_opt);
    
    mpciter = mpciter+1;
end
xmeasureAll = reshape(xmeasureAll,nx,mpciterations);

end
