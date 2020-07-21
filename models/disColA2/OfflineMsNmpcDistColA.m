function [x_nlp_opt, y_nlp_opt, lbwi, ubwi] = OfflineMsNmpcDistColA(optProblem, N, u0, tmeasure, xmeasure, scr)
%OFFLINEMSNMPCDISTCOLA Summary of this function goes here
% 
% [OUTPUTARGS] = OFFLINEMSNMPCDISTCOLA(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2019/01/16 14:07:03 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2019

import casadi.*

%% ITERATE ON THE SCENARIO TREES
% prepare optimization variables
J   = 0;
g   = {};
w0  = [];
w   = {};
lbg = [];
ubg = [];
lbw = [];
ubw = [];
nac = [];
for i=1:scr.numScr
    
    %% GATHER ALL VARIABLES FROM EACH SCENARIO
    % sNmpcCstr(optProblem, system, mpciterations, N, T, tmeasure, xmeasure, u0, scr, varargin)
    %[Tall, xmeasureAll, uAll, ObjVal, primalNLP, params, runtime] = CollectVariablesForEachScenario(optProblem, system, mpciterations, N, T, u0, tmeasure, xmeasure, i);
    [Ji,gi,w0i,wi,lbgi,ubgi,lbwi,ubwi,t0,x0_measure] = CollectVariablesForEachScenarioNew(optProblem, N, u0, tmeasure, xmeasure, i);
    
    %% collection optimization variables for all realization
    J   = J + Ji;
    %g   = [g;gi];
    g   = {g{:},gi{:}};
    w0  = [w0;w0i];
    w   = {w{:},wi{:}};
    lbg = [lbg;lbgi];
    ubg = [ubg;ubgi];
    lbw = [lbw;lbwi];
    ubw = [ubw;ubwi];
    %nac = [nac;wi{3}];  %control input
    nac = [nac wi{3}];
end

%% BUILD NON-ANTICIPATIVITY CONSTRAINTS (NAC)
[numCol, numRow] = size(nac);
countNacSize = 0;
for i=1:(numRow-1)
    temp = nac(:,i+1)-nac(:,i);
    g    = {g{:}, temp};
    lbg  = [lbg; zeros(numCol,1)];
    ubg  = [ubg; zeros(numCol,1)];
    countNacSize = countNacSize + size(temp,1);
end
numNac   = countNacSize;

%% SOLVE NLP FOR ALL REALIZATIONS
prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
options = struct;
options.ipopt.tol                = 1e-12;
options.ipopt.constr_viol_tol    = 1e-12;
solver = nlpsol('solver', 'ipopt', prob, options);

% Solve the NLP
startnlp   = tic;
sol        = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, 'lbg', lbg, 'ubg', ubg);
elapsednlp = toc(startnlp);
fprintf('IPOPT solver runtime = %f\n',elapsednlp);
success = strcmp(solver.stats.return_status,'Infeasible_Problem_Detected');
if (success)
    keyboard;
end

u       = full(sol.x);
y.lam_x = full(sol.lam_x);
y.lam_g = full(sol.lam_g);
% remove the NAC (check how many size of the NAC constraint)
rem                  = numNac-1;
y.lam_g(end-rem:end) = [];
%% INJECT TO PLANT
% reshape optimized variables
numColX    = size(u,1)/scr.numScr;
numColY    = size(y.lam_g,1)/scr.numScr;
x_nlp_opt = reshape(u,numColX,scr.numScr);
y_nlp_opt.lam_x = reshape(y.lam_x,numColX,scr.numScr);
y_nlp_opt.lam_g = reshape(y.lam_g,numColY,scr.numScr);


end
