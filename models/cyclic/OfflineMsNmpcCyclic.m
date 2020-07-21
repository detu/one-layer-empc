function [x_nlp_opt, y_nlp_opt, lbwi, ubwi] = OfflineMsNmpcCyclic(optProblem, N, u0, tmeasure, xmeasure, scr)
%OFFLINEMSNMPCCYCLIC Summary of this function goes here
% 
% [OUTPUTARGS] = OFFLINEMSNMPCCYCLIC(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2018/11/30 13:46:10 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2018

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
    %for i=1:4
    % call NMPC for each scenario
    % divide into 2 realizations: +- 10%
    
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
    nac = [nac;wi{3}];  %control input
end

%% BUILD NON-ANTICIPATIVITY CONSTRAINTS (NAC)
% ingat ada steady-state constraint JANGAN DIPAKE sebagai NAC!
% let try robust horizon = 1
numNac = size(nac,1);
for i=1:(numNac-1)
    g   = {g{:}, nac(i+1,1)-nac(i,1)};
    lbg = [lbg; 0];
    ubg = [ubg; 0];
end


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
%lamda  = full(sol.lam_g);
%objVal = full(sol.f);
% remove the NAC
rem                  = numNac-2;
y.lam_g(end-rem:end) = [];
%% INJECT TO PLANT
% reshape optimized variables
numColX    = size(u,1)/scr.numScr;
numColY    = size(y.lam_g,1)/scr.numScr;
x_nlp_opt = reshape(u,numColX,scr.numScr);
y_nlp_opt.lam_x = reshape(y.lam_x,numColX,scr.numScr);
y_nlp_opt.lam_g = reshape(y.lam_g,numColY,scr.numScr);

% % loop for all scenarios for plotting purposes
% u_nlp_opt     = zeros(scr.numScr, N);
% x_nlp_opt_x1  = zeros(scr.numScr, 801); % 801 is number of state's collocation point
% x_nlp_opt_x2  = zeros(scr.numScr, 801);
% x_nlp_opt_x3  = zeros(scr.numScr, 801);
% 
% for i= 1:scr.numScr
%     uopt   = u(:,i);
%     [u_nlp_opt(i,:), x_nlp_opt, ~] = arrangeOptResults(uopt, lbw(1:numCol,1), ubw(1:numCol,1), N);
%     x_nlp_opt_x1(i,:) = x_nlp_opt(1,:);
%     x_nlp_opt_x2(i,:) = x_nlp_opt(2,:);
%     x_nlp_opt_x3(i,:) = x_nlp_opt(3,:);
% end


end
