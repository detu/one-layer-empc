function [J,g,w0,w,lbg,ubg,lbw,ubw,t0,x0_measure] = CollectVariablesForEachScenarioNew(optProblem, N, u0, tmeasure, xmeasure, scrNo)
%COLLECTVARIABLESFOREACHSCENARIONEW Summary of this function goes here
% 
% [OUTPUTARGS] = COLLECTVARIABLESFOREACHSCENARIONEW(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2019/01/16 14:11:07 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2019


global scenario;

    
%% scenario tree
% apply different price values for different realizations
switch scrNo
    case 1
        scenario = 1;
        
    case 2
        scenario = 2;
        
    case 3
        scenario = 3;
        
    case 4
        scenario = 4;
end

%% building an NLP problem for each realization
% obtain new initial value
[t0, x0] = measurementInitialValue ( tmeasure, xmeasure );

%x0         = max(min(x0,1.0),0); % restrict to boundaries
x0_measure =  x0;    % without noise
%x0_measure = max(min(x0_measure,1.0),0); % restrict to boundaries

[J,g,w0,w,lbg,ubg,lbw,ubw] = buildOptimalControlProblem(optProblem, N, x0, u0, x0_measure);

end
