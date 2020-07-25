function [u_nlp_opt,plotState, xu_opt] = arrangeOptResultsNew(data, lb, ub, N)
%ARRANGEOPTRESULTSNEW Summary of this function goes here
% 
% [OUTPUTARGS] = ARRANGEOPTRESULTSNEW(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2019/01/16 15:25:47 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2019

global nk nx nu d ns;

% optimized initial state
x0_opt       = data(1:nx); 
data(1:nx) = [];
lb(1:nx)   = [];
ub(1:nx)   = []; 

% optimized steady-state
xu_opt        = data(1:nx+nu); 
xu_opt        = xu_opt';
data(1:nx+nu) = [];
lb(1:nx+nu)   = [];
ub(1:nx+nu)   = []; 

data         = reshape(data, (nu + (nx+ns)*d + (nx+ns)), N*nk);
u_nlp_opt    = data(1:nu,1:N*nk);
data(1:nu,:) = [];     % remove optimized controls

lb0          = lb(1:nx+ns);
%lb(1:nx)     = [];
lb           = reshape(lb, (nu + (nx+ns)*d + (nx+ns)), N*nk);
lbU          = lb(1:nu,1:N*nk);
lb(1:nu,:)   = [];

ub0          = ub(1:nx+ns);
%ub(1:nx)     = [];
ub           = reshape(ub, (nu + (nx+ns)*d + (nx+ns)), N*nk);
ubU          = ub(1:nu,1:N*nk);
ub(1:nu,:)   = [];

% prepare a matrix for plotting
nState    = (nx+ns) + N*nk*(d+1)*(nx+ns);
nPoint    = nState/(nx+ns);
plotState = zeros(nx+ns,nPoint);
plotState(1:nx,1) = x0_opt(1:nx);

plotLb         = zeros(nx+ns,nPoint);
plotLb(:,1)    = lb0;

plotUb         = zeros(nx+ns,nPoint);
plotUb(:,1)    = ub0;

% extract states from each collocation point and each time horizon
sInd    = 2; % initial index row
for i=1:N*nk
    temp    = data(:,i);
    numCol  = size(temp,1);
    numRow  = numCol/(nx+ns);
    temp    = reshape(temp,nx+ns,numRow);
    plotState(:,sInd:(numRow+sInd-1)) = temp;
    
    tempLb = lb(:,i);
    tempLb = reshape(tempLb,nx+ns,numRow);
    plotLb(:,sInd:(numRow+sInd-1)) = tempLb;
    
    tempUb = ub(:,i);
    tempUb = reshape(tempUb,nx+ns,numRow);
    plotUb(:,sInd:(numRow+sInd-1)) = tempUb;
    
    sInd    = numRow + sInd;
end

end
