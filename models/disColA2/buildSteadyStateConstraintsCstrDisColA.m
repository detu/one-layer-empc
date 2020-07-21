function [ceq,x] = buildSteadyStateConstraintsCstrDisColA(scenario)
%BUILDSTEADYSTATECONSTRAINTSCSTRDISCOLA Summary of this function goes here
% 
% [OUTPUTARGS] = BUILDSTEADYSTATECONSTRAINTSCSTRDISCOLA(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2019/01/09 15:12:40 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2019


import casadi.*

%% parameter values
NT  = 41;                            % number of trays
LT  = 2.827;                         % Reflux
VB  = 3.454;                         % Boilup
F   = 1.0;                           % Feedrate
zF  = 1.0;                           % Feed composition at CSTR 
D   = 0.5;                           % Distillate flow
B   = 0.5;                           % Bottoms flow 
qF  = 1.0;                           % Feed liquid fraction

% create a struct to make them nice
%dist.F_0  = 0.3;
dist.F_0  = 0.32;
dist.NT   = NT;
dist.zF   = zF;
dist.qF   = qF;

% price setting
price.pf = 1; 
price.pV = 0.02;
price.pB = 2; 
price.pD = 0;

% Symbolic primitives
x = {};
l = {};
for i=1:2*NT+2       % +2 due to CSTR
   x{i} = SX.sym(['x_' num2str(i) num2str(scenario)]);
   l{i} = SX.sym(['l_' num2str(i) num2str(scenario)]);
end
u1  = SX.sym(['u1' num2str(scenario)]);   % LT
u2  = SX.sym(['u2' num2str(scenario)]);   % VB
u3  = SX.sym(['u3' num2str(scenario)]);   % F
u4  = SX.sym(['u4' num2str(scenario)]);   % D
u5  = SX.sym(['u5' num2str(scenario)]);   % B
% concatenate states and controls 

%xo = {x{:}, u1, u2, u3, u4, u5};

x   = vertcat(x{:});
x   = [x;u1;u2;u3;u4;u5];

% define the dynamics as equality constraints and additional inequality
% constraints (lbx, ubx, lbg, and ubg)
ceq = buildModelEq(x,dist);

end

function ceq = buildModelEq(u,dist)
import casadi.* 
NT = dist.NT;
% Location of feed stage (stages are counted from the bottom):
NF = 21;
% Relative volatility
alpha = 1.5;
% Nominal liquid holdups
Muw   = 0.5;
% Data for linearized liquid flow dynamics (does not apply to reboiler and condenser):
taul = 0.063;     	% time constant for liquid dynamics (min)
F0   = 1;	 	    % Nominal feed rate (kmol/min) 
qF0  = 1; 		    % Nominal fraction of liquid in feed 
L0   = 2.70629;     % Nominal reflux flow (from steady-state data)
L0b  = L0 + qF0*F0;	% Nominal liquid flow below feed (kmol/min)
%lambda = 0;		    % Effect of vapor flow on liquid flow ("K2-effect")

% Inputs and disturbances
LT  = u(2*NT+3);                       % Reflux
VB  = u(2*NT+4);                       % Boilup
D   = u(2*NT+6);                       % Distillate
B   = u(2*NT+7);                       % Bottoms
F   = u(2*NT+5);                       % Feedrate
F_0 = dist.F_0;                     
zF  = dist.zF;                         % Feed composition         
qF  = dist.qF;                         % Feed liquid fraction



% THE MODEL

% Vapor-liquid equilibria
for i=1:NT-1
   y{i}  = SX.sym(['y_' num2str(i)]);
   V{i}  = SX.sym(['V_' num2str(i)]);
   L{i}  = SX.sym(['L_' num2str(i)]);
   dMdt{i}  = SX.sym(['dMdt_' num2str(i)]);
   dMxdt{i} = SX.sym(['dMxdt_' num2str(i)]);
end
L{NT}     = SX.sym(['L_' num2str(NT)]);
dMdt{NT}  = SX.sym(['dMdt_' num2str(NT)]);
dMxdt{NT} = SX.sym(['dMxdt_' num2str(NT)]);
dMdt{NT+1}  = SX.sym(['dMdt_' num2str(NT+1)]);
dMxdt{NT+1} = SX.sym(['dMxdt_' num2str(NT+1)]);

y  = vertcat(y{:});
for i=1:NT-1
    y(i)=alpha*u(i)/(1+(alpha-1)*u(i));
end

% Vapor Flows assuming constant molar flows
V = vertcat(V{:});
for i=1:NT-1
    if i >= NF
        V(i) = VB + (1-qF)*F;
    else
        V(i) = VB;
    end
end

L = vertcat(L{:});
L(NT)=LT;
for i=2:NT-1
    if i <=NF
        L(i) = L0b + (u(NT+1+i)-Muw)./taul;
    else
        L(i) = L0  + (u(NT+1+i)-Muw)./taul;
    end
end


% Time derivatives from  material balances for 
% 1) total holdup and 2) component holdup

% Column
dMdt  = vertcat(dMdt{:});
dMxdt = vertcat(dMxdt{:});
for i=2:NT-1
    dMdt(i) = L(i+1)         - L(i)       + V(i-1)         - V(i);
    dMxdt(i)= L(i+1).*u(i+1,1) - L(i).*u(i,1) + V(i-1).*y(i-1) - V(i).*y(i);
end

% Correction for feed at the feed stage
% The feed is assumed to be mixed into the feed stage
dMdt(NF) = dMdt(NF)  + F;
dMxdt(NF)= dMxdt(NF) + F*u(NT+1);

% Reboiler (assumed to be an equilibrium stage)
dMdt(1) = L(2)      - V(1)      - B;
dMxdt(1)= L(2)*u(2) - V(1)*y(1) - B*u(1);

% Total condenser (no equilibrium stage)
dMdt(NT) = V(NT-1)         - LT       - D;
dMxdt(NT)= V(NT-1)*y(NT-1) - LT*u(NT) - D*u(NT);

% Compute the derivative for the mole fractions from d(Mx) = x dM + M dx
for i=1:(2*NT+2)
    ceq{i}  = SX.sym(['ceq_' num2str(i)]);
end


% CSTR model
k1          = 34.1/60.0;
dMdt(NT+1)  = F_0 + D - F;
dMxdt(NT+1) = F_0*zF + D*u(NT) - F*u(NT+1) - k1*u(2*NT+2)*u(NT+1);

ceq  = vertcat(ceq{:});
for i=1:NT+1
    ceq(i) = dMxdt(i);
end

for i=1:NT+1
    ceq(NT+1+i) = dMdt(i);
end

% % bound constraints
% lb_u = [0.1; 0.1; 0.1; 0.1; 0.1];
% ub_u = [10; 4.008; 10; 1.0; 1.0];
% 
% % State bounds and initial guess
% x_min     = zeros(84,1);
% x_max     = ones(84,1);
% xB_max    = 0.1;
% x_max(1)  = xB_max;
% x_min(84) = 0.3;
% x_max(84) = 0.7;
% lbx  = [x_min;lb_u];
% ubx  = [x_max;ub_u];
% lbg  = zeros(2*NT+2,1);
% ubg  = zeros(2*NT+2,1);

end
