function cyclicSsOpt
%CYCLICSSOPT Summary of this function goes here
%
% Steady-state optimization for cyclic parallel reactor.
% This process is nondissipative.
% parallel reaction:
% R --> P1  (desired product)
% R --> P2  (waste product)
%
% [OUTPUTARGS] = CYCLICSSOPT(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2018/01/17 18:58:53 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2018

format long;
import casadi.*

% parameter values
%delta   = 0.55;      % dimensionless
%delta  = 1.1;
delta   = 2.2;

% symbolic primitives
u  = SX.sym('u');     % u  - heat flux through cooling jacket
x1 = SX.sym('x1');    % x1 - concentration of R
x2 = SX.sym('x2');    % x2 - concentration of P1
x3 = SX.sym('x3');    % x3 - concentration of P2 

% concatenate states and controls 
%x   = [x1;x2;x3;u];   % all are dimensionless
x   = [u;x1;x2;x3];

% initial guess
%Uinit = [0.1;0.1;0.1;0.05];
Uinit = [0.1;0.1;0.1;0.1];

% define the dynamics as equality constraints and additional inequality constraints
alpha = 4;
[obj,eq, lbx, ubx, lbg, ubg] = buildModelEq(x,delta,alpha);

prob    = struct('f', obj, 'x', x, 'g', eq);
options = struct;
options.ipopt.tol                   = 1e-12;
% options.acceptable_compl_inf_tol    = 1e-12;
solver  = nlpsol('solver', 'ipopt', prob, options);

% Solve the NLP
startnlp = tic;
sol   = solver('x0', Uinit, 'lbx', lbx, 'ubx', ubx, 'lbg', lbg, 'ubg', ubg);
elapsednlp = toc(startnlp);
fprintf('IPOPT solver runtime = %f\n',elapsednlp);

u      = full(sol.x);
lamda  = full(sol.lam_g);
%Xinit  = u;

%% Compute Hessian and perform Greshgorin convexification

% symbolic variable for dual variables
l1 = SX.sym('l1');  
l2 = SX.sym('l2');  
l3 = SX.sym('l3');  

xsol = u;
lambda.eqnonlin = lamda;
% extract Lagrange multiplier
%l   = vertcat(l{:});
l  = [l1;l2;l3];
L  = obj + l'*eq;

Lagr = Function('Lagr', {x,l}, {L}, char('x','l'), char('Lagr'));
%Jac  = Function(Lagr.jacobian('x','Lagr'));
H    = Function(Lagr.hessian('x','Lagr'));
cons = Function('Const', {x}, {eq}, char('x'), char('cons'));
Jcon = Function(cons.jacobian('x','cons'));

eqVal = cons(xsol);
eqVal = full(eqVal);
Hx   = H(xsol,lambda.eqnonlin);
Hx   = full(Hx);
Jac  = Jcon(xsol);
Jac  = full(Jac);
rH   = null(Jac)'*Hx*null(Jac)
erH  = eig(rH)


[Hxxl,Qmax]   = Greshgorin(Hx);

keyboard;


end

function [J, ceq, lbx, ubx, lbg, ubg] = buildModelEq(u,delta,alpha)
import casadi.* 

% x1 = u(1);
% x2 = u(2);
% x3 = u(3);
% u1 = u(4);

u1 = u(1);
x1 = u(2);
x2 = u(3);
x3 = u(4);

% objective function 
J = -x2;
%J = -x2 + alpha*(u1 - 0.149096766874856);


% CSTR model
% define xdot
x1dot  = SX.sym('x1dot');
x2dot  = SX.sym('x2dot');
x3dot  = SX.sym('x3dot');

x1dot  = 1 - 1e4*x1^2*exp(-1/x3) - 400*x1*exp(-delta/x3) - x1;
x2dot  = 1e4*x1^2*exp(-1/x3) - x2;
x3dot  = u1 - x3;
ceq    = [x1dot;x2dot;x3dot];

% bound constraints
% lbx  = [0;0;0;0.049];
% ubx  = [1;1;1;0.449];
lbx  = [0.049;0;0;0];
ubx  = [0.449;1;1;1];
lbg  = [0;0;0];
ubg  = [0;0;0];

end

function [H,Q] = Greshgorin(H)
numH    = size(H,1);
Q       = zeros(numH,numH);
%Q       = eye(numH);
%delta   = 1e-1;
%delta   = 1e-2;
%delta   = 2.5e-1;  % normal case
%delta   = 0.5;
%delta   = 1;
%delta   = 1.5;
%delta   = 2;
%delta   = 2.5;      % with measurement noise 1 percent
%delta   = 5;
%delta   = 100;
delta   = 1e-6;
for i=1:numH  % iterate all row of Hessian
    sumRow = 0;
    for j=1:numH
        if j ~= i
            sumRow = sumRow + abs(H(i,j));
        end
    end
    %sumRow = sumRow - abs(H(i,i));
    
    if H(i,i) <= sumRow   % include equality 
    %if H(i,i) < sumRow 
        %Q(i,i) = sumRow - H(i,i) + delta;
        Q(i,i) = sumRow - H(i,i) + delta;
    end
end

% % loop over Qm to obtain maximum number
% for i=1:numH
%     if Q(i,i) > Qm(i,1)
%         Qm(i,1) = Q(i,i);
%     end
% end
Q = diag(Q);
end
