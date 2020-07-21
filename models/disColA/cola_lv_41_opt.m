function xprime=cola_lv_41_opt(X) 
% sample usage:   [t,x]=ode15s('cola_lv',[0 5000],0.5*ones(1,82));
%
% cola_lv - Subroutine for simulation with LV-configuration.
%           It calls the model colamod, and 
%           includes control of condenser and reboiler level 
%           using two P-controllers with the LV-configuration. 
%
%            Inputs are reflux (LT) and boilup (VB). Disturbances
%            are feedrate and feed composition. These are set by directly
%            altering 'cola_lv.m'. Outputs are liquid composition and
%            liquid hold up for stages 1 through NT, given in x. 

% Number of stages in the column
global u NT %Make perturbed inputs/disturbances available to model

% Inputs and disturbances
% LT=2.70629;                          % Reflux
% VB=3.20629;                          % Boilup
% F=1.0 + 0.00;                        % Feedrate
% zF=0.5;                              % Feed composition
% qF=1.0;                              % Feed liquid fraction

% P-Controllers for control of reboiler and condenser hold up.
% KcB=10;  KcD=10;         % controller gains
% MDs=0.5; MBs=0.5;        % Nominal holdups - these are rather small  
% Ds=0.5; Bs=0.5;          % Nominal flows
% MB=X(NT+1);  MD=X(2*NT); % Actual reboiler and condenser holdup
% D=Ds+(MD-MDs)*KcD;       % Distillate flow
% B=Bs+(MB-MBs)*KcB;       % Bottoms flow     
% 
% fprintf('-----------------\n');
% fprintf('D value: %f\n',D);
% fprintf('B value: %f\n',B);

D = 0.5;
B = 0.5;

% Store all inputs and disturbances
% Store all inputs and disturbances
u_all(1:2) = u(1:2) ;
u_all(3) = D; 
u_all(4) = B;
u_all(5:7) = u(3:5) ;



t= [];

xprime=colamod_41(t,X,u_all);

