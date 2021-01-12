% Design_manufacturing_rules_PZT_actuator
% The maximum strain (?) and stress (��) at the beam's surface 
% are calculated from the applied deflection (��) and measured force (F) using the standard equations:
% where L_0 and L_i are support span and load span
% epsilon=delta*(6*d/((L_0-L_i)*(L_0-2*L_i))); %  Strain
% sigama=F*1.5*((L_0-L_i)/(b*d^2));          % Stress
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;clc;
%% 
t_pzt=135e-6;   % um % thickness of each PZT plate in the actuator  % t_pzt is the thickness of the PZT
t_cf=50e-6;       %um % the thickness of the actuator's central carbon fiber layer % t_cf is the thickness the carbon fiber
E_cf=340e9;    % Gpa % modulus of the carbon fiber  % E_cf is the modulus of the carbon fiber
L_pzt=9e-3;      % mm % the length of the active portion of the actuator
L_r=0.25;     % the ratio of the length of the actuator's rigid extension to the length of the PZT   % L_r is the length ratio (L_ext/L_act)
w_nom=1.125e-3;   % the mean width of the actuator  % w_nom is the average width of the PZT 
w_r=1.556;   % w_r is the width ratio (the width of the PZT's base divided by w_nom)
m_actuator=40e-6;  % mg����>kg
%%%%%%%%%%%%%
% t_pzt=127e-6;
% t_cf=40e-6; 
% E_cf=350e9;
% L_pzt=6.021e-3; 
% L_r=5.979e-3/L_pzt;       % L_ext=5.979e-3;
% w_nom=1.1755e-3;        % w_nom=(w0+w1)/2; % w0=1.569e-3; w1=0.782e-3;
% w_r=1.569e-3/w_nom;   % w_r =1.3348;
% m_actuator=20e-6;  % mg����>kg
%%%%%%%%%%%%%%%%%%%%%%%%%
L_act=L_pzt;       % L_act is the length of the PZT
L_ext=L_r*L_act;  % L_ext is the length of the alumina extension 
d_31=320e-12;          % pm/V  % d_31 is the piezo coefficient
E_pzt=62e9;  % pa  % E_pzt is the modulus of the PZT
%% 
xi_0=0.4e-6;            % V*um^-1;  
U=300;               % v
xi=U/t_pzt;               % xi is the applied electric field % v
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (1) old
% 
% where ��_pp is the peak-to-peak deflection at the actuator tip
delta_pp=0.5*d_31*E_pzt*t_pzt*(t_pzt+t_cf)*L_act^2*xi*(1+2*(L_ext/L_act))/...
                  ((1/3)*E_pzt*t_pzt*(1.5*t_cf^2+3*t_cf*t_pzt+2*t_pzt^2)+(E_cf*t_cf^3)/12)
%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% F_bexp=F_bmeas/(1-K_b/K_f)=150-170; % mN  % the expected blocked force 
%  % G=��;               % a gain factor (defined in [9]) determined by the shape of the actuator
GF=8*(1-w_r)^3*(1+2*(L_ext/L_act))/(-6*(w_r-1)*(-3+4*L_r*(w_r-1)+2*w_r)+3*(-2+2*L_r*(w_r-1)+w_r)^2*abs(log((2-w_r)/w_r))); 
% F_b is the blocked force (one-way)   
F_b=-0.75*d_31*E_pzt*w_nom*t_pzt*(t_pzt+t_cf)*xi*GF/L_act
F_b_pp=2*F_b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
ED_m=0.5*F_b_pp*delta_pp/m_actuator

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (2) new
%% f_31 
% f_31 is defined as the internal stress produced per unit electric field:  f_31=sigama_1,int/xi_3;
b=-230;                  % recall that we have defined compressive strain to be negative
c=10^(-5);
d=69*10^(-9);       % m*V^-1;
f_31min=14;           % pa*m*V^-1;  
f_31max=29;           % pa*m*V^-1;
% f_31=(1+b*epsilon)*(f_31min+(f_31max*(1-d*xi)-f_31min)*((exp(c*(xi-xi_0)))/(1+exp(c*(xi-xi_0)))));
% f_31_zerostrain=(1+b*epsilon)*(f_31min+(f_31max*(1-d*xi)-f_31min)*((exp(c*(xi-xi_0)))/(1+exp(c*(xi-xi_0))))); % epsilon=0
f_31_zerostrain=f_31min+(f_31max*(1-d*xi)-f_31min)*((exp(c*(xi-xi_0)))/(1+exp(c*(xi-xi_0)))); % epsilon=0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
% Blocked force
% F_bnew=0.75*(f_31min+(f_31max*(1-d*xi)-f_31min)*(exp(c*(xi-xi_0)))/(1+exp(c*(xi-xi_0))))*w_nom*t_pzt*(t_pzt+t_cf)*xi*GF/L_act;
%%%%%%%%%%%%%%%%%%
F_bnew=-0.75*f_31_zerostrain*w_nom*t_pzt*(t_pzt+t_cf)*xi*GF/L_act
F_bnew_pp=2*F_bnew
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Strain
epsilon_0=-0.00047;  % Strain at which the transition occurs
epsilon=-delta_pp*((t_pzt+t_cf)/(2*(1+2*L_ext/L_act)*L_act^2))  %  Strain
%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stress
% epsilon=��;         %  Strain
E_min=38.5*10^9;  % pa;  % the minimum moduli of the PZT under varying strain
E_max=81*10^9;    % pa;  % the maximum moduli of the PZT under varying strain
a=8000;                 % Steepness of the transition
sigama=epsilon*E_min-((E_max-E_min)/a)*(log(1+exp(a*(epsilon_0-epsilon))-log(1+exp(a*epsilon_0))))
%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
E=E_min+(E_max-E_min)*(exp(a*(epsilon_0-epsilon))/(1+exp(a*(epsilon_0-epsilon))));
%%%%%%%%%%%%%%%%%%%%%%%%%
% we can approximate the modulus of the inactive plate as E_min independent of strain.
E_pzt1=E_min; % this is the PZT plate which is inactive and in tension.
E_pzt2=E_min-((E_max-E_min)/(a*epsilon))*(log(1+exp(a*(epsilon_0-epsilon)))-log(1+exp(a*epsilon_0)));
%%%%%%%%%%%%%%%%%%%%%%%%%
% the free deflection
% delta_ppnew=0.5*((1+b*epsilon)*(f_31min+(f_31max*(1-d*xi)-f_31min)*(e^(c*(xi-xi_0))/(1+e^(c*(xi-xi_0)))))*...
%                        t_pzt*(t_pzt+t_cf)*L_act^2*xi)*(1+2*L_ext/L_act)/...
%                        ((1/3)*(E_pzt1+E_pzt2)/2*t_pzt*(1.5*t_cf^2+3*t_cf*t_pzt+2*t_pzt^2)+(E_cf*t_cf^3/12));
%%%%%%%%%%%%%%%%%%
% delta_ppnew=0.5*f_31*t_pzt*(t_pzt+t_cf)*L_act^2*xi*(1+2*L_ext/L_act)/...
%                        ((1/3)*(E_pzt1+E_pzt2)/2*t_pzt*(1.5*t_cf^2+3*t_cf*t_pzt+2*t_pzt^2)+(E_cf*t_cf^3/12));
%%%%%%%%%%%%%%%%%%
delta_ppnew=0.5*((1+b*epsilon)*f_31_zerostrain*t_pzt*(t_pzt+t_cf)*L_act^2*xi)*(1+2*L_ext/L_act)/...
                       ((1/3)*(E_pzt1+E_pzt2)/2*t_pzt*(1.5*t_cf^2+3*t_cf*t_pzt+2*t_pzt^2)+(E_cf*t_cf^3/12))                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
ED_m_new=0.5*F_bnew_pp*delta_ppnew/m_actuator
































