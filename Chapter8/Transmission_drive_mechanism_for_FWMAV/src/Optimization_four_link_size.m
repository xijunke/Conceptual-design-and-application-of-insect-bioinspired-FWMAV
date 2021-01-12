%% Optimization_four_link_size
clear all;clc;
%%%%%%%%%%%%%%%%%%%%%
% L_c=300e-6;   % um！！to！！m:  e-6
% L_s1=300e-6;   % um！！to！！m:  e-6
% L_s2=400e-6;   % um！！to！！m:  e-6
% L_s3=700e-6;   % um！！to！！m:  e-6
%%%%%%%%%%%%%%%%%%%%%
% 歌深��2011-ICIRS-System identification, modeling, and optimization of ...
% an insect-sized flapping-wing micro air vehicle-Finio_IROS11_c-8p
%% 
L_c=300e-6;   % um！！to！！m:  e-6
L_s1=500e-6;   % um！！to！！m:  e-6
L_s2=300e-6;   % um！！to！！m:  e-6
L_s3=630e-6;   % um！！to！！m:  e-6
%%%%%%%%%%%%%%%%%%%%%
%% 
% L_c=312e-6;   % um！！to！！m:  e-6
% L_s1=400e-6;   % um！！to！！m:  e-6
% L_s2=291e-6;   % um！！to！！m:  e-6
% L_s3=498e-6;   % um！！to！！m:  e-6
% L_s3=L_s1;
%%%%%%%%%%%%%%%%%%%%%
x0= [L_c,L_s1,L_s2,L_s3];     
%%%%%%%%%%%%%%%%%%%%%
objfun=@objTdiff_linear;  
% objfun=@objTdiff_linear2;  
% objfun=@objTdiff_nonlinear;  %
%%%%%%%%%%%%%%%%%%%%%
% lb = [0e-6,0e-6,0e-6,0e-6];          % Set lower bounds 
lb = [50e-6,50e-6,50e-6,50e-6];          % Set lower bounds 
% lb = [100e-6,100e-6,100e-6,100e-6];          % Set lower bounds 
% lb = [200e-6,200e-6,200e-6,200e-6];         % Set lower bounds
% ub = [1000e-6, 1000e-6,300e-6,1000e-6];  % Set upper bounds
ub = [1000e-6, 1000e-6,300e-6,1000e-6];  % Set upper bounds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
options=optimset('Algorithm','active-set','Display','iter-detailed','MaxFunEvals',160000,'MaxIter',160000,'TolX',1e-6,'TolFun',1e-6);
[x,fval] =fmincon(objfun,x0,[],[],[],[],lb,ub,[],options);
optimal_four_link_size=x*10^6
%%%%%%%%%%%%%%%%%%%%%%%%%
optimal_four_link_size =[445.8109  564.1168  297.1740  439.1168]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


