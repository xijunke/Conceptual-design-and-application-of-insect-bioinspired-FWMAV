%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% lift and drag coefficients with Alpha
% lift_drag_coefficients_with_alpha
% clear all; clc;
% alpha=atan2(-w_y,w_z);    %Aerodynamic angle of attack  
beta=pi/4; f=0.25; %\Hz
delta=pi/4;    
t=linspace(0,10,200);
alpha=(beta*sin(2*pi*f*t+delta))*180/pi;
% figure(1)
% plot(t,alpha,'r.')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
options=optimset('Display','off');
t1_0=fzero(inline('45*sin(pi/4 + (pi*t)/2)','t'),[3,4],options); 
t2_0=fzero(inline('45*sin(pi/4 + (pi*t)/2)','t'),[5,6],options);
% Result:   t1_0=3.5000;   t2_0=5.5000;
indx1=find(t>3.45 & t<=3.5);  
indx2=find(t>5.45 & t<=5.5);
% Result:   indx1=70;   indx2=110;
alpha_2=alpha(1,indx1:indx2);      
% figure(2)
% plot(t(1,indx1:indx2),alpha_2,'ro')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
[alpha_2_max, locat_max]=max(alpha_2);   
m=length(alpha_2);                                                
alpha_3=alpha_2(1,locat_max:m)+45;  
[alpha_new,index]=sort(alpha_3);       
alpha_4=[0, alpha_2(1,2:locat_max-1), alpha_new];  
% figure(3)
% plot(t(1,indx1:indx2),alpha_4,'m-d')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
%% (1) 1999-Science-MH Dickinson
% C_rot_theo=pi*(0.75-x_0nd(r));　  % 转动环量气动力系数
C_L_dickinson=0.225+1.58*sin((2.13*alpha_4-7.2)*pi/180);  
C_D_dickinson=1.92-1.55*cos((2.04*alpha_4-9.82)*pi/180); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (8) 
% C_L_polhamus=K_p*sin(alpha)*(cos(alpha)).^2+K_v*(cos(alpha))*sin(alpha).^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (9) 2014-JRSI-Nabawy
C_La2d=5.15662;  
lambda=0.70; 
x=[3.004,0.8854,0.3289,0.356737];
% x =[3.6837,1.7620,2.0000,0.0055];
R_wing=x(1);
C_aver=x(2);
xr0=x(3);
C_maxyaxis=x(4);
C_periEdge=C_perimeter_to_Edge_correction(R_wing,C_aver,xr0,C_maxyaxis);
E=lambda*C_periEdge;   
k=1.51;  
% for the span of a single  wing 
AR=3.40158;   
% AR =3.7643;
C_La_Ny=C_La2d/(E+k*C_La2d/(pi*AR));    
C_L_nabawy=0.5*C_La_Ny*sin(2*alpha_4*pi/180);   
C_D_nabawy=C_La_Ny*(sin(alpha_4*pi/180)).^2;  
C_L_polhamus1=(0.5*C_La_Ny*sin(2*alpha_4*pi/180)).*(cos(alpha_4*pi/180)+(1-k*C_La_Ny/(pi*AR)).*sin(alpha_4*pi/180)); 
% figure(11)
% plot(alpha_4,C_L_dickinson,'r*',alpha_4,C_D_dickinson,'b*',alpha_4,C_L_nabawy,'r-',alpha_4,C_D_nabawy,'b-',alpha_4,C_L_polhamus1,'g-')
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (10) 
a_0=2*pi;
AR=2.8;  
C_La_Ta=pi*AR/(1+sqrt((pi*AR/a_0)^2+1));    
% C_L_taha=(0.5*pi*AR/(1+sqrt((pi*AR/a_0)^2+1))).*sin(2*alpha_4*pi/180);
C_L_taha=0.5*C_La_Ta*sin(2*alpha_4*pi/180);
C_D_taha=C_L_taha.*tan(alpha_4*pi/180);
C_L_polhamus2=(0.5*C_La_Ta*sin(2*alpha_4*pi/180)).*(cos(alpha_4*pi/180)+(1-C_La_Ta/(pi*AR)).*sin(alpha_4*pi/180));
figure(12)
C_LD=plot(alpha_4,C_L_dickinson,'r*',alpha_4,C_D_dickinson,'b*',alpha_4,C_L_taha,'r-',alpha_4,C_D_taha,'b',...
       alpha_4,C_L_nabawy,'k-',alpha_4,C_D_nabawy,'k-',alpha_4,C_L_polhamus1,'c-',alpha_4,C_L_polhamus2,'m-');
xlabel('\it\alpha (deg.)')
ylabel('\it C_L (\alpha ) and \it C_D (\alpha )')
title('Aerodynamic coefficients of lift \itvs. \alpha \rm for flapping wing')
legend('\itC_{L,dickinson}','\itC_{D,dickinson}','\itC_{L,taha}','\itC_{D,taha}','\itC_{L,nabawy}',...
            '\itC_{D,nabawy}','\itC_{L,polhamus1}','\itC_{L,polhamus2}')
set(C_LD,'LineWidth',2)
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
figure(13)    %  subplot(224)
C_L=plot(alpha_4,C_L_dickinson,'r*',alpha_4,C_L_nabawy,'k-',alpha_4,C_L_taha,'r:',...
                alpha_4,C_L_polhamus1,'c-',alpha_4,C_L_polhamus2,'m-');
xlabel('\it\alpha (deg.)')
ylabel('\it C_L (\alpha )')
title('Aerodynamic coefficients of lift \itvs. \alpha \rm for flapping wing')
legend('\itC_{L,dickinson}','\itC_{L,nabawy}','\itC_{L,polhamus1}','\itC_{L,taha}','\itC_{L,polhamus2}')
set(C_L,'LineWidth',2) 
grid on
%%
figure(14)
C_D=plot(alpha_4,C_D_dickinson,'b*',alpha_4,C_D_nabawy,'k-',alpha_4,C_D_taha,'b:');
xlabel('\it\alpha (deg.)')
ylabel('\it C_D (\alpha )')
title('Aerodynamic coefficients of drag \itvs. \alpha \rm for flapping wing')
legend('\itC_{D,dickinson}','\itC_{D,nabawy}','\itC_{D,taha}')
set(C_D,'LineWidth',2)   
grid on
%% 
v_axis=axis;     %axis([xmin,xmax,ymin,ymax])
v_axis(3)=-0.5;      
v_axis(4)=3.5;       
axis(v_axis);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%