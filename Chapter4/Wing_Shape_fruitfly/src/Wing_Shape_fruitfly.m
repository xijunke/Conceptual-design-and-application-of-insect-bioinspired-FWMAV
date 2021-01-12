% function wing_para_output=wing_shape_fruitfly_sixteen_good_2()
% clear all;clc;
x =[3.1267,1.1578,0.7896,0.1624]; 
%%%%%%%%%%%%%%%%%%%%%%%%%
R_wing=x(1);
C_aver=x(2);
xr0=x(3);               % x-root offset  \mm
C_maxyaxis=x(4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
R_wingeff=3.004;    
C_avereff=0.8854;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (1) 
xr=0.3289;                     % x-root offset  \mm
xr_nd=xr/R_wingeff;      % x-root offset  
% yr=0;                          % y-root offset  \mm
% yr_nd=yr/C_avereff;   % y-root offset  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r2_nd=0.5801;    
r1_nd=1.106*r2_nd^1.366;  
% r3_nd=0.9*r1_nd^0.581; %1984-JEB-Ellingdon;  
% r3_nd3=r3_nd^3;                   % r3_nd3=r3_nd^3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R_wingeff1=3.007;            
% R_wingeff1=R_wingeff;         
xr_nd1=xr/R_wingeff1;      
F_nd1=r2_nd^2+2*xr_nd1*r1_nd+xr_nd1^2;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (2) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R_proximal=xr;                                                 
R_distal=R_wingeff+xr;                       
x=linspace(R_proximal,R_distal,200);
x_mod_Root=0.636;                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C_maxy=C_maxyaxis;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%
yr_lead=-0.08249*x.^6+0.9167*x.^5-4.04*x.^4+8.872*x.^3-10.06*x.^2+5.674*x-0.413-x_mod_Root;  
yr_trail=-0.0333*x.^6+0.504*x.^5-2.795*x.^4+7.258*x.^3-8.769*x.^2+3.739*x+0.1282-x_mod_Root;
C_rx=yr_lead-yr_trail;    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wing_aera=trapz(x,C_rx);       
C_aver0=wing_aera/R_wingeff; 
%%%%%%%%%%%%%%%%%%%%%%%%%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

%% 
% (a) 
r_nd=(x-xr)/R_wingeff;
yr_leadnd0=yr_lead/C_avereff;
P_coeff_lead=polyfit(r_nd,yr_leadnd0,6);
% (b)  
yr_trailnd0=yr_trail/C_avereff;
P_coeff_trail=polyfit(r_nd,yr_trailnd0,6);
% (c)  
Cr_nd=yr_leadnd0-yr_trailnd0;
P_coeff_Cr=polyfit(r_nd,Cr_nd,6);  
% cftool 
syms r_nd   % 
yr_leadnd=vpa(poly2sym(P_coeff_lead,r_nd),6); 
yr_trailnd=vpa(poly2sym(P_coeff_trail,r_nd),6);
Cr_nd=vpa(poly2sym(P_coeff_Cr,r_nd),6);  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
C_nd=Cr_nd;
R2nd2=double(int(r_nd^2*C_nd,r_nd,0,1)); 
R1nd1=double(int(r_nd*C_nd,r_nd,0,1));   
S_nd=double(int(C_nd,r_nd,0,1));     
disp(['二阶面积矩的回转半径的平方: r2_2nd=' num2str(R2nd2)  ' 量纲单位是mm^4'])
disp(['二阶面积矩的回转半径: r_2nd=' num2str(sqrt(R2nd2))  ' 量纲单位是mm^3'])
disp(['一阶面积矩的回转半径: r_1nd=' num2str(R1nd1)  ' 量纲单位是mm^3'])
disp(['无量纲翅面积Swing_nd=' num2str(S_nd)  ' 量纲单位是mm^2'])
fx1=(r_nd+xr_nd)^2*C_nd;   
% fx2=vpa(fx1,5)
fx3=expand(fx1);
F_ndTrans=double(int(fx3,r_nd,0,1));            
disp(['无量纲气动力F_ndTrans=' num2str(F_ndTrans)  ' 量纲单位是mm^4'])
F_nd2=R2nd2+2*xr_nd*R1nd1+xr_nd^2;   
disp(['无量纲气动力F_ndTrans=' num2str(F_nd2)  ' 量纲单位是mm^4'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                
xr_nd1=xr0/R_wing;   
F_ndTrans=R2nd2+2*xr_nd1*R1nd1+xr_nd1^2; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
f_x_lead=C_aver*yr_leadnd;
f_x_trail =C_aver*yr_trailnd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_x_lead1=inline(vectorize(f_x_lead),'r_nd');  
f_x_trail1=inline(vectorize(f_x_trail),'r_nd');   
r_nd1=linspace(0,1,200);   
f_x_lead2=f_x_lead1(r_nd1); 
f_x_trail2=f_x_trail1(r_nd1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
C_lead_ymax=max(f_x_lead2);     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C_trail_ymin=min(f_x_trail2);      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C_max_LtoT=C_lead_ymax-C_trail_ymin;   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rx=r_nd1*R_wing+xr0;  
r_x_max=max(rx);     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_x=linspace(xr0,r_x_max,200);
yr_leadnd2=-2.307e-009*r_x.^6+7.065e-007*r_x.^5-8.38e-005*r_x.^4+0.004742*r_x.^3-0.1288*r_x.^2+1.698 *r_x-0.5166;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
yr_trailnd2= -9.314e-010*r_x.^6+4.099e-007*r_x.^5-6.329e-005*r_x.^4+0.004316*r_x.^3 -0.1158*r_x.^2+0.02573*r_x+0.034;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)  
hold on;
plot(rx,f_x_lead2,'k-',rx,f_x_trail2,'k-','LineWidth',3);  hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
quiver(-0.25,0,...                                          
            4.25,0,1,'k-','LineWidth',2);   hold on;             
quiver(0,-1.5,...                                                
           0,2.5,1,'k-','LineWidth',2);   hold on;                
quiver(rx(1,1), f_x_lead2(1,1),...                                      
           rx(1,length(rx))-0.5,0,1,'r-.','LineWidth',2);   hold on;    
quiver(rx(1,1), f_x_lead2(1,1)-1,...                                            
           0,2*abs(f_x_lead2(1,1))+1.75,1,'r-.','LineWidth',2);   hold on;
quiver(0, C_maxyaxis*C_max_LtoT,...                  
           rx(1,length(rx))+xr0-0.5,0,1,'r-.','LineWidth',2);   hold on;   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
% plot(r_x,yr_leadnd2,'-.r',r_x,yr_trailnd2,'-.r','LineWidth',3); hold on
xlabel('右侧翅膀的展向(y坐标)')
ylabel('右侧翅膀的弦向(x坐标)')
title('采用前后缘拟合函数绘制右侧翅膀形貌, 扭转轴的位置和压心随时间和空间的分布')
grid off
axis([-0.3,4.75,-1.25,1.05])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
XC=1.920243385-xr; 
YC=-0.149785466;    
R_ratio=R_wing/R_wingeff;    
C_ratio=C_aver/C_avereff;    
XC_tran=R_ratio*XC+xr0;   
YC_tran=C_ratio*YC;          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yr_leadnd2=inline('-2.307e-009*r_x.^6+7.065e-007*r_x.^5-8.38e-005*r_x.^4+0.004742*r_x.^3-0.1288*r_x.^2+1.698 *r_x-0.5166','r_x');
yr_trailnd2= inline('-9.314e-010*r_x.^6+4.099e-007*r_x.^5-6.329e-005*r_x.^4+0.004316*r_x.^3 -0.1158*r_x.^2+0.02573*r_x+0.034','r_x');
Area=quadl(yr_leadnd2,0.3289,100.3289)-quadl(yr_trailnd2,0.3289,100.3289);
integrand_xc=inline('r_x.*((-2.307e-009*r_x.^6+7.065e-007*r_x.^5-8.38e-005*r_x.^4+0.004742*r_x.^3-0.1288*r_x.^2+1.698 *r_x-0.5166)-(-9.314e-010*r_x.^6+4.099e-007*r_x.^5-6.329e-005*r_x.^4+0.004316*r_x.^3 -0.1158*r_x.^2+0.02573*r_x+0.034))','r_x'); 
integrand_zc=inline('(-2.307e-009*r_x.^6+7.065e-007*r_x.^5-8.38e-005*r_x.^4+0.004742*r_x.^3-0.1288*r_x.^2+1.698 *r_x-0.5166).^2-(-9.314e-010*r_x.^6+4.099e-007*r_x.^5-6.329e-005*r_x.^4+0.004316*r_x.^3 -0.1158*r_x.^2+0.02573*r_x+0.034).^2','r_x'); 
X_C=quadl(integrand_xc,xr0,r_x_max)/Area+xr0;                
Z_C=quadl(integrand_zc,xr0,r_x_max)/(2*Area);         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot(X_C,Z_C,'rd',XC_tran,YC_tran,'dk','LineWidth',4)
plot(XC_tran,YC_tran,'dk','LineWidth',4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R_wingeff=R_wing;
C_avereff=C_aver;
%%
% (1) 
Rou=1.225*10^(-3);         
Coeff_liftdragF_N=(1/2)*Rou*C_avereff*R_wingeff^3*F_ndTrans;  
% (2) 
M_xaercoeff=(1/2)*Rou*C_avereff^2*R_wingeff^3*F_ndTrans;   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (3) 
fx4=(r_nd+xr_nd)^3*C_nd;       
fr_nd5=expand(fx4);
I1=double(int(fr_nd5,0,1));         
I1y=(1/2)*Rou*C_avereff*R_wingeff^4*I1;      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
C_maxy=0.138924474377504; 
C_maxynd=C_maxy/C_avereff;
yr_leadnd=yr_leadnd-C_maxynd;
yr_trailnd=yr_trailnd-C_maxynd;
%%%%%%%%%%%%%%%%%%%
y0=yr_leadnd-C_nd;   
y1=yr_leadnd;
yr_nd2=(y1^4+y0^4)/4;   
Z_rnd=double(int(yr_nd2,r_nd,0,1));           
disp(['有效力臂的无量纲位置Z_nd=' num2str(Z_rnd)  ' 量纲单位是mm'])  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M_xrdcoeff=(1/2)*Rou*C_avereff^4*R_wingeff*Z_rnd; 
fx16=(r_nd+xr_nd)*C_nd^3;    
fr_nd17=expand(fx16);
I8=double(int(fr_nd17,0,1));       
I8z=I8;   
X_rnd=I8z;
M_zrdcoeff=(1/6)*Rou*C_avereff^3*R_wingeff^2*X_rnd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
fx6=(r_nd+xr_nd)*C_nd^2;                                  
fr_nd7=expand(fx6);
F_ndRot=double(int(fr_nd7,0,1));     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F_yrotcoeff=(1/2)*Rou*C_avereff^2*R_wingeff^2*F_ndRot;  
M_xRotcoeff=(1/2)*Rou*C_avereff^3*R_wingeff^2*F_ndRot;   
fx8=(r_nd+xr_nd)^2*C_nd^2;             
fr_nd9=expand(fx8);
I2=double(int(fr_nd9,0,1));                             
I2y=(1/2)*Rou*C_avereff^2*R_wingeff^3*I2;   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
yr_hnd=C_nd/2-yr_leadnd+C_maxynd; 
%%%%%%%%%%%%%%%%%%%%%%%%
fx10=(r_nd+xr_nd)*C_nd^2*yr_hnd;
fr_nd11=expand(fx10);
I_xzamnd=int(fr_nd11,0,1);
I3=double(I_xzamnd);                    
I_xzam=pi*Rou*C_avereff^3*R_wingeff^2*I3/4; 
fx12=C_nd^2*(yr_hnd^2+C_nd^2/32);
fr_nd13=expand(fx12);
I_xxamnd=int(fr_nd13,0,1);
I4=double(I_xxamnd);                                      
I_xxam=pi*Rou*C_avereff^4*R_wingeff*I4/4;   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I5z=pi*Rou*C_avereff^2*R_wingeff^2*F_ndRot/4;    
I_xzam2=I_xzam+I5z*C_maxynd^2;       
fx14=C_nd^2*yr_hnd;
fr_nd15=expand(fx14);
I6=double(int(fr_nd15,0,1));                            
I6z=pi*Rou*C_avereff^3*R_wingeff*I6/4;         
I_xxam2=I_xxam+I6z*C_maxynd^2;             
 I7y=pi*Rou*C_avereff^2*R_wingeff^3*I2/4;   
%% 
c_zcopnd_tr=I1/F_ndTrans;                         
c_zcopnd_rot=I2/F_ndRot;                          
c_zcopnd_addaver=-0.3359; 
%% 
r_xcopnd_tr=I1/F_ndTrans;   
r_xcopnd_rot=I2/F_ndRot;    
r_xcopnd_addaver=0.7934;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%