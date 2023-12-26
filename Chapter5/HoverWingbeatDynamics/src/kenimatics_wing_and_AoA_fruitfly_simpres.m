function wing_m_output=kenimatics_wing_and_AoA_fruitfly_simpres()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
syms t        
w =1185.6; 
phi_pres =sym('1.149*sin(w*t+1.571)');  
dphi_pres =diff(phi_pres ,t,1);
ddphi_pres =diff(phi_pres ,t,2); 
dphi=inline(vectorize(dphi_pres ),'w','t');          
ddphi=inline(vectorize(ddphi_pres ),'w','t');  
%% 
psi_sim=sym('-(0.005045+8.152*cos(t*w)+65.33*sin(t*w)+0.01177*cos(2*t*w)-0.006058*sin(2*t*w)-8.152*cos(3*t*w)+10.67*sin(3*t*w)+0.03641*cos(4*t*w)-0.008574*sin(4*t*w)-1.666*cos(5*t*w)-0.8538*sin(5*t*w)+0.0234*cos(6*t*w)-0.003305*sin(6*t*w)-0.1157*cos(7*t*w)+0.1072*sin(7*t*w)+0.02172*cos(8*t*w)-0.001624*sin(8*t*w))*pi/180');  %·ûºÅº¯Êý
dpsi_sim=diff(psi_sim,t,1);   
ddpsi_sim=diff(psi_sim,t,2);    
dpsi=inline(vectorize(dpsi_sim),'w','t');    
ddpsi=inline(vectorize(ddpsi_sim),'w','t');    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
f=188.7; T=1/f;  
t=linspace(0.0052824335,0.0052824335+3*T,1000);  % t_steady1   
dphi=dphi(w,t);    
ddphi=ddphi(w,t);    
dpsi=dpsi(w,t);        
ddpsi=ddpsi(w,t);     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              
phi_pres=1.149*sin(w*t+1.571);
psi_sim=-(0.005045+8.152*cos(t*w)+65.33*sin(t*w)+0.01177*cos(2*t*w)-0.006058*sin(2*t*w)-8.152*cos(3*t*w)+10.67*sin(3*t*w)+...
              0.03641*cos(4*t*w)-0.008574*sin(4*t*w)-1.666*cos(5*t*w)-0.8538*sin(5*t*w)+0.0234*cos(6*t*w)...
             -0.003305*sin(6*t*w)-0.1157*cos(7*t*w)+0.1072*sin(7*t*w)+0.02172*cos(8*t*w)-0.001624*sin(8*t*w))*pi/180;   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
alpha1=pi/2-psi_sim.*sign(dphi);        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   
omega_y=dphi.*sin(psi_sim);
omega_z=dphi.*cos(psi_sim);
v_y_nonr=omega_z; 
v_z_nonr=-omega_y;   
alpha2=atan2(v_y_nonr,-v_z_nonr);   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
alpha=alpha2;
C_L =(0.225+1.58*sin((2.13*abs(alpha)*180/pi-7.2)*pi/180)); 
C_D =(1.92-1.55*cos((2.04*abs(alpha)*180/pi-9.82)*pi/180));
C_N1=cos(abs(alpha)).*C_L+sin(abs(alpha)).*C_D; 
C_L2 =0.225+1.58*sin((2.13*abs(alpha)*180/pi-7.2)*pi/180);  
C_D2 =1.92-1.55*cos((2.04*abs(alpha)*180/pi-9.82)*pi/180); 
C_N2=sign(alpha2).*sqrt(C_L2.^2+C_D2.^2); 
C_N3=3.4*sin(abs(alpha));      
C_T=0.4*(cos(2*alpha)).^2.*(alpha>-pi/4 & alpha<pi/4);       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
wing_m_output=[t',phi_pres',psi_sim',alpha',dphi',dpsi',ddphi',ddpsi',C_L',C_D',C_N1',C_T'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

