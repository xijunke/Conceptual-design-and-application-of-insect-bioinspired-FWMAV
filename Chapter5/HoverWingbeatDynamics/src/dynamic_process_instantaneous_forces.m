%% The dynamic process of instantaneous forces for a complete stroke
clear all; clc;
wing_kenimatics=kenimatics_wing_and_AoA_fruitfly_simpres(); 
t=wing_kenimatics(:,1);                        %ms
phi_pres=wing_kenimatics(:,2);            % rad
psi_sim=wing_kenimatics(:,3);              % rad
alpha_sim=wing_kenimatics(:,4);          % alpha=atan2(omega_z,-omega_y);  
dphi_pres=wing_kenimatics(:,5);          % rad/s
dpsi_sim=wing_kenimatics(:,6);            % rad/s
ddphi_pres=wing_kenimatics(:,7);        % rad/s^2
ddpsi_sim=wing_kenimatics(:,8);          % rad/s^2
% C_L=wing_kenimatics(:,9);          
% C_D=wing_kenimatics(:,10);     
% C_N1=wing_kenimatics(:,11);   
% C_T=wing_kenimatics(:,12);  
% C_N=C_N1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
f=188.7;                % Hz % T=1/f;  
[psi_sim_0, locat_0]=min(abs(psi_sim));  
T=1/f;
t_T=t(locat_0)+T;
index=find(t<=t_T);      
num=length(index)+4;     
t1=t((locat_0:num),1);                
phi1=phi_pres((locat_0:num),1);     
psi1=psi_sim((locat_0:num),1);       
alpha1=alpha_sim((locat_0:num),1);   
dphi1=dphi_pres((locat_0:num),1);    
dpsi1=dpsi_sim((locat_0:num),1);    
ddphi1=ddphi_pres((locat_0:num),1);  
ddpsi1=ddpsi_sim((locat_0:num),1);     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F_normal=xlsread('Forcenormal_oneT_simpres.xlsx','A1:C1000');
F_ytran=F_normal((locat_0:num),1);  
F_yadd=F_normal((locat_0:num),2); 
F_yrot=F_normal((locat_0:num),3);   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C_avereff =0.8854; 
C_025=0.25*C_avereff;
C_075=0.75*C_avereff;
C_05=C_avereff/2-C_025; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  
N_skip=12;  
num=338;  
k_index=(1:N_skip:num); 
Nstep=num/2;                        
Y_axisup=linspace(-2.5,2.5,15); 
Z_axisup=zeros(1,length(Y_axisup));
Y_axisdown=fliplr(Y_axisup); 
Z_axisdown=fliplr(Z_axisup);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
k_index1=(1:12:Nstep);
N_down=length(k_index1);    
psi_down=psi1(k_index1);  
alpha_down=alpha1(k_index1);
F_ytrandown=F_ytran(k_index1);  
F_yadddown=F_yadd(k_index1); 
F_yrotdown=F_yrot(k_index1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_step=2*N_down;  
np=2;            
Nframe=np*15;     
Nskip1=floor(N_step/Nframe);
Nspeed=3;  
duration=Nframe/Nspeed;  
animation=1;     
display(['duratin of movie will be ', num2str(duration), ' seconds ']); 
if animation == 1
        vidObj = VideoWriter('video_S2.MP4','MPEG-4'); 
        vidObj.FrameRate = Nspeed;
        open(vidObj);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig=figure;
set(gcf,'position',[40,40, 1400, 700]); 
set(gcf,'color','white'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:1:N_down   
    displacement=[Y_axisdown(1,i);Z_axisdown(1,i)];  
    psi=psi_down(i,1);    
    if (psi<0)
        Pitch = [cos(psi)   sin(abs(psi)); ...
                     -sin(abs(psi))    cos(psi)];   
        R_w2s =Pitch;     
         y1=0;
         z1=C_025;
         p_lead0=[y1;z1];                  
         p_leadrot=R_w2s*p_lead0;   
         pn1=[sign(psi),0;
                    0,           1];
         p_lead=pn1*p_leadrot+displacement;    
         y5=0;
         z5=C_075;
         p_tail0=[y5;z5];    
         p_tailrot=R_w2s*p_tail0;  
         pn2=[1,0;
                  0,sign(psi)];
         p_tail=pn2*p_tailrot+displacement; 
         d_cp=(0.82*abs(alpha_down(i,1))/pi+0.05)*C_avereff;        
         if d_cp>C_025
             y3=0;
             z3=(d_cp-C_025);
             p_cop0=[y3;z3];
             p_coprot=pn2*R_w2s*p_cop0;
         elseif d_cp<=C_025;
             y3=0;
             z3=C_025-d_cp; 
             p_cop0=[y3;z3];
             p_coprot=pn1*R_w2s*p_cop0;
         end
         p_cop=p_coprot+displacement; 
         y4=0;
         z4=C_05;   
         p_mid0=[y4;z4];
         p_midrot=pn2*R_w2s*p_mid0;
         p_mid=p_midrot+displacement; 
    else
        Pitch = [cos(psi)   sin(psi); ...
                     -sin(psi)    cos(psi)];     
        R_w2s =Pitch;     
         y1=0;
         z1=C_025;
         p_lead0=[y1;z1];        
         p_leadrot=R_w2s*p_lead0;  
         p_lead=p_leadrot+displacement;    
         y5=0;
         z5=-C_075;
         p_tail0=[y5;z5];         
         p_tailrot=R_w2s*p_tail0;     
         p_tail=p_tailrot+displacement;       
         d_cp=(0.82*abs(alpha_down(i,1))/pi+0.05)*C_avereff;        
         if d_cp>C_025
             y3=0;
             z3=-(d_cp-C_025);
         elseif d_cp<=C_025;
             y3=0;
             z3=C_025-d_cp; 
         end
         p_cop0=[y3;z3];
         p_coprot=R_w2s*p_cop0;
         p_cop=p_coprot+displacement; 
         y4=0;
         z4=-C_05;   
         p_mid0=[y4;z4];
         p_midrot=R_w2s*p_mid0;
         p_mid=p_midrot+displacement; 
    end
    plot([Y_axisdown(1,1),Y_axisdown(1,N_down)],[Z_axisdown(1,1),Z_axisdown(1,N_down)],'r:','LineWidth',4) 
    hold on
    xlabel('The forward and backward direction of the projection line of the wing torsion axis on the body longitudinal plane ({\ity}_{rr})','FontSize',18,'FontName','Times','FontWeight','Bold')
    ylabel('The vertical direction of the wing stroke plane ({\itz}_{rr})','FontSize',18,'FontName','Times','FontWeight','Bold')  
    grid on
    axis([-3.5,3.5,-1,2]) 
    set(gca,'LineWidth',1.5,'FontSize',14,'FontName','Times','FontWeight','Bold')
    hold on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h3=text(-0.5,2.125,'Xijun Ke @ SJTU');
    set(h3,'Color','k','FontSize',12.0,'FontName','Times','FontWeight','Bold')
    h4=text(0.5,2.125,'01-09-2020');
    set(h4,'Color','k','FontSize',12.0,'FontName','Times','FontWeight','Bold')
    hold on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot([p_lead(1,1),p_tail(1,1)],[p_lead(2,1),p_tail(2,1)],'k-','LineWidth',3)
    hold on
    plot(Y_axisdown(1,i),Z_axisdown(1,i),'ro','LineWidth',1.5,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',4) 
    hold on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot(p_lead(1,1),p_lead(2,1),'ko-','LineWidth',1.5,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)
    hold on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot(p_cop(1,1),p_cop(2,1),'bo','LineWidth',1.5,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',6);  
    hold on   
    y_cop=(F_ytrandown(i,1)+F_yrotdown(i,1))+Y_axisdown(1,i); 
    z_cop=z3+Z_axisdown(1,i);
    F_ytran_cop0=[y_cop;z_cop]; 
    pn2=[1,0;
             0,sign(psi)];
    Fytran_coprot=pn2*R_w2s*F_ytran_cop0;     
    Fytran_cop=Fytran_coprot; 
    quiver(p_cop(1,1),p_cop(2,1),...     
               Fytran_cop(1,1),Fytran_cop(2,1),0.1,'b-','LineWidth',3); 
    hold on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y_mid=F_yadddown(i,1)+Y_axisdown(1,i);  
    z_mid=z4+Z_axisdown(1,i);
    F_yadd_mid0=[y_mid;z_mid];  
    Fyadd_midrot=pn2*R_w2s*F_yadd_mid0;  
    Fyadd_mid=Fyadd_midrot; 
    quiver(p_mid(1,1),p_mid(2,1),...                   
               Fyadd_mid(1,1),Fyadd_mid(2,1),0.1,'g-','LineWidth',3); 
    hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if animation == 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%
       currFrame = getframe(fig);
       writeVideo(vidObj,currFrame);
    else
        pause(0.5)
    end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
k_index2=(Nstep:12:num);  
N_up=length(k_index2);  
psi_up=psi1(k_index2);  
alpha_up=alpha1(k_index2);
F_ytranup=F_ytran(k_index2);  
F_yaddup=F_yadd(k_index2);   
F_yrotup=F_yrot(k_index2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:1:N_up  
    displacement=[Y_axisup(1,i);Z_axisup(1,i)]; 
    psi=psi_up(i,1);  
     if (psi>0)
        Pitch = [cos(psi)   sin(psi); ...
                     -sin(psi)    cos(psi)];    
        R_w2s =Pitch;     
         y1=0;
         z1=C_025;
         p_lead0=[y1;z1];         
         p_leadrot=R_w2s*p_lead0; 
         p_lead=p_leadrot+displacement;   
         y5=0;
         z5=-C_075;
         p_tail0=[y5;z5];                 
         p_tailrot=R_w2s*p_tail0;     
         p_tail=p_tailrot+displacement;          
         d_cp=(0.82*abs(alpha_down(i,1))/pi+0.05)*C_avereff;        
         if d_cp>C_025
             y3=0;
             z3=-(d_cp-C_025);
         elseif d_cp<=C_025;
             y3=0;
             z3=C_025-d_cp; 
         end
         p_cop0=[y3;z3];
         p_coprot=R_w2s*p_cop0;
         p_cop=p_coprot+displacement; 
         y4=0;
         z4=-C_05;   
         p_mid0=[y4;z4];
         p_midrot=R_w2s*p_mid0;
         p_mid=p_midrot+displacement; 
    else  
        Pitch = [cos(psi)   sin(abs(psi)); ...
                     -sin(abs(psi))    cos(psi)];    
        R_w2s =Pitch;    
         y1=0;
         z1=C_025;
         p_lead0=[y1;z1];                  
         p_leadrot=R_w2s*p_lead0;    
         pn1=[sign(psi),0;
                    0,           1];
         p_lead=pn1*p_leadrot+displacement;       
          y5=0;
         z5=C_075;
         p_tail0=[y5;z5];                  
         p_tailrot=R_w2s*p_tail0;    
         pn2=[1,0;
                  0,sign(psi)];
         p_tail=pn2*p_tailrot+displacement;      
         d_cp=(0.82*abs(alpha_down(i,1))/pi+0.05)*C_avereff;        
         if d_cp>C_025
             y3=0;
             z3=(d_cp-C_025);
             p_cop0=[y3;z3];
             p_coprot=pn2*R_w2s*p_cop0;
         elseif d_cp<=C_025;
             y3=0;
             z3=C_025-d_cp; 
             p_cop0=[y3;z3];
             p_coprot=pn1*R_w2s*p_cop0;
         end
         p_cop=p_coprot+displacement; 
         y4=0;
         z4=C_05;   
         p_mid0=[y4;z4];
         p_midrot=pn2*R_w2s*p_mid0;
         p_mid=p_midrot+displacement; 
    end   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot([Y_axisup(1,1),Y_axisup(1,N_up)],[Z_axisup(1,1),Z_axisup(1,N_up)],'r:','LineWidth',4) 
    hold on
    xlabel('The forward and backward direction of the projection line of the wing torsion axis on the body longitudinal plane ({\ity}_{rr})','FontSize',20,'FontName','Times','FontWeight','Bold')   
    ylabel('The vertical direction of the wing stroke plane ({\itz}_{rr})','FontSize',20,'FontName','Times','FontWeight','Bold')  
    grid on
    axis([-3.5,3.5,-1,2])
    set(gca,'LineWidth',1.5,'FontSize',14,'FontName','Times','FontWeight','Bold')
    hold on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h3=text(-0.5,2.125,'Xijun Ke @ SJTU');
    set(h3,'Color','k','FontSize',12.0,'FontName','Times','FontWeight','Bold')
    h4=text(0.5,2.125,'01-09-2020');
    set(h4,'Color','k','FontSize',12.0,'FontName','Times','FontWeight','Bold')
    hold on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot([p_lead(1,1),p_tail(1,1)],[p_lead(2,1),p_tail(2,1)],'k-','LineWidth',3) 
    hold on
    plot(Y_axisup(1,i),Z_axisup(1,i),'ro','LineWidth',1.5,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',4) 
    hold on
    plot(p_lead(1,1),p_lead(2,1),'ko-','LineWidth',1.5,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)
    hold on
    plot(p_cop(1,1),p_cop(2,1),'bo','LineWidth',1.7,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',6);  
    hold on
    y_cop=(F_ytranup(i,1)+F_yrotup(i,1))+Y_axisup(1,i); 
    z_cop=z3+Z_axisup(1,i);
    F_ytran_cop0=[y_cop;z_cop];                  
    Fytran_coprot= R_w2s*F_ytran_cop0;       
    Fytran_cop=Fytran_coprot; 
    quiver(p_cop(1,1),p_cop(2,1),...         
               Fytran_cop(1,1),Fytran_cop(2,1),0.1,'b-','LineWidth',3);    
    hold on
    y_mid=F_yaddup(i,1)+Y_axisup(1,i); 
    z_mid=z4+Z_axisup(1,i);
    F_yadd_mid0=[y_mid;z_mid];  
    Fyadd_midrot=R_w2s*F_yadd_mid0;  
    Fyadd_mid=Fyadd_midrot; 
    quiver(p_mid(1,1),p_mid(2,1),...          
               Fyadd_mid(1,1),Fyadd_mid(2,1),0.1,'g-','LineWidth',3);   
    hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if animation == 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%
       currFrame = getframe(fig);
       writeVideo(vidObj,currFrame);        
    else
        pause(0.5)
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if animation == 1
      close(vidObj);    % Close the file.
end
PROFILES = VideoWriter.getProfiles()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%