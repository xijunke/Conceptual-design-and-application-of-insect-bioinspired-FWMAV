%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The temporal and spacial variable process of chordwise position of center of pressure (CoP)
clear all; clc;
tic                             
load('FRUITFLY_LOOMINGRESPONSE_DATABASE.mat')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
z_cop=xlsread('z_cop2.xlsx','A1:E2000');
alpha2=z_cop(:,2)*pi/180;  % rad
[alpha2_min,local_min]=min(abs(alpha2(1:350,1)));
[alpha2_max,local_max]=max(alpha2(1:350,1)); 
alpha1=linspace(alpha2_min,alpha2_max,20); 
alpha=[alpha1';fliplr(alpha1)']; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COP_distribution_timespace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% animation
N_step=20; 
np=3;   
Nframe=np*20; 
Nspeed=2;  % number of frames per second ------------------------ for change
duration=Nframe/Nspeed;  % duration of movie  % 5 seconds @Nframe=1*15
animation=1;         %  1=save animation file --------------------------- for change
display(['duratin of movie will be ', num2str(duration), ' seconds ']); 
if animation == 1
        vidObj = VideoWriter('video_S1.MP4','MPEG-4');        
        vidObj.FrameRate = Nspeed;
        open(vidObj);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig=figure;
set(gcf,'position',[40,40, 1400, 700]); 
set(gcf,'color','white'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:1:length(alpha)
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
y_lead=linspace(0.3289351765,3.332792062,1000); 
y_trail=linspace(0.3289351765,3.332792062,1000);  
%%%%%%%%%%%%%%%%%%%%%%%%%
x_mod_Root =0.636;   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C_maxylb =0.464385778290230;
C_maxy25 =0.138924474377504; 
C_maxyub =-0.186536829535222; 
C_maxy=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_x_lead=-0.08249*y_lead.^6 +0.9167*y_lead.^5-4.04*y_lead.^4+8.872*y_lead.^3-10.06*y_lead.^2+5.674*y_lead-0.413-x_mod_Root-C_maxy;
f_x_trail=-0.0333*y_trail.^6+0.504*y_trail.^5-2.795*y_trail.^4+7.258*y_trail.^3-8.769*y_trail.^2+3.739*y_trail+0.1282-x_mod_Root-C_maxy;
C_x=f_x_lead-f_x_trail;
plot(y_lead,f_x_lead,'k-',y_trail,f_x_trail,'k-','LineWidth',3) 
xlabel('The spanwise {\itx}-axis of the right wing ({\itx}_{rw})')  
ylabel('The chordwise {\itz}-axis of the right wing ({\itz}_{rw})') 
axis([0,3.5,-1,0.5])
% set(gca,'LineStyle','-','LineWidth',1.5,'FontSize',18,'FontName','Times','FontWeight','Bold')
set(gca,'LineWidth',1.5,'FontSize',18,'FontName','Times','FontWeight','Bold')
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
C_avereff=0.8854;  % mm
[C_lead_ymax,k_leadmax]=max(f_x_lead);   
yr_leadnd_max=C_lead_ymax/C_avereff; 
C_max_x=y_lead(1,k_leadmax);             
C_trail_y=f_x_trail(1,k_leadmax);            
C_leadmax=C_lead_ymax-C_trail_y;           
[C_trail_ymin,k_trailmin]=min(f_x_trail);    
yr_trailnd_min=C_trail_ymin/C_avereff;  
C_min_x=y_lead(1,k_trailmin);   
C_lead_y=f_x_lead(1,k_trailmin);     
C_trailmin=C_lead_y-C_trail_ymin;      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
C_max_LtoT=C_lead_ymax-C_trail_ymin; 
x_0lb=0*C_max_LtoT;     
x_025=0.25*C_max_LtoT;    
x_0ub=0.5*C_max_LtoT;  
x_0=x_025;
zero_x=C_lead_ymax-x_0;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%
x_mod_Root1=0;  
plot([0,3.5],[x_mod_Root1,x_mod_Root1],'k-.','LineWidth',3)
hold on
plot(3.332792062,x_mod_Root1,'k*','LineWidth',5) 
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h3=text(1.35,0.563,'Xijun Ke @ SJTU');
set(h3,'Color','k','FontSize',12.0,'FontName','Times','FontWeight','Bold')
h4=text(1.85,0.563,'01-09-2020');
set(h4,'Color','k','FontSize',12.0,'FontName','Times','FontWeight','Bold')
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
XC=1.920243385; 
YC=-0.149785466-C_maxy; 
plot(XC,YC,'dk','LineWidth',4)
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
C_rx=f_x_lead-f_x_trail;
dealta_x0=0.356737-0.25;
d_cprnd1=0.82*abs(alpha(i,1))/pi+0.05+dealta_x0; 
Cop_r1=f_x_lead-d_cprnd1*C_rx;
A=Cop_r1'-x_mod_Root;  
plot(y_lead,Cop_r1,'r-.','LineWidth',3) 
hold on
%% 
Cop_add=f_x_lead-(f_x_lead-f_x_trail)/2;
plot(y_lead,Cop_add,'g-.','LineWidth',4) 
hold off
%%%%%%%%%%%%%%%%%%%%%%%%
    if animation == 1
       currFrame = getframe(fig);
       writeVideo(vidObj,currFrame);           
    else
        pause(0.5)
    end
end 
hold on
% pause(5)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:1:length(alpha)/2
%% 
y_lead=linspace(0.3289351765,3.332792062,1000);  
y_trail=linspace(0.3289351765,3.332792062,1000);   
x_mod_Root =0.636;   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C_maxylb =0.464385778290230;
C_maxy25 =0.138924474377504;  
C_maxyub =-0.186536829535222; 
C_maxy=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_x_lead=-0.08249*y_lead.^6 +0.9167*y_lead.^5-4.04*y_lead.^4+8.872*y_lead.^3-10.06*y_lead.^2+5.674*y_lead-0.413-x_mod_Root-C_maxy;
f_x_trail=-0.0333*y_trail.^6+0.504*y_trail.^5-2.795*y_trail.^4+7.258*y_trail.^3-8.769*y_trail.^2+3.739*y_trail+0.1282-x_mod_Root-C_maxy;
C_x=f_x_lead-f_x_trail;
plot(y_lead,f_x_lead,'k-',y_trail,f_x_trail,'k-','LineWidth',3)  
xlabel('The spanwise {\itx}-axis of the right wing ({\itx}_{rw})')  
ylabel('The chordwise {\itz}-axis of the right wing ({\itz}_{rw})') 
axis([0,3.5,-1,0.5])
% set(gca,'LineStyle','-','LineWidth',1.5,'FontSize',18,'FontName','Times','FontWeight','Bold')
set(gca,'LineWidth',1.5,'FontSize',18,'FontName','Times','FontWeight','Bold')
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 
C_avereff=0.8854;  % mm
[C_lead_ymax,k_leadmax]=max(f_x_lead);   
yr_leadnd_max=C_lead_ymax/C_avereff;    
C_max_x=y_lead(1,k_leadmax);                   
C_trail_y=f_x_trail(1,k_leadmax);               
C_leadmax=C_lead_ymax-C_trail_y;           
[C_trail_ymin,k_trailmin]=min(f_x_trail);    
yr_trailnd_min=C_trail_ymin/C_avereff;    
C_min_x=y_lead(1,k_trailmin);              
C_lead_y=f_x_lead(1,k_trailmin);            
C_trailmin=C_lead_y-C_trail_ymin;            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
C_max_LtoT=C_lead_ymax-C_trail_ymin;    
x_0lb=0*C_max_LtoT;   
x_025=0.25*C_max_LtoT;     
x_0ub=0.5*C_max_LtoT;     
x_0=x_025;
zero_x=C_lead_ymax-x_0;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
x_mod_Root1=0;   
plot([0,3.5],[x_mod_Root1,x_mod_Root1],'k-.','LineWidth',3)
hold on
plot(3.332792062,x_mod_Root1,'k*','LineWidth',5)  
hold on
h3=text(1.35,0.563,'Xijun Ke @ SJTU');
set(h3,'Color','k','FontSize',12.0,'FontName','Times','FontWeight','Bold')
h4=text(1.85,0.563,'01-09-2020');
set(h4,'Color','k','FontSize',12.0,'FontName','Times','FontWeight','Bold')
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
XC=1.920243385; 
YC=-0.149785466-C_maxy; 
plot(XC,YC,'dk','LineWidth',4)
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
C_rx=f_x_lead-f_x_trail;
dealta_x0=0.356737-0.25;
d_cprnd1=0.82*abs(alpha(i,1))/pi+0.05+dealta_x0; 
Cop_r1=f_x_lead-d_cprnd1*C_rx;
A=Cop_r1'-x_mod_Root; 
plot(y_lead,Cop_r1,'r-.','LineWidth',3) 
hold on

%% 
Cop_add=f_x_lead-(f_x_lead-f_x_trail)/2;
plot(y_lead,Cop_add,'g-.','LineWidth',4) 
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if animation == 1
       currFrame = getframe(fig);
       writeVideo(vidObj,currFrame);           
    else
        pause(0.5)
    end
end 
hold on
% pause(5)

for i=1:1:5
R_wingeff=3.004; 
xr=0.3289;    
r_xcopnd_tr_vari=-z_cop(:,3);
r_xcopnd_rot_vari=-z_cop(:,4);  
r_xcopnd_vari=-z_cop(:,5);
r_xcopnd_tr =0.788691874094779; 
R_tr=R_wingeff*r_xcopnd_tr;  
R_tr_vari=R_wingeff*r_xcopnd_tr_vari;  
r_xcopnd_rot =0.7126;  
R_rot=R_wingeff*r_xcopnd_rot;   
R_rot_vari=R_wingeff*abs(r_xcopnd_rot_vari);       
r_xcopnd_addaver=0.7320;  
R_add=R_wingeff*r_xcopnd_addaver; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
syms y_lead  y_trail
f_x_lead=-0.08249*y_lead.^6 +0.9167*y_lead.^5-4.04*y_lead.^4+8.872*y_lead.^3-10.06*y_lead.^2+5.674*y_lead-0.413-x_mod_Root-C_maxy;
f_x_trail=-0.0333*y_trail.^6+0.504*y_trail.^5-2.795*y_trail.^4+7.258*y_trail.^3-8.769*y_trail.^2+3.739*y_trail+0.1282-x_mod_Root-C_maxy;
f_x_lead1=inline(vectorize(f_x_lead),'y_lead');
f_x_trail1=inline(vectorize(f_x_trail),'y_trail');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
xr=0.3289;  
R_ref=xr+0.7*R_wingeff;                
f_x_lead_ref=f_x_lead1(R_ref);     
f_x_trail_ref=f_x_trail1(R_ref);        
C_ref=f_x_lead_ref-f_x_trail_ref;     
%%%%%%%%%%%%%%%%
yr_leadbem=1.0958;         
x_mod_Root =0.6360;                    
L_arm1=yr_leadbem-x_mod_Root;  
C_bem =1.1257;                   
L_arm=C_bem/2-L_arm1;   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
f_x_lead_cop_tr=f_x_lead1(R_tr);         
f_x_trail_cop_tr=f_x_trail1(R_tr);                     
C_cop_tr=f_x_lead_cop_tr-f_x_trail_cop_tr;       
Ratio_axistr=f_x_lead_cop_tr/C_cop_tr;               
f_x_lead_cop_rot=f_x_lead1(R_rot);    
f_x_trail_cop_rot=f_x_trail1(R_rot);                   
C_cop_rot=f_x_lead_cop_rot-f_x_trail_cop_rot;  
Ratio_axisrot=f_x_lead_cop_rot/C_cop_rot;         
f_x_lead_cop_add=f_x_lead1(R_add);       
f_x_trail_cop_add=f_x_trail1(R_add);               
C_cop_add=f_x_lead_cop_add-f_x_trail_cop_add;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
z_add=-(C_cop_add/2-f_x_lead_cop_add);  
plot(R_add,z_add,'dm','LineWidth',4)
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
x_cop=xlsread('x_cop2.xlsx','A1:I2000');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y_rcpnd=x_cop(:,3);  
Y_rcpndaver=mean(Y_rcpnd);    
%%%%%%%%%%%%%%%%%%%%%%%%%
Z_trans1=x_cop(:,4);    
Z_transaver1=mean(Z_trans1);    
Y_rcpnd_trans=x_cop(:,5);  
Y_rcpnd_transaver=mean(Y_rcpnd_trans);  
Z_trans=x_cop(:,5);
Z_transaver=mean(Z_trans); 
Y_rcpnd_rot=x_cop(:,7); 
Y_rcpnd_rotaver=mean(Y_rcpnd_rot); 
Z_rot=x_cop(:,7); 
Z_rotaver=mean(Z_rot);  
Z_circ=x_cop(:,9);  
Z_circaver=mean(Z_circ);   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
z_tr=1.125*C_cop_tr*Z_trans;   
for i=1:1:length(Y_rcpnd_trans)
    plot(R_tr_vari(i,1),z_tr(i,1),'c.','LineWidth',1) 
    hold on
end
z_traver=1.125*C_cop_tr*Z_transaver; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
z_rot=1.04*C_cop_rot*Z_rot;  
for i=1:1:length(Y_rcpnd_trans)
    plot(R_rot_vari(i,1),z_rot(i,1),'b.','LineWidth',1) 
    hold on
end
z_rotaver=1.04*C_cop_rot*Z_rotaver; 
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if animation == 1
       currFrame = getframe(fig);
       writeVideo(vidObj,currFrame);   
    else
        pause(0.5)
    end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end 

pause(10)
if animation == 1
      close(vidObj);    % Close the file.
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PROFILES = VideoWriter.getProfiles()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc
