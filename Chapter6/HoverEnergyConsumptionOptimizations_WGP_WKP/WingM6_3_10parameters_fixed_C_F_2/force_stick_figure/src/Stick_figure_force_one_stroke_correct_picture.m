%% Stick figure force
clear all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 第一部分——调用翅膀运动学-翅运动规律和几何攻角(AOA)等数据
% wing_m_output=[t',phi_pres',psi_sim',alpha',dphi',dpsi',ddphi',ddpsi',C_L',C_D',C_N1',C_T'];
% wing_kenimatics=kenimatics_wing_and_AoA_fruitfly_exp();
% wing_kenimatics=kenimatics_wing_and_AoA_fruitfly_simpres(); 
%%%%%%%%%%%%%%%%%%%%%%%%%% %调用翅膀运动学-翅运动规律和几何攻角(AOA)等数据
% x =[46.6104,1.3076,0.1891,1.2715,2.5103,-1.5169,3.9813,1.96,1.7845,0.0001]; % P_asterisk =4.7867;L =1.0000; delta = -2.9524e-008; AR =2.9417; Re =186.1461;% 下一行为精确值
x =[4.6610447e+001,1.3075639e+0,1.8911207e-001,1.2714657e+0,2.5103035e+0, -1.5168925e+0,3.9812634e+0,1.9600006e+0,1.7844515e+0,7.9074669e-005];
f=x(1);  T=1/f;   % Hz 
phi_m=x(2);
K=x(3);                  %  epsilon=x(3);    
eta_m=x(4);          % phi_0=x(4);       
C_eta=x(5);           %  psi_m=x(5);
Phi_eta=x(6);        % zeta=x(6); 
% eta_0=x(7);       % psi_0=x(7);   
wing_kenimatics=kenimatics_wing_and_AoA_fruitfly_sim(f,phi_m,K,eta_m,C_eta,Phi_eta); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%调用函数kenimatics_wing_and_AoA; % size(wing_kenimatics)    % (1000,12)
% wing_kenimatics=kenimatics_wing_and_AoA_fruitfly_exp();
t=wing_kenimatics(:,1);                        % 单位是ms
phi_pres=wing_kenimatics(:,2);            % 拍打角——单位是rad
psi_sim=wing_kenimatics(:,3);              % 拍打角——单位是rad
alpha_sim=wing_kenimatics(:,4);          % alpha=atan2(omega_z,-omega_y);  % 几何攻角——弧度制   %输出——有正有负
dphi_pres=wing_kenimatics(:,5);          % 单位是rad/s
dpsi_sim=wing_kenimatics(:,6);            % 单位是rad/s
ddphi_pres=wing_kenimatics(:,7);        % 单位是rad/s^2
ddpsi_sim=wing_kenimatics(:,8);          % 单位是rad/s^2
% C_N1=wing_kenimatics(:,9);   
%%%%%%%%%%%%%%%%%%%%%%%%%
% C_L=wing_kenimatics(:,9);          
% C_D=wing_kenimatics(:,10);     
% C_N1=wing_kenimatics(:,11);   
% C_T=wing_kenimatics(:,12);  
% C_N=C_N1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
% f=188.7;                % Hz——翅拍频率 % T=1/f;  
% figure(1)
% plot(t*f,psi_sim*180/pi,'r-',t*f,phi_pres*180/pi,'g-',t*f,alpha_sim*180/pi,'b:','LineWidth',2)
% xlabel('\itNormalized time');  
% ylabel('\psi_{sim}(t) and \phi_{pres}(t) and \alpha_{sim}(t)');  
% legend('\psi_{sim}(t)','\phi_{pres}(t)','\alpha_{sim}(t)');  
% title('被动转动角\psi_{sim}(t)和拍打角phi_{pres}(t)随时间的变化')  
% grid on  % 被动转动角psi_sim(t)和扭转角phi_pres(t)随时间的变化
% axis([0.9,4.05,-110,110])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) 截取单个稳定的周期——以攻角变化周期为区间
% 采用攻角为90度对应时刻为上下冲程的结束时间点
[psi_sim_0, locat_0]=min(abs(psi_sim));   % psi_sim_0 =6.2058e-004; locat_0 =3;
T=1/f;
t_T=t(locat_0)+T;
index=find(t<=t_T);          % 336
num=length(index)+4;     % num=336+2;
% [alpha_sim_max, locat_0]=min(abs(alpha_sim));   % alpha_sim_max =0.3581; locat_max =219;
% T=1/f;
% t_T=t(locat_0)+T;
% index=find(t<=t_T);      % 552
% num=length(index);     % num=552;
% 刚好104个数据点=275-172+1;
t1=t((locat_0:num),1);                             % 单位是ms
phi1=phi_pres((locat_0:num),1);             % 拍打角——单位是rad
psi1=psi_sim((locat_0:num),1);               % 扭转角——单位是rad
alpha1=alpha_sim((locat_0:num),1);       % 扭转角——单位是rad%  alpha=atan2(omega_z,-omega_y);  %输出——有正有负
dphi1=dphi_pres((locat_0:num),1);        % 拍打角——单位是rad/s
dpsi1=dpsi_sim((locat_0:num),1);           % 扭转角——单位是rad/s
ddphi1=ddphi_pres((locat_0:num),1);     % 拍打角——单位是rad/s^2
ddpsi1=ddpsi_sim((locat_0:num),1);       % 扭转角——单位是rad/s^2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 问题(1)———————已经修改,由其是虚质量力臂L_arm=C_bem/2-L_arm1; 
% (3)将单个稳定的周期的总(合成)法向气动力数据读入(mN)—这里输入的数据应该是单个片元的受力
% % 以备thr_dim_chord2程序的调用用于绘制球棍图——翅坐标系下的法向力失终点y坐标
% B=[F_ytran,F_yadd1,F_yrot];    % size(B)  % (2000*1)     %——mN
% F_normal=xlsread('D:\KXJ\PassiveRot_dynamic_Science_fruitfly\optimal\optimal_wing_para\Forcenormal_oneT_simpres.xlsx','A1:C1000');
% F_normal=xlsread('D:\KXJ\PassiveRot_dynamic_Science_fruitfly\optimal\optimal_wing_para\Forcenormal_3T.xlsx','A1:D1000');
F_normal=xlsread('Forcenormal_3T.xlsx',1,'A1:D1000'); % 读入数据
% F_norm=F_normal((locat_0:num),1);    %  被缩小了0.1倍显示quiver3
F_ytran=F_normal((locat_0:num),1);         %  被缩小了0.1倍显示quiver3 
F_yadd=F_normal((locat_0:num),2);         %  被缩小了0.1倍显示quiver3
F_yrot=F_normal((locat_0:num),3);          %  被缩小了0.1倍显示quiver3
F_inert_y=0*F_normal((locat_0:num),4); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(2)
% plot(t1*f,psi1*180/pi,'r*',t1*f,phi1*180/pi,'g>',t1*f,alpha1*180/pi,'k-',...
%        t1*f,3*F_ytran,'b-',t1*f,3*F_yadd,'c-',t1*f,3*F_yrot,'m-','LineWidth',1.5)  % 转换为度数
% xlabel('\itNormalized time');  
% ylabel('\psi_{sim}(t) and \phi_{pres}(t) and \alpha_{sim}(t) and 3*F_{norm}(t))');  
% legend('\psi_{sim}(t)','\phi_{pres}(t)','\alpha_{sim}(t)','3*F_{y,tran}(t))','3*F_{y,add}(t))','3*F_{y,rot}(t))');  
% title('被动转动角\psi_{sim}(t), 拍打角phi_{pres}(t) and 3*F_{norm}(t))随时间的变化')  
% grid on  % 被动转动角psi_sim(t)和扭转角phi_pres(t)和法向气动力分量随时间的变化
% % axis([1.6,2.7,-110,110])
% axis([0.95,2.05,-inf,inf])
% set(gca,'XTick',(0.95:0.05:2.05))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 第二部分—— 绘制单片翅弦微元的二维空间位置
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % R_wingeff =3.0040;        % 单位是 mm   
% % xr=0.3289;                     % x-root offset  \mm
% % R_eff=R_wingeff+xr;
% 问题(2)—————建议不去平均弦长作为片条的长度，应该采用特征展向长度对应的片条弦长，比如下面的% 特征弦长C_bem=1.1257mm;   
C_avereff =0.8854;          % 单位是 mm
C_025=0.25*C_avereff;
C_075=0.75*C_avereff;
% C_max_LtoT= 1.3018;
% C_avereff =C_max_LtoT;  
% C_025=0.36*C_avereff;
% C_075=0.64*C_avereff;
C_05=C_avereff/2-C_025; 
%%%%%%%%%%%%%%%%%%%%%%%%%
%% 针对环量最大处所对应的翅片条的展向长度R_ref和弦长C_bem——来自wing_model_88_translation_yaxis程序
% % R_wingeff =3.0040;        % 单位是 mm
% % xr=0.3289;  % x-root offset  \mm      % R_wingeff=wing_para(1,1);  % R_wingeff =3.0040;   % 单位是 mm
% % R_ref=(xr+0.7*R_wingeff);                  % mm——特征展向长度R_ref      
% f_x_lead_ref =0.4598; %  mm % f_x_lead_ref=f_x_lead1(R_ref);       
% f_x_trail_ref =-0.6658; %  mm %  f_x_trail_ref=f_x_trail1(R_ref); 
% C_ref=f_x_lead_ref-f_x_trail_ref;     % C_ref =1.1257; %  mm
% %%%%%%%%%%%%%%%%
% % yr_leadbem=1.0958;                         % mm
% % x_mod_Root =0.6360;                       % mm
% % L_arm1=yr_leadbem-x_mod_Root;   % L_arm1=0.4598;     % 前缘至扭转轴的距离:  mm
% % C_bem =1.1257;                                % C_bem =1.1257;    % 特征弦长C_bem:  mm
% % L_arm=C_bem/2-L_arm1;                  % % L_arm =0.1030; % 虚质量片元力臂被修正了(/2):  mm
% C_avereff =C_ref;          % 单位是 mm
% C_025=f_x_lead_ref;
% C_075=-f_x_trail_ref;  % C_075取为正号
% C_bem=C_ref;
% C_05=C_bem/2-f_x_lead_ref;   % 中弦点到扭转轴的距离——修改了针对实际的特征弦长的中弦点，而非0.25*C_avereff的中弦点
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  分解上下冲程
N_skip=12;  % 一个周期取28个点，半个周期取14个点
num=338;  % num=336+2;
k_index=(1:N_skip:num); %  length=length(k_index)  % (1*28)  
Nstep=num/2;                        % Nstep=169;
% k_index1=(1:12:Nstep);       % 1...13..25......169  % k_index1(1,1:15);  
% N_down=length(k_index1);      % (15)
% k_index2=(Nstep:12:num);  % 169...181...193......337  % k_index2(1,15:29);
% N_up=length(k_index2); % (15)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 扭转轴点坐标——人为设定等时间距离——上冲程
Y_axisup=linspace(-2.5,2.5,15); 
Z_axisup=zeros(1,length(Y_axisup));
%  Z_axisup=zeros(1,length(Y_axisup))+2;
% 扭转轴点坐标——人为设定等时间距离——下冲程
Y_axisdown=fliplr(Y_axisup); 
Z_axisdown=fliplr(Z_axisup);
% Z_axisdown=zeros(1,length(Y_axisdown))+2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (1) 下冲程
k_index1=(1:12:Nstep);         % 
N_down=length(k_index1);        % (15)
% t_down=t1(k_index1);    % size(t_down)  %  (15*1)
psi_down=psi1(k_index1);  % 针对一个周期之内的数据，每隔8个点取一个数据形成向量
alpha_down=alpha1(k_index1);
% figure(3)
% plot(t_down*f,psi_down*180/pi,'rd','LineWidth',2)      % 转换为度数
% Forcenormaldown=Forcenormal(k_index1);
F_ytrandown=F_ytran(k_index1);  % size(F_ytrandown)   % (15*1)     %  10^(-5)*
F_yadddown=F_yadd(k_index1);   % size(F_yadddown)   % (15*1)
F_yrotdown=F_yrot(k_index1);
F_inert_ydown=F_inert_y(k_index1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 翅拍上冲程——不考虑扭转
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i = 1:1:N_down   % N_down=15
%     % (1) 前缘点坐标——实心圆点绘制
%    y1=Y_axisdown(1,i);
%    z1=C_025+Z_axisdown(1,i);
%     % (2) 弦向压心点————实心圆点绘制 
%     alpha=pi/2;             % 扭转角——单位是rad
%     d_cp=(0.82*abs(alpha)/pi+0.05)*C_avereff;        
%     if d_cp>C_025
%         y3=Y_axisdown(1,i);
%         z3=-(d_cp-C_025)+Z_axisdown(1,i);
%     elseif d_cp<=C_025;
%         y3=Y_axisdown(1,i);
%         z3=(C_025-d_cp)+Z_axisdown(1,i);
%     end
%     % (3) 中弦点坐标—力失的起点—实心圆绘制
%     y4=Y_axisdown(1,i);
%     z4=-C_025+Z_axisdown(1,i);
%     % (4) 后缘点坐标——实心圆点绘制
%     y5=Y_axisdown(1,i);
%     z5=-C_075+Z_axisdown(1,i);
%     % 绘制片元弦向的特殊点阵的分布
%     figure(4)
%     plot([y1,y5],[z1,z5],'k-','LineWidth',1.5)          %  绘制直线——连接前后缘
%     hold on
%     plot(Y_axisdown(1,i),Z_axisdown(1,i),'ro','LineWidth',1.5,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',3) % 扭转轴心——实心圆点绘制
%     hold on
%     plot(y1,z1,'ko','LineWidth',1.5,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',4);   % 前缘——实心圆点绘制
%     hold on
%     plot(y4,z4,'go','LineWidth',1.5,'MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',3);  % 中弦点—力失的起点—实心圆绘制
%     hold on
%     % plot(y5,z5,'ko','LineWidth',1.5,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',1.5)   % 后缘——实心圆点绘制
%     % hold on
%     % (4) 弦向压心点—绘制力失       
%        % 弦向压心点力失的起点
%     plot(y3,z3,'bo','LineWidth',1.5,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',3);  % 弦向压心点————实心圆点绘制
%     hold on
% %     % p_cop=p_coprot+displacement;  
% %     % displacement=[Y_axisdown(1,i);Z_axisdown(1,i)];  % (2*1)
% %     y_cop=-F_ytrandown(i,1)+Y_axisdown(1,i);    %——注意符号
% % %     z_cop=z3+Z_axisdown(1,i);
% %      z_cop=0;
% %     F_ytran_cop0=[y_cop;z_cop];  % size(wing_chord)   % (2*2)    % 弦向压心点力失的终点
% % %     Fytran_coprot= R_w2s*F_ytran_cop0;        % wing_chord in body frame  % (2*2)     %—R_w2s'—注意符号
% %     Fytran_coprot= F_ytran_cop0;        % wing_chord in body frame  % (2*2)     %—R_w2s'—注意符号
% %     Fytran_cop=Fytran_coprot; 
% %     quiver(y3,z3,...                                                                 % 起点
% %                Fytran_cop(1,1),Fytran_cop(2,1),0.5,'b-','LineWidth',2.0); % 比例因子0.5        % 终点
% %     hold on
% end
% xlabel('纵向(y)')
% ylabel('垂直方向向(z)')
% grid on
% axis([-4,4,-3,3])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % animation———————— 动画片段程序开始段——(1)
% % N_step=N_down;         % N_down=15; 
% % np=1;                    %  调节总的帧数和时间步长
% N_step=2*N_down;         % N_down=15; 
% np=2;                    %  调节总的帧数和时间步长
% Nframe=np*15;     % number of movie frames -------------------------- for change
% Nskip1=floor(N_step/Nframe);
% Nspeed=3;  % number of frames per second ------------------------ for change
% duration=Nframe/Nspeed;  % duration of movie  % 5 seconds @Nframe=1*15
% animation=1;         %  1=save animation file --------------------------- for change
% display(['duratin of movie will be ', num2str(duration), ' seconds ']); %输出显示动画持续的时间(s)
% %输出如下: duratin of movie will be 5 seconds 
% if animation == 1
% %     aviobj = avifile('Stick_figure_force_downstroke.avi','fps',Nspeed); % filename --- for change
%         aviobj = avifile('Stick_figure_force_one_stroke.avi','fps',Nspeed); % filename --- for change
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fig=figure;
fig=figure(1);
for i = 1:1:N_down   % N_down=15
    displacement=[Y_axisdown(1,i);Z_axisdown(1,i)];  % (2*1)
    % 下面一段是坐标系之间的变换
    psi=psi_down(i,1);               % 扭转角——单位是rad    %——注意为负号 
    % psi_deg=psi_down(i,1)*180/pi  
    if (psi<0)
         % 下面一段是坐标系之间的变换
        Pitch = [cos(psi)   sin(abs(psi)); ...
                     -sin(abs(psi))    cos(psi)];        % pitch matrix——符号可能需要方向
        R_w2s =Pitch;     % from wing frame to spar frame——从翅坐标系到体坐标系
        % (1) 前缘点坐标——实心圆点绘制
         y1=0;
         z1=C_025;
         p_lead0=[y1;z1];                   % (2*1)
         p_leadrot=R_w2s*p_lead0;    % (2*2)*(2*1)=(2*1) 
         pn1=[sign(psi),0;
                    0,           1];
         p_lead=pn1*p_leadrot+displacement;          % (2*1)
        % (2) 后缘点坐标——实心圆点绘制
         y5=0;
         z5=C_075;
         p_tail0=[y5;z5];                   % (2*1)
         p_tailrot=R_w2s*p_tail0;      % (2*2)*(2*1)=(2*1)
         pn2=[1,0;
                  0,sign(psi)];
         p_tail=pn2*p_tailrot+displacement;          % (2*1)   
         % (3) 弦向压心点————实心圆点绘制 
         % alpha_up=alpha1(i,1);             % 扭转角——单位是rad
         d_cp=(0.82*abs(alpha_down(i,1))/pi+0.05)*C_avereff;        
         if d_cp>C_025
             y3=0;
             z3=(d_cp-C_025);% 负值
             p_cop0=[y3;z3];
             p_coprot=pn2*R_w2s*p_cop0;
         elseif d_cp<=C_025;
             y3=0;
             z3=C_025-d_cp; % 正值
             p_cop0=[y3;z3];
             p_coprot=pn1*R_w2s*p_cop0;
         end
         p_cop=p_coprot+displacement; 
         % (4) 中弦点坐标—力失的起点—实心圆绘制
         y4=0;
         % z4=-C_025; % 中弦点到扭转轴的距离
         z4=C_05;   % 中弦点到扭转轴的距离——修改了针对实际的特征弦长的中弦点，而非0.25*C_avereff的中弦点
         p_mid0=[y4;z4];
         p_midrot=pn2*R_w2s*p_mid0;
         p_mid=p_midrot+displacement; 
    else  % (psi>0)
        Pitch = [cos(psi)   sin(psi); ...
                     -sin(psi)    cos(psi)];        % pitch matrix——符号可能需要方向
        R_w2s =Pitch;     % from wing frame to spar frame——从翅坐标系到体坐标系
        % (1) 前缘点坐标——实心圆点绘制
         y1=0;
         z1=C_025;
         p_lead0=[y1;z1];                   % (2*1)
         p_leadrot=R_w2s*p_lead0;    % (2*2)*(2*1)=(2*1)   % [-sign(psi),-sign(psi)]*
         p_lead=p_leadrot+displacement;          % (2*1)
        % (2) 后缘点坐标——实心圆点绘制
         y5=0;
         z5=-C_075;
         p_tail0=[y5;z5];                   % (2*1)
         p_tailrot=R_w2s*p_tail0;      % (2*2)*(2*1)=(2*1)
         p_tail=p_tailrot+displacement;          % (2*1)   
         % (3) 弦向压心点————实心圆点绘制 
         % alpha_up=alpha1(i,1);             % 扭转角——单位是rad
         d_cp=(0.82*abs(alpha_down(i,1))/pi+0.05)*C_avereff;        
         if d_cp>C_025
             y3=0;
             z3=-(d_cp-C_025);% 负值
         elseif d_cp<=C_025;
             y3=0;
             z3=C_025-d_cp; % 正值
         end
         p_cop0=[y3;z3];
         p_coprot=R_w2s*p_cop0;
         p_cop=p_coprot+displacement; 
         % (4) 中弦点坐标—力失的起点—实心圆绘制
         y4=0;
         % z4=-C_025; % 中弦点到扭转轴的距离
         z4=-C_05;   % 中弦点到扭转轴的距离——修改了针对实际的特征弦长的中弦点，而非0.25*C_avereff的中弦点
         p_mid0=[y4;z4];
         p_midrot=R_w2s*p_mid0;
         p_mid=p_midrot+displacement; 
    end
    % 片元弦向的特殊点阵
    %  y=[p_lead(1,1),Y_axisdown(1,i),p_cop(1,1),p_mid(1,1),p_tail(1,1)];
    %  z=[p_lead(2,1),Z_axisdown(1,i),p_cop(2,1),p_mid(2,1),p_tail(2,1)];
     %  wing_C=[y;z];  % size(wing_chord)   % (3*5)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % figure(5)
    % figure(i)
    plot([Y_axisdown(1,1),Y_axisdown(1,N_down)],[Z_axisdown(1,1),Z_axisdown(1,N_down)],'r:','LineWidth',1.5) % 扭转轴心连线——虚线绘制
    hold on
%     % h1=text(-2.4,2.32,'Transient normal force vectors generated during one flapping stroke');
%     h1=text(-2.4,3.42,'Transient normal force vectors generated during one flapping stroke');
%     set(h1,'Color','k','FontSize',8,'FontName','Times','FontWeight','Bold')
%    % h2=text(-1.0,2.2,'for optimized wing geometry');
%     h2=text(-3.0,3.25,'for combined optimal WGP and WKP with 2D translational aerodynamic coefficients');
%     set(h2,'Color','k','FontSize',8,'FontName','Times','FontWeight','Bold')
%     h3=text(-0.25,-0.8,'Xijun Ke @ SJTU');
%     set(h3,'Color','k','FontSize',7.0,'FontName','Times','FontWeight','Bold')
%     h4=text(-0.15,-0.915,'17.03.2016');
%     set(h4,'Color','k','FontSize',7.0,'FontName','Times','FontWeight','Bold')
%     hold on
    % xlabel('纵向(y)')
    % ylabel('垂直方向向(z)')
    % grid on
    % % axis([-4,4,-1,3.1])     % 比例不合适于法向气动力的表达或者显示
    %     axis([-4,4,-1.9,2.1])   % axis([-4,4,-0.6808,1.4952])
    % axis([-3.5,3.5,-1,2]) 
    % axis([-3.5,4.25,-1,2.5])
    axis([-3.5,4.25,-1,3.5])
    hold on
    % (1) 绘制直线——连接前后缘
    plot([p_lead(1,1),p_tail(1,1)],[p_lead(2,1),p_tail(2,1)],'k-','LineWidth',2) % 绘制前缘后缘实线
    hold on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (2) 扭转轴心——实心圆点绘制
    plot(Y_axisdown(1,i),Z_axisdown(1,i),'ro','LineWidth',1.5,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',3) % 扭转轴心——实心圆点绘制
    hold on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (3) 绘制前缘实心圆点
    plot(p_lead(1,1),p_lead(2,1),'ko-','LineWidth',1.5,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5)
    hold on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % (4) 弦向压心点—绘制力失
    plot(p_cop(1,1),p_cop(2,1),'bo','LineWidth',1.5,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',5);  % 弦向压心点————实心圆点绘制
    hold on   
    % p_cop=p_coprot+displacement;                 % 弦向压心点力失的起点
    % displacement=[Y_axisdown(1,i);Z_axisdown(1,i)];          % (2*1)
    % y_cop=-F_ytrandown(i,1)+Y_axisdown(1,i);                     %——注意符号——仅含有平动法向力分量
    y_cop=(F_ytrandown(i,1)+F_yrotdown(i,1))+Y_axisdown(1,i); %——注意符号——含有平动法向力分量和转动环量法向力分量    
    z_cop=z3+Z_axisdown(1,i);
    F_ytran_cop0=[y_cop;z_cop];                           % size(wing_chord)   % (2*2)    % 弦向压心点力失的终点
    %   pn1=[sign(psi),0;
    %            0,           1];
    pn2=[1,0;
             0,sign(psi)];
    Fytran_coprot=pn2*R_w2s*F_ytran_cop0;              % wing_chord in body frame  % (2*2)     %—R_w2s'—注意符号
    Fytran_cop=Fytran_coprot; 
    quiver(p_cop(1,1),p_cop(2,1),...                        % 起点
               Fytran_cop(1,1),Fytran_cop(2,1),0.1,'b-','LineWidth',2);    % 终点   % 比例因子0.05——————注意力失被缩小了0.14倍 
    hold on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % (5) 中弦点—绘制力失
   % plot(y4,z4,'go','LineWidth',1.5,'MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',3);  % 中弦点—力失的起点—实心圆绘制
   % hold on  
    y_mid=F_yadddown(i,1)+F_inert_ydown(i,1)+Y_axisdown(1,i);  % 中弦点力失的起点——>终点   % ——第一方案
    % y_mid=F_yadddown(i,1)+Y_axisdown(1,i);  % 中弦点力失的起点——>终点        % ——第二方案
    z_mid=z4+Z_axisdown(1,i);
    F_yadd_mid0=[y_mid;z_mid];  % size(wing_chord)   % (2*2)
    Fyadd_midrot=pn2*R_w2s*F_yadd_mid0;        % wing_chord in body frame  % (2*2)
    Fyadd_mid=Fyadd_midrot; 
    quiver(p_mid(1,1),p_mid(2,1),...                                                      % 起点
               Fyadd_mid(1,1),Fyadd_mid(2,1),0.1,'g-','LineWidth',2);   % 终点   % 比例因子0.05——————注意力失被缩小了0.135倍
    %%%%%%%%%%%%%%%%%%%%%%%%%      
    % set(gca,'Color','w')   %XXX
    % axes('Color','w');    %XXX
    % axes('Ycolor','w','Xcolor','w');  %XXX
    % axes('YTicklabel',[],'XTicklabel',[]); %XXX
    set(gca,'YTick',[],'XTick',[],'Ycolor','w','Xcolor','w'); 
    % box off
    box on
    % hold off
    hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 动画片段程序结束段——(2)
%     if animation == 1
%         frame = getframe(fig);
%         aviobj = addframe(aviobj,frame);
%     else
%         pause(0.5)
%     end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
% close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  动画片段程序结束段——(3)
% if animation == 1
%     aviobj = close(aviobj);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (2) 上冲程
k_index2=(Nstep:12:num);     % 169...181...193......337
N_up=length(k_index2);    % (15)
% t_up=t1(k_index2);    %  size(t_up)  %  (15*1)
psi_up=psi1(k_index2);  % 针对一个周期之内的数据，每隔个8个点取一个数据形成向量
alpha_up=alpha1(k_index2);
% figure(6)
% plot(t_up*f,psi_up*180/pi,'rd','LineWidth',2)      % 转换为度数
% Forcenormalup=Forcenormal(k_index2);
F_ytranup=F_ytran(k_index2);  % size(F_ytranup)   % (15*1)     %  10^(-5)*
F_yaddup=F_yadd(k_index2);   % size(F_yaddup)   % (15*1)
F_yrotup=F_yrot(k_index2);
F_inert_yup=F_inert_y(k_index2); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % animation———————— 动画片段程序开始段——(1)
% Nstep=N_up;    % N_up=15; 
% np=1;                     %  调节总的帧数和时间步长
% Nframe=np*15;   % number of movie frames -------------------------- for change
% Nskip2=floor(Nstep/Nframe);
% Nspeed=3;  % number of frames per second ------------------------ for change
% animation=1;       %  1=save animation file --------------------------- for change
% duration=Nframe/Nspeed;  % duration of movie  % 80 seconds @Nframe=4*40
% display(['duratin of movie will be ', num2str(duration), ' seconds ']); %输出显示动画持续的时间(s)
% %输出如下: duratin of movie will be 10 seconds 
% if animation == 1
%     aviobj = avifile('Stick_figure_force_upstroke.avi','fps',Nspeed); % filename --- for change
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % fig=figure;
fig=figure(2);
for i = 1:1:N_up  % N_up=15
    displacement=[Y_axisup(1,i);Z_axisup(1,i)];  % (2*1)
    psi=psi_up(i,1);             % 扭转角——单位是rad    %——注意为正号
    % psi_deg=psi_up(i,1)*180/pi  
     if (psi>0)
         % 下面一段是坐标系之间的变换
        Pitch = [cos(psi)   sin(psi); ...
                     -sin(psi)    cos(psi)];        % pitch matrix——符号可能需要方向
        R_w2s =Pitch;     % from wing frame to spar frame——从翅坐标系到体坐标系
        % (1) 前缘点坐标——实心圆点绘制
         y1=0;
         z1=C_025;
         p_lead0=[y1;z1];                   % (2*1)
         p_leadrot=R_w2s*p_lead0;    % (2*2)*(2*1)=(2*1)   % [-sign(psi),-sign(psi)]*
         p_lead=p_leadrot+displacement;          % (2*1)
        % (2) 后缘点坐标——实心圆点绘制
         y5=0;
         z5=-C_075;
         p_tail0=[y5;z5];                   % (2*1)
         p_tailrot=R_w2s*p_tail0;      % (2*2)*(2*1)=(2*1)
         p_tail=p_tailrot+displacement;          % (2*1)   
         % (3) 弦向压心点————实心圆点绘制 
         % alpha_up=alpha1(i,1);             % 扭转角——单位是rad
         d_cp=(0.82*abs(alpha_down(i,1))/pi+0.05)*C_avereff;        
         if d_cp>C_025
             y3=0;
             z3=-(d_cp-C_025);% 负值
         elseif d_cp<=C_025;
             y3=0;
             z3=C_025-d_cp; % 正值
         end
         p_cop0=[y3;z3];
         p_coprot=R_w2s*p_cop0;
         p_cop=p_coprot+displacement; 
         % (4) 中弦点坐标—力失的起点—实心圆绘制
         y4=0;
         % z4=-C_025; % 中弦点到扭转轴的距离
         z4=-C_05;   % 中弦点到扭转轴的距离——修改了针对实际的特征弦长的中弦点，而非0.25*C_avereff的中弦点
         p_mid0=[y4;z4];
         p_midrot=R_w2s*p_mid0;
         p_mid=p_midrot+displacement; 
    else  % (psi<0)
       % 下面一段是坐标系之间的变换
        Pitch = [cos(psi)   sin(abs(psi)); ...
                     -sin(abs(psi))    cos(psi)];        % pitch matrix——符号可能需要方向
        R_w2s =Pitch;     % from wing frame to spar frame——从翅坐标系到体坐标系
        % (1) 前缘点坐标——实心圆点绘制
         y1=0;
         z1=C_025;
         p_lead0=[y1;z1];                   % (2*1)
         p_leadrot=R_w2s*p_lead0;    % (2*2)*(2*1)=(2*1) 
         pn1=[sign(psi),0;
                    0,           1];
         p_lead=pn1*p_leadrot+displacement;          % (2*1)
        % (2) 后缘点坐标——实心圆点绘制
         y5=0;
         z5=C_075;
         p_tail0=[y5;z5];                   % (2*1)
         p_tailrot=R_w2s*p_tail0;      % (2*2)*(2*1)=(2*1)
         pn2=[1,0;
                  0,sign(psi)];
         p_tail=pn2*p_tailrot+displacement;          % (2*1)   
         % (3) 弦向压心点————实心圆点绘制 
         % alpha_up=alpha1(i,1);             % 扭转角——单位是rad
         d_cp=(0.82*abs(alpha_down(i,1))/pi+0.05)*C_avereff;        
         if d_cp>C_025
             y3=0;
             z3=(d_cp-C_025);% 负值
             p_cop0=[y3;z3];
             p_coprot=pn2*R_w2s*p_cop0;
         elseif d_cp<=C_025;
             y3=0;
             z3=C_025-d_cp; % 正值
             p_cop0=[y3;z3];
             p_coprot=pn1*R_w2s*p_cop0;
         end
         p_cop=p_coprot+displacement; 
         % (4) 中弦点坐标—力失的起点—实心圆绘制
         y4=0;
         % z4=-C_025; % 中弦点到扭转轴的距离
         z4=C_05;   % 中弦点到扭转轴的距离——修改了针对实际的特征弦长的中弦点，而非0.25*C_avereff的中弦点
         p_mid0=[y4;z4];
         p_midrot=pn2*R_w2s*p_mid0;
         p_mid=p_midrot+displacement; 
    end   
    % 片元弦向的特殊点阵
%     y=[p_lead(1,1),Y_axisup(1,i),p_cop(1,1),p_mid(1,1),p_tail(1,1)];
%     z=[p_lead(2,1),Z_axisup(1,i),p_cop(2,1),p_mid(2,1),p_tail(2,1)];
%     wing_C=[y;z];  % size(wing_chord)   % (3*5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % figure(7)
    %  figure(i)
    plot([Y_axisup(1,1),Y_axisup(1,N_up)],[Z_axisup(1,1),Z_axisup(1,N_up)],'r:','LineWidth',1.5) % 扭转轴心——实心圆点绘制
    hold on
%     % h1=text(-2.4,2.32,'Transient normal force vectors generated during one flapping stroke');
%     h1=text(-2.4,3.42,'Transient normal force vectors generated during one flapping stroke');
%     set(h1,'Color','k','FontSize',8,'FontName','Times','FontWeight','Bold')
%     % h2=text(-1.0,2.2,'for optimized wing geometry');
%     h2=text(-3.0,3.25,'for combined optimal WGP and WKP with 2D translational aerodynamic coefficients');
%     set(h2,'Color','k','FontSize',8,'FontName','Times','FontWeight','Bold')
%     h3=text(-0.25,-0.8,'Xijun Ke @ SJTU');
%     set(h3,'Color','k','FontSize',7.0,'FontName','Times','FontWeight','Bold')
%     h4=text(-0.15,-0.915,'17.03.2016');
%     set(h4,'Color','k','FontSize',7.0,'FontName','Times','FontWeight','Bold')
%     hold on
    % xlabel('纵向(y)')
    % ylabel('垂直方向向(z)')
    % grid on
    % % axis([-4,4,-1,3.1])   % 比例不合适于法向气动力的表达或者显示
    % axis([-4,4,-1.9,2.1])    % axis([-4,4,-0.6808,1.4952]) 
    % axis([-3.5,3.5,-1,2]) 
    % axis([-3.5,4.25,-1,2.5])
    axis([-3.5,4.25,-1,3.5])
    hold on
    % (1) 绘制直线——连接前后缘
    plot([p_lead(1,1),p_tail(1,1)],[p_lead(2,1),p_tail(2,1)],'k-','LineWidth',2) % 绘制前缘后缘实线
    hold on
    % (2) 扭转轴心——实心圆点绘制
    plot(Y_axisup(1,i),Z_axisup(1,i),'ro','LineWidth',1.5,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',3) % 扭转轴心——实心圆点绘制
    hold on
     % (3) 绘制前缘实心圆点
    plot(p_lead(1,1),p_lead(2,1),'ko-','LineWidth',1.5,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5)
    hold on
   % (4) 弦向压心点—绘制力失
    plot(p_cop(1,1),p_cop(2,1),'bo','LineWidth',1.7,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',5);  % 弦向压心点————实心圆点绘制
    hold on
    % p_cop=p_coprot+displacement;                 % 弦向压心点力失的起点
    % displacement=[Y_axisup(1,i);Z_axisup(1,i)];          % (2*1)
    % y_cop=F_ytranup(i,1)+Y_axisup(1,i);                     %——注意符号——仅含有平动法向力分量
    y_cop=(F_ytranup(i,1)+F_yrotup(i,1))+Y_axisup(1,i); %——注意符号——含有平动法向力分量和转动环量法向力分量
    % y_cop=(F_ytranup(i,1)+F_yrotup(i,1))+Y_axisup(1,i); %——注意符号——含有平动法向力分量和转动环量法向力分量
    z_cop=z3+Z_axisup(1,i);
    F_ytran_cop0=[y_cop;z_cop];                           % size(wing_chord)   % (2*2)    % 弦向压心点力失的终点
    Fytran_coprot= R_w2s*F_ytran_cop0;              % wing_chord in body frame  % (2*2)     %—R_w2s'—注意符号
    Fytran_cop=Fytran_coprot; 
    quiver(p_cop(1,1),p_cop(2,1),...                        % 起点
               Fytran_cop(1,1),Fytran_cop(2,1),0.1,'b-','LineWidth',2);    % 终点   % 比例因子0.05——————注意力失被缩小了0.14倍 
    hold on
   % (5) 中弦点—绘制力失
   % plot(y4,z4,'go','LineWidth',1.5,'MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',3);  % 中弦点—力失的起点—实心圆绘制
   % hold on  
    y_mid=F_yaddup(i,1)+F_inert_yup(i,1)+Y_axisup(1,i);  % 中弦点力失的起点——>终点   % ——第一方案
    % y_mid=F_yaddup(i,1)+Y_axisup(1,i);  % 中弦点力失的起点——>终点        % ——第二方案
    z_mid=z4+Z_axisup(1,i);
    F_yadd_mid0=[y_mid;z_mid];  % size(wing_chord)   % (2*2)
    Fyadd_midrot=R_w2s*F_yadd_mid0;        % wing_chord in body frame  % (2*2)
    Fyadd_mid=Fyadd_midrot; 
    quiver(p_mid(1,1),p_mid(2,1),...                                                      % 起点
               Fyadd_mid(1,1),Fyadd_mid(2,1),0.1,'g-','LineWidth',2);   % 终点   % 比例因子0.05——————注意力失被缩小了0.135倍
    set(gca,'YTick',[],'XTick',[],'Ycolor','w','Xcolor','w');       
    % box off
    box on    
    % hold off
    hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % 动画片段程序结束段——(2)
%     if animation == 1
%         frame = getframe(fig);
%         aviobj = addframe(aviobj,frame);
%     else
%         pause(0.5)
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
% close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %  动画片段程序结束段——(3)
% if animation == 1
%     aviobj = close(aviobj);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
