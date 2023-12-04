% function obj_function=objfunction_aeroF_stdev(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 原始文件名——Aero_F3_fruitfly_exp
% objfunction_aeroF_stdev
% (1)涉及求解准稳态理论预测而得的虫体坐标系的气动力和实验测试机械果蝇模型获得的数据之差的标准方差
% (2) 可以考虑翅膀惯性力和不考虑翅膀惯性力，分别进行优化
% (3) 涉及针对四个气动力系数的惩罚项的约束;——F_N=F_ytran+F_yrot+F_yadd1+F_inert_y;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tic
% x=[1,1,0.356737,1,1];
% x=[1,1,1,1,0];  % obj_function =40.1316;
% k_C_N=x(1);
% k_C_R=x(2);
% x0_nd=x(3);  % 扭转轴的位置 x0_nd=0.356737;
% k_add=x(4);
% k_inert=x(5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;clc;
% % x=[1,1.281,1,1];
x=[1,1.8,0.35,1]; % 比x=[1,1.8,0.15,1]; 要好点哦
% x=[0.85,0.75,0.05,1]; 
% x =[1,0.5,0.571,0.01];  % obj_function =36.6581;
% x =[0.5,0.5,0.2000];  %fval =18.7559;
% x =[0.5, 0.5,1.0889]; % fval = 12.6197;
% x =[1,0.5,0.5,0.01];
% x =[1,0.5,0.5,0.5,0];
% x =[1.0000,1.0000,0.7500,0.5643,1.5000];  % fval =35.5047;  % Elapsed time is 583.787292 seconds. %F_N没有考虑惯性力，k_inert的随机结果1.5，与实测结果不符合
% x =[1.0000,1.6003,0.7500,0.5574,0.0100]; % fval =35.5099;
% x= [1.0000,0.5000,0.3634, 0.5000];  % fval =36.8167;
% x =[1.0000,0.5000,0.3358,0.5000];       % fval=36.5981; 
% x =[1.0000,1.6003,0.7500,1]; 
k_C_N=x(1);
k_C_R=x(2);
% x0_nd=x(3);  % 扭转轴的位置 
x0_nd=0.356737; % 扭转轴的位置 
k_add=x(3);
k_inert=x(4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 气动力: Solution of the aerodynamic force
% 修改时间——2014年12月20日,11:54——转动轴偏离弦向中点的偏移量坐标,在扭转轴向上偏移C_maxy之后
% clear all;clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 调用含翅形貌参数化和气动力系数的函数
% wing_para_output=zeros(1,16);
% wing_para_output=[F_ndTrans,Coeff_liftdragF_N,M_xaercoeff,I1y,Z_rnd, M_xrdcoeff,...
%     F_ndRot,F_yrotcoeff,M_xRotcoeff, I2y,...
%     I_xzam,I_xxam, I5z,I6z,I7y,M_zrdcoeff];
wing_para=wing_shape_fruitfly_sixteen_good();   %调用函数wing_shape_fruitfly;  % size(wing_para)
% F_ndTrans=wing_para(1,1);             % F_ndTrans=0.46392;  无量纲, 量纲化单位为mm^4
Coeff_liftdragF_N=wing_para(1,2);      % Coeff_liftdragF_N=0.00682;  %单位是mg*mm
% M_xaercoeff=wing_para(1,3);          % M_xaercoeff=0.006038;   %单位是: mg.mm^2
% I1z=wing_para(1,4);                         % I1y=0.016158   % 单位是 mg.mm^2
% Z_rnd=wing_para(1,5);                     % Z_rnd=0.16265;  无量纲, 量纲化单位为mm
% M_xrdcoeff=wing_para(1,6);           % M_xrdcoeff=0.0001839; % 单位是mg.mm^2
% F_ndRot=wing_para(1,7);                % F_ndRot=0.74851;  无量纲, 量纲化单位为mm^4
F_yrotcoeff=wing_para(1,8);               % F_yrotcoeff =0.0032433;  % 单位是 mg.mm
% M_xRotcoeff=wing_para(1,9);         % M_xRotcoeff=0.002871;   % 单位是 mg.mm^2
% I2z=wing_para(1,10);                      % I2y=0.006943;        % 单位是 mg.mm^2
% I_xzam=wing_para(1,11);                % I_xzam=0.001424  % 单位是 mg.mm^2
% I_xxam=wing_para(1,12);                % I_xxam=0.000338  % 单位是 mg.mm^2
I5y=wing_para(1,13);                          % I5z=0.0050926   % 单位是 mg.mm
I6y=wing_para(1,14);                          % I6z=0.00077164  % 单位是 mg.mm
% I7z=wing_para(1,15);                      % I7y=0.0109056;        % 单位是 mg.mm^2
% M_zrdcoeff=wing_para(1,16);         % M_zrdcoeff=0.0011169; % 单位是 mg.mm % 转动气动阻尼力矩参数—绕翅平面下的弦向轴
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I5y=wing_para(1,13);                       % I5z=0.0050945   % 单位是 mg.mm        % 下文中I3y应该改为I5y
% I6y=wing_para(1,14);                       % I6z=0.0011         % 单位是 mg.mm         % 下文中I4y应该改为I6y
% I7z=wing_para(1,15);                       % I7y=0.0109;        % 单位是 mg.mm^2    % 下文中I5z应该改为I7z
% I1z=wing_para(1,4);                         % I1y=0.0162;        % 单位是 mg.mm^2    % 下文中I7z应该改为I1z
% I2z=wing_para(1,10);                       % I2y=0.0069;        % 单位是 mg.mm^2    % 下文中I6z应该改为I2z 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 调用翅膀运动学-翅运动规律和几何攻角(AOA)等数据
% wing_m_output=[t',phi',psi',alpha',dphi',dpsi',ddphi',ddpsi',C_L',C_D',C_N1',C_T'];
% wing_kenimatics=kenimatics_wing_and_AoA_fruitfly_sim(); %调用函数kenimatics_wing_and_AoA;  % size(wing_kenimatics)  % (1000,11)
wing_kenimatics=kenimatics_wing_and_AoA_fruitfly_exp();
t=wing_kenimatics(:,1);                % 单位是ms
phi=wing_kenimatics(:,2);        % 拍打角——单位是rad
psi=wing_kenimatics(:,3);            % 拍打角——单位是rad
alpha=wing_kenimatics(:,4);      % alpha=atan2(omega_z,-omega_y);  % 几何攻角——弧度制   %输出——有正有负
dphi=wing_kenimatics(:,5);          % 单位是rad/s
dpsi=wing_kenimatics(:,6);          % 单位是rad/s
ddphi=wing_kenimatics(:,7);       % 单位是rad/s^2
ddpsi=wing_kenimatics(:,8);     % 单位是rad/s^2
C_L=wing_kenimatics(:,9);          
C_D=wing_kenimatics(:,10);     
C_N1=wing_kenimatics(:,11);   
C_T=wing_kenimatics(:,12);   
C_N=C_N1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 翅膀坐标系下的角速率和角加速率——————这组数据来自翅2DOF运动
% 翅膀坐标系下的角速度
omega_x=dpsi;                     % 展向
omega_y=dphi.*sin(psi);       % 法向(初始向左)
omega_z=dphi.*cos(psi);      % 弦向朝上
omega_h=dphi;       % 铰链的角速度% omega_h=-sqrt(omega_y^2+omega_z^2);  % 右手法则顺时针
% 翅膀坐标系下的角加速度——用于虚质量力的计算
domega_x=ddpsi;
domega_y=ddphi.*sin(psi)+dphi.*dpsi.*cos(psi);
domega_z=ddphi.*cos(psi)-dphi.*dpsi.*sin(psi); 
%% 气动攻角和当地流场速度的计算——取流产速度相对于刚体好了速度，所以整体加负号；
v_y_nonr=-omega_z;   % v_y=r*dphi*cos(psi)
v_z_nonr=omega_y;    % v_z=-r*dphi*sin(psi)
% alpha2=atan2(-v_y_nonr,v_z_nonr);   % 正确——注意与下文的alpha=atan2(omega_z,-omega_y)*180/pi; 不同
% % 由于alpha=atan2(cot(psi))=atan2(cot(pi/2-alpha))=atan2(tan(alpha)); %这里atan2给出象限正负值，尤其是alpha>pi/2时
V_nonr=sqrt(v_y_nonr.^2+v_z_nonr.^2); % 当地来流速度V_nonr=omega_h=dphi;   % 单位是 rad/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xr=0.3289;  % x-root offset  \mm      % R_wingeff=wing_para(1,1);  % R_wingeff =3.0040;   % 单位是 mm
% R_ref=(xr+R_wingeff)*10^(-3);  
% R_ref_non=3.3293*10^(-3);
R_ref_non=1;   % r2_nd=0.5801;   
v_y=R_ref_non*omega_z;
v_z=-R_ref_non*omega_y;
V_ref=sqrt(v_y.^2+v_z.^2);
f=188.7;  T=1/f;          % Hz
V_ref_aver=trapz(t,V_ref)/(3*T);    % 翼尖参考速度: V_ref_aver=867.0901rad/s;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 第一种气动力——翅坐标系下：平动气动环量力——升阻力分量；作用点：片元弦长中点
%% 输入参数为：4个
% C_avereff;  R_wingeff;    F_nd;                ——————这组数据来自翅形貌参数化;% F_nd=0.50236;
% omega_h;  C_L(alpha);  C_D(alpha);       ——————这组数据来自翅2DOF运动铰链角速率, 攻角
%%  Coeff_liftdragF_N   %单位  mg.mm=*10^(-9)kg.m
% g=9.821*10^(-6);   % 这里重力加速度:g=9.821N/kg=9.821*10^(-6)=9.821e-006   N/mg  ——g的国际单位是m*s^-2或N*kg^-1
g=9.821;   % 这里重力加速度:g=9.821N/kg=9.821*10^(-6)=9.821e-006   N/mg  ——g的国际单位是m*s^-2或N*kg^-1
% lift_inst=V_nonr.^2.*C_L*Coeff_liftdragF_N*10^(-3)/g;      % 单位是rad^2*s^-2*kg*m / (m*s^-2)=mg
% drag_inst=V_nonr.^2.*C_D*Coeff_liftdragF_N*10^(-3)/g;  % 单位是rad^2*s^-2*kg*m / (m*s^-2)=mg
lift_inst=V_nonr.^2.*C_L*Coeff_liftdragF_N*10^(-3);      % 单位是rad^2*s^-2*kg*m / (m*s^-2)=mg
drag_inst=V_nonr.^2.*C_D*Coeff_liftdragF_N*10^(-3);  % 单位是rad^2*s^-2*kg*m / (m*s^-2)=mg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%对瞬时气动升阻力使用梯形积分函数trapz进行数值积分，除以周期, 求解平均升阻力
lift_aver=trapz(t,lift_inst)/(3*T);                            % lift_aver =1.0913 
drag_aver=trapz(t,drag_inst)/(3*T);                      % drag_aver =0.9887
drag_instabs_aver=trapz(t,abs(drag_inst))/(3*T);        
lift2drag_aver=lift_aver/drag_instabs_aver;         % 升力和阻力比值: lift2drag_aver =1.1037——不是升阻系数的比值，没有实际意义？
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(10)
% F_LD=plot(t/T,lift_inst,'r-',t/T,drag_inst,'b:','LineWidth',2);   
% xlabel('\itNormalized time')
% ylabel('\itAerodynamic force (mg)')                    %单位换算到mg了哦
% title('Aerodynamic lift and drag force \itvs. \itt \rm for flapping wing')
% % legend('\itinstantaneous lift (\itF_L)','\itinstantaneous drag (\it|F_D|)')
% legend('\itinstantaneous lift (\itF_L)','\itinstantaneous drag (\itF_D)')
% set(F_LD,'LineWidth',2)
% grid on
% axis([0.9,3.0,-inf,inf])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 翅坐标系下：平动气动环量力(1)——法向力；作用点：片元弦长中点
% C_N=C_N1=cos(abs(alpha)).*C_L+sin(abs(alpha)).*C_D; %由升阻力系数合成——2010-JFM-RJ Wood
% Coeff_liftdragF_N=6.8201e-012——单位是:kg/m^3*mm*mm^3= 10^(-12) kg.m
% F_ytran=-sign(alpha).*V_nonr.^2.*C_N*Coeff_liftdragF_N/g;  %单位是rad^2*s^-2*kg*m / (m*s^-2)=mg——仅用于显示
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(11)
% plot(t/T,F_ytran,'k-','LineWidth',2)
% xlabel('\itNormalized time')
% ylabel('\itF_{y,tran} (mg)')
% legend('F_{y,tran}')
% title('平动气动环量力(法向)随时间的变化规律') 
% grid on
% axis([0.9,4.05,-inf,inf])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 单位是rad^2*s^-2*kg*m / (m*s^-2)=mg
F_ytran=-k_C_N*sign(alpha).*V_nonr.^2.*abs(C_N)*Coeff_liftdragF_N*10^(-3);   % 单位是rad^2*s^-2*mg*mm=10^(-9)N=10^(-9)*10^6uN
F_ytran_aver=trapz(t,abs(F_ytran))/(3*T);   % F_ytran_aver =14.5887uN;
F_ztran=-sign(alpha).*V_nonr.^2.*C_T*Coeff_liftdragF_N*10^(-3);  
% F_ytran=V_nonr.^2.*C_N*Coeff_liftdragF_N*10^(-3); 
% F_ztran=V_nonr.^2.*C_T*Coeff_liftdragF_N*10^(-3);  
% figure(12)
% plot(t/T,F_ytran,'k-',t/T,F_ztran,'r-','LineWidth',2)
% xlabel('\itNormalized time')
% ylabel('\itF_{y,tran} & F_{z,tran} (uN)')
% legend('F_{y,tran}','F_{z,tran}')
% title('平动气动环量力(法向和切向)随时间的变化规律') 
% grid on
% axis([0.9,4.05,-inf,inf])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F_vert_tran=-sign(alpha).*F_ytran.*cos(alpha);              % 垂直方向vertical, %单位是uN  %  
F_horiz_tran=F_ytran.*sin(abs(alpha));                           % 水平方向horizontal,%单位是uN
% figure(13)                  % 该气动力有些小  
% hold on
% plot(t/T,F_vert_tran,'k-',t/T,F_horiz_tran,'r-','LineWidth',2)
% xlabel('\itNormalized time')
% ylabel('\itF_{vert,tran} &  F_{horiz,tran} (uN)')
% legend('\itF_{vert,tran}','\itF_{horiz,tran}')
% title('旋转气动环量力的垂直方向和水平方向的分量随时间的变化规律')   
% grid on
% axis([0.9,4.05,-inf,inf])
% set(gca,'XTick',(0.9:0.1:4.05))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 第二种气动力—— 翅坐标系下： 旋转气动环量力：法向力；作用点：片元弦长中点
%% 输入参数为：3个
% C_avereff;   R_wingeff;   F_ndRot;  ——————这组数据来自翅形貌参数化
% omega_h;  omega_x;                    ——————这组数据来自翅2DOF运动铰链角速率, 翅角速率
% F_zrot=(1/2)*Rou*R_wingeff^2*C_avereff^2*C_R*omega_x.*omega_h*F_ndRot*(10^(-12)*10^3);%单位是mN——未考虑方向
% F_y_rot=C_R*sign(alpha2).*abs(omega_x).*abs(omega_h)*F_y_rotcoeff*10^6;    % 单位是uN
% x0_nd=0.356737;
C_R=pi*(0.75-x0_nd);
% C_R=1;  
% C_R=1.55;
% F_yrot=-C_R*sign(alpha).*omega_x.*V_nonr*F_yrotcoeff*10^(-3); % 
F_yrot=k_C_R*C_R*omega_x.*V_nonr*F_yrotcoeff*10^(-3);      % 单位是rad^2*s^-2*kg*m=10^(-9)N=10^(-9)*10^6uN
F_yrot_aver=trapz(t,abs(F_yrot))/(3*T);                             % F_yrot_aver =2.8944uN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(14)                  % 该气动力有些小  
% hold on
% plot(t/T,F_ytran,'k-',t/T,F_yrot,'r-','LineWidth',2)
% xlabel('\itNormalized time')
% ylabel('\itF_{y,tran} &  F_{y,rot} (uN)')
% legend('\itF_{y,tran}','\itF_{y,rot}')
% title('平动气动环量力(法向)和旋转气动环量力(法向)随时间的变化规律')   
% grid on
% axis([0.9,4.05,-inf,inf])
% set(gca,'XTick',(0.9:0.1:4.05))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F_vert_rot=-sign(alpha).*F_yrot.*cos(alpha);              % 垂直方向vertical, %单位是uN  %  
F_horiz_rot=F_yrot.*sin(abs(alpha));                          % 水平方向horizontal,%单位是uN
% figure(15)                  % 该气动力有些小  
% hold on
% plot(t/T,F_vert_rot,'k-',t/T,F_horiz_rot,'r-','LineWidth',2)
% xlabel('\itNormalized time')
% ylabel('\itF_{vert,rot} &  F_{horiz,rot} (uN)')
% legend('\itF_{vert,rot}','\itF_{horiz,rot}')
% title('旋转气动环量力的垂直方向和水平方向的分量随时间的变化规律')   
% grid on
% axis([0.9,4.05,-inf,inf])
% set(gca,'XTick',(0.9:0.1:4.05))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 第三种气动力——翅坐标系下：虚质量气动力：法向力；作用点：片元弦长中点
%%% 拍打和被动扭转的协同作用？——虚质量力——参考SK Agrawal的文献推导过程，以便理解其物理意义
%% 输入参数为：10个
% I5y;   I6y;     % 单位是 mg.mm=10^(-9) kg.m ——————这组数据来自翅形貌参数化
% domega_z; omega_x; omega_y;                      ——————这组数据来自翅2DOF运动翅角加速率, 角速率
% dphi; dalpha; ddphi; ddalpha; alpha;               ——————这组数据来自翅2DOF运动角速率, 角加速率, 攻角
%%  (1) 2010_JFM_Aeromechanics-虚质量气动力公式——% I5y和I6y单位是 mg.mm=*10^(-9) kg.m
% F_yadd=sign(phi).*(I5y*(domega_z-omega_x.*omega_y)+I6y*domega_x)*(10^(-9)*10^6);   %单位是uN
% F_yadd1=(I5y*(domega_z-omega_x.*omega_y)+I6y*domega_x)*10^(-3); % 原文推出的公式——更为合理些哦
% 注意公式前(-)符——根据果蝇翅膀坐标系修改了角速度和角加速度
% F_yadd1=-(I5y*(domega_z+omega_x.*omega_y)-I6y*domega_x/4)*10^(-3); %单位是rad^2*s^-2*kg*m=N=10^6uN % 与2001-JEB匹配
% F_yadd1=-(I5y*(domega_z+omega_x.*omega_y)+I6y*domega_x)*10^(-3); % 修改了原文推出的公式的正负号,且不含/4  % ——第一方案
% F_yadd1=k_add*(-I5y*(domega_z+omega_x.*omega_y)+I6y*domega_x)*10^(-3);    % 修改了原文推出的公式的正负号,且不含/4  % ——第二方案
F_yadd1=k_add*(I5y*(domega_z+omega_x.*omega_y)+I6y*domega_x)*10^(-3);  
% domega_z=ddphi.*cos(psi)-dphi.*dpsi.*sin(psi); 
% F_yadd1=k_add*(I5y*(ddphi.*cos(psi)+dphi.*dpsi.*sin(psi))+I6y*domega_x)*10^(-3);  
F_yadd1_aver=trapz(t,abs(F_yadd1))/(3*T);       % F_yadd1_aver=6.0929uN(old);  F_yadd1_aver=4.5472;
% (2)下面是2001-JEB-Sanjay SP-虚质量气动力公式
% F_yadd2=(I5y*(ddphi.*sin(alpha)+dphi.*dalpha.*cos(alpha))-(I6y/4)*ddalpha/4)*10^(-3);% 原文推出的公式——不动
% 下面的alpha1为全正几何攻角，有正有负的dalpha1和ddalpha1
% F_yadd2=-(I5y*(ddphi.*sin(alpha1)+dphi.*dalpha1.*cos(alpha1))-(I6y/4)*ddalpha1/4)*10^(-3); 
% F_yadd2=-(I5y*(ddphi.*sin(alpha1)+dphi.*dalpha1.*cos(alpha1))+(I6y/4)*ddalpha1/4)*10^(-3); 
% F_yadd2=-(I5y*(ddphi.*sin(alpha1)+dphi.*dpsi.*cos(alpha1))-(I6y/4)*ddpsi/4)*10^(-3); %这里的dpsi和ddpsi方向变化需要确定
% F_yadd2=-(I5y*(ddphi.*sin(alpha1)+dphi.*dpsi.*cos(alpha1))+(I6y/4)*ddpsi/4)*10^(-3); %这里的dpsi和ddpsi方向变化需要确定
% 下面的alpha2有正有负, dalpha和ddalpha需要求导
% F_yadd2=-(I5y*(ddphi.*sin(alpha)+dphi.*dalpha.*cos(alpha))-(I6y/4)*ddalpha/4)*10^(-3); 
i=length(alpha);   % size(alpha)=(2000*1)
dalpha=[diff(alpha);alpha(i,1)];    
j=length(dalpha);
ddalpha=[diff(dalpha);dalpha(j,1)];  % size(ddalpha)
% F_yadd2=-(I5y*(ddphi.*sin(abs(alpha))+dphi.*dalpha.*cos(abs(alpha)))-(I6y/4)*ddalpha/4)*10^(-3);% ——第一方案
F_yadd2=(I5y*(ddphi.*sin(abs(alpha))+dphi.*dalpha.*cos(abs(alpha)))-(I6y/4)*ddalpha/4)*10^(-3);
%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(16)  % 单位是 mg.mm*rad^2.*s^-2=10^(-9) kg.m*s^-2=10^(-9)N=10^(6)uN;   %单位是uN
% F_yadd1_tran=-(I5y*(domega_z+omega_x.*omega_y))*10^(-3);  % 法向虚质量气动力平动分量——翅膀
% F_yadd1_rot=(I6y*domega_x)*10^(-3);                                         % 法向虚质量气动力转动分量——翅膀
% plot(t/T,F_yadd1_tran,'r-.',t/T,F_yadd1_rot,'g:',t/T,F_yadd1,'k-','LineWidth',2) 
% ylabel('法向气动力分量F_{norm}(t) (uN)'); 
% legend('F_{y,add1,tran}(t)','F_{y,add1,rot}(t)','F_{y,add1}(t)');  
% title('法向虚质量气动力分量F_{y,add1}(t))随时间的变化')  
% grid on
% axis([0.9,4.05,-inf,inf])
% set(gca,'XTick',(0.9:0.1:4.05))
%%%%%%%%%%%%%%%%%%%%%%%%%
F_vert_add=-sign(alpha).*F_yadd1.*cos(alpha);              % 垂直方向vertical, %单位是uN  %  
F_horiz_add=F_yadd1.*sin(abs(alpha));                          % 水平方向horizontal,%单位是uN
% figure(17) 
% plot(t/T,F_vert_add,'r-.',t/T,F_horiz_add,'g:','LineWidth',2) 
% ylabel('法向气动力分量F_{norm}(t) (uN)'); 
% legend('F_{vert,add}(t)','F_{horiz,add}(t)');  
% title('法向虚质量气动力的垂直方向和水平方向的分量随时间的变化')  
% grid on
% axis([0.9,4.05,-inf,inf])
% set(gca,'XTick',(0.9:0.1:4.05))
%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(18)
% hold on
% plot(t/T,F_yadd1,'r-',t/T,F_yadd2,'b-','LineWidth',2)    % 注意这里的domega_x=ddpsi与ddalpha的符号存在不同哦
% xlabel('\itNormalized time')
% ylabel('\itF_{y,add} (uN)')
% legend('F_{y,add1}','F_{y,add2}')
% title('虚质量气动力(法向)随时间的变化规律')  
% grid on
% axis([0.9,4.05,-inf,inf])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 翅膀自身惯性力分量
% % 翅膀坐标系下的角速度
% omega_x=dpsi;                     % 展向
% omega_y=dphi.*sin(psi);       % 法向(初始向左)
% omega_z=dphi.*cos(psi);      % 弦向朝上
% omega_h=dphi;       % 铰链的角速度% omega_h=-sqrt(omega_y^2+omega_z^2);  % 右手法则顺时针
% % 翅膀坐标系下的角加速度——用于虚质量力的计算
% domega_x=ddpsi;
% % domega_y=ddphi.*sin(psi)+dphi.*dpsi.*cos(psi);
% domega_z=ddphi.*cos(psi)-dphi.*dpsi.*sin(psi); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% m_wing=2.4*10^-9;                    % kg
% x_com=1.920243385*10^-3;      %  m                           % 质心的展向坐标
% % z_com=-0.149785466+0.636=0.486215*10^-3;      % 到扭转轴的弦向距离
% z_com=0.149785466*10^-3;       %  m                          % 到扭转轴的弦向距离
% F_inert_y=-k_inert*m_wing*(domega_z*x_com-domega_x*z_com+omega_y.*(omega_x*x_com+omega_z*z_com))*10^6;  % 翅膀自身惯性力—法向
% % F_inert_y=-m_wing*(domega_z*x_com-domega_x*z_com+omega_y.*(omega_x*x_com+omega_z*z_com))*10^6; 
% F_inert_z=-k_inert*m_wing*(-domega_y*x_com-omega_y.^2*z_com+omega_x.*(omega_z*x_com-omega_x*z_com))*10^6; % 翅膀自身惯性力—弦向
% % F_inert_y=-m_wing*(ddphi.*cos(psi)*x_com-(ddpsi-dphi.^2.*sin(2*psi)/2)*z_com)*10^6;       % 翅膀自身惯性力—法向
% % F_inert_z=-m_wing*(-ddphi.*sin(psi)*x_com-(dpsi.^2+dphi.^2.*(sin(psi)).^2)*z_com)*10^6;  % 翅膀自身惯性力—弦向
% figure(19) 
% plot(t/T,F_inert_y,'r-',t/T,F_inert_z,'g-','LineWidth',2) 
% ylabel('翅膀自身惯性气动力分量F_{inert,y}(t) & F_{inert,z}(t) (uN)'); 
% legend('F_{inert,y}(t)','F_{inert,z}(t)');  
% title('翅膀自身惯性气动力分量F_{inert,y}(t) & F_{inert,z}(t)随时间的变化')  
% grid on
% axis([0.9,4.05,-inf,inf])
% set(gca,'XTick',(0.9:0.1:4.05))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F_vert_inert=-sign(alpha).*F_inert_y.*cos(alpha);              % 垂直方向vertical, %单位是uN  %  
% F_horiz_inert=F_inert_y.*sin(abs(alpha));                          % 水平方向horizontal,%单位是uN
% figure(20) 
% plot(t/T,F_vert_inert,'r-',t/T,F_horiz_inert,'g-','LineWidth',2) 
% ylabel('翅膀法向自身惯性气动力的垂直方向F_{vert,inert}(t)和水平方向F_{horiz,inert}(t)的分量  (uN)'); 
% legend('F_{vert,inert}(t)','F_{horiz,inert}(t)');  
% title('翅膀自身法向惯性气动力的垂直方向和水平方向的分量随时间的变化')  
% grid on
% axis([0.9,4.05,-inf,inf])
% set(gca,'XTick',(0.9:0.1:4.05))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 总气动力—— 翅坐标系下：总(合成)法向气动力：平动气动环量力+旋转气动环量力+虚质量气动力；作用点：片元弦长中点
F_N=F_ytran+F_yrot+F_yadd1;    % 单位是rad^2*s^-2*kg*m=N=10^6uN  ——不考虑翅膀惯性力
% F_N=F_ytran+F_yrot+F_yadd1+F_inert_y;    % 单位是rad^2*s^-2*kg*m=N=10^6uN——考虑翅膀惯性力
% figure(21)
% plot(t/T,F_N,'k-','LineWidth',2)
% xlabel('\itNormalized time')
% ylabel('\itF_N (uN)')
% legend('F_N')
% title('平动气动环量力+旋转气动环量力+虚质量气动力(法向)随时间的变化规律')
% grid on
% axis([0.9,4.05,-inf,inf])
% hold on
% L=length(t);
% plot([0,t(L)/T],[0,0],'k-','LineWidth',2);     %画x-axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 翅坐标系下：片元法向合力分解到片元垂直方向和水平方向；作用点：片元弦长中点
F_vertical=-sign(alpha).*cos(alpha).*F_N;                % 垂直方向vertical, %单位是uN  %——实际上该式计算更快，而且与实验结果对比性方便
% F_vertical=sign(alpha).*sin(psi).*F_N;                   % 垂直方向vertical, %单位是uN  % ——根据公式推导该采用本式
F_horizontal=sin(abs(alpha)).*F_N;                          % 水平方向horizontal,%单位是uN
% % F_horizontal=-sign(alpha).*F_N.*sin(alpha);     % 水平方向horizontal,%单位是uN
% F_vertical=abs(F_N.*cos(alpha));                          % 垂直方向vertical, %单位是uN
% F_horizontal=abs(sign(alpha).*F_N.*sin(alpha));  % 水平方向horizontal,%单位是uN
%% 虫体坐标系下： 侧边力(side force)F_horiz_X和水平推阻力(drag)=F_horiz_y
% F_Z=F_vertical;
% F_horiz_X=-sin(phi).*F_horizontal;
% F_horiz_Y=cos(phi).*F_horizontal;  
F_horiz_X=-sin(phi).*cos(psi).*F_N; % 侧边力(side force)F_horiz_X
F_horiz_Y=cos(phi).*cos(psi).*F_N;  % 水平推阻力(drag)=F_horiz_y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g_uN=9.821;   % 这里重力加速度:g=9.821N/kg=9.821*10^6/10^6=9.821 uN/mg  ——g的国际单位是m*s^-2或N*kg^-1
F_vaver=trapz(t,F_vertical)/(3*T)/g_uN;               % F_vaver =1.2532; %单位是mg      1.2532(含惯性力)
F_haver=trapz(t,F_horiz_Y)/(3*T)/g_uN;          % F_haver =-0.1109; %单位是mg     -0.1109(含惯性力)
F_haverabs=trapz(t,abs(F_horiz_Y))/(3*T);    
F_v2haver=F_vaver/F_haverabs;                % 升力和阻力比值: F_v2haver =0.1150;  0.1045(含惯性力)——不是升阻系数的比值，没有实际意义？
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=9.81;         % 这里重力加速度:g=9.821N/kg=9.821*10^6/10^6=9.821 uN/mg  ——g的国际单位是m*s^-2或N*kg^-1
M_body =1.8;  %果蝇虫体的质量(mg)—Science——M_body =1.8e-06;(kg)
W=M_body*g;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 读入science机械果蝇翅膀测得气动力数据——归一化之后的数据
% %%下面的读入数据有误%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % output_1=[t_NOfreq,Fx_norm_stroke_butterfilt_mean,Fy_norm_stroke_butterfilt_mean,Fz_norm_stroke_butterfilt_mean];
% force_science=xlsread('D:\KXJ\PassiveRot_dynamic_Science_fruitfly\aerodynamic_moment\aerodynamic_moment_model2_2\force_science.xlsx','A1:D1145'); % 读入数据
%下面读入数据正确%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% output_1=[t_NOfreq,Fx_norm_all_butterfilt_steady,Fy_norm_all_butterfilt_steady,Fz_norm_all_butterfilt_steady];
% force_science=xlsread('D:\KXJ\PassiveRot_dynamic_Science_fruitfly\wing_parameter\datanalysis_science_fruitfly\robotForcesTorques\ForceModulations\aeroforce_for_steady_wingbeat.xlsx','A1:D1145');
force_science=xlsread('aeroforce_for_steady_wingbeat.xlsx',1,'A1:D1145'); % 读入数据
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_NOfreq1=force_science(:,1);  % t_NOfreq(1,1)=0;
Fx_norm1=-force_science(:,2);
Fy_norm1=force_science(:,3);
Fz_norm1=-force_science(:,4); 
t_NOfreq=[t_NOfreq1+0.0052824335;t_NOfreq1+T+0.0052824335;t_NOfreq1+2*T+0.0052824335];
Fx_norm=[Fx_norm1;Fx_norm1;Fx_norm1];
Fy_norm=[Fy_norm1;Fy_norm1;Fy_norm1];
Fz_norm=[Fz_norm1;Fz_norm1;Fz_norm1];       % size(Fz_norm)  % (3435*1)
% figure(22)
% plot(t_NOfreq/T,Fz_norm,'r--',t_NOfreq/T,Fx_norm,'b--','LineWidth',2)
% xlabel('\itNormalized time')
% ylabel('\itF_{z,exp} & F_{x,exp} (uN)')
% legend('F_{z,exp}','F_{x,exp}')
% title('虫体坐标下实验测得的垂直方向升力和水平方向的推力——2014-Science-MH Dickinson')   
% grid on
% axis([0.9,4.05,-inf,inf])
% set(gca,'XTick',(0.9:0.1:4.05))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 法向合力分解到垂直方向和水平方向的分量力随时间的变化规律与实验测试结果的对比
figure(44)
plot(t_NOfreq/T,Fz_norm,'r-',t_NOfreq/T,Fx_norm,'b-','LineWidth',2.5)  % 实验测试结果
hold on
% 1.525——考虑翅膀惯性力  % 1.65——不考虑翅膀惯性力——————严格意义上讲应该是两个翅膀产生的总作用力
% plot(t/T,2*F_vertical/W,'r-.',t/T,2*F_horiz_Y/W,'b-.','LineWidth',2)       % 理论预测结果
% plot(t/T,1.65*F_vertical/W,'r-.',t/T,1.65*F_horiz_Y/W,'b-.','LineWidth',2)  
% plot(t/T,1.5*F_vertical/W,'r-.',t/T,2.4*F_horiz_Y/W,'b-.','LineWidth',2)    % 用于与实验结果对比哦
% plot(t/T,1.45*F_vertical/W,'r-.',t/T,2.3*F_horiz_Y/W,'b-.','LineWidth',2.5)  
% plot(t/T,1.65*F_vertical/W,'r-.',t/T,2.3*F_horiz_Y/W,'b-.','LineWidth',2.5)  % 含惯性力 % 用于与实验结果对比哦
plot(t/T,1.45*F_vertical/W,'r-.',t/T,2.0*F_horiz_Y/W,'b-.','LineWidth',2.5)     % 不含惯性力 % 用于与实验结果对比哦 x=[1,1.8,0.35,1];
% plot(t/T,2*F_vertical/W,'r-.',t/T,2.0*F_horiz_Y/W,'b-.','LineWidth',2.5) 
% plot(t/T,1.425*F_vertical/W,'r-.',t/T,2.225*F_horiz_Y/W,'b-.','LineWidth',2.5)  
xlabel('\rmNormalized time','Fontsize',24,'FontName','Times','FontWeight','Bold')
% ylabel('\itF_{vertical} & F_{horiz,Y} (uN)')
% legend('F_{vertical}','F_{horiz,Y}')
ylabel('\rmF / m_{body}g','Fontsize',24,'FontName','Times','FontWeight','Bold')
legend('F_{vert,z,exp}','F_{horiz,y,exp}','F_{vert,z,cal}','F_{horiz,y,cal}')
box on
set(gca,'LineStyle','-','LineWidth',1.5,'FontSize',20,'FontName','Times','FontWeight','Bold')
% title('法向合力分解到垂直方向和水平方向的分量力随时间的变化规律与实验测试结果的对比')   
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%
axis([0.9,4.05,-4,4])
% set(gca,'XTick',(0.9:0.1:4.05))
hold on
L=length(t);
plot([0,t(L)/T],[0,0],'k-','LineWidth',2);     %画x-axis
axis([min(t/T),min(t/T)+1,-4,3])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% figure(45)
% plot(t_NOfreq/T,Fz_norm,'r-',t_NOfreq/T,Fx_norm,'b-',t_NOfreq/T,Fy_norm,'g-','LineWidth',2)  % 实验测试结果
% hold on
% % 1.525——考虑翅膀惯性力  % 1.65——不考虑翅膀惯性力——————严格意义上讲应该是两个翅膀产生的总作用力
% plot(t/T,2*F_vertical/W,'r-.',t/T,2*F_horiz_Y/W,'b-.',t/T,0*F_horiz_X/W,'g-.','LineWidth',2)       % 理论预测结果
% xlabel('\rmNormalized time','Fontsize',20,'FontName','Times','FontWeight','Bold')
% ylabel('\rmF / m_{body}g','Fontsize',20,'FontName','Times','FontWeight','Bold')
% legend('F_{vert,z,exp}','F_{horiz,y,exp}','F_{lateral,x,exp}','F_{vert,z,cal}','F_{horiz,y,cal}','F_{lateral,x,cal}')
% set(gca,'FontSize',12,'FontName','Times','FontWeight','Bold')
% axis([0.9,4.05,-4,4])
% hold on
% L=length(t);
% plot([0,t(L)/T],[0,0],'k-','LineWidth',2);     %画x-axis



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% 翅坐标系下： 由法向气动力系数求解: 垂直方向和水平方向气动力系数
% % Coeff_liftdragF_N=0.0068201——单位是:mg*mm^-3*mm*mm^3= 10^(-9) kg*m
% % Coeff_liftdragF_N=(1/2)*Rou*C_avereff*R_wingeff^3*F_nd*10^(-12);  % kg*m
% C_N_total=F_N*10^(-6)./(V_ref_aver.^2*Coeff_liftdragF_N*10^(-9));          % 单位: N/(rad^2*s^-2*kg*m)=一无量纲
% % C_v=F_vertical*10^(-6)./(V_ref_aver.^2*Coeff_liftdragF_N*10^(-9)); 
% % C_h=F_horizontal*10^(-6)./(V_ref_aver.^2*Coeff_liftdragF_N*10^(-9));    % *10^3); 
% C_v=-sign(alpha).*C_N_total.*cos(alpha); 
% C_h=C_N_total.*sin(abs(alpha)); 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %对瞬时气动力系数使用梯形积分函数trapz进行数值积分，除以周期, 求解平均垂直方向和水平方向气动力系数
% C_vaver=trapz(t,C_v)/(3*T)                 % C_vaver =2.3852;       2.4003(含惯性力)   
% C_haver=trapz(t,C_h)/(3*T)                 % C_haver =-0.2655;   -0.2125(含惯性力)  
% C_habsaver=trapz(t,abs(C_h))/(3*T);     
% C_v2haver=C_vaver/C_habsaver           % 垂直升力系数与水平阻力系数的比值: C_v2haver = 1.1298;  1.026(含惯性力)
% figure(24)
% plot(t/T,C_N_total,'g-','LineWidth',2)
% xlabel('\itNormalized time')
% ylabel('\itC_{N,total}')
% legend('C_{N,total}')
% title('法向合气动力系数随时间的变化规律')
% grid on
% axis([0.9,4.05,-inf,inf])
% figure(25)
% plot(t/T,C_v,'r-',t/T,C_h,'b-','LineWidth',2)
% xlabel('\itNormalized time')
% ylabel('\itC_v & C_h')
% legend('C_v','C_h')
% title('垂直方向和水平方向气动力系数随时间的变化规律')
% grid on
% axis([0.9,4.05,-inf,inf])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 平均值的对比
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F_v_norm_aver=trapz(t,1.5*F_vertical/W)/(3*T)
% Fz_norm_aver=trapz(t,Fz_norm)/(3*T)
% F_haver_norm=trapz(t,2.4*F_horiz_Y/W)/(3*T)
% Fx_norm_aver=trapz(t,Fx_norm)/(3*T)
% F_vert_relaerror=abs(F_v_norm_aver-Fz_norm_aver)/Fz_norm_aver*100
% F_horiz_relaerror=abs(abs(F_haver_norm)-abs(Fx_norm_aver))/abs(Fx_norm_aver)*100
% F_horiz_relaerror=abs(abs(F_haver_norm)-abs(Fx_norm_aver))/abs(F_haver_norm)*100
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 针对四个气动力系数的惩罚项的约束
coeff_con=aeroF_coeff_constraint(x);
% 构建气动力矩系数的惩罚项——penaltyfun
s=2000;
penaltyfun=s*coeff_con;           % penaltyfun =;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stdev=std((2*F_vertical-Fz_norm),0,1);   % 按照S1(0)求解各列(1)元素的标准方差
stdev_v=std((abs(2*F_vertical)-abs(Fz_norm)),0,1);   % 按照S1(0)求解各列(1)元素的标准方差
stdev_h=std((abs(2*F_horiz_Y)-abs(Fx_norm)),0,1);   % 按照S1(0)求解各列(1)元素的标准方差
obj_function=stdev_v+stdev_h+penaltyfun;  % obj_function =40.1316;
% toc  % Elapsed time is 3.401964 seconds.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
