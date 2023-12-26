%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 绘制三维空间球棍模型图——受力图
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (1) 读取翅膀运动学-翅运动规律和几何攻角(AOA)等数据
clear all;clc;
% (1)下面的扭转角是数值求解获得的被动扭转角和角速度以及攻角
pitch_AOA=xlsread('PassiveRot_angle_Alpha.xlsx','A1150:D3170');
t=pitch_AOA(:,1);  % tspan=linspace(0.0001,6*T,1200);
% psi_sim=-pitch_AOA(:,2);     % 注意加了负号——被动扭转角
psi_sim=pitch_AOA(:,2);           % 被动扭转角——度数
alpha_sim=pitch_AOA(:,3);       % 被动扭转攻角AoA——度数
% dpsi_sim=pitch_AOA(:,4);     % 被动扭转角速度——rad/s
w =1185.6;                               %  角频率  
f=188.7;                                    % Hz——翅拍频率 % T=1/f;  
% (2)下面的拍打角是2014-Science-Dickinson等人实验测试拟合而得的拍打角
phi_exp=(3.9008+65.0445*cos(w*t)+4.2642*sin(w*t)+3.5806*cos(2*w*t)-2.9492*sin(2*w*t)+...
               0.1319*cos(3*w*t)+0.3639*sin(3*w*t)+0.7844*cos(4*w*t)+0.2098*sin(4*w*t))*pi/180;
dphi_exp=(4.2642.*w.*cos(t.*w)-5.8984.*w.*cos(2*t.*w)+1.0917.*w.*cos(3*t.*w)+0.8392.*w.*cos(4*t.*w)...
                -65.0445.*w.*sin(t.*w)-7.1612.*w.*sin(2*t.*w)-0.3957.*w.*sin(3*t.*w)-3.1376.*w.*sin(4*t.*w))*pi/180;
ddphi_exp=(11.7968.*w.^2.*sin(2*t.*w)-14.3224.*w.^2.*cos(2*t.*w)-1.1871.*w.^2.*cos(3*t.*w)-12.5504.*w.^2.*cos(4*t.*w)...
                  -4.2642.*w.^2.*sin(t.*w)-65.0445.*w.^2.*cos(t.*w)-3.2751.*w.^2.*sin(3*t.*w)-3.3568.*w.^2.*sin(4*t.*w))*pi/180;  
% figure(20)
% plot(t*f,psi_sim,'r-',t*f,phi_exp*180/pi,'g-',t*f,alpha_sim,'k-')
% xlabel('t');  
% ylabel('\psi_{sim}(t) and \phi_{exp}(t) and \alpha_{sim}(t)');  
% legend('\psi_{sim}(t)','\phi_{exp}(t)','\alpha_{sim}(t)');  
% title('被动转动角\psi_{sim}(t)和拍打角phi_{exp}(t)随时间的变化')  
% grid on  % 被动转动角psi_sim(t)和扭转角phi_exp(t)随时间的变化
%% (2) 截取单个稳定的周期——以攻角变化周期为区间
[alpha_sim_max, locat_max]=max(alpha_sim); % alpha_sim_max =89.8917; locat_max =6;
T=1/f;
t_T=t(locat_max)+T;
indxx=find(t<=t_T);     % 2005
% 刚好2000个数据点
t1=t([6:2005],1);               % 单位是ms
psi1=psi_sim([6:2005],1)*pi/180;           % 扭转角——单位是rad
phi1=phi_exp([6:2005],1);                        % 拍打角——单位是rad
alpha1=alpha_sim([6:2005],1)*pi/180;
% figure(21)
% plot(t1*f,psi1*180/pi,'r-',t1*f,phi1*180/pi,'g-',t1*f,alpha1*180/pi,'k-')  % 转换为度数
% xlabel('t');  
% ylabel('\psi_{sim}(t) and \phi_{exp}(t) and \alpha_{sim}(t)');  
% legend('\psi_{sim}(t)','\phi_{exp}(t)','\alpha_{sim}(t)');  
% title('被动转动角\psi_{sim}(t)和拍打角phi_{exp}(t)随时间的变化')  
% grid on  % 被动转动角psi_sim(t)和扭转角phi_exp(t)随时间的变化
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (3) 绘制单片翅弦微元的三维空间位置——冲程中点时刻
x_lead=2.928631633;
y_lead=0;
z_lead=1.02437628;
x_tail=2.928631633;
y_tail=0;
z_tail=0.14037623;
x=[x_lead,x_tail];
y=[y_lead,y_tail];
z=[z_lead,z_tail];
wing_chord=[x;y;z];  % size(wing_chord)   % (3*2)
figure(22)
plot3(x_lead,y_lead,z_lead,'ko','LineWidth',1.5,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)   % 前缘实心圆点绘制
hold on

mid_chord=[(x_lead+x_tail)/2,(y_lead+y_tail)/2,(z_lead+z_tail)/2];  %中弦点实心圆绘制——力失的起点

plot3(mid_chord(1,1),mid_chord(1,2),mid_chord(1,3),'bo','LineWidth',1.5,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',6) %中弦点实心圆绘制
hold on
% plot3(x,y,z,'k-','LineWidth',1.5)
plot3(wing_chord(1,:),wing_chord(2,:),wing_chord(3,:),'k-','LineWidth',1.5)
xlabel('侧向(x)')
ylabel('纵向(y)')
zlabel('垂直方向(z)')
grid on
% % axis square              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (4) 绘制单片翅弦微元的三维空间位置
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% animation——动画
% 动画片段程序开始段
Nstep=2000;
animation = 1;       %  1=save animation file --------------------------- for change
np=4;                     %  调节总的帧数和时间步长
Nframe = np*40;   % number of movie frames -------------------------- for change
Nskip = floor(Nstep/Nframe);
Nspeed = 2;  % number of frames per second ------------------------ for change
duration = Nframe/Nspeed;  % duration of movie  % 80 seconds @Nframe=4*40
display(['duratin of movie will be ', num2str(duration), ' seconds ']); %输出显示动画持续的时间(s)
%输出如下：duratin of movie will be 10 seconds 
if animation == 1
    aviobj = avifile('ElmentMotion.avi','fps',Nspeed); % filename --- for change
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nskip=50;
% Nstep=2000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:Nskip:Nstep   % L=[1:Nskip:Nstep]; % size(L)=1*40;
    %% 下面一段是坐标系之间的变换
    phi=phi1(i,1);            % 拍打角——单位是rad
    psi=-psi1(i,1);             % 扭转角——单位是rad
    Yaw = [ cos(phi)   -sin(phi)   0; ...
                 sin(phi)   cos(phi)   0; ...
                 0              0            1];         % yaw matrix——符号可能需要方向
    Pitch = [1    0             0; ...
                 0    cos(psi)   -sin(psi); ...
                 0    sin(psi)    cos(psi)];        % pitch matrix——符号可能需要方向

    R_lb = Yaw*Pitch;                                            % from wing frame to body frame
    % R_bl=R_lb';                                                  % from body frame to wing frame
    % wing_chord=[x;y;z];  % size(wing_chord)    % (3*2)
    wing_C = R_lb*wing_chord;                              % wing_chord in body frame  % (3*2)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%    theta=(-90:10:90)'*pi/180;   n=length(theta);     %  size(theta)=(19*1)
%    p1=[1:n;zeros(1,n)]';   % size(p1)=(19*2)            % y=0时x轴上分布的19个起始点的x坐标
%    r=2*ones(n,1);            % size(r)=(19*1)
%    % 极坐标转换为笛卡尔坐标
%    [u,v]=pol2cart(theta,r)  % 调用函数pol2cart:Transform polar or cylindrical coordinates to Cartesian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     for j = 1:4
%         wing_C(:,j) = xl_B(i,:)'+bodyt(:,j);   % body shape in lab frame——可用于添加额外位矢
%     end
figure(23)
% 绘制前缘实心圆点
plot3(wing_C(1,1),wing_C(2,1),wing_C(3,1),'ro-','LineWidth',1.5,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',6)  
hold on
plot3(wing_C(1,1:2),wing_C(2,1:2),wing_C(3,1:2),'r-','LineWidth',1.5) % 绘制前缘后缘实线
hold on
xlabel('侧向(x)')
ylabel('纵向(y)')
zlabel('垂直方向(z)')
grid on
%  axis([-5 5 -5 5 -20 5]);
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 动画片段程序结束段
    if animation == 1
        frame = getframe;
        aviobj = addframe(aviobj,frame);
    else
        pause(0.5)
    end
end
if animation == 1
    aviobj = close(aviobj);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%