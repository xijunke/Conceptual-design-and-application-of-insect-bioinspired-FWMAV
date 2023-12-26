%% Stick figure
clear all; clc;
%% 调用翅膀运动学-翅运动规律和几何攻角(AOA)等数据
% wing_m_output=[t',phi',theta',eta',alpha',dphi',dtheta',deta',ddphi',ddtheta',ddeta',v_x',v_y',a_x',a_y'];
wing_kenimatics=kenimatics_wing_and_AoA_Stickfig();      %调用函数kenimatics_wing_and_AoA;  % size(wing_kenimatics)
% size(wing_kenimatics);                 % (32*15)
t=wing_kenimatics(:,1);                    % 单位是ms
phi=wing_kenimatics(:,2);                % 拍打角——单位是rad
theta=wing_kenimatics(:,3);             % 拍打角——单位是rad
eta=wing_kenimatics(:,4);                % 单位是rad
alpha=wing_kenimatics(:,5);            % 单位是rad/s
dphi=wing_kenimatics(:,6);             % 单位是rad/s
dtheta=wing_kenimatics(:,7);          % 单位是rad/s
deta=wing_kenimatics(:,8);             % 单位是rad/s^2
ddphi=wing_kenimatics(:,9);           % 单位是rad/s^2
ddtheta=wing_kenimatics(:,10);      % 单位是rad/s^2
ddeta=wing_kenimatics(:,11);         % 单位是rad/s^2
v_x=wing_kenimatics(:,12);             % 单位是m/s
v_y=wing_kenimatics(:,13);             % 单位是m/s
a_x=wing_kenimatics(:,14);             % 单位是m/s^2
a_y=wing_kenimatics(:,15);             % 单位是m/s^2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 绘制球棍图
figure
C_aver=300*4.02;  % 平均弦长
x_1=0;
y_1=0;
% dx=10;  dy=0;
plot([0,0],[-5000,5000],'c-',[-5000,5000],[0,0],'c-');  hold on;  % 绘制对称线
%% (1) 前半个冲程
%%%%%%%%%%
eta1=eta(1:16,1);                           %  扭转角
v_x1=v_x(1:16,1);                          %  线速度
v_y1=v_y(1:16,1);                          %  线速度
%%%%%%%%%%
% X_1=[linspace(-1000,x_1,8),linspace(x_1,1000,8)]'; 
X_1=[linspace(-4000,4000,16)]';   % 扭转轴点的x坐标
Y_1=zeros(length(X_1),1)+2000;    % 扭转轴点的y坐标
% plot(X_1,Y_1,'ko')
X_2=X_1+C_aver*sin(eta1)./4;       %  前缘点x坐标
Y_2=Y_1+C_aver*cos(eta1)./4;      %  前缘点y坐标
X_3=X_1-3*C_aver*sin(eta1)./4;    %  后缘点x坐标
Y_3=Y_1-3*C_aver*cos(eta1)./4;    % 后缘点y坐标
X_mid=(X_2+X_3)./2;                    % 生成中弦点的x位置坐标
Y_mid=(Y_2+Y_3)./2;                    % 生成中弦点的x位置坐标
%% 弦向速度方向失和法向速度方向失
U_chord=v_x1.*sin(eta1);             % 翅坐标系下弦向速度——U方向——红箭头
V_chord=v_x1.*cos(eta1);             % 翅坐标系下弦向速度——V方向——红箭头
U_normal=v_y1.*cos(eta1);          % 翅坐标系下法向速度——U方向——绿箭头
V_normal=v_y1.*sin(eta1);           % 翅坐标系下法向速度——V方向——绿箭头
%%
% F_normalU=F_x.*sin(eta1)+F_y.*cos(eta1);
% F_normalV=F_x.*cos(eta1)+F_y.*sin(eta1);
for i=1:16
line([X_2(i),X_3(i)],[Y_2(i),Y_3(i)],'color','k','LineStyle','-','LineWidth',1.5)   % 弦长连线
hold on
plot(X_2(i),Y_2(i),'ko-','LineWidth',1.5,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)   % 前缘实心圆点绘制
end
% %% 弦向——从[X_mid,Y_mid]指向[U_normal, V_normal]的带箭头线段
% 一是采用annotation命令，可以控制箭头特性，语法参见帮助help annotation
% 二是采用text命令，也可以画一些简单箭头；采用第一个命令更容易控制箭头达到你想要的模式。
% h=annotation('arrow',[.9 .5],[.9,.5],'Color','r');
% quiver(X_mid, Y_mid,...                                            % 绘制x,y坐标
%            U_chord, V_chord, 'r','Linewidth',2,'MarkerSize',3);             % 翅坐标系下弦向速度——U,V方向——红箭头
% hold on
%% 法向——从[X_mid,Y_mid]指向[U_normal, V_normal]的带箭头线段
h=quiver(X_mid, Y_mid,...                                             % 绘制x,y坐标
           U_normal, V_normal,'g','Linewidth',2,'MarkerSize',3);   % 翅坐标系下法向速度——U,V方向——绿箭头
% set(h,'maxheadsize',8);
% quiver(X_mid, Y_mid,....                                        % 绘制x,y坐标
%            F_normalU, F_normalV,'b','Linewidth',1);   % 翅坐标系下法向气动力——U,V方向——蓝箭头
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (2) 后半个冲程
%%%%%%%%%%
eta2=eta(17:32,1);                          %  扭转角
v_x2=v_x(17:32,1);                         %  线速度
v_y2=v_y(17:32,1);
%%%%%%%%%%
X_11=sort(X_1,'descend');               % 扭转轴点的x位置坐标
Y_11=zeros(length(X_11),1)-2000;    % 扭转轴点的y位置坐标
X_2=X_11+C_aver*sin(eta2)./4;       %  前缘点x坐标
Y_2=Y_11+C_aver*cos(eta2)./4;      %  前缘点y坐标
X_3=X_11-3*C_aver*sin(eta2)./4;    %  后缘点x坐标
Y_3=Y_11-3*C_aver*cos(eta2)./4;    % 后缘点y坐标
X_mid=(X_2+X_3)./2;                      % 生成中弦点的x位置坐标
Y_mid=(Y_2+Y_3)./2;                      % 生成中弦点的x位置坐标
%% 弦向速度方向失和法向速度方向失
U_chord=v_x2.*sin(eta2);              % 翅坐标系下弦向速度——U方向——红箭头
V_chord=v_x2.*cos(eta2);             % 翅坐标系下弦向速度——V方向——红箭头
U_normal=v_y2.*cos(eta2);           % 翅坐标系下法向速度——U方向——绿箭头
V_normal=v_y2.*sin(eta2);            % 翅坐标系下法向速度——V方向——绿箭头
% F_normalU=F_x.*sin(eta2)+F_y.*cos(eta2);
% F_normalV=F_x.*cos(eta2)+F_y.*sin(eta2);
for i=1:16
   line([X_2(i),X_3(i)],[Y_2(i),Y_3(i)],'color','k','LineStyle','-','LineWidth',1.5)   % 弦长连线
   hold on
   plot(X_2(i),Y_2(i),'ko-','LineWidth',1.5,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)   % 前缘实心圆点绘制
end
%% 弦向——从[X_mid,Y_mid]指向[U_normal, V_normal]的带箭头线段
% quiver(X_mid, Y_mid,...                                        % 绘制x,y坐标
%            U_chord, V_chord, 'r','Linewidth',2,'MarkerSize',3);         % 翅坐标系下弦向速度——U,V方向——红箭头
hold on
%% 法向——从[X_mid,Y_mid]指向[U_normal, V_normal]的带箭头线段
quiver(X_mid, Y_mid,...                                        % 绘制x,y坐标
           U_normal, V_normal,'g','Linewidth',2,'MarkerSize',3);     % 翅坐标系下法向速度——U,V方向——绿箭头
% quiver(X_mid, Y_mid,....                                       % 绘制x,y坐标
%            F_normalU, F_normalV,'b','Linewidth',1); % 翅坐标系下法向气动力——U,V方向——蓝箭头
title('参数化翅运动上叠加的瞬时弦向速度分量和法向速度分量')
axis([-5000,5000,-5000,5000])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





