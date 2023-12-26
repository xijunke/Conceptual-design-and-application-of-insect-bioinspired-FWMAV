%% Stick figure
clear all; clc;
%% ���ó���˶�ѧ-���˶����ɺͼ��ι���(AOA)������
% wing_m_output=[t',phi',theta',eta',alpha',dphi',dtheta',deta',ddphi',ddtheta',ddeta',v_x',v_y',a_x',a_y'];
wing_kenimatics=kenimatics_wing_and_AoA_Stickfig();      %���ú���kenimatics_wing_and_AoA;  % size(wing_kenimatics)
% size(wing_kenimatics);                 % (32*15)
t=wing_kenimatics(:,1);                    % ��λ��ms
phi=wing_kenimatics(:,2);                % �Ĵ�ǡ�����λ��rad
theta=wing_kenimatics(:,3);             % �Ĵ�ǡ�����λ��rad
eta=wing_kenimatics(:,4);                % ��λ��rad
alpha=wing_kenimatics(:,5);            % ��λ��rad/s
dphi=wing_kenimatics(:,6);             % ��λ��rad/s
dtheta=wing_kenimatics(:,7);          % ��λ��rad/s
deta=wing_kenimatics(:,8);             % ��λ��rad/s^2
ddphi=wing_kenimatics(:,9);           % ��λ��rad/s^2
ddtheta=wing_kenimatics(:,10);      % ��λ��rad/s^2
ddeta=wing_kenimatics(:,11);         % ��λ��rad/s^2
v_x=wing_kenimatics(:,12);             % ��λ��m/s
v_y=wing_kenimatics(:,13);             % ��λ��m/s
a_x=wing_kenimatics(:,14);             % ��λ��m/s^2
a_y=wing_kenimatics(:,15);             % ��λ��m/s^2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �������ͼ
figure
C_aver=300*4.02;  % ƽ���ҳ�
x_1=0;
y_1=0;
% dx=10;  dy=0;
plot([0,0],[-5000,5000],'c-',[-5000,5000],[0,0],'c-');  hold on;  % ���ƶԳ���
%% (1) ǰ������
%%%%%%%%%%
eta1=eta(1:16,1);                           %  Ťת��
v_x1=v_x(1:16,1);                          %  ���ٶ�
v_y1=v_y(1:16,1);                          %  ���ٶ�
%%%%%%%%%%
% X_1=[linspace(-1000,x_1,8),linspace(x_1,1000,8)]'; 
X_1=[linspace(-4000,4000,16)]';   % Ťת����x����
Y_1=zeros(length(X_1),1)+2000;    % Ťת����y����
% plot(X_1,Y_1,'ko')
X_2=X_1+C_aver*sin(eta1)./4;       %  ǰԵ��x����
Y_2=Y_1+C_aver*cos(eta1)./4;      %  ǰԵ��y����
X_3=X_1-3*C_aver*sin(eta1)./4;    %  ��Ե��x����
Y_3=Y_1-3*C_aver*cos(eta1)./4;    % ��Ե��y����
X_mid=(X_2+X_3)./2;                    % �������ҵ��xλ������
Y_mid=(Y_2+Y_3)./2;                    % �������ҵ��xλ������
%% �����ٶȷ���ʧ�ͷ����ٶȷ���ʧ
U_chord=v_x1.*sin(eta1);             % ������ϵ�������ٶȡ���U���򡪡����ͷ
V_chord=v_x1.*cos(eta1);             % ������ϵ�������ٶȡ���V���򡪡����ͷ
U_normal=v_y1.*cos(eta1);          % ������ϵ�·����ٶȡ���U���򡪡��̼�ͷ
V_normal=v_y1.*sin(eta1);           % ������ϵ�·����ٶȡ���V���򡪡��̼�ͷ
%%
% F_normalU=F_x.*sin(eta1)+F_y.*cos(eta1);
% F_normalV=F_x.*cos(eta1)+F_y.*sin(eta1);
for i=1:16
line([X_2(i),X_3(i)],[Y_2(i),Y_3(i)],'color','k','LineStyle','-','LineWidth',1.5)   % �ҳ�����
hold on
plot(X_2(i),Y_2(i),'ko-','LineWidth',1.5,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)   % ǰԵʵ��Բ�����
end
% %% ���򡪡���[X_mid,Y_mid]ָ��[U_normal, V_normal]�Ĵ���ͷ�߶�
% һ�ǲ���annotation������Կ��Ƽ�ͷ���ԣ��﷨�μ�����help annotation
% ���ǲ���text���Ҳ���Ի�һЩ�򵥼�ͷ�����õ�һ����������׿��Ƽ�ͷ�ﵽ����Ҫ��ģʽ��
% h=annotation('arrow',[.9 .5],[.9,.5],'Color','r');
% quiver(X_mid, Y_mid,...                                            % ����x,y����
%            U_chord, V_chord, 'r','Linewidth',2,'MarkerSize',3);             % ������ϵ�������ٶȡ���U,V���򡪡����ͷ
% hold on
%% ���򡪡���[X_mid,Y_mid]ָ��[U_normal, V_normal]�Ĵ���ͷ�߶�
h=quiver(X_mid, Y_mid,...                                             % ����x,y����
           U_normal, V_normal,'g','Linewidth',2,'MarkerSize',3);   % ������ϵ�·����ٶȡ���U,V���򡪡��̼�ͷ
% set(h,'maxheadsize',8);
% quiver(X_mid, Y_mid,....                                        % ����x,y����
%            F_normalU, F_normalV,'b','Linewidth',1);   % ������ϵ�·�������������U,V���򡪡�����ͷ
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (2) �������
%%%%%%%%%%
eta2=eta(17:32,1);                          %  Ťת��
v_x2=v_x(17:32,1);                         %  ���ٶ�
v_y2=v_y(17:32,1);
%%%%%%%%%%
X_11=sort(X_1,'descend');               % Ťת����xλ������
Y_11=zeros(length(X_11),1)-2000;    % Ťת����yλ������
X_2=X_11+C_aver*sin(eta2)./4;       %  ǰԵ��x����
Y_2=Y_11+C_aver*cos(eta2)./4;      %  ǰԵ��y����
X_3=X_11-3*C_aver*sin(eta2)./4;    %  ��Ե��x����
Y_3=Y_11-3*C_aver*cos(eta2)./4;    % ��Ե��y����
X_mid=(X_2+X_3)./2;                      % �������ҵ��xλ������
Y_mid=(Y_2+Y_3)./2;                      % �������ҵ��xλ������
%% �����ٶȷ���ʧ�ͷ����ٶȷ���ʧ
U_chord=v_x2.*sin(eta2);              % ������ϵ�������ٶȡ���U���򡪡����ͷ
V_chord=v_x2.*cos(eta2);             % ������ϵ�������ٶȡ���V���򡪡����ͷ
U_normal=v_y2.*cos(eta2);           % ������ϵ�·����ٶȡ���U���򡪡��̼�ͷ
V_normal=v_y2.*sin(eta2);            % ������ϵ�·����ٶȡ���V���򡪡��̼�ͷ
% F_normalU=F_x.*sin(eta2)+F_y.*cos(eta2);
% F_normalV=F_x.*cos(eta2)+F_y.*sin(eta2);
for i=1:16
   line([X_2(i),X_3(i)],[Y_2(i),Y_3(i)],'color','k','LineStyle','-','LineWidth',1.5)   % �ҳ�����
   hold on
   plot(X_2(i),Y_2(i),'ko-','LineWidth',1.5,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)   % ǰԵʵ��Բ�����
end
%% ���򡪡���[X_mid,Y_mid]ָ��[U_normal, V_normal]�Ĵ���ͷ�߶�
% quiver(X_mid, Y_mid,...                                        % ����x,y����
%            U_chord, V_chord, 'r','Linewidth',2,'MarkerSize',3);         % ������ϵ�������ٶȡ���U,V���򡪡����ͷ
hold on
%% ���򡪡���[X_mid,Y_mid]ָ��[U_normal, V_normal]�Ĵ���ͷ�߶�
quiver(X_mid, Y_mid,...                                        % ����x,y����
           U_normal, V_normal,'g','Linewidth',2,'MarkerSize',3);     % ������ϵ�·����ٶȡ���U,V���򡪡��̼�ͷ
% quiver(X_mid, Y_mid,....                                       % ����x,y����
%            F_normalU, F_normalV,'b','Linewidth',1); % ������ϵ�·�������������U,V���򡪡�����ͷ
title('���������˶��ϵ��ӵ�˲ʱ�����ٶȷ����ͷ����ٶȷ���')
axis([-5000,5000,-5000,5000])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





