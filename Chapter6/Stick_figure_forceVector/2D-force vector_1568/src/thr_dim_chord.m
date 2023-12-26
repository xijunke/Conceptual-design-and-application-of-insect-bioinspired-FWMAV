%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ������ά�ռ����ģ��ͼ��������ͼ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (1) ��ȡ����˶�ѧ-���˶����ɺͼ��ι���(AOA)������
clear all;clc;
% (1)�����Ťת������ֵ����õı���Ťת�Ǻͽ��ٶ��Լ�����
pitch_AOA=xlsread('PassiveRot_angle_Alpha.xlsx','A1150:D3170');
t=pitch_AOA(:,1);  % tspan=linspace(0.0001,6*T,1200);
% psi_sim=-pitch_AOA(:,2);     % ע����˸��š�������Ťת��
psi_sim=pitch_AOA(:,2);           % ����Ťת�ǡ�������
alpha_sim=pitch_AOA(:,3);       % ����Ťת����AoA��������
% dpsi_sim=pitch_AOA(:,4);     % ����Ťת���ٶȡ���rad/s
w =1185.6;                               %  ��Ƶ��  
f=188.7;                                    % Hz��������Ƶ�� % T=1/f;  
% (2)������Ĵ����2014-Science-Dickinson����ʵ�������϶��õ��Ĵ��
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
% title('����ת����\psi_{sim}(t)���Ĵ��phi_{exp}(t)��ʱ��ı仯')  
% grid on  % ����ת����psi_sim(t)��Ťת��phi_exp(t)��ʱ��ı仯
%% (2) ��ȡ�����ȶ������ڡ����Թ��Ǳ仯����Ϊ����
[alpha_sim_max, locat_max]=max(alpha_sim); % alpha_sim_max =89.8917; locat_max =6;
T=1/f;
t_T=t(locat_max)+T;
indxx=find(t<=t_T);     % 2005
% �պ�2000�����ݵ�
t1=t([6:2005],1);               % ��λ��ms
psi1=psi_sim([6:2005],1)*pi/180;           % Ťת�ǡ�����λ��rad
phi1=phi_exp([6:2005],1);                        % �Ĵ�ǡ�����λ��rad
alpha1=alpha_sim([6:2005],1)*pi/180;
% figure(21)
% plot(t1*f,psi1*180/pi,'r-',t1*f,phi1*180/pi,'g-',t1*f,alpha1*180/pi,'k-')  % ת��Ϊ����
% xlabel('t');  
% ylabel('\psi_{sim}(t) and \phi_{exp}(t) and \alpha_{sim}(t)');  
% legend('\psi_{sim}(t)','\phi_{exp}(t)','\alpha_{sim}(t)');  
% title('����ת����\psi_{sim}(t)���Ĵ��phi_{exp}(t)��ʱ��ı仯')  
% grid on  % ����ת����psi_sim(t)��Ťת��phi_exp(t)��ʱ��ı仯
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (3) ���Ƶ�Ƭ����΢Ԫ����ά�ռ�λ�á�������е�ʱ��
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
plot3(x_lead,y_lead,z_lead,'ko','LineWidth',1.5,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)   % ǰԵʵ��Բ�����
hold on

mid_chord=[(x_lead+x_tail)/2,(y_lead+y_tail)/2,(z_lead+z_tail)/2];  %���ҵ�ʵ��Բ���ơ�����ʧ�����

plot3(mid_chord(1,1),mid_chord(1,2),mid_chord(1,3),'bo','LineWidth',1.5,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',6) %���ҵ�ʵ��Բ����
hold on
% plot3(x,y,z,'k-','LineWidth',1.5)
plot3(wing_chord(1,:),wing_chord(2,:),wing_chord(3,:),'k-','LineWidth',1.5)
xlabel('����(x)')
ylabel('����(y)')
zlabel('��ֱ����(z)')
grid on
% % axis square              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (4) ���Ƶ�Ƭ����΢Ԫ����ά�ռ�λ��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% animation��������
% ����Ƭ�γ���ʼ��
Nstep=2000;
animation = 1;       %  1=save animation file --------------------------- for change
np=4;                     %  �����ܵ�֡����ʱ�䲽��
Nframe = np*40;   % number of movie frames -------------------------- for change
Nskip = floor(Nstep/Nframe);
Nspeed = 2;  % number of frames per second ------------------------ for change
duration = Nframe/Nspeed;  % duration of movie  % 80 seconds @Nframe=4*40
display(['duratin of movie will be ', num2str(duration), ' seconds ']); %�����ʾ����������ʱ��(s)
%������£�duratin of movie will be 10 seconds 
if animation == 1
    aviobj = avifile('ElmentMotion.avi','fps',Nspeed); % filename --- for change
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nskip=50;
% Nstep=2000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:Nskip:Nstep   % L=[1:Nskip:Nstep]; % size(L)=1*40;
    %% ����һ��������ϵ֮��ı任
    phi=phi1(i,1);            % �Ĵ�ǡ�����λ��rad
    psi=-psi1(i,1);             % Ťת�ǡ�����λ��rad
    Yaw = [ cos(phi)   -sin(phi)   0; ...
                 sin(phi)   cos(phi)   0; ...
                 0              0            1];         % yaw matrix�������ſ�����Ҫ����
    Pitch = [1    0             0; ...
                 0    cos(psi)   -sin(psi); ...
                 0    sin(psi)    cos(psi)];        % pitch matrix�������ſ�����Ҫ����

    R_lb = Yaw*Pitch;                                            % from wing frame to body frame
    % R_bl=R_lb';                                                  % from body frame to wing frame
    % wing_chord=[x;y;z];  % size(wing_chord)    % (3*2)
    wing_C = R_lb*wing_chord;                              % wing_chord in body frame  % (3*2)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%    theta=(-90:10:90)'*pi/180;   n=length(theta);     %  size(theta)=(19*1)
%    p1=[1:n;zeros(1,n)]';   % size(p1)=(19*2)            % y=0ʱx���Ϸֲ���19����ʼ���x����
%    r=2*ones(n,1);            % size(r)=(19*1)
%    % ������ת��Ϊ�ѿ�������
%    [u,v]=pol2cart(theta,r)  % ���ú���pol2cart:Transform polar or cylindrical coordinates to Cartesian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     for j = 1:4
%         wing_C(:,j) = xl_B(i,:)'+bodyt(:,j);   % body shape in lab frame������������Ӷ���λʸ
%     end
figure(23)
% ����ǰԵʵ��Բ��
plot3(wing_C(1,1),wing_C(2,1),wing_C(3,1),'ro-','LineWidth',1.5,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',6)  
hold on
plot3(wing_C(1,1:2),wing_C(2,1:2),wing_C(3,1:2),'r-','LineWidth',1.5) % ����ǰԵ��Եʵ��
hold on
xlabel('����(x)')
ylabel('����(y)')
zlabel('��ֱ����(z)')
grid on
%  axis([-5 5 -5 5 -20 5]);
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ����Ƭ�γ��������
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