%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ������ά�ռ����ģ��ͼ��������ͼ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (1) ��ȡ����˶�ѧ-���˶����ɺͼ��ι���(AOA)������
clear all;clc;
% %%  (A) Euler_Motion_Eq6_fruitfly�������ֵ���㱻��Ťת�Ǻͽ��ٶ�
% % A=[t,psi,alpha,dpsi];   % size(B)  % (1200*4)
% % xlswrite('PassiveRot_angle_Alpha.xlsx',B,'sheet1','A1:D1200');
% % pitch_AOA=xlsread('PassiveRot_angle_Alpha.xlsx','A1:D12000'); % ��������
% pitch_AOA=xlsread('PassiveRot_angle_Alpha.xlsx','A1:D12000'); % �������ݡ�������������Ҫ�޸�
% t=pitch_AOA(:,1);
% % psi_sim=pitch_AOA(:,2);          % ����Ťת�� ��������Ҳ���
% % alpha_sim=pitch_AOA(:,3);      % ����Ťת����AoA��������Ҳ���
% psi_sim=-pitch_AOA(:,2);             % ����Ťת�ǡ�����������ע����˸��š���Ŀ����Ԥʵ����Խ���Ա���ʾ
% alpha_sim=pitch_AOA(:,3);          % ����Ťת����AoA��������
% dpsi_sim=pitch_AOA(:,4);            % ����Ťת���ٶȡ���rad/s 
% %%  (B) ��Ϊ�趨�Ĵ��: prescribed=�涨��
% % B=[t,phi_pres,dphi,ddphi];
% stroke=xlsread('stroke_pres_angle_Alpha.xlsx','A1:D12000');% �������ݡ�������������Ҫ�޸�
% t_stroke=stroke(:,1);
% phi_pres=stroke(:,2);
% dphi_pres=stroke(:,3);
% ddphi_pres=stroke(:,4);
%% ���ó���˶�ѧ-���˶����ɺͼ��ι���(AOA)������
% wing_m_output=[t',phi',psi',alpha',dphi',dpsi',ddphi',ddpsi',C_L',C_D',C_N1'];
wing_kenimatics=kenimatics_wing_and_AoA_fruitfly_simpres();      %���ú���kenimatics_wing_and_AoA;  % size(wing_kenimatics)  % (1000,11)
% wing_kenimatics=kenimatics_wing_and_AoA_fruitfly_exp();
t=wing_kenimatics(:,1);                % ��λ��ms
phi_pres=wing_kenimatics(:,2);        % �Ĵ�ǡ�����λ��rad
psi_sim=wing_kenimatics(:,3);            % �Ĵ�ǡ�����λ��rad
alpha_sim=wing_kenimatics(:,4);      % alpha=atan2(omega_z,-omega_y);  % ���ι��ǡ���������   %������������и�
dphi_pres=wing_kenimatics(:,5);          % ��λ��rad/s
dpsi_sim=wing_kenimatics(:,6);          % ��λ��rad/s
ddphi_pres=wing_kenimatics(:,7);       % ��λ��rad/s^2
ddpsi_sim=wing_kenimatics(:,8);     % ��λ��rad/s^2
% C_L=wing_kenimatics(:,9);          
% C_D=wing_kenimatics(:,10);     
% C_N1=wing_kenimatics(:,11);   
% C_N=C_N1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
f=188.7;                % Hz��������Ƶ�� % T=1/f;  
% figure(1)
% plot(t*f,psi_sim*180/pi,'r-',t*f,phi_pres*180/pi,'g-',t*f,alpha_sim*180/pi,'b:','LineWidth',2)
% xlabel('\itNormalized time');  
% ylabel('\psi_{sim}(t) and \phi_{pres}(t) and \alpha_{sim}(t)');  
% legend('\psi_{sim}(t)','\phi_{pres}(t)','\alpha_{sim}(t)');  
% title('����ת����\psi_{sim}(t)���Ĵ��phi_{pres}(t)��ʱ��ı仯')  
% grid on  % ����ת����psi_sim(t)��Ťת��phi_pres(t)��ʱ��ı仯
% axis([0.9,4.05,-110,110])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% (4) ��ȡ�����ȶ������ڡ����Թ��Ǳ仯����Ϊ����
% ���ù���Ϊ90�ȶ�Ӧʱ��Ϊ���³�̵Ľ���ʱ���
% [alpha_sim_max, locat_max]=max(alpha_sim)   % alpha_sim_max =1.5913; locat_max =835;
% T=1/f;
% t_T=t(locat_max)+T;
% indxx=find(t<=t_T)      % 2005
% num=length(indxx)     % num=275;
[alpha_sim_max, locat_max]=min(abs(alpha_sim));   % alpha_sim_max =0.3581; locat_max =219;
T=1/f;
t_T=t(locat_max)+T;
indxx=find(t<=t_T);      % 552
num=length(indxx);     % num=552;
% �պ�104�����ݵ�=275-172+1;
t1=t([locat_max:num],1);                     % ��λ��ms
psi1=psi_sim([locat_max:num],1);           % Ťת�ǡ�����λ��rad
dpsi1=dpsi_sim([locat_max:num],1);           % Ťת�ǡ�����λ��rad
% ddpsi1=ddpsi_sim([locat_max:num],1);           % Ťת�ǡ�����λ��rad
ddpsi1=dpsi_sim([locat_max:num],1);          % Ťת�ǡ�����λ��rad��ע�������Ťת�Ǽ��ٶ�ȡ��Ťת���ٶ�ddpsi1=dpsi_sim, ��ΪODE����û�����ddpsi
phi1=phi_pres([locat_max:num],1);                        % �Ĵ�ǡ�����λ��rad
dphi1=dphi_pres([locat_max:num],1);                        % �Ĵ�ǡ�����λ��rad
ddphi1=ddphi_pres([locat_max:num],1);                        % �Ĵ�ǡ�����λ��rad
alpha1=alpha_sim([locat_max:num],1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% �������ȶ������ڵ���(�ϳ�)�������������ݶ��롪��mN
% �Ա�thr_dim_chord2����ĵ������ڻ������ͼ����������ϵ�µķ�����ʧ�յ�y����
F_normal1=xlsread('Forcenormal_oneT_simpres.xlsx','A1:A2000');
Forcenormal=0.2*F_normal1([locat_max:num],1);   %  ����С��1/5����ʾquiver3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(2)
% plot(t1*f,psi1*180/pi,'r-',t1*f,phi1*180/pi,'g-',t1*f,alpha1*180/pi,'k-',t1*f,2*abs(Forcenormal/0.2),'b-','LineWidth',2)  % ת��Ϊ����
% xlabel('\itNormalized time');  
% ylabel('\psi_{sim}(t) and \phi_{pres}(t) and \alpha_{sim}(t) and 2*abs(Force_{normal}(t))');  
% legend('\psi_{sim}(t)','\phi_{pres}(t)','\alpha_{sim}(t)','2*abs(Force_{normal}(t))');  
% title('����ת����\psi_{sim}(t), �Ĵ��phi_{pres}(t) and 2*abs(Force_{normal}(t))��ʱ��ı仯')  
% grid on  % ����ת����psi_sim(t)��Ťת��phi_pres(t)��ʱ��ı仯
% axis([1.6,2.7,-110,110])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (5) ����ȡ�ĵ����ȶ������ڡ����Թ��Ǳ仯����Ϊ������˶�ѧ���������Excel����������
% �Ա�kenimatics_wing_and_AoA_fruitfly�����ĵ������ڼ���������ϵ����
A=[t1,psi1,dpsi1,ddpsi1,phi1,dphi1,ddphi1];     % size(A)  % (2000*7)     % ����������
xlswrite('wingmotion_oneT_simpres.xlsx',A,'sheet1','A1:G2000');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (7) ���Ƶ�Ƭ����΢Ԫ����ά�ռ�λ�á�������е�ʱ��
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
% figure(3)
% mid_chord=[(x_lead+x_tail)/2,(y_lead+y_tail)/2,(z_lead+z_tail)/2];  %���ҵ�ʵ��Բ���ơ�����ʧ�����
% plot3(x_lead,y_lead,z_lead,'ko','LineWidth',1.5,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)   % ǰԵʵ��Բ�����
% hold on
% plot3(mid_chord(1,1),mid_chord(1,2),mid_chord(1,3),'bo','LineWidth',1.5,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',6) %���ҵ�ʵ��Բ����
% hold on
% % plot3(x,y,z,'k-','LineWidth',1.5)
% plot3(wing_chord(1,:),wing_chord(2,:),wing_chord(3,:),'k-','LineWidth',1.5)
% xlabel('����(x)')
% ylabel('����(y)')
% zlabel('��ֱ����(z)')
% grid on
% % axis square              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (8) ���Ƶ�Ƭ����΢Ԫ����ά�ռ�λ��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% animation��������
% ����Ƭ�γ���ʼ��
Nstep=num- locat_max+1;    % Nstep=334; 
animation = 1;       %  1=save animation file --------------------------- for change
np=2;                     %  �����ܵ�֡����ʱ�䲽��
Nframe = np*16.7;   % number of movie frames -------------------------- for change
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
    psi=psi1(i,1);             % Ťת�ǡ�����λ��rad
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(23)
% ����ǰԵʵ��Բ��
plot3(wing_C(1,1),wing_C(2,1),wing_C(3,1),'ro-','LineWidth',1.5,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',6)
hold on
grid on
axis([0 4 -3 3 -0.5 3.3]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot3(Fnormal (1,1:2),Fnormal (2,1:2),Fnormal (3,1:2),'b-','LineWidth',1.5) % �������������ҵ�ķ�����ʧ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) ʹ��arrow3������ά��ͷ
% % mid_chord=[(x_lead+x_tail)/2,(y_lead+y_tail)/2,(z_lead+z_tail)/2];  %���������ҵ�ķ�����ʧ�������
% % F_normal_termini=[(x_lead+x_tail)/2,Forcenormal(i,1),(z_lead+z_tail)/2];  % ���������ҵ�ķ�����ʧ�����յ�
% x1=[(x_lead+x_tail)/2,(x_lead+x_tail)/2];
% y1=[(y_lead+y_tail)/2,-Forcenormal(i,1)];
% z1=[(z_lead+z_tail)/2,(z_lead+z_tail)/2];
% F_normal=[x1;y1;z1];  % size(wing_chord)   % (3*2)
% Fnormal = R_lb*F_normal;        % wing_chord in body frame  % (3*2)
% F_mid=[Fnormal(1,1),Fnormal(2,1),Fnormal(3,1)];
% F_termini=[Fnormal(1,2),Fnormal(2,2),Fnormal(3,2)];
% % % pbaspect([4 4 1])    % ���ú���pbaspect: Set or query plot box aspect ratio(ͼ���ݺ��)
% % % pbaspect('auto'); % daspect('auto')
% arrow3(F_mid,F_termini,'b',0.9)     % ���ú���arrow3���ƴ����������ҵ�ķ�����ʧ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) ע��quiver��ʹ��: �յ��Ӧ���������
% A three-dimensional quiver plot displays vectors with components (u,v,w) at the points (x,y,z), 
% where u,v,w,x,y, and z all have real (non-complex) values.
x1=[(x_lead+x_tail)/2,0];
y1=[(y_lead+y_tail)/2,Forcenormal(i,1)];
z1=[(z_lead+z_tail)/2,0];
F_normal=[x1;y1;z1];  % size(wing_chord)   % (3*2)
Fnormal = R_lb*F_normal;        % wing_chord in body frame  % (3*2)
quiver3(Fnormal(1,1),Fnormal(2,1),Fnormal(3,1),Fnormal(1,2),Fnormal(2,2),Fnormal(3,2),0.5,'b-','LineWidth',2.0); % ��������0.5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot3(wing_C(1,1:2),wing_C(2,1:2),wing_C(3,1:2),'r-','LineWidth',1.5) % ����ǰԵ��Եʵ��
hold on
xlabel('����(x)')
ylabel('����(y)')
zlabel('��ֱ����(z)')
% hold off
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