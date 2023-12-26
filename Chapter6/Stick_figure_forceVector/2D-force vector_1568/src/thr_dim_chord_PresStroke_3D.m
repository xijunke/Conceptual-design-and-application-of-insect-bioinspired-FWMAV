%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ������ά�ռ����ģ��ͼ��������ͼ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (1) ��ȡ����˶�ѧ-���˶����ɺͼ��ι���(AOA)������
clear all;clc;
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

%% (2) ��ȡ�����ȶ������ڡ����Թ��Ǳ仯����Ϊ����
% ���ù���Ϊ90�ȶ�Ӧʱ��Ϊ���³�̵Ľ���ʱ���
[phi_pres_max, locat_max]=max(phi_pres);   % phi_pres_max =1.149; locat_max =2;
T=1/f;
t_T=t(locat_max)+T;
index=find(t<=t_T);      % 2005
num=length(index);     % num=335;
% [alpha_sim_max, locat_max]=min(abs(alpha_sim));   % alpha_sim_max =0.3581; locat_max =219;
% T=1/f;
% t_T=t(locat_max)+T;
% index=find(t<=t_T);      % 552
% num=length(index);     % num=552;
% �պ�104�����ݵ�=275-172+1;
t1=t((locat_max:num),1);                     % ��λ��ms
psi1=psi_sim((locat_max:num),1);           % Ťת�ǡ�����λ��rad
alpha1=alpha_sim((locat_max:num),1);         % Ťת�ǡ�����λ��rad
dpsi1=dpsi_sim((locat_max:num),1);           % Ťת�ǡ�����λ��rad
% ddpsi1=ddpsi_sim((locat_max:num),1);           % Ťת�ǡ�����λ��rad
ddpsi1=dpsi_sim((locat_max:num),1);          % Ťת�ǡ�����λ��rad��ע�������Ťת�Ǽ��ٶ�ȡ��Ťת���ٶ�ddpsi1=dpsi_sim, ��ΪODE����û�����ddpsi
phi1=phi_pres((locat_max:num),1);                        % �Ĵ�ǡ�����λ��rad
dphi1=dphi_pres((locat_max:num),1);                        % �Ĵ�ǡ�����λ��rad
ddphi1=ddphi_pres((locat_max:num),1);                        % �Ĵ�ǡ�����λ��rad
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (3) �������ȶ������ڵ���(�ϳ�)�������������ݶ��롪��mN
% �Ա�thr_dim_chord2����ĵ������ڻ������ͼ����������ϵ�µķ�����ʧ�յ�y����
% B=[F_ytran,F_yadd1,F_yrot];    % size(B)  % (2000*1)     %����mN
F_normal=xlsread('Forcenormal_oneT_simpres.xlsx','A1:C1000');
% Forcenormal=0.2*F_normal((locat_max:num),1);   %  ����С��1/5����ʾquiver3
F_ytran=0.1*F_normal((locat_max:num),1);
F_yadd=0.5*F_normal((locat_max:num),2);
F_yrot=0.2*F_normal((locat_max:num),3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(2)
% plot(t1*f,psi1*180/pi,'r-',t1*f,phi1*180/pi,'g-',t1*f,alpha1*180/pi,'k-',t1*f,2*abs(F_ytran/0.2),'b-','LineWidth',2)  % ת��Ϊ����
% xlabel('\itNormalized time');  
% ylabel('\psi_{sim}(t) and \phi_{pres}(t) and \alpha_{sim}(t) and 2*abs(Force_{normal}(t))');  
% legend('\psi_{sim}(t)','\phi_{pres}(t)','\alpha_{sim}(t)','2*abs(Force_{normal}(t))');  
% title('����ת����\psi_{sim}(t), �Ĵ��phi_{pres}(t) and 2*abs(Force_{normal}(t))��ʱ��ı仯')  
% grid on  % ����ת����psi_sim(t)��Ťת��phi_pres(t)��ʱ��ı仯
% % axis([1.6,2.7,-110,110])
% axis([0.9,2.1,-110,110])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (4) ����ȡ�ĵ����ȶ������ڡ����Թ��Ǳ仯����Ϊ������˶�ѧ���������Excel����������
% �Ա�kenimatics_wing_and_AoA_fruitfly�����ĵ������ڼ���������ϵ����
A=[t1,psi1,dpsi1,ddpsi1,phi1,dphi1,ddphi1];     % size(A)  % (2000*7)     % ����������
xlswrite('wingmotion_oneT_simpres.xlsx',A,'sheet1','A1:G2000');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (5) ���Ƶ�Ƭ����΢Ԫ����ά�ռ�λ�á�������е�ʱ��
R_wingeff =3.0040;        % ��λ�� mm   
xr=0.3289;                     % x-root offset  \mm
R_eff=R_wingeff+xr;
C_avereff =0.8854;          % ��λ�� mm
C_025=0.25*C_avereff;
C_075=0.75*C_avereff;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % psi_sim=0;
% % alpha_sim=0;
% psi_sim=pi/4;
% alpha_sim=pi/4;
% % % ǰԵ����ʵ��Բ�����
% x1=R_eff;
% y1=C_025*sin(psi_sim);
% z1=C_025*cos(psi_sim);
% % Ťת���ġ���ʵ��Բ�����
% x2=R_eff;
% y2=0;
% z2=0;
% % ����ѹ�ĵ㡪������ʵ��Բ����� 
% d_cp=(0.82*abs(alpha_sim)/pi+0.05)*C_avereff;
% if d_cp>C_025
%     x3=R_eff;
%     y3=-(d_cp-C_025)*cos(psi_sim);
%     z3=-(d_cp-C_025)*sin(psi_sim);
% elseif d_cp<=C_025;
%     x3=R_eff;
%     y3=(C_025-d_cp)*sin(psi_sim);
%     z3=(C_025-d_cp)*cos(psi_sim);
% end
% % ���ҵ㡪��ʧ����㡪ʵ��Բ����
% x4=R_eff;
% y4=-C_025*sin(psi_sim);
% z4=-C_025*cos(psi_sim);
% % ��Ե����ʵ��Բ�����
% x5=R_eff;
% y5=-C_075*sin(psi_sim);
% z5=-C_075*cos(psi_sim);
% % ƬԪ������������
% x=[x1,x2,x3,x4,x5];
% y=[y1,y2,y3,y4,y5];
% z=[z1,z2,z3,z4,z5];
% wing_chord=[x;y;z];  % size(wing_chord)   % (3*5)
% % ����ƬԪ������������ķֲ�
% figure(3)
% plot3(x1,y1,z1,'ko','LineWidth',1.5,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',8);   % ǰԵ����ʵ��Բ�����
% hold on
% plot3(x2,y2,z2,'ro','LineWidth',1.5,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',5);   % Ťת���ġ���ʵ��Բ�����
% hold on
% plot3(x3,y3,z3,'bo','LineWidth',1.5,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',4);  % ����ѹ�ĵ㡪������ʵ��Բ�����
% hold on
% plot3(x4,y4,z4,'go','LineWidth',1.5,'MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',6);  % ���ҵ㡪��ʧ����㡪ʵ��Բ����
% hold on
% plot3(x5,y5,z5,'ko','LineWidth',1.5,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',3)   % ��Ե����ʵ��Բ�����
% hold on
% plot3(x,y,z,'k-','LineWidth',2)
% xlabel('����(x)')
% ylabel('����(y)')
% zlabel('��ֱ����(z)')
% grid on
% % axis square              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (6) ���Ƶ�Ƭ����΢Ԫ����ά�ռ�λ��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% animation��������
% ����Ƭ�γ���ʼ��
Nstep=num- locat_max+1;    % Nstep=334; 
np=2;                     %  �����ܵ�֡����ʱ�䲽��
% Nframe = np*10.4375;   % number of movie frames -------------------------- for change
Nframe = np*16.7;   % number of movie frames -------------------------- for change
Nskip = floor(Nstep/Nframe);
Nspeed = 2;  % number of frames per second ------------------------ for change
animation = 1;       %  1=save animation file --------------------------- for change
duration = Nframe/Nspeed;  % duration of movie  % 80 seconds @Nframe=4*40
display(['duratin of movie will be ', num2str(duration), ' seconds ']); %�����ʾ����������ʱ��(s)
%������£�duratin of movie will be 10 seconds 
if animation == 1
    aviobj = avifile('ElmentMotion_3D.avi','fps',Nspeed); % filename --- for change
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
    R_w2b = Yaw*Pitch;     % from wing frame to body frame�����ӳ�����ϵ��������ϵ
    %     R_b2w=R_w2b';                    % from body frame to wing frame������������ϵ��������ϵ
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % ǰԵ����ʵ��Բ�����
    x1=R_eff;
    y1=0;
    z1=C_025;
    p_lead=[x1;y1;z1];
    p_lead=R_w2b*p_lead;
    % Ťת���ġ���ʵ��Բ�����
    x2=R_eff;
    y2=0;
    z2=0;
    p_rotaxis=[x2;y2;z2];
    p_rotaxis=R_w2b*p_rotaxis;
    % ����ѹ�ĵ㡪������ʵ��Բ����� 
    alpha=alpha1(i,1);             % Ťת�ǡ�����λ��rad
    d_cp=(0.82*abs(alpha)/pi+0.05)*C_avereff;
    y_dcp=0;
    z_dcp=-d_cp;
%     if d_cp>C_025
%         x3=R_eff;
%         y3=y_dcp-y1;
%         z3=z_dcp-z1;
%     elseif d_cp<=C_025;
%         x3=R_eff;
%         y3=y_dcp+y1;
%         z3=z_dcp+z1;
%     end
    if d_cp>C_025
        x3=R_eff;
        y3=0;
        z3=-(d_cp-C_025);
    elseif d_cp<=C_025;
        x3=R_eff;
        y3=0;
        z3=(C_025-d_cp);
    end
   p_cop=[x3;y3;z3];
   p_cop=R_w2b*p_cop;
    % ���ҵ㡪��ʧ����㡪ʵ��Բ����
    x4=R_eff;
    y4=0;
    z4=-C_025;
    p_mid=[x4;y4;z4];
    p_mid=R_w2b*p_mid;
    % ��Ե����ʵ��Բ�����
    x5=R_eff;
    y5=0;
    z5=-C_075;
    p_tail=[x5;y5;z5];
    p_tail=R_w2b*p_tail;
    % ƬԪ������������
    x=[p_lead(1,1),p_rotaxis(1,1),p_cop(1,1),p_mid(1,1),p_tail(1,1)];
    y=[p_lead(2,1),p_rotaxis(2,1),p_cop(2,1),p_mid(2,1),p_tail(2,1)];
    z=[p_lead(3,1),p_rotaxis(3,1),p_cop(3,1),p_mid(3,1),p_tail(3,1)];
    wing_C=[x;y;z];  % size(wing_chord)   % (3*5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
figure(4)
xlabel('����(x)')
ylabel('����(y)')
zlabel('��ֱ����(z)')
grid on
axis([-0.5 3.5 -4 4 -1 1.5]);
% ����ǰԵʵ��Բ��
plot3(p_lead(1,1),p_lead(2,1),p_lead(3,1),'ro-','LineWidth',1.5,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',6)
hold on
plot3([p_lead(1,1),p_tail(1,1)],[p_lead(2,1),p_tail(2,1)],[p_lead(3,1),p_tail(3,1)],'r-','LineWidth',1.5) % ����ǰԵ��Եʵ��
hold on
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
% Fnormal = R_w2b*F_normal;        % wing_chord in body frame  % (3*2)
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
%% (a) ����ѹ�ĵ㡪������ʧ
x_cop=[x3,0];
y_cop=[y3,F_ytran(i,1)];  % ����ѹ�ĵ���ʧ����㡪��>�յ�
z_cop=[z3,0];
F_ytran_cop=[x_cop;y_cop;z_cop];  % size(wing_chord)   % (3*2)
Fytran_cop= R_w2b*F_ytran_cop;        % wing_chord in body frame  % (3*2)
quiver3(Fytran_cop(1,1),Fytran_cop(2,1),Fytran_cop(3,1),...                                                               % ���
             Fytran_cop(1,2),Fytran_cop(2,2),Fytran_cop(3,2),0.5,'b-','LineWidth',2.0); % ��������0.5      % �յ�
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (b) ���ҵ㡪������ʧ
x_mid=[x4,0];
y_mid=[y4,F_yadd(i,1)];  % ���ҵ���ʧ����㡪��>�յ�
z_mid=[z4,0];
F_yadd_mid=[x_mid;y_mid;z_mid];  % size(wing_chord)   % (3*2)
Fyadd_mid= R_w2b*F_yadd_mid;        % wing_chord in body frame  % (3*2)
quiver3(Fyadd_mid(1,1),Fyadd_mid(2,1),Fyadd_mid(3,1),...                                                               % ���
             Fyadd_mid(1,2),Fyadd_mid(2,2),Fyadd_mid(3,2),0.5,'g-','LineWidth',2.0); % ��������0.5      % �յ�
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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