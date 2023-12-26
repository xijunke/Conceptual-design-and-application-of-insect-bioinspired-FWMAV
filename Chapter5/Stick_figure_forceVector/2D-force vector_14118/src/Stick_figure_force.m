%% Stick figure force
clear all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ��һ���֡������ó���˶�ѧ-���˶����ɺͼ��ι���(AOA)������
% wing_m_output=[t',phi',psi',alpha',dphi',dpsi',ddphi',ddpsi',C_L',C_D',C_N1'];
wing_kenimatics=kenimatics_wing_and_AoA_fruitfly_simpres(); %���ú���kenimatics_wing_and_AoA; % size(wing_kenimatics) % (1000,11)
% wing_kenimatics=kenimatics_wing_and_AoA_fruitfly_exp();
t=wing_kenimatics(:,1);                        % ��λ��ms
phi_pres=wing_kenimatics(:,2);            % �Ĵ�ǡ�����λ��rad
psi_sim=wing_kenimatics(:,3);              % �Ĵ�ǡ�����λ��rad
alpha_sim=wing_kenimatics(:,4);          % alpha=atan2(omega_z,-omega_y);  % ���ι��ǡ���������   %������������и�
dphi_pres=wing_kenimatics(:,5);          % ��λ��rad/s
dpsi_sim=wing_kenimatics(:,6);            % ��λ��rad/s
ddphi_pres=wing_kenimatics(:,7);        % ��λ��rad/s^2
ddpsi_sim=wing_kenimatics(:,8);          % ��λ��rad/s^2
% C_L=wing_kenimatics(:,9);          
% C_D=wing_kenimatics(:,10);     
% C_N1=wing_kenimatics(:,11);   
% C_N=C_N1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
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
% (2) ��ȡ�����ȶ������ڡ����Թ��Ǳ仯����Ϊ����
% ���ù���Ϊ90�ȶ�Ӧʱ��Ϊ���³�̵Ľ���ʱ���
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
% �պ�104�����ݵ�=275-172+1;
t1=t((locat_0:num),1);                             % ��λ��ms
psi1=psi_sim((locat_0:num),1);               % Ťת�ǡ�����λ��rad
alpha1=alpha_sim((locat_0:num),1);       % Ťת�ǡ�����λ��rad
dpsi1=dpsi_sim((locat_0:num),1);           % Ťת�ǡ�����λ��rad
% ddpsi1=ddpsi_sim((locat_0:num),1);  % Ťת�ǡ�����λ��rad
ddpsi1=dpsi_sim((locat_0:num),1);     % Ťת�ǡ�����λ��rad��ע�������Ťת�Ǽ��ٶ�ȡ��Ťת���ٶ�ddpsi1=dpsi_sim, ��ΪODE����û�����ddpsi
phi1=phi_pres((locat_0:num),1);         % �Ĵ�ǡ�����λ��rad
dphi1=dphi_pres((locat_0:num),1);      % �Ĵ�ǡ�����λ��rad
ddphi1=ddphi_pres((locat_0:num),1);  % �Ĵ�ǡ�����λ��rad
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (3) �������ȶ������ڵ���(�ϳ�)�������������ݶ��롪��mN
% �Ա�thr_dim_chord2����ĵ������ڻ������ͼ����������ϵ�µķ�����ʧ�յ�y����
% B=[F_ytran,F_yadd1,F_yrot];    % size(B)  % (2000*1)     %����mN
F_normal=xlsread('Forcenormal_oneT_simpres.xlsx','A1:C1000');
% Forcenormal=F_normal((locat_0:num),1);   %  ����С��1/5=0.2����ʾquiver3
F_ytran=F_normal((locat_0:num),1);                % 10^(-1)*
F_yadd=F_normal((locat_0:num),2);
% F_yrot=F_normal((locat_0:num),3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(2)
% plot(t1*f,psi1*180/pi,'r-',t1*f,phi1*180/pi,'g-',t1*f,alpha1*180/pi,'k-',t1*f,2*abs(F_ytran/0.2),'b-','LineWidth',2)  % ת��Ϊ����
% xlabel('\itNormalized time');  
% ylabel('\psi_{sim}(t) and \phi_{pres}(t) and \alpha_{sim}(t) and 2*abs(Force_{normal}(t))');  
% legend('\psi_{sim}(t)','\phi_{pres}(t)','\alpha_{sim}(t)','2*abs(Force_{normal}(t))');  
% title('����ת����\psi_{sim}(t), �Ĵ��phi_{pres}(t) and 2*abs(Force_{normal}(t))��ʱ��ı仯')  
% grid on  % ����ת����psi_sim(t)��Ťת��phi_pres(t)��ʱ��ı仯
% % axis([1.6,2.7,-110,110])
% axis([0.9,2.1,-inf,inf])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% �ڶ����֡��� ���Ƶ�Ƭ����΢Ԫ�Ķ�ά�ռ�λ��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % R_wingeff =3.0040;        % ��λ�� mm   
% % xr=0.3289;                     % x-root offset  \mm
% % R_eff=R_wingeff+xr;
C_avereff =0.8854;          % ��λ�� mm
C_025=0.25*C_avereff;
C_075=0.75*C_avereff;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  �ֽ����³��
N_skip=12;  % һ������ȡ28���㣬�������ȡ14����
num=338;  % num=336+2;
k_index=(1:N_skip:num); %  length=length(k_index)  % (1*28)  
Nstep=num/2;                        % Nstep=169;
% k_index1=(1:12:Nstep);       % 1...13..25......169  % k_index1(1,1:15);  
% N_up=length(k_index1);      % (15)
% k_index2=(Nstep:12:num);  % 169...181...193......337  % k_index2(1,15:29);
% N_down=length(k_index2); % (15)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (1) �ϳ��
k_index1=(1:12:Nstep);         % 
N_up=length(k_index1);        % (15)
% t_up=t1(k_index1);    % size(t_up)  %  (15*1)
psi_up=psi1(k_index1);  % ���һ������֮�ڵ����ݣ�ÿ��8����ȡһ�������γ�����
alpha_up=alpha1(k_index1);
% figure(3)
% plot(t_up*f,psi_up*180/pi,'rd','LineWidth',2)      % ת��Ϊ����
% Forcenormalup=Forcenormal(k_index1);
F_ytranup=F_ytran(k_index1);  % size(F_ytranup)   % (15*1)     %  10^(-5)*
F_yaddup=F_yadd(k_index1);   % size(F_yaddup)   % (15*1)
% Ťת������ꡪ����Ϊ�趨��ʱ�����
Y_axisup=linspace(-2.5,2.5,15); 
Z_axisup=zeros(1,length(Y_axisup));
% Z_axisup=zeros(1,length(Y_axisup))+2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% �����ϳ�̡���������Ťת
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i = 1:1:N_up   % N_up=15
%     % (1) ǰԵ�����ꡪ��ʵ��Բ�����
%    y1=Y_axisup(1,i);
%    z1=C_025+Z_axisup(1,i);
%     % (2) ����ѹ�ĵ㡪������ʵ��Բ����� 
%     alpha=pi/2;             % Ťת�ǡ�����λ��rad
%     d_cp=(0.82*abs(alpha)/pi+0.05)*C_avereff;        
%     if d_cp>C_025
%         y3=Y_axisup(1,i);
%         z3=-(d_cp-C_025)+Z_axisup(1,i);
%     elseif d_cp<=C_025;
%         y3=Y_axisup(1,i);
%         z3=(C_025-d_cp)+Z_axisup(1,i);
%     end
%     % (3) ���ҵ����ꡪ��ʧ����㡪ʵ��Բ����
%     y4=Y_axisup(1,i);
%     z4=-C_025+Z_axisup(1,i);
%     % (4) ��Ե�����ꡪ��ʵ��Բ�����
%     y5=Y_axisup(1,i);
%     z5=-C_075+Z_axisup(1,i);
%     % ����ƬԪ������������ķֲ�
%     figure(4)
%     plot([y1,y5],[z1,z5],'k-','LineWidth',1.5)          %  ����ֱ�ߡ�������ǰ��Ե
%     hold on
%     plot(Y_axisup(1,i),Z_axisup(1,i),'ro','LineWidth',1.5,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',3) % Ťת���ġ���ʵ��Բ�����
%     hold on
%     plot(y1,z1,'ko','LineWidth',1.5,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',4);   % ǰԵ����ʵ��Բ�����
%     hold on
%     plot(y4,z4,'go','LineWidth',1.5,'MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',3);  % ���ҵ㡪��ʧ����㡪ʵ��Բ����
%     hold on
%     % plot(y5,z5,'ko','LineWidth',1.5,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',1.5)   % ��Ե����ʵ��Բ�����
%     % hold on
%     % (4) ����ѹ�ĵ㡪������ʧ       
%        % ����ѹ�ĵ���ʧ�����
%     plot(y3,z3,'bo','LineWidth',1.5,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',3);  % ����ѹ�ĵ㡪������ʵ��Բ�����
%     hold on
% %     % p_cop=p_coprot+displacement;  
% %     % displacement=[Y_axisup(1,i);Z_axisup(1,i)];  % (2*1)
% %     y_cop=-F_ytranup(i,1)+Y_axisup(1,i);    %����ע�����
% % %     z_cop=z3+Z_axisup(1,i);
% %      z_cop=0;
% %     F_ytran_cop0=[y_cop;z_cop];  % size(wing_chord)   % (2*2)    % ����ѹ�ĵ���ʧ���յ�
% % %     Fytran_coprot= R_w2s*F_ytran_cop0;        % wing_chord in body frame  % (2*2)     %��R_w2s'��ע�����
% %     Fytran_coprot= F_ytran_cop0;        % wing_chord in body frame  % (2*2)     %��R_w2s'��ע�����
% %     Fytran_cop=Fytran_coprot; 
% %     quiver(y3,z3,...                                                                 % ���
% %                Fytran_cop(1,1),Fytran_cop(2,1),0.5,'b-','LineWidth',2.0); % ��������0.5        % �յ�
% %     hold on
% end
% xlabel('����(y)')
% ylabel('��ֱ������(z)')
% grid on
% axis([-4,4,-3,3])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% animation���������������� ����Ƭ�γ���ʼ�Ρ���(1)
N_step=N_up;         % N_up=15; 
animation=1;         %  1=save animation file --------------------------- for change
np=1;                    %  �����ܵ�֡����ʱ�䲽��
Nframe=np*15;     % number of movie frames -------------------------- for change
Nskip1=floor(N_step/Nframe);
Nspeed=3;  % number of frames per second ------------------------ for change
duration=Nframe/Nspeed;  % duration of movie  % 5 seconds @Nframe=1*15
display(['duratin of movie will be ', num2str(duration), ' seconds ']); %�����ʾ����������ʱ��(s)
%�������: duratin of movie will be 5 seconds 
if animation == 1
    aviobj = avifile('Stick_figure_force_upstroke.avi','fps',Nspeed); % filename --- for change
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:1:N_up   % N_up=15
    % ����һ��������ϵ֮��ı任
    psi=-psi_up(i,1);             % Ťת�ǡ�����λ��rad    %����ע�����
    Pitch = [cos(psi)   -sin(psi); ...
                 sin(psi)    cos(psi)];        % pitch matrix�������ſ�����Ҫ����
    R_w2s =Pitch;     % from wing frame to spar frame�����ӳ�����ϵ��������ϵ
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    displacement=[Y_axisup(1,i);Z_axisup(1,i)];  % (2*1)
   % (1) ǰԵ�����ꡪ��ʵ��Բ�����
    y1=0;
    z1=C_025;
    p_lead0=[y1;z1];                   % (2*1)
    p_leadrot=R_w2s*p_lead0;    % (2*2)*(2*1)=(2*1)
    p_lead=p_leadrot+displacement;          % (2*1)
   % (4) ��Ե�����ꡪ��ʵ��Բ�����
    y5=0;
    z5=-C_075;
    p_tail0=[y5;z5];                   % (2*1)
    p_tailrot=R_w2s*p_tail0;      % (2*2)*(2*1)=(2*1)
    p_tail=p_tailrot+displacement;          % (2*1)   
    % (2) ����ѹ�ĵ㡪������ʵ��Բ����� 
    % alpha_up=alpha1(i,1);             % Ťת�ǡ�����λ��rad
    d_cp=(0.82*abs(alpha_up(i,1))/pi+0.05)*C_avereff;        
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
    % (3) ���ҵ����ꡪ��ʧ����㡪ʵ��Բ����
    y4=0;
    z4=-C_025;
    p_mid0=[y4;z4];
    p_midrot=R_w2s*p_mid0;
    p_mid=p_midrot+displacement; 
    % ƬԪ������������
%     y=[p_lead(1,1),Y_axisup(1,i),p_cop(1,1),p_mid(1,1),p_tail(1,1)];
%     z=[p_lead(2,1),Z_axisup(1,i),p_cop(2,1),p_mid(2,1),p_tail(2,1)];
%     wing_C=[y;z];  % size(wing_chord)   % (3*5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(5)
    plot([Y_axisup(1,1),Y_axisup(1,N_up)],[Z_axisup(1,1),Z_axisup(1,N_up)],'r:','LineWidth',1.5) % Ťת���ġ���ʵ��Բ�����
    hold on
    xlabel('����(y)')
    ylabel('��ֱ������(z)')
    grid on
    axis([-4,4,-1,3.1])
    hold on
    % (1) ����ֱ�ߡ�������ǰ��Ե
    plot([p_lead(1,1),p_tail(1,1)],[p_lead(2,1),p_tail(2,1)],'k-','LineWidth',2) % ����ǰԵ��Եʵ��
    hold on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (2) Ťת���ġ���ʵ��Բ�����
    plot(Y_axisup(1,i),Z_axisup(1,i),'ro','LineWidth',1.5,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',3) % Ťת���ġ���ʵ��Բ�����
    hold on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (3) ����ǰԵʵ��Բ��
    plot(p_lead(1,1),p_lead(2,1),'ko-','LineWidth',1.5,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5)
    hold on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % (4) ����ѹ�ĵ㡪������ʧ
    plot(p_cop(1,1),p_cop(2,1),'bo','LineWidth',1.7,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',5);  % ����ѹ�ĵ㡪������ʵ��Բ�����
    hold on   
    % p_cop=p_coprot+displacement;                 % ����ѹ�ĵ���ʧ�����
    % displacement=[Y_axisup(1,i);Z_axisup(1,i)];          % (2*1)
    y_cop=-F_ytranup(i,1)+Y_axisup(1,i);                     %����ע�����
    z_cop=z3+Z_axisup(1,i);
    F_ytran_cop0=[y_cop;z_cop];                           % size(wing_chord)   % (2*2)    % ����ѹ�ĵ���ʧ���յ�
    Fytran_coprot= R_w2s*F_ytran_cop0;              % wing_chord in body frame  % (2*2)     %��R_w2s'��ע�����
    Fytran_cop=Fytran_coprot; 
    quiver(p_cop(1,1),p_cop(2,1),...                        % ���
               Fytran_cop(1,1),Fytran_cop(2,1),0.14,'b-','LineWidth',1.5);    % �յ�   % ��������0.05������������ע����ʧ����С��0.14�� 
    hold on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % (5) ���ҵ㡪������ʧ
   % plot(y4,z4,'go','LineWidth',1.5,'MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',3);  % ���ҵ㡪��ʧ����㡪ʵ��Բ����
   % hold on  
    y_mid=-F_yaddup(i,1)+Y_axisup(1,i);  % ���ҵ���ʧ����㡪��>�յ�
    z_mid=z4+Z_axisup(1,i);
    F_yadd_mid0=[y_mid;z_mid];  % size(wing_chord)   % (2*2)
    Fyadd_midrot=R_w2s*F_yadd_mid0;        % wing_chord in body frame  % (2*2)
    Fyadd_mid=Fyadd_midrot; 
    quiver(p_mid(1,1),p_mid(2,1),...                                                      % ���
               Fyadd_mid(1,1),Fyadd_mid(2,1),0.135,'g-','LineWidth',1.5);   % �յ�   % ��������0.05������������ע����ʧ����С��0.135�� 
    hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ����Ƭ�γ�������Ρ���(2)
    if animation == 1
        frame = getframe;
        aviobj = addframe(aviobj,frame);
    else
        pause(0.5)
    end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ����Ƭ�γ�������Ρ���(3)
if animation == 1
    aviobj = close(aviobj);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (2) �³��
k_index2=(Nstep:12:num);     % 169...181...193......337
N_down=length(k_index2);    % (15)
% t_down=t1(k_index2);    %  size(t_down)  %  (15*1)
psi_down=psi1(k_index2);  % ���һ������֮�ڵ����ݣ�ÿ����8����ȡһ�������γ�����
alpha_down=alpha1(k_index2);
% figure(6)
% plot(t_down*f,psi_down*180/pi,'rd','LineWidth',2)      % ת��Ϊ����
% Forcenormaldown=Forcenormal(k_index2);
F_ytrandown=F_ytran(k_index2);  % size(F_ytrandown)   % (15*1)     %  10^(-5)*
F_yadddown=F_yadd(k_index2);   % size(F_yadddown)   % (15*1)
% Ťת������ꡪ����Ϊ�趨��ʱ�����
% Y_axisdown=linspace(-2.5,2.5,15); 
% Z_axisdown=zeros(1,length(Y_axisdown));
Y_axisdown=fliplr(Y_axisup); 
Z_axisdown=fliplr(Z_axisup);
%%  Z_axisdown=zeros(1,length(Y_axisdown))+2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % animation���������������� ����Ƭ�γ���ʼ�Ρ���(1)
Nstep=N_down;    % N_down=15; 
animation=1;       %  1=save animation file --------------------------- for change
np=1;                     %  �����ܵ�֡����ʱ�䲽��
Nframe=np*15;   % number of movie frames -------------------------- for change
Nskip2=floor(Nstep/Nframe);
Nspeed=3;  % number of frames per second ------------------------ for change
duration=Nframe/Nspeed;  % duration of movie  % 80 seconds @Nframe=4*40
display(['duratin of movie will be ', num2str(duration), ' seconds ']); %�����ʾ����������ʱ��(s)
%�������: duratin of movie will be 10 seconds 
if animation == 1
    aviobj = avifile('Stick_figure_force_downstroke.avi','fps',Nspeed); % filename --- for change
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:1:N_down  % N_down=15
    % ����һ��������ϵ֮��ı任
    psi=psi_down(i,1);             % Ťת�ǡ�����λ��rad
    Pitch = [cos(psi)   sin(psi); ...
                 -sin(psi)    cos(psi)];        % pitch matrix�������ſ�����Ҫ����
    R_w2s =Pitch;     % from wing frame to spar frame�����ӳ�����ϵ��������ϵ
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   displacement=[Y_axisdown(1,i);Z_axisdown(1,i)];  % (2*1)
    % (1) ǰԵ�����ꡪ��ʵ��Բ�����
    y1=0;
    z1=C_025;
    p_lead0=[y1;z1];                   % (2*1)
    p_leadrot=R_w2s*p_lead0;    % (2*2)*(2*1)=(2*1)
    p_lead=p_leadrot+displacement;          % (2*1)
    % (4) ��Ե�����ꡪ��ʵ��Բ�����
    y5=0;
    z5=-C_075;
    p_tail0=[y5;z5];                   % (2*1)
    p_tailrot=R_w2s*p_tail0;      % (2*2)*(2*1)=(2*1)
    p_tail=p_tailrot+displacement;          % (2*1)
    % (2) ����ѹ�ĵ㡪������ʵ��Բ����� 
    %  alpha_down=alpha1(i,1);             % Ťת�ǡ�����λ��rad
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
    % (3) ���ҵ����ꡪ��ʧ����㡪ʵ��Բ����
    y4=0;
    z4=-C_025;
    p_mid0=[y4;z4];
    p_midrot=R_w2s*p_mid0;
    p_mid=p_midrot+displacement; 
    % ƬԪ������������
%     y=[p_lead(1,1),Y_axisdown(1,i),p_cop(1,1),p_mid(1,1),p_tail(1,1)];
%     z=[p_lead(2,1),Z_axisdown(1,i),p_cop(2,1),p_mid(2,1),p_tail(2,1)];
%     wing_C=[y;z];  % size(wing_chord)   % (3*5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(7)
    plot([Y_axisdown(1,1),Y_axisdown(1,N_down-1)],[Z_axisdown(1,1),Z_axisdown(1,N_down-1)],'r:','LineWidth',1.5) % Ťת���ġ���ʵ��Բ�����
    hold on
    xlabel('����(y)')
    ylabel('��ֱ������(z)')
    grid on
    axis([-4,4,-1,3.1])
    hold on
    % (1) ����ֱ�ߡ�������ǰ��Ե
    plot([p_lead(1,1),p_tail(1,1)],[p_lead(2,1),p_tail(2,1)],'k-','LineWidth',1.5) % ����ǰԵ��Եʵ��
    hold on
    % (2) Ťת���ġ���ʵ��Բ�����
    plot(Y_axisdown(1,i),Z_axisdown(1,i),'ro','LineWidth',1.5,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',3) % Ťת���ġ���ʵ��Բ�����
    hold on
     % (3) ����ǰԵʵ��Բ��
    plot(p_lead(1,1),p_lead(2,1),'ko-','LineWidth',1.5,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)
    hold on
   % (4) ����ѹ�ĵ㡪������ʧ
    plot(p_cop(1,1),p_cop(2,1),'bo','LineWidth',1.7,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',5);  % ����ѹ�ĵ㡪������ʵ��Բ�����
    hold on
    % p_cop=p_coprot+displacement;                 % ����ѹ�ĵ���ʧ�����
    % displacement=[Y_axisdown(1,i);Z_axisdown(1,i)];          % (2*1)
    y_cop=-F_ytrandown(i,1)+Y_axisdown(1,i);                     %����ע�����
    z_cop=z3+Z_axisdown(1,i);
    F_ytran_cop0=[y_cop;z_cop];                           % size(wing_chord)   % (2*2)    % ����ѹ�ĵ���ʧ���յ�
    Fytran_coprot= R_w2s*F_ytran_cop0;              % wing_chord in body frame  % (2*2)     %��R_w2s'��ע�����
    Fytran_cop=Fytran_coprot; 
    quiver(p_cop(1,1),p_cop(2,1),...                        % ���
               Fytran_cop(1,1),Fytran_cop(2,1),0.14,'b-','LineWidth',1.5);    % �յ�   % ��������0.05������������ע����ʧ����С��0.14�� 
    hold on
   % (5) ���ҵ㡪������ʧ
   % plot(y4,z4,'go','LineWidth',1.5,'MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',3);  % ���ҵ㡪��ʧ����㡪ʵ��Բ����
   % hold on  
    y_mid=-F_yadddown(i,1)+Y_axisdown(1,i);  % ���ҵ���ʧ����㡪��>�յ�
    z_mid=z4+Z_axisdown(1,i);
    F_yadd_mid0=[y_mid;z_mid];  % size(wing_chord)   % (2*2)
    Fyadd_midrot=R_w2s*F_yadd_mid0;        % wing_chord in body frame  % (2*2)
    Fyadd_mid=Fyadd_midrot; 
    quiver(p_mid(1,1),p_mid(2,1),...                                                      % ���
               Fyadd_mid(1,1),Fyadd_mid(2,1),0.135,'g-','LineWidth',1.5);   % �յ�   % ��������0.05������������ע����ʧ����С��0.135�� 
    hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ����Ƭ�γ�������Ρ���(2)
    if animation == 1
        frame = getframe;
        aviobj = addframe(aviobj,frame);
    else
        pause(0.5)
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  ����Ƭ�γ�������Ρ���(3)
if animation == 1
    aviobj = close(aviobj);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
