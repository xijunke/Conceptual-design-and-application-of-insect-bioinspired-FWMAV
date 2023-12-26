%% Stick figure force
clear all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ��һ���֡������ó���˶�ѧ-���˶����ɺͼ��ι���(AOA)������
% wing_m_output=[t',phi_pres',psi_sim',alpha',dphi',dpsi',ddphi',ddpsi',C_L',C_D',C_N1',C_T'];
wing_kenimatics=kenimatics_wing_and_AoA_fruitfly_simpres(); 
%���ú���kenimatics_wing_and_AoA; % size(wing_kenimatics)    % (1000,12)
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
% C_T=wing_kenimatics(:,12);  
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
phi1=phi_pres((locat_0:num),1);             % �Ĵ�ǡ�����λ��rad
psi1=psi_sim((locat_0:num),1);               % Ťת�ǡ�����λ��rad
alpha1=alpha_sim((locat_0:num),1);       % Ťת�ǡ�����λ��rad%  alpha=atan2(omega_z,-omega_y);  %������������и�
dphi1=dphi_pres((locat_0:num),1);        % �Ĵ�ǡ�����λ��rad/s
dpsi1=dpsi_sim((locat_0:num),1);           % Ťת�ǡ�����λ��rad/s
ddphi1=ddphi_pres((locat_0:num),1);     % �Ĵ�ǡ�����λ��rad/s^2
ddpsi1=ddpsi_sim((locat_0:num),1);       % Ťת�ǡ�����λ��rad/s^2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ����(1)���������������Ѿ��޸�,����������������L_arm=C_bem/2-L_arm1; 
% (3)�������ȶ������ڵ���(�ϳ�)�������������ݶ���(mN)���������������Ӧ���ǵ���ƬԪ������
% % �Ա�thr_dim_chord2����ĵ������ڻ������ͼ����������ϵ�µķ�����ʧ�յ�y����
% B=[F_ytran,F_yadd1,F_yrot];    % size(B)  % (2000*1)     %����mN
F_normal=xlsread('D:\KXJ\PassiveRot_dynamic_Science_fruitfly\wing_parameter\datanalysis_science_fruitfly\wing_model\Stick_figure\Forcenormal_oneT_simpres.xlsx','A1:C1000');
% F_norm=F_normal((locat_0:num),1);    %  ����С��0.1����ʾquiver3
F_ytran=F_normal((locat_0:num),1);         %  ����С��0.1����ʾquiver3 
F_yadd=F_normal((locat_0:num),2);         %  ����С��0.1����ʾquiver3
F_yrot=F_normal((locat_0:num),3);          %  ����С��0.1����ʾquiver3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(2)
% plot(t1*f,psi1*180/pi,'r*',t1*f,phi1*180/pi,'g>',t1*f,alpha1*180/pi,'k-',...
%        t1*f,3*F_ytran,'b-',t1*f,3*F_yadd,'c-',t1*f,3*F_yrot,'m-','LineWidth',1.5)  % ת��Ϊ����
% xlabel('\itNormalized time');  
% ylabel('\psi_{sim}(t) and \phi_{pres}(t) and \alpha_{sim}(t) and 3*F_{norm}(t))');  
% legend('\psi_{sim}(t)','\phi_{pres}(t)','\alpha_{sim}(t)','3*F_{y,tran}(t))','3*F_{y,add}(t))','3*F_{y,rot}(t))');  
% title('����ת����\psi_{sim}(t), �Ĵ��phi_{pres}(t) and 3*F_{norm}(t))��ʱ��ı仯')  
% grid on  % ����ת����psi_sim(t)��Ťת��phi_pres(t)�ͷ���������������ʱ��ı仯
% % axis([1.6,2.7,-110,110])
% axis([0.95,2.05,-inf,inf])
% set(gca,'XTick',(0.95:0.05:2.05))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% �ڶ����֡��� ���Ƶ�Ƭ����΢Ԫ�Ķ�ά�ռ�λ��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % R_wingeff =3.0040;        % ��λ�� mm   
% % xr=0.3289;                     % x-root offset  \mm
% % R_eff=R_wingeff+xr;
% ����(2)�������������鲻ȥƽ���ҳ���ΪƬ���ĳ��ȣ�Ӧ�ò�������չ�򳤶ȶ�Ӧ��Ƭ���ҳ������������% �����ҳ�C_bem=1.1257mm;   
% C_avereff =0.8854;          % ��λ�� mm
% C_025=0.25*C_avereff;
% C_075=0.75*C_avereff;
%%%%%%%%%%%%%%%%%%%%%%%%%
%% ��Ի����������Ӧ�ĳ�Ƭ����չ�򳤶�R_ref���ҳ�C_bem��������wing_model_88_translation_yaxis����
R_wingeff =3.0040;        % ��λ�� mm
xr=0.3289;  % x-root offset  \mm      % R_wingeff=wing_para(1,1);  % R_wingeff =3.0040;   % ��λ�� mm
R_ref=(xr+0.7*R_wingeff);                  % mm��������չ�򳤶�R_ref      
f_x_lead_ref =0.4598; %  mm % f_x_lead_ref=f_x_lead1(R_ref);       
f_x_trail_ref =-0.6658; %  mm %  f_x_trail_ref=f_x_trail1(R_ref); 
C_ref=f_x_lead_ref-f_x_trail_ref;     % C_ref =1.1257; %  mm
%%%%%%%%%%%%%%%%
% yr_leadbem=1.0958;                         % mm
% x_mod_Root =0.6360;                       % mm
% L_arm1=yr_leadbem-x_mod_Root;   % L_arm1=0.4598;     % ǰԵ��Ťת��ľ���:  mm
% C_bem =1.1257;                                % C_bem =1.1257;    % �����ҳ�C_bem:  mm
% L_arm=C_bem/2-L_arm1;                  % % L_arm =0.1030; % ������ƬԪ���۱�������(/2):  mm
C_avereff =C_ref;          % ��λ�� mm
C_025=f_x_lead_ref;
C_075=-f_x_trail_ref;  % C_075ȡΪ����
C_bem=C_ref;
C_05=C_bem/2-f_x_lead_ref;   % ���ҵ㵽Ťת��ľ��롪���޸������ʵ�ʵ������ҳ������ҵ㣬����0.25*C_avereff�����ҵ�
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  �ֽ����³��
N_skip=12;  % һ������ȡ28���㣬�������ȡ14����
num=338;  % num=336+2;
k_index=(1:N_skip:num); %  length=length(k_index)  % (1*28)  
Nstep=num/2;                        % Nstep=169;
% k_index1=(1:12:Nstep);       % 1...13..25......169  % k_index1(1,1:15);  
% N_down=length(k_index1);      % (15)
% k_index2=(Nstep:12:num);  % 169...181...193......337  % k_index2(1,15:29);
% N_up=length(k_index2); % (15)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ťת������ꡪ����Ϊ�趨��ʱ����롪���ϳ��
Y_axisup=linspace(-2.5,2.5,15); 
Z_axisup=zeros(1,length(Y_axisup));
%  Z_axisup=zeros(1,length(Y_axisup))+2;
% Ťת������ꡪ����Ϊ�趨��ʱ����롪���³��
Y_axisdown=fliplr(Y_axisup); 
Z_axisdown=fliplr(Z_axisup);
% Z_axisdown=zeros(1,length(Y_axisdown))+2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (1) �³��
k_index1=(1:12:Nstep);         % 
N_down=length(k_index1);        % (15)
% t_down=t1(k_index1);    % size(t_down)  %  (15*1)
psi_down=psi1(k_index1);  % ���һ������֮�ڵ����ݣ�ÿ��8����ȡһ�������γ�����
alpha_down=alpha1(k_index1);
% figure(3)
% plot(t_down*f,psi_down*180/pi,'rd','LineWidth',2)      % ת��Ϊ����
% Forcenormaldown=Forcenormal(k_index1);
F_ytrandown=F_ytran(k_index1);  % size(F_ytrandown)   % (15*1)     %  10^(-5)*
F_yadddown=F_yadd(k_index1);   % size(F_yadddown)   % (15*1)
F_yrotdown=F_yrot(k_index1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% �����ϳ�̡���������Ťת
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i = 1:1:N_down   % N_down=15
%     % (1) ǰԵ�����ꡪ��ʵ��Բ�����
%    y1=Y_axisdown(1,i);
%    z1=C_025+Z_axisdown(1,i);
%     % (2) ����ѹ�ĵ㡪������ʵ��Բ����� 
%     alpha=pi/2;             % Ťת�ǡ�����λ��rad
%     d_cp=(0.82*abs(alpha)/pi+0.05)*C_avereff;        
%     if d_cp>C_025
%         y3=Y_axisdown(1,i);
%         z3=-(d_cp-C_025)+Z_axisdown(1,i);
%     elseif d_cp<=C_025;
%         y3=Y_axisdown(1,i);
%         z3=(C_025-d_cp)+Z_axisdown(1,i);
%     end
%     % (3) ���ҵ����ꡪ��ʧ����㡪ʵ��Բ����
%     y4=Y_axisdown(1,i);
%     z4=-C_025+Z_axisdown(1,i);
%     % (4) ��Ե�����ꡪ��ʵ��Բ�����
%     y5=Y_axisdown(1,i);
%     z5=-C_075+Z_axisdown(1,i);
%     % ����ƬԪ������������ķֲ�
%     figure(4)
%     plot([y1,y5],[z1,z5],'k-','LineWidth',1.5)          %  ����ֱ�ߡ�������ǰ��Ե
%     hold on
%     plot(Y_axisdown(1,i),Z_axisdown(1,i),'ro','LineWidth',1.5,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',3) % Ťת���ġ���ʵ��Բ�����
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
% %     % displacement=[Y_axisdown(1,i);Z_axisdown(1,i)];  % (2*1)
% %     y_cop=-F_ytrandown(i,1)+Y_axisdown(1,i);    %����ע�����
% % %     z_cop=z3+Z_axisdown(1,i);
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
% N_step=N_down;         % N_down=15; 
% np=1;                    %  �����ܵ�֡����ʱ�䲽��
N_step=2*N_down;         % N_down=15; 
np=2;                    %  �����ܵ�֡����ʱ�䲽��
Nframe=np*15;     % number of movie frames -------------------------- for change
Nskip1=floor(N_step/Nframe);
Nspeed=3;  % number of frames per second ------------------------ for change
duration=Nframe/Nspeed;  % duration of movie  % 5 seconds @Nframe=1*15
animation=1;         %  1=save animation file --------------------------- for change
display(['duratin of movie will be ', num2str(duration), ' seconds ']); %�����ʾ����������ʱ��(s)
%�������: duratin of movie will be 5 seconds 
if animation == 1
%     aviobj = avifile('Stick_figure_force_downstroke.avi','fps',Nspeed); % filename --- for change
        aviobj = avifile('Stick_figure_force_one_stroke.avi','fps',Nspeed); % filename --- for change
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig=figure;
for i = 1:1:N_down   % N_down=15
    % ����һ��������ϵ֮��ı任
    psi=psi_down(i,1);               % Ťת�ǡ�����λ��rad    %����ע��Ϊ����
%     psi_deg=psi_down(i,1)*180/pi  
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
    % alpha_up=alpha1(i,1);             % Ťת�ǡ�����λ��rad
    d_cp=(0.82*abs(alpha_down(i,1))/pi+0.05)*C_avereff;        
    if d_cp>C_025
        y3=0;
        z3=-(d_cp-C_025);% ��ֵ
    elseif d_cp<=C_025;
        y3=0;
        z3=C_025-d_cp; % ��ֵ
    end
    p_cop0=[y3;z3];
    p_coprot=R_w2s*p_cop0;
    p_cop=p_coprot+displacement; 
    % (3) ���ҵ����ꡪ��ʧ����㡪ʵ��Բ����
    y4=0;
    % z4=-C_025; % ���ҵ㵽Ťת��ľ���
    z4=-C_05;   % ���ҵ㵽Ťת��ľ��롪���޸������ʵ�ʵ������ҳ������ҵ㣬����0.25*C_avereff�����ҵ�
    p_mid0=[y4;z4];
    p_midrot=R_w2s*p_mid0;
    p_mid=p_midrot+displacement; 
    % ƬԪ������������
    %  y=[p_lead(1,1),Y_axisdown(1,i),p_cop(1,1),p_mid(1,1),p_tail(1,1)];
    %  z=[p_lead(2,1),Z_axisdown(1,i),p_cop(2,1),p_mid(2,1),p_tail(2,1)];
     %  wing_C=[y;z];  % size(wing_chord)   % (3*5)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % figure(5)
    % figure(i)
    plot([Y_axisdown(1,1),Y_axisdown(1,N_down)],[Z_axisdown(1,1),Z_axisdown(1,N_down)],'r:','LineWidth',1.5) % Ťת�������ߡ������߻���
    hold on
    xlabel('����(y)')
    ylabel('��ֱ������(z)')
    grid on
    % % axis([-4,4,-1,3.1])     % �����������ڷ����������ı�������ʾ
    axis([-4,4,-1.9,2.1])   % axis([-4,4,-0.6808,1.4952]) 
    hold on
    % (1) ����ֱ�ߡ�������ǰ��Ե
    plot([p_lead(1,1),p_tail(1,1)],[p_lead(2,1),p_tail(2,1)],'k-','LineWidth',2) % ����ǰԵ��Եʵ��
    hold on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (2) Ťת���ġ���ʵ��Բ�����
    plot(Y_axisdown(1,i),Z_axisdown(1,i),'ro','LineWidth',1.5,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',3) % Ťת���ġ���ʵ��Բ�����
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
    % displacement=[Y_axisdown(1,i);Z_axisdown(1,i)];          % (2*1)
    % y_cop=-F_ytrandown(i,1)+Y_axisdown(1,i);                     %����ע����š���������ƽ������������
    y_cop=(F_ytrandown(i,1)+F_yrotdown(i,1))+Y_axisdown(1,i); %����ע����š�������ƽ��������������ת����������������    
    z_cop=z3+Z_axisdown(1,i);
    F_ytran_cop0=[y_cop;z_cop];                           % size(wing_chord)   % (2*2)    % ����ѹ�ĵ���ʧ���յ�
    Fytran_coprot= R_w2s*F_ytran_cop0;              % wing_chord in body frame  % (2*2)     %��R_w2s'��ע�����
    Fytran_cop=Fytran_coprot; 
    quiver(p_cop(1,1),p_cop(2,1),...                        % ���
               Fytran_cop(1,1),Fytran_cop(2,1),0.1,'b-','LineWidth',2);    % �յ�   % ��������0.05������������ע����ʧ����С��0.14�� 
    hold on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % (5) ���ҵ㡪������ʧ
   % plot(y4,z4,'go','LineWidth',1.5,'MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',3);  % ���ҵ㡪��ʧ����㡪ʵ��Բ����
   % hold on  
    y_mid=F_yadddown(i,1)+Y_axisdown(1,i);  % ���ҵ���ʧ����㡪��>�յ�   % ������һ����
    % y_mid=F_yadddown(i,1)+Y_axisdown(1,i);  % ���ҵ���ʧ����㡪��>�յ�        % �����ڶ�����
    z_mid=z4+Z_axisdown(1,i);
    F_yadd_mid0=[y_mid;z_mid];  % size(wing_chord)   % (2*2)
    Fyadd_midrot=R_w2s*F_yadd_mid0;        % wing_chord in body frame  % (2*2)
    Fyadd_mid=Fyadd_midrot; 
    quiver(p_mid(1,1),p_mid(2,1),...                                                      % ���
               Fyadd_mid(1,1),Fyadd_mid(2,1),0.1,'g-','LineWidth',2);   % �յ�   % ��������0.05������������ע����ʧ����С��0.135�� 
    hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ����Ƭ�γ�������Ρ���(2)
    if animation == 1
        frame = getframe(fig);
        aviobj = addframe(aviobj,frame);
    else
        pause(0.5)
    end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
% close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  ����Ƭ�γ�������Ρ���(3)
% if animation == 1
%     aviobj = close(aviobj);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (2) �ϳ��
k_index2=(Nstep:12:num);     % 169...181...193......337
N_up=length(k_index2);    % (15)
% t_up=t1(k_index2);    %  size(t_up)  %  (15*1)
psi_up=psi1(k_index2);  % ���һ������֮�ڵ����ݣ�ÿ����8����ȡһ�������γ�����
alpha_up=alpha1(k_index2);
% figure(6)
% plot(t_up*f,psi_up*180/pi,'rd','LineWidth',2)      % ת��Ϊ����
% Forcenormalup=Forcenormal(k_index2);
F_ytranup=F_ytran(k_index2);  % size(F_ytranup)   % (15*1)     %  10^(-5)*
F_yaddup=F_yadd(k_index2);   % size(F_yaddup)   % (15*1)
F_yrotup=F_yrot(k_index2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % animation���������������� ����Ƭ�γ���ʼ�Ρ���(1)
% Nstep=N_up;    % N_up=15; 
% np=1;                     %  �����ܵ�֡����ʱ�䲽��
% Nframe=np*15;   % number of movie frames -------------------------- for change
% Nskip2=floor(Nstep/Nframe);
% Nspeed=3;  % number of frames per second ------------------------ for change
% animation=1;       %  1=save animation file --------------------------- for change
% duration=Nframe/Nspeed;  % duration of movie  % 80 seconds @Nframe=4*40
% display(['duratin of movie will be ', num2str(duration), ' seconds ']); %�����ʾ����������ʱ��(s)
% %�������: duratin of movie will be 10 seconds 
% if animation == 1
%     aviobj = avifile('Stick_figure_force_upstroke.avi','fps',Nspeed); % filename --- for change
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fig=figure;
for i = 1:1:N_up  % N_up=15
    % ����һ��������ϵ֮��ı任
    psi=psi_up(i,1);             % Ťת�ǡ�����λ��rad    %����ע��Ϊ����
    Pitch = [cos(psi)   sin(psi); ...
                 -sin(psi)    cos(psi)];        % pitch matrix�������ſ�����Ҫ����
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
    %  alpha_up=alpha1(i,1);             % Ťת�ǡ�����λ��rad
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
    % z4=-C_025;
    z4=-C_05;  % �޸������ʵ�ʵ������ҳ������ҵ㣬����0.25*C_avereff�����ҵ�
    p_mid0=[y4;z4];
    p_midrot=R_w2s*p_mid0;
    p_mid=p_midrot+displacement; 
    % ƬԪ������������
%     y=[p_lead(1,1),Y_axisup(1,i),p_cop(1,1),p_mid(1,1),p_tail(1,1)];
%     z=[p_lead(2,1),Z_axisup(1,i),p_cop(2,1),p_mid(2,1),p_tail(2,1)];
%     wing_C=[y;z];  % size(wing_chord)   % (3*5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % figure(7)
    %  figure(i)
    plot([Y_axisup(1,1),Y_axisup(1,N_up)],[Z_axisup(1,1),Z_axisup(1,N_up)],'r:','LineWidth',1.5) % Ťת���ġ���ʵ��Բ�����
    hold on
    xlabel('����(y)')
    ylabel('��ֱ������(z)')
    grid on
    % % axis([-4,4,-1,3.1])   % �����������ڷ����������ı�������ʾ
    axis([-4,4,-1.9,2.1])    % axis([-4,4,-0.6808,1.4952]) 
    hold on
    % (1) ����ֱ�ߡ�������ǰ��Ե
    plot([p_lead(1,1),p_tail(1,1)],[p_lead(2,1),p_tail(2,1)],'k-','LineWidth',2) % ����ǰԵ��Եʵ��
    hold on
    % (2) Ťת���ġ���ʵ��Բ�����
    plot(Y_axisup(1,i),Z_axisup(1,i),'ro','LineWidth',1.5,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',3) % Ťת���ġ���ʵ��Բ�����
    hold on
     % (3) ����ǰԵʵ��Բ��
    plot(p_lead(1,1),p_lead(2,1),'ko-','LineWidth',1.5,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5)
    hold on
   % (4) ����ѹ�ĵ㡪������ʧ
    plot(p_cop(1,1),p_cop(2,1),'bo','LineWidth',1.7,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',5);  % ����ѹ�ĵ㡪������ʵ��Բ�����
    hold on
    % p_cop=p_coprot+displacement;                 % ����ѹ�ĵ���ʧ�����
    % displacement=[Y_axisup(1,i);Z_axisup(1,i)];          % (2*1)
    % y_cop=F_ytranup(i,1)+Y_axisup(1,i);                     %����ע����š���������ƽ������������
    y_cop=(F_ytranup(i,1)+F_yrotup(i,1))+Y_axisup(1,i); %����ע����š�������ƽ��������������ת����������������
    % y_cop=(F_ytranup(i,1)+F_yrotup(i,1))+Y_axisup(1,i); %����ע����š�������ƽ��������������ת����������������
    z_cop=z3+Z_axisup(1,i);
    F_ytran_cop0=[y_cop;z_cop];                           % size(wing_chord)   % (2*2)    % ����ѹ�ĵ���ʧ���յ�
    Fytran_coprot= R_w2s*F_ytran_cop0;              % wing_chord in body frame  % (2*2)     %��R_w2s'��ע�����
    Fytran_cop=Fytran_coprot; 
    quiver(p_cop(1,1),p_cop(2,1),...                        % ���
               Fytran_cop(1,1),Fytran_cop(2,1),0.1,'b-','LineWidth',2);    % �յ�   % ��������0.05������������ע����ʧ����С��0.14�� 
    hold on
   % (5) ���ҵ㡪������ʧ
   % plot(y4,z4,'go','LineWidth',1.5,'MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',3);  % ���ҵ㡪��ʧ����㡪ʵ��Բ����
   % hold on  
    y_mid=F_yaddup(i,1)+Y_axisup(1,i);  % ���ҵ���ʧ����㡪��>�յ�   % ������һ����
    % y_mid=F_yaddup(i,1)+Y_axisup(1,i);  % ���ҵ���ʧ����㡪��>�յ�        % �����ڶ�����
    z_mid=z4+Z_axisup(1,i);
    F_yadd_mid0=[y_mid;z_mid];  % size(wing_chord)   % (2*2)
    Fyadd_midrot=R_w2s*F_yadd_mid0;        % wing_chord in body frame  % (2*2)
    Fyadd_mid=Fyadd_midrot; 
    quiver(p_mid(1,1),p_mid(2,1),...                                                      % ���
               Fyadd_mid(1,1),Fyadd_mid(2,1),0.1,'g-','LineWidth',2);   % �յ�   % ��������0.05������������ע����ʧ����С��0.135�� 
    hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ����Ƭ�γ�������Ρ���(2)
    if animation == 1
        frame = getframe(fig);
        aviobj = addframe(aviobj,frame);
    else
        pause(0.5)
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
% close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  ����Ƭ�γ�������Ρ���(3)
if animation == 1
    aviobj = close(aviobj);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
