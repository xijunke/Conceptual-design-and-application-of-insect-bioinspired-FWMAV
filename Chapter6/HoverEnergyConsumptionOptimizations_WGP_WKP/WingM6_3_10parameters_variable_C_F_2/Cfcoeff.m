function C_periEdge=Cfcoeff(R_wing,C_aver)
% C_periEdge=C_perimeter_to_Edge_correction(R_wing,C_aver,xr0,C_maxyaxis);
% C_perimeter_to_Edge_correction�������ǳ�򼸺���ò����Ӱ���������ϵ��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x=[3.004,0.8854,0.3289,0.356737];
% R_wing=x(1);
% C_aver=x(2);
% xr0=x(3);
% C_maxyaxis=x(4);
%%%%%%%%%%%%%%%%%%%%%%%%%
R_wingeff=3.004;    %��Ч��򳤶�(mm)  
C_avereff=0.8854;  % ��λ:mm
xr=0.3289;                % x-root offset  \mm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ���-ƫ��-����ϵԭ��ľ���
R_proximal=xr;                                                    % xr=3.19;     %RJ Wood��Ƶĳ��\mm
R_distal=R_wingeff+xr;                                        % yr=0.73;    %RJ Wood��Ƶĳ��\mm
x=linspace(R_proximal,R_distal,200);
x_mod_Root=0.636;                                            % mm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ��Գ���ͳ�����ߵ�Ťת�ᶨ����ǰԵ     
% C_maxy=C_lead_ymax-C_maxy*C_max_LtoT; %ת����ƫ�������е��ƫ��������,��Ťת������ƫ�� -C_maxy֮��  %XXX
yr_lead=-0.08249*x.^6+0.9167*x.^5-4.04*x.^4+8.872*x.^3-10.06*x.^2+5.674*x-0.413-x_mod_Root;  
yr_trail=-0.0333*x.^6+0.504*x.^5-2.795*x.^4+7.258*x.^3-8.769*x.^2+3.739*x+0.1282-x_mod_Root;
%% ����1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_coeff_le=polyfit(x,yr_lead,6);
P_coeff_tr=polyfit(x,yr_trail,6);
%% ����2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% ���ö���ʽ������ϻ�á���������ǰ��Ե�ֲ��������ҳ��ֲ�����
% % (a) ������ǰԵ�ֲ�����
% r_nd=(x-xr)/R_wingeff;
% yr_leadnd0=yr_lead/C_avereff;
% P_coeff_lead=polyfit(r_nd,yr_leadnd0,6);
% % (b)  �����ٺ�Ե�ֲ�����
% yr_trailnd0=yr_trail/C_avereff;
% P_coeff_trail=polyfit(r_nd,yr_trailnd0,6);
% % % (c)  �������ҳ��ֲ�����
% % Cr_nd=yr_leadnd0-yr_trailnd0;
% % P_coeff_Cr=polyfit(r_nd,Cr_nd,6);  % ����ʽϵ��  % Cr_nd2=polyval(Coeff,r_nd1);
% syms  r_nd   % �������ҳ��ֲ�Ϊ6�׶���ʽ����ת������������ָ��
% yr_leadnd=vpa(poly2sym(P_coeff_lead,r_nd),6); 
% yr_trailnd=vpa(poly2sym(P_coeff_trail,r_nd),6);
% % Cr_nd=vpa(poly2sym(P_coeff_Cr,r_nd),6);  
% % yr_leadnd =-68.4639*r_nd^6+208.297*r_nd^5-245.231*r_nd^4+137.467*r_nd^3-36.8594*r_nd^2+4.79235*r_nd-0.156071;
% % yr_trailnd =-27.6379*r_nd^6+121.093*r_nd^5-185.804*r_nd^4+125.603*r_nd^3-33.1053*r_nd^2-0.145527*r_nd-0.156013;
% % Cr_nd =-40.826*r_nd^6+87.204*r_nd^5-59.4267*r_nd^4+11.8645*r_nd^3-3.75408*r_nd^2+4.93788*r_nd-0.0000578215;
% % C_nd=Cr_nd;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ��ʼ��Ӭ���
% % x_mod_Root=0.636;    
% % y_r_lead=-0.08249*r.^6+0.9167*r.^5-4.04*r.^4+8.872*r.^3-10.06*r.^2+5.674*r-0.413-x_mod_Root;  
% % y_r_trail=-0.0333*r.^6+0.504*r.^5-2.795*r.^4+7.258*r.^3-8.769*r.^2+3.739*r+0.1282-x_mod_Root;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%
% f_x_lead=C_avereff*yr_leadnd;  % ���ź���
% f_x_trail =C_avereff*yr_trailnd;  % ���ź���
% %%%%%%%%%%%%%%%%%%%%%%%%%
% f_x_lead1=inline(vectorize(f_x_lead),'r_nd');  %��ֵ����
% f_x_trail1=inline(vectorize(f_x_trail),'r_nd');   %��ֵ����
% r_nd1=linspace(0,1,200);     % x=linspace(R_proximal,R_distal,200);  % x=r_nd*R_wingeff+xr;
% f_x_lead2=f_x_lead1(r_nd1); % ��ֵ
% f_x_trail2=f_x_trail1(r_nd1);   % ��ֵ
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rx=r_nd1*R_wing+xr0;  % rx=linspace(0,1,200)*R_wingeff+xr; 
% r_x_min=min(rx);  
% r_x_max=max(rx);   
% cftool
% warning('off')
% P_coeff_le=polyfit(rx,f_x_lead2,6);     % ����ʽϵ��
% P_coeff_tr=polyfit(rx,f_x_trail2,6);      % ����ʽϵ��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_x_min=min(x);  
r_x_max=max(x);   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms  rx
y_r_lead=vpa(poly2sym(P_coeff_le,rx),6);    % ���ź���  % y_r_lead=poly2sym(P_coeff_le,rx);    % ���ź���
y_r_trail=vpa(poly2sym(P_coeff_tr,rx),6);     % ���ź���  % y_r_trail=poly2sym(P_coeff_tr,rx);     % ���ź���
% %%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%
fun_le=expand((1+(diff(y_r_lead,rx,1))^2)^(1/2));               % fun_le=simplify(expand((1+(diff(y_r_lead,rx,1))^2)^(1/2)))
warning('off')
C_perimeter_le=double((int(fun_le,rx,r_x_min,r_x_max)));     % C_perimeter_le =3.2791; % change double to vpa
fun_tr=expand((1+(diff(y_r_trail,rx,1))^2)^(1/2));                 % fun_tr=simplify(expand((1+(diff(y_r_trail,rx,1))^2)^(1/2)))
C_perimeter_tr=double((int(fun_tr,rx,r_x_min,r_x_max)));      % C_perimeter_le =3.6067; % change double to vpa
E_C_peri=((C_perimeter_le+C_perimeter_tr)/2)/R_wingeff;    % E_C_peri =1.1461;
% ��̬���ű������
R_ratio=R_wing/R_wingeff;
C_ratio=C_aver/C_avereff;
C_periEdge=R_ratio*C_ratio*E_C_peri;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%