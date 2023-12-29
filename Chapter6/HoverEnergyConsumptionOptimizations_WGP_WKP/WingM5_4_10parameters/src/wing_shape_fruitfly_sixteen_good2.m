function wing_para_output=wing_shape_fruitfly_sixteen_good2(R_wing,C_aver,xr0,C_maxyaxis)
%%%%%%%%%%%%%%%%%%%%%%%%%
% ����ʱ�䡪��2014��6��19��,0:16:04
% �޸�ʱ�䡪��2014��6��14��,23:35:03
% �޸�ʱ�䡪��2014��12��21��,23:36
% �޸�ʱ�䡪��2015��01��20��,11:10
% �޸�ʱ�䡪��2015��05��20��,11:10
%%%%%%%%%%%%%%%%%%%%%%%%%
% (1)��ȷ���������ҳ����㡪��������ȷ;
% (2)��ȷ����Ч���۵�������λ��Z_nd���㡪��������ȷ;
% (3)��Ťת������ƫ��C_maxy֮�󡪡������ǰԵ�����С��Ե������������Ӧ��Ƭ������������Ťת���λ��
% ���������ǰԵ�����0.25*[���ǰԵ��y�������С��Ե��y����(�����������)=C_max]λ��ʱ
% (4) ��ƽ������ϵ��ԭ����Ťת��Ľ��˵��غ�;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all;clc;
% x=[3.004,0.8854,0.3289,0.25]; 
% % x=[3.004,0.8854,0.3289,0.356737];
% % x=[2.5604,0.9805,0.8083,0.0946];
% % x=[2.5223,1.0679,0.7204,0.2377];
% R_wing=x(1);
% C_aver=x(2);
% xr0=x(3);
% C_maxyaxis=x(4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ����ò������������Hedrick����������
%��ϵõ� R_wingeff=3.3328-0.3289=3.004;  Hedrick����������: R_wingeff=3.007; Science������: R_wingeff=2.99; 
R_wingeff=3.004;    %��Ч��򳤶�(mm)  
% Hedrick����������: C_avereff=0.884;     % Science������C_avereff=0.9468mm;  
C_avereff=0.8854;  % ��λ:mm---�����ġ�C_avereff=C_aver =0.8854;����ǰ��Եʵ��������ߺ����������־�ֵ����
% ���Կ�����֪չ�ұȣ�������������ŵĳ���ò
% AR=R_wingeff/C_avereff;     %aspect ratio: AR=R^2/A_w;  �������Ϊ: AR=3.40158;  % Science������: AR=3.1519;
% C_avereff=R_wingeff/AR;     %mean chord length: C_aver=A_w/R=R^2/AR/R=R/AR; % C_aver=0.884mm;
% A_w=R_wingeff^2/AR;        %Area of wing: �������Ϊ: A_w=2.66mm^2;   %RJ Wood��Ƶĳ��: A_w=2.84345 mm^2 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (1) ��һ�ַ�ʽ���������������������������ò���ӡ�����Ҫ����Ĳ�������
xr=0.3289;                % x-root offset  \mm
% xr_nd=xr/R_wingeff;      % x-root offset  ������չ��ƫ�þ���
% yr=0;                          % y-root offset  \mm
% yr_nd=yr/C_avereff;   % y-root offset  ����������ƫ�þ���;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ��֪r2_nd,���r1_nd, r3_nd, ��r3_nd3
% r2_nd=0.5801;      %��֪�����ٶ��������r2_nd;    %��Hedrick����������non_dimensional_2ndMoment: 0.5801(��׼) 
% r1_nd=1.106*r2_nd^1.366;  %2013-ICRA-Deng XY   % ���: r1_nd =0.5257
% % r3_nd=0.9*r1_nd^0.581; %1984-JEB-Ellingdon;  % ���: r3_nd =0.6194 % ��Hedrick���������� non_dimensional_3rdMoment: 0.6184(��׼)     
% % r3_nd3=r3_nd^3;                   % r3_nd3=r3_nd^3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R_wingeff1=3.007;            % Hedrick����������
% xr_nd1=xr/R_wingeff1;      % x-root offset  ������չ��ƫ�þ���
% F_nd1=r2_nd^2+2*xr_nd1*r1_nd+xr_nd1^2;   %����������������F_nd1, �������: F_nd1 =0.4635
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (2) �ڶ��ַ�ʽ�����������������������wing_shape������ǰ��Ե��Ϻ���
% ��ȡǰ��Ե�Ľ������� 
% yr_lead= -2.607e-015*x.^11.78+0.8139-0.806122+0.73;    % ǰԵ����Ϻ��������˷�����
% yr_lead= -2.607e-015*x.^11.78+0.737778;                          % ������һ��ָ����
% yr_trail=-0.0017*x.^3+0.1073*x.^2-1.3182*x+0.7783;        % ��Ե����Ϻ����������׶���ʽ
% p(x)=-0.0017*x.^3+0.1073*x.^2-1.3182*x+0.7783+2.607e-015*x.^11.78-0.737778;  %����Maple���
% fsolve(p2(x) = 0, x, x = 10 .. 18);     %����Maple���   % x=16.13067567;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ���-ƫ��-����ϵԭ��ľ���
R_proximal=xr;                                                    % xr=3.19;     %RJ Wood��Ƶĳ��\mm
R_distal=R_wingeff+xr;                                        % yr=0.73;    %RJ Wood��Ƶĳ��\mm
x=linspace(R_proximal,R_distal,200);
x_mod_Root=0.636;                                            % mm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C_maxylb=0.464385778290230;
% C_maxy25=0.138924474377504;        % ��Գ���wing_model_88_yaxis��: ��122��; C_maxy =0.1389; 
% C_maxyub=-0.186536829535222; 
% C_maxy=C_maxylb;
% C_maxy=C_maxy25;
% C_maxy=C_maxyub;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ԭʼ��Ӭ��� R_wing=3.004;  C_aver=0.8854;   
C_lead_ymax=0.4644;   % C_trail_ymin =-0.8374;  
C_max_LtoT= 1.3018;    % @C_maxy=0;
Ratio_leadmax=C_lead_ymax/C_max_LtoT; %Ratio_leadmax=0.4644/1.3018=0.356737; %��Գ���ͳ�����ߵ�Ťת�ᶨ����ǰԵ
% x_start=[3.004,0.8854,0.3289,0.356737]; % ��ʼֵ. % δ�Ż��Ĺ�Ӭ�����ò������������ͳ�����ߵ�Ťת�ᶨ����ǰԵ
% C_lead_ymax =0.3255;  C_trail_ymin =-0.9764;  C_max_LtoT =1.3018;
% Ratio_leadmax=C_lead_ymax/C_max_LtoT; %Ratio_leadmax=0.3255/1.3018=0.25;%���Ťת��λ�����ǰԵ������º�Ե��������0.25��ʱ
% x_start=[3.004,0.8854,0.3289,0.25];     % ��ʼֵ. % δ�Ż��Ĺ�Ӭ�����ò��������Ťת��λ�����ǰԵ������º�Ե��������0.25��ʱ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta_pitchaxis=Ratio_leadmax-C_maxyaxis;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ��Գ���ͳ�����ߵ�Ťת�ᶨ����ǰԵ     
% C_maxy=C_lead_ymax-C_maxy*C_max_LtoT; %ת����ƫ�������е��ƫ��������,��Ťת������ƫ�� -C_maxy֮��  %XXX
yr_lead=-0.08249*x.^6+0.9167*x.^5-4.04*x.^4+8.872*x.^3-10.06*x.^2+5.674*x-0.413-x_mod_Root;  
yr_trail=-0.0333*x.^6+0.504*x.^5-2.795*x.^4+7.258*x.^3-8.769*x.^2+3.739*x+0.1282-x_mod_Root;
% C_rx=yr_lead-yr_trail;      % ��ȷ�������ٻ�ʵ���ҳ��ֲ�
% figure(1)  % ����ò��������ǰԵ��Ϻ���
% plot(x,yr_lead,'r-',x,yr_trail,'b-')
% xlabel('չ��r (mm)')
% ylabel('ǰԵ�ͺ�Ե�ҳ� (mm)')
% legend('ǰԵ�������','��Ե�������')
% title('R_wingeff=3.004mm: �Ἰ����ò��ǰԵ�ͺ�Ե�������')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x_bem=xr+0.7*R_wingeff;
% yr_leadbem=-0.08249*x_bem.^6+0.9167*x_bem.^5-4.04*x_bem.^4+8.872*x_bem.^3-10.06*x_bem.^2+5.674*x_bem-0.413-x_mod_Root-C_maxy; % yr_leadbem =0.4598;  
% yr_trailbem=-0.0333*x_bem.^6+0.504*x_bem.^5-2.795*x_bem.^4+7.258*x_bem.^3-8.769*x_bem.^2+3.739*x_bem+0.1282-x_mod_Root-C_maxy; % yr_trailbem =-0.6658; 
% C_bem=yr_leadbem-yr_trailbem   % C_bem =1.1257mm;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ����ǰԵ��Ϻ�����⡪��ʵ��ƽ���ҳ�=���/R_wingeff
% wing_aera=trapz(x,C_rx);             %���: wing_aera =2.6597; % mm^2
% C_aver=wing_aera/R_wingeff;   % ������ٻ�ƽ���ҳ�: C_avereff=C_aver =0.8854; % mm
% % ���õڶ��ֻ��ַ�����⣬ò�Ʋ��ԡ�������XXXX
% yr_lead1=inline('-0.08249*x.^6+0.9167*x.^5-4.04*x.^4+8.872*x.^3-10.06*x.^2+5.674*x-0.413-0.636-0.1389','x'); 
% yr_trail1=inline('-0.0333*x.^6+0.504*x.^5-2.795*x.^4+7.258*x.^3-8.769*x.^2+3.739*x+0.1282-0.636-0.1389','x'); 
% wing_aera1=abs(quadl(yr_lead1,R_proximal,R_distal))+abs(quadl(yr_trail1,R_proximal,R_distal)); %���: wing_aera1 =3.1119; % mm^2
% C_aver1=wing_aera1/R_wingeff;    % ���: C_aver1 =1.0359;

%%%%%%%%%%%%%%%%%%%%%%%%%
%% ����(1)
%% ��⡪���������ҳ��ֲ�����ʽ���������������󡪡�XXX
% r_nd1=linspace(0,1,200);    % r=r_nd*R_wingeff;
% Cr_nd1=C_rx/C_aver;          % �������ҳ��ֲ�  
% % cftool
% % ���r_nd1��Cr_nd1���ú���cftool������Ϲ���������������ҳ��ֲ�������ϣ��������K�׶���ʽ
% % ������
% % Linear model Poly6:  f(x) = p1*x^6 + p2*x^5 + p3*x^4 + p4*x^3 + p5*x^2 + p6*x + p7
% % Coefficients (with 95% confidence bounds):
% % p1 =-40.83  (-40.83, -40.83);p2 =87.21  (87.21, 87.21);p3 = -59.43  (-59.43, -59.43);p4 =11.86  (11.86, 11.86);
% % p5 =-3.754  (-3.754, -3.754);p6 =4.938  (4.938, 4.938);p7 = -5.782e-005  (-5.782e-005, -5.782e-005);
% % Goodness of fit: SSE: 4.697e-026; R-square: 1; Adjusted R-square: 1;  RMSE: 1.56e-014;
% % ����������ҳ��ֲ�����6�׶���ʽ
% Cr_nd2=-40.83*r_nd1.^6 + 87.21*r_nd1.^5-59.43*r_nd1.^4+11.86*r_nd1.^3-3.754*r_nd1.^2+4.938*r_nd1-5.782e-005;
% % figure(2)
% % subplot(211)
% % plot(x,C_rx,'r-')
% % xlabel('չ��r (mm)')
% % ylabel('���ٻ�ʵ���ҳ�(mm) ')
% % legend('���ٻ�ʵ���ҳ�')
% % title('R_wingeff=3.004mm: �Ἰ����ò�����ٻ�ʵ���ҳ��ֲ�')
% % subplot(212)
% % plot(r_nd1,Cr_nd2,'b-')
% % xlabel('չ��\itr_{\rmnd} (I)')
% % ylabel('�������ҳ�(I)')
% % legend('�������ҳ�')
% % title('R_wingeff_nd=1: �Ἰ����ò���������ҳ��ֲ�')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% ����������������(nondimention_aerodynamic_component)����⡪��F_nd
% Coeff=polyfit(r_nd1,Cr_nd1,6);  % ����ʽϵ��  % Cr_nd2=polyval(Coeff,r_nd1);
% syms r_nd
% Cr_nd2=vpa(poly2sym(Coeff,r_nd),6);   % �������ҳ��ֲ�Ϊ6�׶���ʽ��������Cr_nd3=Cr_nd2;ֻ���Ա�����r_nd1�����r_nd
% % Cr_nd2 =-40.827*r_nd^6+87.2061*r_nd^5-59.4281*r_nd^4+11.8648*r_nd^3-3.75417*r_nd^2+4.938*r_nd-0.0000578229;
% C_nd=Cr_nd2;
% S_nd=double(int(C_nd,r_nd,0,1))                 %�����ٳ���� % S_nd =1.0000;
%%%%%%%%%%%%%%%%%%%%%%%%%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

%% ����(2)���������������ö���ʽ������ϻ�á���������ǰ��Ե�ֲ��������ҳ��ֲ�����
% (a) ������ǰԵ�ֲ�����
r_nd=(x-xr)/R_wingeff;
yr_leadnd0=yr_lead/C_avereff;
P_coeff_lead=polyfit(r_nd,yr_leadnd0,6);
% (b)  �����ٺ�Ե�ֲ�����
yr_trailnd0=yr_trail/C_avereff;
P_coeff_trail=polyfit(r_nd,yr_trailnd0,6);
% (c)  �������ҳ��ֲ�����
Cr_nd=yr_leadnd0-yr_trailnd0;
P_coeff_Cr=polyfit(r_nd,Cr_nd,6);  % ����ʽϵ��  % Cr_nd2=polyval(Coeff,r_nd1);
syms r_nd   % �������ҳ��ֲ�Ϊ6�׶���ʽ����ת������������ָ��
yr_leadnd=vpa(poly2sym(P_coeff_lead,r_nd),6); 
yr_trailnd=vpa(poly2sym(P_coeff_trail,r_nd),6);
Cr_nd=vpa(poly2sym(P_coeff_Cr,r_nd),6);  
% yr_leadnd =-68.4639*r_nd^6+208.297*r_nd^5-245.231*r_nd^4+137.467*r_nd^3-36.8594*r_nd^2+4.79235*r_nd-0.156071;
% yr_trailnd =-27.6379*r_nd^6+121.093*r_nd^5-185.804*r_nd^4+125.603*r_nd^3-33.1053*r_nd^2-0.145527*r_nd-0.156013;
% Cr_nd =-40.826*r_nd^6+87.204*r_nd^5-59.4267*r_nd^4+11.8645*r_nd^3-3.75408*r_nd^2+4.93788*r_nd-0.0000578215;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cftool % ���ַ�ʽ���ᳫ
% % (1) ��϶���ʽ����yr_leadnd(r_nd)
% yr_leadnd=-68.46*r_nd.^6+208.3*r_nd.^5-245.2*r_nd.^4+137.5 *r_nd.^3-36.86*r_nd.^2+4.792*r_nd+0.0008352;
% % Linear model Poly6: f(x)=p1*x^6+p2*x^5+p3*x^4+p4*x^3+p5*x^2+p6*x + p7
% % Coefficients (with 95% confidence bounds):
% %        p1= -68.46(-68.46, -68.46);p2=208.3(208.3, 208.3); p3=-245.2(-245.2, -245.2); p4 =137.5(137.5, 137.5);
% %        p5= -36.86(-36.86, -36.86); p6=4.792(4.792, 4.792); p7= 0.0008352(0.0008352, 0.0008352);
% % Goodness of fit:  SSE: 9.451e-026; R-square: 1; Adjusted R-square: 1; RMSE: 2.213e-014
% % (2) ��϶���ʽ����yr_trailnd(r_nd)
% yr_trailnd=-27.64*r_nd.^6+121.1*r_nd.^5-185.8*r_nd.^4+125.6*r_nd.^3-33.11*r_nd.^2-0.1455*r_nd+0.000893;
% % Linear model Poly6:  f(x) = p1*x^6 + p2*x^5 + p3*x^4 + p4*x^3 + p5*x^2 +p6*x + p7
% % Coefficients (with 95% confidence bounds):
% %        p1=-27.64(-27.64, -27.64);p2=121.1(121.1, 121.1);p3=-185.8(-185.8, -185.8);p4 =125.6(125.6, 125.6);
% %        p5=-33.11  (-33.11, -33.11);p6=-0.1455  (-0.1455, -0.1455);p7=0.000893(0.000893, 0.000893);
% % Goodness of fit:  SSE: 3.744e-026;  R-square: 1;  Adjusted R-square: 1; RMSE: 1.393e-014;  
% % (3) ��϶���ʽ����Cr_nd(r_nd)
% Cr_nd=-40.83*r_nd.^6+87.2*r_nd.^5-59.43 *r_nd.^4+11.86*r_nd.^3-3.754*r_nd.^2+4.938*r_nd-5.782e-005; 
% % Cr_nd2=-40.83*r_nd1.^6+87.21*r_nd1.^5-59.43*r_nd1.^4+11.86*r_nd1.^3-3.754*r_nd1.^2+4.938*r_nd1-5.782e-005;%���ĵ����ݴ���ʽ
% % Linear model Poly6:   f(x) = p1*x^6+p2*x^5+p3*x^4+p4*x^3 +p5*x^2+p6*x + p7
% % Coefficients (with 95% confidence bounds):
% %        p1 =-40.83(-40.83, -40.83); p2 =87.2(87.2, 87.2);p3=-59.43(-59.43, -59.43);p4 = 11.86(11.86, 11.86);
% %        p5 =-3.754(-3.754, -3.754);p6 =4.938(4.938, 4.938);p7=-5.782e-005(-5.782e-005, -5.782e-005);
% % Goodness of fit: SSE: 4.898e-026; R-square: 1; Adjusted R-square: 1; RMSE: 1.593e-014;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R_wingeff=R_wing;       % ����������������������������������������������
C_avereff=C_aver;        % ����������������������������������������������  
xr_nd_vari=xr0/R_wing;
xr_nd=xr_nd_vari;
yr_leadnd=yr_leadnd-delta_pitchaxis;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ע�⡪���öγ����мǲ����޸ģ�ǰ��ֻҪ��֤������ȷ���������ҳ��ֲ����ɡ�
%���µĹ�ʽӦʹ�ú���������ٵ��ҳ��ֲ���ʽC_nd
C_nd=Cr_nd;
% R2nd2=double(int(r_nd^2*C_nd,r_nd,0,1)); %��������صĻ�ת�뾶��ƽ��
% R1nd1=double(int(r_nd*C_nd,r_nd,0,1));     %һ������صĻ�ת�뾶
% S_nd=double(int(C_nd,r_nd,0,1));                %�����ٳ���� % S_nd =1.0000;
% disp(['��������صĻ�ת�뾶��ƽ��: r2_2nd=' num2str(R2nd2)  ' ���ٵ�λ��mm^4'])
% disp(['��������صĻ�ת�뾶: r_2nd=' num2str(sqrt(R2nd2))  ' ���ٵ�λ��mm^3'])
% disp(['һ������صĻ�ת�뾶: r_1nd=' num2str(R1nd1)  ' ���ٵ�λ��mm^3'])
% disp(['�����ٳ����Swing_nd=' num2str(S_nd)  ' ���ٵ�λ��mm^2'])
% C_nd=vpa(C_nd,5)  %����vpa�����ű��ʽ�е���ֵ(��ת��Ϊ���������ı�ֵ,������)��ת��Ϊʮ����С����ʾ��
% ����xr_nd=0�����Ա���r_nd��ȡֵ��ΧΪ:r_nd��[0,1], ��ȡ1ʱr_nd=R=11.95 /mm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fx1=(r_nd+xr_nd)^2*C_nd;    % ������������F_nd��ԭʼ��������
% fx2=vpa(fx1,5)
fx3=expand(fx1);
F_ndTrans=double(int(fx3,r_nd,0,1));                    % Result: F_ndTrans =0.46391;
% disp(['������������F_ndTrans=' num2str(F_ndTrans)  ' ���ٵ�λ��mm^4'])
% F_nd2=R2nd2+2*xr_nd*R1nd1+xr_nd^2;    %ʹ����������Ҳ��ȷ; ���:F_nd2 =0.46391;
% disp(['������������F_ndTrans=' num2str(F_nd2)  ' ���ٵ�λ��mm^4'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ���:
% ��������صĻ�ת�뾶��ƽ��: r2_2nd=0.3358 ���ٵ�λ��mm^4
% ��������صĻ�ת�뾶: r_2nd=0.57949 ���ٵ�λ��mm^3  ��������������Ҫ������
% һ������صĻ�ת�뾶: r_1nd=0.53032 ���ٵ�λ��mm^3
% �����ٳ����Swing_nd=1 ���ٵ�λ��mm^2
% ������������F_nd=0.46391 ���ٵ�λ��mm^4    ����������������������Ҫ������
% ������������F_nd=0.46391 ���ٵ�λ��mm^4    ����������������������Ҫ������
%  �Ա������� (1) ��һ�ַ�ʽ����õ�����������������
% F_nd1=r2_nd^2+2*xr_nd1*r1_nd+xr_nd1^2;   %����������������F_nd1, �������: F_nd1 =0.4635
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ��Ӭ�Ŵ��򡪡�fruitgfly_magnified_wing
%%%%%%%%%%%%%%%%%%%%%%%%%
% R_wing=100;             % mm %
% C_aver=33.7;            % mm %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ������Ϻõ�ǰ��Ե
% syms r_nd    % r_nd=(x-xr)/R_wingeff;  % r_nd��(0,1)
% yr_leadnd=-68.4639*r_nd^6+208.297*r_nd^5-245.231*r_nd^4+137.467*r_nd^3-36.8594*r_nd^2+4.79235*r_nd+0.000835197; % @���ǰԵ��������;
% yr_trailnd=-27.6379*r_nd^6+121.093*r_nd^5-185.804*r_nd^4+125.603*r_nd^3-33.1053*r_nd^2-0.145527*r_nd+0.000893019; % @���ǰԵ��������;
f_x_lead=C_aver*yr_leadnd;  % ��Ӭ�������ݱ��Ŵ���
f_x_trail =C_aver*yr_trailnd;  % ��Ӭ�������ݱ��Ŵ���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_x_lead1=inline(vectorize(f_x_lead),'r_nd');        %��ֵ����
f_x_trail1=inline(vectorize(f_x_trail),'r_nd');    %��ֵ����
r_nd1=linspace(0,1,200);     % x=linspace(R_proximal,R_distal,200);  % x=r_nd*R_wingeff+xr;
f_x_lead2=f_x_lead1(r_nd1); 
f_x_trail2=f_x_trail1(r_nd1);
C_lead_ymax=max(f_x_lead2);  % ���: C_lead_ymax=0.4644; k_leadmax=644;%ԭʼ��Ӭ���@���ǰԵ�㵽���������ߵľ���;@C_maxy=0;
C_trail_ymin=min(f_x_trail2);     % ���: C_trail_ymin=-0.8375; k_trailmin=409;%ԭʼ��Ӭ��� @��С��Ե�㵽���������ߵľ���;@C_maxy=0;
C_max_LtoT=C_lead_ymax-C_trail_ymin;  % C_max_LtoT =1.3018;% ԭʼ��Ӭ���@���ǰԵ�㵽��С��Ե��ľ���@C_maxy=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (һ) �ֱ�ƽ���ҳ��ͳ��Ťת������λ�� 
% C_avereff=0.8854;  % mm
% % (1) ���ǰԵ���Ӧ��Ƭ��������ҳ�
% [C_lead_ymax,k_leadmax]=max(f_x_lead);   % ���: C_lead_ymax=0.4644; k_leadmax=644;
% yr_leadnd_max=C_lead_ymax/C_avereff;    % yr_leadnd_max =0.5245;  % ���ǰԵ��������ٳ���
% C_max_x=y_lead(1,k_leadmax);                   % C_max_x= 2.2623;
% C_trail_y=f_x_trail(1,k_leadmax);                 % C_trail_y=-0.7023;
% C_leadmax=C_lead_ymax-C_trail_y;            % C_leadmax=1.1667;
% % (2) ��С��Ե���Ӧ��Ƭ��������ҳ�
% [C_trail_ymin,k_trailmin]=min(f_x_trail);    % ���: C_trail_ymin =-0.8375;  k_trailmin =409;
% yr_trailnd_min=C_trail_ymin/C_avereff;     % yr_trailnd_min =-0.9459;  % ���ǰԵ��������ٳ���
% C_min_x=y_lead(1,k_trailmin);                   % C_min_x =1.5557;
% C_lead_y=f_x_lead(1,k_trailmin);                % C_lead_y =0.3549;
% C_trailmin=C_lead_y-C_trail_ymin;             % C_trailmin =1.1924;
% C_max_LtoT=C_lead_ymax-C_trail_ymin;       % C_max_LtoT =1.3018;
% x_0lb=0*C_max_LtoT;                                    % x_0lb=0;
% x_025=0.25*C_max_LtoT;                              % x_025=0.3255;
% x_0ub=0.5*C_max_LtoT;                                % x_0ub =0.6509;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ��һ���֡���ƽ�����������������ز���
% (1) ƽ���������������������Ƴ�ƽ���µķ���
Rou=1.225*10^(-3);           %��λ��Kg/m^3=10^6/(10^3)^3=10^(-3)mg/mm^3
% Rou=1.225;                         %��λ��Kg/m^3   
% Coeff_liftdragF_N=6.8201e-012������λ��:mg*mm^-3*mm^4= 10^(-9) kg*m
Coeff_liftdragF_N=(1/2)*Rou*C_avereff*R_wingeff^3*F_ndTrans;  % mg*mm
% (2) ƽ�������������ز��������Ƴ�ƽ���µ�չ����
% mg*mm^-3*mm^5=Kg.m^2:[10^(-12)]=mg.mm^2: [10^12*10^(-12)]= mg.mm^2
M_xaercoeff=(1/2)*Rou*C_avereff^2*R_wingeff^3*F_ndTrans;   % ��ת����������ϵ��  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (3) ƽ�������������ز�������I1y�����Ƴ�ƽ������ϵ�µ�������
fx4=(r_nd+xr_nd)^3*C_nd;                                 % ����C_nd��ǰ��Ե����ʽ����֮�����
fr_nd5=expand(fx4);
I1=double(int(fr_nd5,0,1));                                 % ������,���ٻ���λ��mm^5; 
I1y=(1/2)*Rou*C_avereff*R_wingeff^4*I1;         % ��λ�� mg.mm^2  % I1y=0.0162
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �ڶ����֡��� ��������������ء��������Ч���۵�������λ��Y_nd
% % ���-ƫ��-����ϵԭ��ľ���
% R_proximal=xr;                                                    % xr=3.19;     %RJ Wood��Ƶĳ��\mm
% R_distal=R_wingeff+xr;                                        % yr=0.73;    %RJ Wood��Ƶĳ��\mm
% x=linspace(R_proximal,R_distal,200);
% yr_lead=-0.08249*x.^6+0.9167*x.^5-4.04*x.^4+8.872*x.^3-10.06*x.^2+5.674*x-0.413-x_mod_Root-C_maxy; 
% yr_trail=-0.0333*x.^6+0.504*x.^5-2.795*x.^4+7.258*x.^3-8.769*x.^2+3.739*x+0.1282-x_mod_Root-C_maxy;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ����(1)������ǰ��Ե�����������ٻ�yr_leadnd����yr_trailnd�������
% yr_leadnd =-68.4639*r_nd^6+208.297*r_nd^5-245.231*r_nd^4+137.467*r_nd^3-36.8594*r_nd^2+4.79235*r_nd+0.000835197;
% yr_trailnd =-27.6379*r_nd^6+121.093*r_nd^5-185.804*r_nd^4+125.603*r_nd^3-33.1053*r_nd^2-0.145527*r_nd+0.000893019;
%%%%%%%%%%%%%%%%%%%
% (3)��Ťת������ƫ��C_maxy֮�󡪡������ǰԵ�����С��Ե������������Ӧ��Ƭ������������Ťת���λ��
% ���������ǰԵ�����0.25*[���ǰԵ��y�������С��Ե��y����(�����������)=C_max]λ��ʱ
% C_maxy=0.138924474377504;
% C_maxynd=C_maxy/C_avereff;
% yr_leadnd=yr_leadnd-C_maxynd;
% yr_trailnd=yr_trailnd-C_maxynd;
%%%%%%%%%%%%%%%%%%%
% yr_nd1=expand((yr_leadnd^4+yr_trailnd^4)/4);     %��������2
% Z_rnd=double(int(yr_nd1,r_nd,0,1));                         % ���: Y_nd=0.08802(old);  %  mm  % Z_rnd =0.1626;
% ����(2)������ǰԵ�����������ٻ�yr_leadnd�������ٻ��ҳ��ֲ�C_nd���
y0=yr_leadnd-C_nd; % y0=yr_trailnd;  % ��ʱ���Z_nd=0.08802; ����������(1)
y1=yr_leadnd;
yr_nd2=(y1^4+y0^4)/4;   % ע�������Ťת��ֱ��ͨ�����ͳ�������ߣ�ǰԵ�����ĳ��չ������Ψһ
Z_rnd=double(int(yr_nd2,r_nd,0,1));                          % ���: Z_nd=0.08802(old);   %  mm  % Z_rnd =0.1627;
% disp(['��Ч���۵�������λ��Z_nd=' num2str(Z_rnd)  ' ���ٵ�λ��mm'])  % �������Maple���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) ת�������������ز��������Ƴ�ƽ���µ�չ����
% ����Ϊת��������������ϵ��% mg*mm^-3*mm^5=Kg.m^2:[10^(-12)]= mg.mm^2
M_xrdcoeff=(1/2)*Rou*C_avereff^4*R_wingeff*Z_rnd;  %M_xrdcoeff=0.0001839;
% (2) ת�������������ز��������Ƴ�ƽ���µ�������
% ����Ϊת��������������ϵ��% mg*mm^-3*mm^5=Kg.m^2:[10^(-12)]= mg.mm^2
fx16=(r_nd+xr_nd)*C_nd^3;                                 % ����C_nd��ǰ��Ե����ʽ����֮�����
fr_nd17=expand(fx16);
I8=double(int(fr_nd17,0,1));                                  % ������,���ٻ���λ��mm^5;
I8z=I8;           % ��λ�� mg.mm^2
X_rnd=I8z;
M_zrdcoeff=(1/6)*Rou*C_avereff^3*R_wingeff^2*X_rnd; % M_zrdcoeff=0.001169
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �������֡���ת�����������������ز���
fx6=(r_nd+xr_nd)*C_nd^2;                                  % ����C_nd��ǰ��Ե����ʽ����֮�����   % ������������F_ndRot��ԭʼ��������
fr_nd7=expand(fx6);
% F_ndRot1=double(int(fr_nd7,r_nd,0,1))           % �������Ҳ��Ŷ
F_ndRot=double(int(fr_nd7,0,1));                        % ������,���ٻ���λ��mm^4; F_ndRot=0.74851
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) ת��������������������F_yrotcoeff�����Ƴ�ƽ���µķ���   % ��λ��mg*mm^-3*mm^4= kg.m 10^(-9)
F_yrotcoeff=(1/2)*Rou*C_avereff^2*R_wingeff^2*F_ndRot;  % ��λ��mg*mm
% (2) ת�������������ز�������M_xRotcoeff�����Ƴ�ƽ���µ�չ����
% mg*mm^-3*mm^5=Kg.m^2:[10^(-12)]=[10^12*10^(-12)] mg.mm^2
M_xRotcoeff=(1/2)*Rou*C_avereff^3*R_wingeff^2*F_ndRot;   % ת��������������ϵ�������Ƴ�ƽ���µ�Ťת��
% (3) ת�������������ز�������I6y�����Ƴ�ƽ���µ�������
fx8=(r_nd+xr_nd)^2*C_nd^2;                            % ����C_nd��ǰ��Ե����ʽ����֮�����
fr_nd9=expand(fx8);
I2=double(int(fr_nd9,0,1));                                 % ������,���ٻ���λ��mm^5;
I2y=(1/2)*Rou*C_avereff^2*R_wingeff^3*I2;     % ��λ�� mg.mm^2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ���Ĳ��֡���������ЧӦ��ϵ������������ת���������������������������ز���
% C_maxy=0.138924474377504;
% C_maxynd=C_maxy/C_avereff;
% yr_hnd=C_nd/2-yr_leadnd+C_maxynd; %ת����ƫ�������е��ƫ��������,��Ťת������ƫ��C_maxy֮��, ���������������۽�ȫ��Ϊ��ֵ
yr_hnd=C_nd/2-yr_leadnd; %ת����ƫ�������е��ƫ��������,��Ťת������ƫ��C_maxy֮��, ���������������۽�ȫ��Ϊ��ֵ
% ����������ϵ��
% a=C_nd/2;
% lambda_z=pi*Rou*a^2;
% lambda_zw=-pi*Rou*a^2*yr_hnd;
% lambda_w=pi*Rou*a^2*yr_hnd^2+pi*Rou*a^4/8;
% dOmega0=-r*(domega_y-omega_x*omega_z);           %������ٶ�
% Z0=-lambda_z*dOmega0-lambda_zw*domega_x;      %������������
% Y0=0;                                                                            %������������
% M0=-lambda_zw*dOmega0-lambda_w*domega_x;    %��������������
%%%%%%%%%%%%%%%%%%%%%%%%
% (1) �������������ز�������������ת������I_xzam�����Ƴ�ƽ���µ�Ťת��
fx10=(r_nd+xr_nd)*C_nd^2*yr_hnd;
fr_nd11=expand(fx10);
I_xzamnd=int(fr_nd11,0,1);
I3=double(I_xzamnd);                                         % ������,���ٻ���λ��mm^5;  I3 =0.2362; 
I_xzam=pi*Rou*C_avereff^3*R_wingeff^2*I3/4;  % ��λ�� mg.mm^2;   % I_xzam =0.001424
% (2) �������������ز�������������ת������I_xxam�����Ƴ�ƽ���µ�Ťת��
fx12=C_nd^2*(yr_hnd^2+C_nd^2/32);
fr_nd13=expand(fx12);
I_xxamnd=int(fr_nd13,0,1);
I4=double(I_xxamnd);                                         % ������,���ٻ���λ��mm^5;  I4 =0.1903;
I_xxam=pi*Rou*C_avereff^4*R_wingeff*I4/4;      % ��λ�� mg.mm^2;   % I_xxam =0.0003380;
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (3) ��������������������I5z
I5z=pi*Rou*C_avereff^2*R_wingeff^2*F_ndRot/4;    % ��λ�� mg.mm  I5z =0.005094;
% (4) ��������������������I6z
fx14=C_nd^2*yr_hnd;
fr_nd15=expand(fx14);
I6=double(int(fr_nd15,0,1));                               % ������,���ٻ���λ��mm^4; 
I6z=pi*Rou*C_avereff^3*R_wingeff*I6/4;           % ��λ�� mg.mm  % I6z =0.000771644; 
% (5) �������������ز�������������ת������I7y�����Ƴ�ƽ���µ�������
 I7y=pi*Rou*C_avereff^2*R_wingeff^3*I2/4;     % ��λ�� mg.mm^2
%% (a) ����ѹ�ġ���չ������������/����������% ���ٻ���Ҫ����*C_avereff  or c(r)@r=(R+x_rnd)*r_xcopnd_tr;...
% % c_zcopnd_tr=I1/F_ndTrans;  % XXX       % Y_rcpnd_transaver=0.1715; % *C_avereff  or C(r) @r=R*r_xcopnd_tr;...
% % c_zcopnd_rot=I2/F_ndRot;    % XXX       % Y_rcpnd_rotaver=0.1677;   % *C_avereff  or C(r) @r=R*r_xcopnd_rot;...
% % ���������������ɵ�ƽ��ѹ�ġ�����Ӧ����ÿ��Ƭ�������ҵ�
% c_zcopnd_addaver=-0.3058;        % c_zcopnd_add=M_xam./F_yadd1;%c_zcopnd_addaver=mean(c_zcopnd_add); 
% % c_zcopnd_addtr=I3/F_ndRot;   % c_zcopnd_addtr=0.1587;    % *C_avereff or C(r) @r=R*r_xcopnd_addtr;...
% % c_zcopnd_addrot=I4/I6;           % c_zcopnd_addrot=0.4816;  % *C_avereff or C(r) @r=R*r_xcopnd_addrot;...
%% (b) չ��ѹ�ġ�����������������/���������� % ���ٻ���Ҫ����*R_wingeff ������  *(R_wingeff+xr)  
% % format long
% r_xcopnd_tr=I1/F_ndTrans;   % r_xcopnd_tr =0.788691874094779; % r_xcopnd_tr= 0.7887;          % *R_wingeff
% r_xcopnd_rot=I2/F_ndRot;   % r_xcopnd_rot =0.712638285096407;  % r_xcopnd_rot =0.7126;        % *R_wingeff  
% % չ�������������ɵ�ƽ��ѹ�ġ�����Ӧһ������չ��λ��
% r_xcopnd_addaver=0.7320;  % r_xcopnd_add=M_zadd./F_yadd1;  % r_xcopnd_addaver=mean(r_xcopnd_add); 
% % r_xcopnd_addtr=I2/F_ndRot;   % r_xcopnd_addtr=0.7126;     % *R_wingeff
% % r_xcopnd_addrot=I3/I6;           % r_xcopnd_addrot=0.5837;   % *R_wingeff
% % format short
%% (c) չ��ѹ�ġ�����������������/����������
% r_ycopnd_tr1=I1/F_ndTrans;                          % r_ycopnd_tr= 0.7887;
% r_ycopnd_rot1=I2/F_ndRot;                           % r_ycopnd_rot =0.7126;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (1) ����Aero_F3_fruitfly_simpres &  Aero_M_fruitfly1 & Aero_M_fruitfly3������Ϊ��Ƶ�����
% ���ڵ����
% % I_xyam =-0.001245;  % I_xxam =2.5080e-004;   % I4y =-5.3461e-004;       % XXX
% % I_xyam =-0.0045;      % I_xxam =0.0020;             % I4y =-0.0022;                % XXX
% wing_para_output=zeros(1,9);
% wing_para_output=[R_wingeff,C_avereff,F_ndTrans,Z_rnd,I_xzam,I_xxam,F_ndRot,I5z,I6z];
% % I_xyam =0.002;      % I_xxam =6.2892e-004;   % I4y =0.0011;                    % ��ȷ���
%% (2) ����Aero_F3_fruitfly_exp & Aero_M_fruitfly2����ʵ����Ե�����
% 20141122-�޸ĺ�����
% % wing_para_output=zeros(1,17);
wing_para_output=[F_ndTrans,Coeff_liftdragF_N,M_xaercoeff,I1y,Z_rnd, M_xrdcoeff,...
    F_ndRot,F_yrotcoeff,M_xRotcoeff, I2y,...
    I_xzam,I_xxam, I5z,I6z,I7y,M_zrdcoeff,C_max_LtoT];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
