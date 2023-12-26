function wing_para_output=wing_shape_fruitfly()
% ����ʱ�䡪��2014��6��19��,0:16:04
% �޸�ʱ�䡪��2014��6��14��,23:35:03
% clear all;clc;
%% ����ò������������Hedrick����������
%��ϵõ� R_wingeff=3.3328-0.3289=3.004;  Hedrick����������: R_wingeff=3.007; Science������: R_wingeff=2.99; 
R_wingeff=3.004;    %��Ч��򳤶�(mm)  
% Hedrick����������: C_avereff=0.884;     % Science������C_avereff=0.9468mm;  
C_avereff=0.8854;  %�����ġ�C_avereff=C_aver =0.8854;����ǰ��Եʵ��������ߺ����������־�ֵ����
% ���Կ�����֪չ�ұȣ�������������ŵĳ���ò
% AR=R_wingeff/C_avereff;     %aspect ratio: AR=R^2/A_w;  �������Ϊ: AR=3.40158;  % Science������: AR=3.1519;
% C_avereff=R_wingeff/AR;     %mean chord length: C_aver=A_w/R=R^2/AR/R=R/AR; % C_aver=0.884mm;
% A_w=R_wingeff^2/AR;        %Area of wing: �������Ϊ: A_w=2.66mm^2;   %RJ Wood��Ƶĳ��: A_w=2.84345 mm^2 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (1) ��һ�ַ�ʽ���������������������������ò���ӡ�����Ҫ����Ĳ�������
xr=0.3289;                     % x-root offset  \mm
xr_nd=xr/R_wingeff;      % x-root offset  ������չ��ƫ�þ���
% yr=0;                          % y-root offset  \mm
% yr_nd=yr/C_avereff;   % y-root offset  ����������ƫ�þ���;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ��֪r2_nd,���r1_nd, r3_nd, ��r3_nd3
r2_nd=0.5801;      %��֪�����ٶ��������r2_nd;    %��Hedrick����������non_dimensional_2ndMoment: 0.5801(��׼) 
r1_nd=1.106*r2_nd^1.366;  %2013-ICRA-Deng XY   % ���: r1_nd =0.5257
% r3_nd=0.9*r1_nd^0.581; %1984-JEB-Ellingdon;  % ���: r3_nd =0.6194 % ��Hedrick���������� non_dimensional_3rdMoment: 0.6184(��׼)     
% r3_nd3=r3_nd^3;                   % r3_nd3=r3_nd^3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R_wingeff1=3.007;            % Hedrick����������
xr_nd1=xr/R_wingeff1;      % x-root offset  ������չ��ƫ�þ���
F_nd1=r2_nd^2+2*xr_nd1*r1_nd+xr_nd1^2;   %����������������F_nd1, �������: F_nd1 =0.4635
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
yr_lead=-0.08249*x.^6+0.9167*x.^5-4.04*x.^4+8.872*x.^3-10.06*x.^2+5.674*x-0.413;  
yr_trail=-0.0333*x.^6+0.504*x.^5-2.795*x.^4+7.258*x.^3-8.769*x.^2+3.739*x+0.1282;
% figure(1)  % ����ò��������ǰԵ��Ϻ���
% plot(x,yr_lead,'r-',x,yr_trail,'b-')
% xlabel('չ��r (mm)')
% ylabel('ǰԵ�ͺ�Ե�ҳ� (mm)')
% legend('ǰԵ�������','��Ե�������')
% title('R_wingeff=3.004mm: �Ἰ����ò��ǰԵ�ͺ�Ե�������')
C_rx=yr_lead-yr_trail;              % ��ȷ�������ٻ�ʵ���ҳ��ֲ�
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ����ǰԵ��Ϻ�����⡪��ʵ��ƽ���ҳ�=���/R_wingeff
wing_aera=trapz(x,C_rx);             %���: wing_aera =2.6597; % mm^2
C_aver=wing_aera/R_wingeff;   % ������ٻ�ƽ���ҳ�: C_avereff=C_aver =0.8854;
% % ���õڶ��ֻ��ַ�����⣬ò�Ʋ��ԡ�������XXXX
% yr_lead1=inline('-0.08249*x.^6+0.9167*x.^5-4.04*x.^4+8.872*x.^3-10.06*x.^2+5.674*x-0.413','x'); 
% yr_trail1=inline('-0.0333*x.^6+0.504*x.^5-2.795*x.^4+7.258*x.^3-8.769*x.^2+3.739*x+0.1282','x'); 
% wing_aera1=abs(quadl(yr_lead1,R_proximal,R_distal))+abs(quadl(yr_trail1,R_proximal,R_distal)); %���: wing_aera1 =3.1119; % mm^2
% C_aver1=wing_aera1/R_wingeff;    % ���: C_aver1 =1.0359;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ��⡪���������ҳ��ֲ�����ʽ����
r_nd1=linspace(0,1,200);    % r=r_nd*R_wingeff;
Cr_nd1=C_rx/C_aver;          % �������ҳ��ֲ�  
% cftool
% ���r_nd1��Cr_nd1���ú���cftool������Ϲ���������������ҳ��ֲ�������ϣ��������K�׶���ʽ
% ������
% Linear model Poly6:
%      f(x) = p1*x^6 + p2*x^5 + p3*x^4 + p4*x^3 + p5*x^2 + 
%                     p6*x + p7
% Coefficients (with 95% confidence bounds):
%        p1 =      -40.83  (-40.83, -40.83)
%        p2 =       87.21  (87.21, 87.21)
%        p3 =      -59.43  (-59.43, -59.43)
%        p4 =       11.86  (11.86, 11.86)
%        p5 =      -3.754  (-3.754, -3.754)
%        p6 =       4.938  (4.938, 4.938)
%        p7 = -5.782e-005  (-5.782e-005, -5.782e-005)
% Goodness of fit:
%   SSE: 4.697e-026
%   R-square: 1
%   Adjusted R-square: 1
%   RMSE: 1.56e-014
% ����������ҳ��ֲ�����6�׶���ʽ
% Cr_nd2=-40.83*r_nd1.^6 + 87.21*r_nd1.^5-59.43*r_nd1.^4+11.86*r_nd1.^3-3.754*r_nd1.^2+4.938*r_nd1-5.782e-005;
% figure(2)
% subplot(211)
% plot(x,C_rx,'r-')
% xlabel('չ��r (mm)')
% ylabel('���ٻ�ʵ���ҳ�(mm) ')
% legend('���ٻ�ʵ���ҳ�')
% title('R_wingeff=3.004mm: �Ἰ����ò�����ٻ�ʵ���ҳ��ֲ�')
% subplot(212)
% plot(r_nd,Cr_nd2,'b-')
% xlabel('չ��\itr_{\rmnd} (I)')
% ylabel('�������ҳ�(I)')
% legend('�������ҳ�')
% title('R_wingeff_nd=1: �Ἰ����ò���������ҳ��ֲ�')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ����������������(nondimention_aerodynamic_component)����⡪��F_nd
Coeff=polyfit(r_nd1,Cr_nd1,6);  % ����ʽϵ��  % Cr_nd2=polyval(Coeff,r_nd1);
syms r_nd
Cr_nd2=vpa(poly2sym(Coeff,r_nd),6);   % �������ҳ��ֲ�Ϊ6�׶���ʽ��������Cr_nd3=Cr_nd2;ֻ���Ա�����r_nd1�����r_nd
% Cr_nd2 =-40.827*r_nd^6+87.2061*r_nd^5-59.4281*r_nd^4+11.8648*r_nd^3-3.75417*r_nd^2+4.938*r_nd-0.0000578229;
C_nd=Cr_nd2;
% ע�⡪���öγ����мǲ����޸ģ�ǰ��ֻҪ��֤������ȷ���������ҳ��ֲ����ɡ�
%���µĹ�ʽӦʹ�ú���������ٵ��ҳ��ֲ���ʽC_nd
R2nd2=double(int(r_nd^2*C_nd,r_nd,0,1)); %��������صĻ�ת�뾶��ƽ��
R1nd1=double(int(r_nd*C_nd,r_nd,0,1));     %һ������صĻ�ת�뾶
S_nd=double(int(C_nd,r_nd,0,1));                %�����ٳ����
disp(['��������صĻ�ת�뾶��ƽ��: r2_2nd=' num2str(R2nd2)  ' ���ٵ�λ��mm^4'])
disp(['��������صĻ�ת�뾶: r_2nd=' num2str(sqrt(R2nd2))  ' ���ٵ�λ��mm^3'])
disp(['һ������صĻ�ת�뾶: r_1nd=' num2str(R1nd1)  ' ���ٵ�λ��mm^3'])
disp(['�����ٳ����Swing_nd=' num2str(S_nd)  ' ���ٵ�λ��mm^2'])
% C_nd=vpa(C_nd,5)  %����vpa�����ű��ʽ�е���ֵ(��ת��Ϊ���������ı�ֵ,������)��ת��Ϊʮ����С����ʾ��
% ����xr_nd=0�����Ա���r_nd��ȡֵ��ΧΪ:r_nd��[0,1], ��ȡ1ʱr_nd=R=11.95 /mm
fx2=(r_nd+xr_nd)^2*C_nd;    % ������������F_nd��ԭʼ��������
% fx3=vpa(fx2,5)
fx4=expand(fx2);
F_nd=double(int(fx4,r_nd,0,1));                    % Result: F_nd =0.46392;
disp(['������������F_nd=' num2str(F_nd)  ' ���ٵ�λ��mm^4'])
F_nd2=R2nd2+2*xr_nd*R1nd1+xr_nd^2;    %ʹ����������Ҳ��ȷ; ���:F_nd2 =0.46392;
disp(['������������F_nd=' num2str(F_nd2)  ' ���ٵ�λ��mm^4'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ���:
% ��������صĻ�ת�뾶��ƽ��: r2_2nd=0.3358 ���ٵ�λ��mm^4
% ��������صĻ�ת�뾶: r_2nd=0.57949 ���ٵ�λ��mm^3  ����������������������Ҫ������
% һ������صĻ�ת�뾶: r_1nd=0.53033 ���ٵ�λ��mm^3
% �����ٳ����Swing_nd=1 ���ٵ�λ��mm^2
% ������������F_nd=0.46392 ���ٵ�λ��mm^4    ����������������������Ҫ������
% ������������F_nd=0.46392 ���ٵ�λ��mm^4    ����������������������Ҫ������
%  �Ա������� (1) ��һ�ַ�ʽ����õ�����������������
% F_nd1=r2_nd^2+2*xr_nd1*r1_nd+xr_nd1^2;   %����������������F_nd1, �������: F_nd1 =0.4635
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (2) ��������������ء��������Ч���۵�������λ��Y_nd
% % ���-ƫ��-����ϵԭ��ľ���
% R_proximal=xr;                                                    % xr=3.19;     %RJ Wood��Ƶĳ��\mm
% R_distal=R_wingeff+xr;                                        % yr=0.73;    %RJ Wood��Ƶĳ��\mm
% x=linspace(R_proximal,R_distal,200);
% yr_lead=-0.08249*x.^6+0.9167*x.^5-4.04*x.^4+8.872*x.^3-10.06*x.^2+5.674*x-0.413;  
% yr_trail=-0.0333*x.^6+0.504*x.^5-2.795*x.^4+7.258*x.^3-8.769*x.^2+3.739*x+0.1282;
syms r_nd
yr=0;                              % Ťת��ͨ�������������ߡ�\mm
yr_nd=yr/C_avereff;       % y-root offset  ����������ƫ�þ��� yr_nd = 0.0214;
yr_lead1=-0.08249*r_nd^6+0.9167*r_nd^5-4.04*r_nd^4+8.872*r_nd^3-10.06*r_nd^2+5.674*r_nd-0.413;
yr_leadnd=yr_lead1/C_avereff;
% % ����(1)������ǰ��Ե�����������ٻ�yr_leadnd����yr_trailnd�������
% yr_trail1=-0.0333*r_nd^6+0.504*r_nd^5-2.795*r_nd^4+7.258*r_nd^3-8.769*r_nd^2+3.739*r_nd+0.1282;
% yr_trailnd=yr_trail1/C_avereff;
% yr_nd1=expand((yr_leadnd^4+yr_trailnd^4)/4);     %��������2
% Y_rnd=double(int(yr_nd1,r_nd,0,1));                         % ���: Y_nd=0.15341;  %  mm
% ����(2)������ǰԵ�����������ٻ�yr_leadnd�������ٻ��ҳ��ֲ�C_nd���
y0=yr_nd+yr_leadnd-C_nd;
y1=yr_nd+yr_leadnd;
yr_nd2=(y1^4+y0^4)/4;   % ע�������Ťת��ֱ��ͨ�����ͳ�������ߣ�ǰԵ�����ĳ��չ������Ψһ
Y_rnd=double(int(yr_nd2,r_nd,0,1));                           % ���: Y_rnd2=0.1402;
disp(['��Ч���۵�������λ��Y_nd=' num2str(Y_rnd)  ' ���ٵ�λ��mm'])  % �������Maple���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (3) ������ЧӦ��ϵ������������ת����������������������������
% yr_lead1=-0.08249*r_nd^6+0.9167*r_nd^5-4.04*r_nd^4+8.872*r_nd^3-10.06*r_nd^2+5.674*r_nd-0.413;
% yr_leadnd=yr_lead1/C_avereff;
% yr=0;                              % Ťת��ͨ�������������ߡ�\mm
% yr_nd=yr/C_avereff;       % y-root offset  ����������ƫ�þ��� yr_nd = 0.0214;
yr_hnd=C_nd/2-yr_leadnd-yr_nd;    %ת����ƫ�������е��ƫ��������
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
% �������������ز�������������ת������
fx5=(r_nd+xr_nd)*C_nd^2*yr_hnd;
fr_nd1=expand(fx5);
fx6=C_nd^2*(yr_hnd^2+C_nd^2/32);
fr_nd2=expand(fx6);
I_xzamnd=int(fr_nd1,0,1);
I_xxamnd=int(fr_nd2,0,1);
I1=double(I_xzamnd);                                          % ������,���ٻ���λ��mm^5;          % I1 =0.2666;  
I2=double(I_xxamnd);                                          % ������,���ٻ���λ��mm^5;          % I2 =0.1637;  
Rou=1.225*10^(-3);           %��λ��Kg/m^3=10^6/(10^3)^3=10^(-3)mg/mm^3
I_xzam=pi*Rou*C_avereff^3*R_wingeff^2*I1/4;  % ��λ�� mg.mm^2;   % I_xzam =2.3951;
I_xxam=pi*Rou*C_avereff^4*R_wingeff*I2/4;      % ��λ�� mg.mm^2;   % I_xxam = 0.3279;
%%%%%%%%%%%%%%%%%%%%%%%%
% ����������������
fx7=(r_nd+xr_nd)*C_nd^2;                                  % ����C_nd��ǰ��Ե����ʽ����֮�����
fr_nd3=expand(fx7);
fx8=C_nd^2*yr_hnd;
fr_nd4=expand(fx8);
I3=double(int(fr_nd3,0,1));                                  % ������,���ٻ���λ��mm^4; ��ת����������Ҳ�õ���I3
I4=double(int(fr_nd4,0,1));                                  % ������,���ٻ���λ��mm^4; 
I3z=pi*Rou*C_avereff^2*R_wingeff^2*I3/4;       % ��λ�� mg.mm
I4z=pi*Rou*C_avereff^3*R_wingeff*I4/4;            % ��λ�� mg.mm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wing_para_output=zeros(1,9);
wing_para_output=[R_wingeff,C_avereff,F_nd,Y_rnd,I_xzam,I_xxam,I3,I3z,I4z];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
