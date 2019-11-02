function C_periEdge=C_perimeter_to_Edge_correction(R_wing,C_aver,xr0,C_maxyaxis)
% C_perimeter_to_Edge_correction
%%%%%%%%%%%%%%%%%%%%%%%%%
R_wingeff=3.004;    %有效翅膀长度(mm)  
C_avereff=0.8854;  % 单位:mm
xr=0.3289;                % x-root offset  \mm
%%%%%%%%%%%%%%%%%%%%%%%%%
% R_wing=R_wingeff;
% C_aver=C_avereff;
% xr0=xr;
x_0_vari=C_maxyaxis;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 翅根-偏离-坐标系原点的距离
R_proximal=xr;                                                    % xr=3.19;     %RJ Wood设计的翅膀—\mm
R_distal=R_wingeff+xr;                                        % yr=0.73;    %RJ Wood设计的翅膀—\mm
x=linspace(R_proximal,R_distal,200);
x_mod_Root=0.636;                                            % mm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 针对翅根和翅尖连线的扭转轴定出的前缘     
% C_maxy=C_lead_ymax-C_maxy*C_max_LtoT; %转动轴偏离弦向中点的偏移量坐标,在扭转轴向上偏移 -C_maxy之后  %XXX
yr_lead=-0.08249*x.^6+0.9167*x.^5-4.04*x.^4+8.872*x.^3-10.06*x.^2+5.674*x-0.413-x_mod_Root;  
yr_trail=-0.0333*x.^6+0.504*x.^5-2.795*x.^4+7.258*x.^3-8.769*x.^2+3.739*x+0.1282-x_mod_Root;
%% 采用多项式函数拟合获得——无量纲前后缘分布函数和弦长分布函数
% (a) 无量纲前缘分布函数
r_nd=(x-xr)/R_wingeff;
yr_leadnd0=yr_lead/C_avereff;
P_coeff_lead=polyfit(r_nd,yr_leadnd0,6);
% (b)  无量纲后缘分布函数
yr_trailnd0=yr_trail/C_avereff;
P_coeff_trail=polyfit(r_nd,yr_trailnd0,6);
% (c)  无量纲弦长分布函数
Cr_nd=yr_leadnd0-yr_trailnd0;
P_coeff_Cr=polyfit(r_nd,Cr_nd,6);  % 多项式系数  % Cr_nd2=polyval(Coeff,r_nd1);
syms r   r_nd   % 无量纲弦长分布为6阶多项式——转换必须有这条指令
yr_leadnd=vpa(poly2sym(P_coeff_lead,r_nd),6); 
yr_trailnd=vpa(poly2sym(P_coeff_trail,r_nd),6);
Cr_nd=vpa(poly2sym(P_coeff_Cr,r_nd),6);  
C_nd=Cr_nd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 初始果蝇翅膀
% x_mod_Root=0.636;    
% y_r_lead=-0.08249*r.^6+0.9167*r.^5-4.04*r.^4+8.872*r.^3-10.06*r.^2+5.674*r-0.413-x_mod_Root;  
% y_r_trail=-0.0333*r.^6+0.504*r.^5-2.795*r.^4+7.258*r.^3-8.769*r.^2+3.739*r+0.1282-x_mod_Root;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 动态缩放比例翅膀
%%%%%%%%%%%%%%%%%%%%%%%%%
f_x_lead=C_aver*yr_leadnd;  % 符号函数
f_x_trail =C_aver*yr_trailnd;  % 符号函数
%%%%%%%%%%%%%%%%%%%%%%%%%
f_x_lead1=inline(vectorize(f_x_lead),'r_nd');  %数值函数
f_x_trail1=inline(vectorize(f_x_trail),'r_nd');   %数值函数
r_nd1=linspace(0,1,200);     % x=linspace(R_proximal,R_distal,200);  % x=r_nd*R_wingeff+xr;
f_x_lead2=f_x_lead1(r_nd1); 
f_x_trail2=f_x_trail1(r_nd1);
rx=r_nd1*R_wing+xr0;  % rx=linspace(0,1,200)*R_wingeff+xr; 
r_x_min=min(rx);  
r_x_max=max(rx);          % r_x_max=100.3289;
% cftool
P_coeff_le=polyfit(rx,f_x_lead2,6);     % 多项式系数
P_coeff_tr=polyfit(rx,f_x_trail2,6);      % 多项式系数
y_r_lead=vpa(poly2sym(P_coeff_le,r),6);    % 符号函数
y_r_trail=vpa(poly2sym(P_coeff_tr,r),6);     % 符号函数
% %%%%%%%%%%%%%%%%%%%%%%%%%
% fun_le=simplify(expand((1+(diff(y_r_lead,r,1))^2)^(1/2)))
fun_le=expand((1+(diff(y_r_lead,r,1))^2)^(1/2));
C_perimeter_le=double(int(fun_le,r,r_x_min,r_x_max));     % C_perimeter_le =3.2791;
% fun_tr=simplify(expand((1+(diff(y_r_trail,r,1))^2)^(1/2)))
fun_tr=expand((1+(diff(y_r_trail,r,1))^2)^(1/2));
C_perimeter_tr=double(int(fun_tr,r,r_x_min,r_x_max));     % C_perimeter_le =3.6067;
E_C_peri=((C_perimeter_le+C_perimeter_tr)/2)/R_wing;               % E_C_peri =1.1461;
C_periEdge=E_C_peri;