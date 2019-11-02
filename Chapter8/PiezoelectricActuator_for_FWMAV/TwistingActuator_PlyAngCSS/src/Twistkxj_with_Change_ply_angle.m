clear all;
clc;
syms E v E1 E2 v12 v21 G12  z1 z2 z3 z0 E3 L W d31 d32 y p1 p2;
E3=200;                                    %外加驱动的电场强度v
%t=0:0.0001:0.025;
%E3=100*sin(2*pi*120*t)+100;     
%plot(t,E3)
%PZT力学性能参数
p1=7800;   %密度kg/m
E=6.2e10;   v=0.31; %弹性和剪切模量Gpa &泊松比 其中G12=24e9
d31=-320e-12; d32=-320e-12; %压电常数pm/v
d=[d31;d32;0];
%CF力学性能参数
p2=1500;   %密度kg/m
E1=350e9;  E2=7e9; G12=5e9;   %弹性和剪切模量Gpa
v12=0.33; v21=E2*v12/E1;         %泊松比
%扭转压电驱动的几何参数
W=2e-3; L=10e-3;    %宽度和长度
z0=-143.5e-6; z1=-63.5e-6; z2=63.5e-6; z3=143.5e-6;     %各层厚度m

%刚度矩阵的元素   各同异性材料PZT
Q112=E/(1-v^2);  Q222=Q112;
Q122=E*v/(1-v^2);
Q662=E/(2*(1+v)); G=Q662;
%刚度矩阵的元素   各向异性材料CF
Q11=E1/(1-v12*v21);
Q12=v12*E2/(1-v12*v21);
Q22=E2/(1-v12*v21);
Q66=G12;
%赋常数
i=1;j=1;k=1;
y=zeros(91,1);
Q=zeros(91,1); 
theta=zeros(91,1); 
Tm=zeros(91,1);   
Ud=zeros(91,1);

% PZT layer刚度矩阵
Q0=[Q112 Q122 0
      Q122 Q222 0
      0      0    Q662];
% CF layer刚度矩阵
Q1=[Q11 Q12 0
      Q12 Q22 0
      0      0    Q66];
%ply angle
angle=0:1:90;  angle=angle';
y=angle*(pi/180);

% 需加循环
for i=1:91
% T变换矩阵
m=cos(y(i,1));n=sin(y(i,1));
% PZT mid layer变换矩阵
 Tmid=[m.^2,n.^2,2*m.*n;
      n.^2 m.^2 -2*m.*n;
      -m.*n,m.*n,m.^2-n.^2];
 Tinv=inv(Tmid);Tinvz=Tinv';
 % CF1 layer变换矩阵
 Tfirst=[m.^2,n.^2,-2*m.*n;
      n.^2 m.^2 2*m.*n;
     m.*n,-m.*n,m.^2-n.^2];
 Tinvfirst=inv(Tfirst);Tinvzfirst=Tinvfirst';
 % CF3 layer变换矩阵
 Tthird=[m.^2,n.^2,2*m.*n;
      n.^2,m.^2,-2*m.*n;
      -m.*n,m.*n,m.^2-n.^2];
 Tinvthird=inv(Tthird);Tinvzthird=Tinvthird';
% 调整刚度矩阵
Qmid=Tinv*Q0*Tinvz;
Qfirst=Tinvfirst*Q1*Tinvzfirst;
Qthird=Tinvthird*Q1*Tinvzthird;

% 多层调整刚度矩阵
A=Qfirst*(z1-z0)+Qmid*(z2-z1)+Qthird*(z3-z2);
B=(Qfirst*(z1^2-z0^2)+Qmid*(z2^2-z1^2)+Qthird*(z3^2-z2^2))*(1/2);
D=(Qfirst*(z1^3-z0^3)+Qmid*(z2^3-z1^3)+Qthird*(z3^3-z2^3))*(1/3);
F=[A,B;B,D];
F=inv(F);
%只考虑单压电层的厚度,电场驱动时的力内力和力矩
Fp=E3*(z2-z1)*Qmid*d;     
Mp=E3*1/2*(z2^2-z1^2)*Qmid*d;
Np=[Fp;Mp];
%外力和外力矩
Fext=[0;0;0];
Mext=[0;0;0];   %含Mx,My,Mxy;
Next=[Fext;Mext];
kxy=F*(Next+Np);
Q(i,1)=atan(kxy(6,1)*L);          %L为长度
theta(i,1)=Q(i,1)*(180/pi);       % 角位移
Tm(i,1)=W*(F(6,1)*Fp(1,1)+F(6,2)*Fp(2,1)+F(6,3)*Fp(3,1))./F(6,6);    % 力矩
m=L*W*(p1*(z1-z0+z3-z2)+p2*(z2-z1));   % 质量
Ud(i,1)=Q(i,1)*Tm(i,1)/m/2;   % 能量密度
end    

figure(1);
plot(angle,theta)
xlabel('Ply Angle (±°)');
ylabel('Output twist (deg)')
hold on

figure(2);
plot(angle,Tm)
xlabel('Ply Angle (±°)');
ylabel('Moment (N/m)');
hold on

figure(3);
plot(angle,Ud)
xlabel('Ply Angle (±°)');
ylabel('Energy density (J/kg)');
