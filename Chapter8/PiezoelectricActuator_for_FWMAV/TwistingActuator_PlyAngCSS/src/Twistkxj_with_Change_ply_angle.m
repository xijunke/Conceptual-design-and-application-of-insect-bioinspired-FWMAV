clear all;
clc;
syms E v E1 E2 v12 v21 G12  z1 z2 z3 z0 E3 L W d31 d32 y p1 p2;
E3=200;                                    %v
%t=0:0.0001:0.025;
%E3=100*sin(2*pi*120*t)+100;     
%plot(t,E3)
%PZT
p1=7800;   %kg/m
E=6.2e10;   v=0.31; %Gpa &G12=24e9
d31=-320e-12; d32=-320e-12; %pm/v
d=[d31;d32;0];
%CF
p2=1500;   %kg/m
E1=350e9;  E2=7e9; G12=5e9;   %Gpa
v12=0.33; v21=E2*v12/E1;         %
%
W=2e-3; L=10e-3;    %
z0=-143.5e-6; z1=-63.5e-6; z2=63.5e-6; z3=143.5e-6;     %m

%PZT
Q112=E/(1-v^2);  Q222=Q112;
Q122=E*v/(1-v^2);
Q662=E/(2*(1+v)); G=Q662;
%CF
Q11=E1/(1-v12*v21);
Q12=v12*E2/(1-v12*v21);
Q22=E2/(1-v12*v21);
Q66=G12;
%
i=1;j=1;k=1;
y=zeros(91,1);
Q=zeros(91,1); 
theta=zeros(91,1); 
Tm=zeros(91,1);   
Ud=zeros(91,1);

% PZT layer
Q0=[Q112 Q122 0
      Q122 Q222 0
      0      0    Q662];
% CF layer
Q1=[Q11 Q12 0
      Q12 Q22 0
      0      0    Q66];
%ply angle
angle=0:1:90;  angle=angle';
y=angle*(pi/180);

% 
for i=1:91
% T
m=cos(y(i,1));n=sin(y(i,1));
% PZT mid layer
 Tmid=[m.^2,n.^2,2*m.*n;
      n.^2 m.^2 -2*m.*n;
      -m.*n,m.*n,m.^2-n.^2];
 Tinv=inv(Tmid);Tinvz=Tinv';
 % CF1 layer
 Tfirst=[m.^2,n.^2,-2*m.*n;
      n.^2 m.^2 2*m.*n;
     m.*n,-m.*n,m.^2-n.^2];
 Tinvfirst=inv(Tfirst);Tinvzfirst=Tinvfirst';
 % CF3 layer
 Tthird=[m.^2,n.^2,2*m.*n;
      n.^2,m.^2,-2*m.*n;
      -m.*n,m.*n,m.^2-n.^2];
 Tinvthird=inv(Tthird);Tinvzthird=Tinvthird';
% 
Qmid=Tinv*Q0*Tinvz;
Qfirst=Tinvfirst*Q1*Tinvzfirst;
Qthird=Tinvthird*Q1*Tinvzthird;

% 
A=Qfirst*(z1-z0)+Qmid*(z2-z1)+Qthird*(z3-z2);
B=(Qfirst*(z1^2-z0^2)+Qmid*(z2^2-z1^2)+Qthird*(z3^2-z2^2))*(1/2);
D=(Qfirst*(z1^3-z0^3)+Qmid*(z2^3-z1^3)+Qthird*(z3^3-z2^3))*(1/3);
F=[A,B;B,D];
F=inv(F);
%
Fp=E3*(z2-z1)*Qmid*d;     
Mp=E3*1/2*(z2^2-z1^2)*Qmid*d;
Np=[Fp;Mp];
%
Fext=[0;0;0];
Mext=[0;0;0];   %Mx,My,Mxy;
Next=[Fext;Mext];
kxy=F*(Next+Np);
Q(i,1)=atan(kxy(6,1)*L);          %L
theta(i,1)=Q(i,1)*(180/pi);       % 
Tm(i,1)=W*(F(6,1)*Fp(1,1)+F(6,2)*Fp(2,1)+F(6,3)*Fp(3,1))./F(6,6);    % 
m=L*W*(p1*(z1-z0+z3-z2)+p2*(z2-z1));   %
Ud(i,1)=Q(i,1)*Tm(i,1)/m/2;   % 
end    

figure(1);
plot(angle,theta)
xlabel('Ply Angle (б└бу)');
ylabel('Output twist (deg)')
hold on

figure(2);
plot(angle,Tm)
xlabel('Ply Angle (б└бу)');
ylabel('Moment (N/m)');
hold on

figure(3);
plot(angle,Ud)
xlabel('Ply Angle (б└бу)');
ylabel('Energy density (J/kg)');
