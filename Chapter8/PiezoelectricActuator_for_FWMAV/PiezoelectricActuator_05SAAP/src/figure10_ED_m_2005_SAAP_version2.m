% 201305272005
clear all; clc;
%PZT
%aa=input('Enter the A:');bb=input('Enter the B:');
p1=7800;   %kg/m
t=127e-6;
U=300;
E3=U/t;            %v
E=6.2e10;   v=0.31;%Gpa &G12=24e9
%E=E/(1-v^2);
d31=-320e-12; d32=d31; %pm/v
d=[d31;d32;0];
%CF
p2=1500;   %kg/m
E1=350e9;  E2=7e9; G12=5e9;   %Gpa
v12=0.35; v21=E2*v12/E1;%
%E1=E1/(1-v12^2);
%E2=E2/(1-v21^2);
%% PZT
Q112=E/(1-v^2);  Q222=Q112;
Q122=E*v/(1-v^2);
% Q662=24E9; G=Q662; % version1
Q662=E/(2*(1+v)); G=Q662; % version2����good
%CF
Q11=E1/(1-v12*v21);
Q12=v12*E2/(1-v12*v21);
Q22=E2/(1-v12*v21);
Q66=G12;
% PZT layer
Q0=[Q112 Q122 0
      Q122 Q222 0
      0      0    Q662];
% CF laye
Q1=[Q11 Q12 0
      Q12 Q22 0
      0      0    Q66];
% 
Qmid=Q1;
Qfirst=Q0;
Qthird=Q0;
%% 
% wnom=4.5/2*1e-3; l=8e-3;   % version1
% l_r=0.5;w_r=7/4.5; % version1
wnom=1e-3; l=5e-3;    % version2����good
l_r=1; w_r=1.5;  % version2����good
%% 
g_delta=(1+2*l_r);
g_c=8*(1-w_r)^3;
% g_d=6*(w_r-1)*(3+4*l_r-2*w_r-4*l_r*w_r); % version1
g_d=-6*(w_r-1)*(3+4*l_r-2*w_r-4*l_r*w_r);  % version2����good
g_e=3*(-2-2*l_r+w_r+2*l_r*w_r)^2*log((2-w_r)/w_r);
G_le=(g_d+g_e)/g_c;
G_F=g_delta/G_le;
G_U=g_delta*G_F
%%
Du=[];%
TD=[];%tip displacement
BF=[];%block force
for tr=0:0.01:2;  % 
   tc=tr*t ;
   %m
   z0=-tc/2-t; 
   z1=-tc/2; 
   z2=tc/2; 
   z3=tc/2+t; 
  
   A=Qfirst*(z1-z0)+Qmid*(z2-z1)-Qthird*(z3-z2);
   B=(Qfirst*(z1^2-z0^2)+Qmid*(z2^2-z1^2)+Qthird*(z3^2-z2^2))*(1/2);
   D=(Qfirst*(z1^3-z0^3)+Qmid*(z2^3-z1^3)+Qthird*(z3^3-z2^3))*(1/3);
   C=pinv([A B;B D]); % pinv_Moore-Penrose pseudoinverse of matrix

   Fp=E3*(z1-z0)*Qfirst*d+E3*(z3-z2)*Qthird*d;  
   Mp=E3*1/2*(z1^2-z0^2)*Qfirst*d-E3*1/2*(z3^2-z2^2)*Qthird*d;
   Np=[Fp;Mp];
   P=C(4,1)*Np(1)+C(4,2)*Np(2)+C(4,4)*Np(4)+C(4,5)*Np(5);
   % du=3/8*P*P/C(4,4)/(2*p1*t+p2*tc)*G_U; % version1
   du=-3/8*P*P/C(4,4)/(2*p1*t+p2*tc)*G_U;  % version2����good
  
   Du=[Du du];
   td=P*l^2/2*g_delta;
   td=td*10^6; % um
   TD=[TD td];
   bf=3*P*wnom/(2*C(4,4)*l)*G_F;
   bf=bf*10^3;  % mN
   BF=[BF bf];
end
%%
tr=0:0.01:2;
tc=127*tr ;
%figure(1);
%plot(tc,Du)
[Dumax ind]=max(Du)  
trDumax = tr(ind);        
td=trDumax*t
m=(2*p1*t+p2*tc)*l*wnom; % 
disp(['td(����CF���)=',num2str(td)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
figure(1);
% title('Predicted effect of passive layer thickness on energy density for a bimorph')
plot(tr,Du)
xlabel('Passive Layer Thickness Ratio(tr)');
ylabel('Energy density (J/kg)');
grid on
axis([0 2 0.8 2.6]);
%%%%%%%%%%%%%%%%%
% plot(tc,Du)
% xlabel('��������tc��um��');
% ylabel('�����ܶ�(J/kg)');
% grid on
% axis([0 250 0.8 2.6]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
figure(2);
plot(tr,TD)
xlabel('Passive Layer Thickness Ratio(tr)');
ylabel('Tip Displacement (um)');
% title('Predicted effect of passive layer thickness on displacement for a bimorph')
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
figure(3);
plot(tr,BF)
xlabel('Passive Layer Thickness Ratio(tr)');
ylabel('Block Force (mN)');
% title('Predicted effect of passive layer thickness on block force for a bimorph')
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%