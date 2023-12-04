function coeff_con=aeroF_coeff_constraint(x)
% ������������ϵ��Լ��֮��ĳͷ�����
% k_C_N=x(1);  k_C_R=x(2);  k_add=x(3);  k_inert=x(4);
% fcn_con=aeroF_coeff_constraint(x);
LB=[0.5,0.5,0.2,0.1];     % Lower bound  % 4������ 
UB=[2,3,3,2];             % Upper bound  % 4������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=length(x);
con_min=LB;
con_max=UB;
zeta=zeros(N,1);
% y=zeros(N,1);
parfor i=1:N
    if x(i) < con_min(i)
        zeta(i)=abs(con_min(i)-x(i))/(con_max(i)-con_min(i));
    elseif x(i) > con_max(i)  % 
        zeta(i)=abs(x(i)-con_max(i))/(con_max(i)-con_min(i));
    else % x(i)>=con_min(i) && x(i)<=con_max(i);  % ���������ʽ
        zeta(i)=0;  % disp('����δ�����߽�Լ��');  
    end
end
coeff_con=sum(zeta);
end

