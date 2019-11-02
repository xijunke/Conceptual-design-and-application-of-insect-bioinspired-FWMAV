%% resonant_frequency_cal
m_actuator=38.5182e-6;  % mg
M=0.0566;
m_eff=M*m_actuator;
k_a=581.82; % N/m
zeta=0.1;
f_d=(1/(2*pi))*sqrt(1-zeta^2)*sqrt(k_a/m_eff)
% f_d =2.5870e+003; % Hz