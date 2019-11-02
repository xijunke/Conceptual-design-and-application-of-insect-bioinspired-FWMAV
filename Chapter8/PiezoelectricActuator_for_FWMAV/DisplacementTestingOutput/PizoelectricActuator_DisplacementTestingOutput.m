% pizoelectric_actuator_displacement_output experimentally mearured by KEYENCE close range laser displacement sensor
delta_t=xlsread('1dis.xls','A4:B1003');                 %∂¡»Î ˝æ›
t=delta_t(:,1);
delta=delta_t(:,2);
plot(t,delta,'r-o','LineWidth',2)
xlabel('t')
ylabel('\delta_{actuator}')
