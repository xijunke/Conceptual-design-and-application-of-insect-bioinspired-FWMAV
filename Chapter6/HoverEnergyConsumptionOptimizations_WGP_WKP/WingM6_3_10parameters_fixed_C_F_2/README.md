# Here,the optimal post-processed resuts are for fixed aerodynamic lift and drag coefficients (C_F) 

Notes: Aerodynamic lift and drag coefficients (C_F) do not vary with wing morphology parameters

# 说明注释区别

本文件夹的程序原始名称为: ***WingM6_3_10variable_group_Wang_wingbeatM_fixed_C_F_2***

早期改编自WingM5_4_10variable_group，20160106日之前该论文的一些气动力和力矩以及功率的计算程序有小错误

目前改程序已经修改，几近正确

0.hybrid_GA_fminsearch_WingM4_4_2由hybrid_GA_fminsearch_WingM4_4_1进化修改而来

1.该文件夹下为10变量GA混合优化之程序—fmincon—搜索算法

2.含翅膀形貌学和运动学参数(采用Wang ZJ wingbeat pattern 运动学规律)共计10个变量

3.功率调用函数和气动力调用函数一起调用同一个函数Aero_F_M_fruitfly——程序较简洁

4.变量约束区间改小——注意数据的上下限修改——直指频率约束f∈[0,300];

5.kenimatics_wing_and_AoA_fruitfly_sim输出(1000*9)矩阵

## For the combined optimal wing shape parameters and kinematic parameters, aerodynamic force, inertial force, aerodynamic power, inertial power, wing pitch and flapping motion power, and total power output are listed here.

## Aerodynamic forces and inertia force
```
Aerodynamic force inbcludes normal force to wing plane, which are decomposed into aerodynamic lift and thrust again. And inertial force is also given out.
```

## Aerodynamic power and inertial power for pitch and flapping axis and total power output

# Figures for above mentioned aerodynamic force and power, inertial force and power, and wing pitch and flapping motion power, and total power output

## aerodynamic_forces_inertia_force_lift_and_thrust
![calculated results](https://github.com/xijunke/HoverEnergyConsumptionOptimizations_WGP_WKP/blob/main/WingM6_3_10parameters_fixed_C_F_2/force_moment_power_20160316/pic_png/aerodynamic_forces_inertia_force_lift_and_thrust.png)

## 沿着扭转轴x-axis的气动力矩
![calculated results](https://github.com/xijunke/HoverEnergyConsumptionOptimizations_WGP_WKP/blob/main/WingM6_3_10parameters_fixed_C_F_2/force_moment_power_20160316/pic_png/%E6%B2%BF%E7%9D%80%E6%89%AD%E8%BD%AC%E8%BD%B4x-axis%E7%9A%84%E6%B0%94%E5%8A%A8%E5%8A%9B%E7%9F%A9.png)

## 沿着拍打轴z-axis的气动力矩
![calculated results](https://github.com/xijunke/HoverEnergyConsumptionOptimizations_WGP_WKP/blob/main/WingM6_3_10parameters_fixed_C_F_2/force_moment_power_20160316/pic_png/%E6%B2%BF%E7%9D%80%E6%8B%8D%E6%89%93%E8%BD%B4z-axis%E7%9A%84%E6%B0%94%E5%8A%A8%E5%8A%9B%E7%9F%A9.png)

## 沿着扭转轴x_{rw}-axis的气动功率和惯性功率
![calculated results](https://github.com/xijunke/HoverEnergyConsumptionOptimizations_WGP_WKP/blob/main/WingM6_3_10parameters_fixed_C_F_2/force_moment_power_20160316/pic_png/%E6%B2%BF%E7%9D%80%E6%89%AD%E8%BD%AC%E8%BD%B4x_%7Brw%7D-axis%E7%9A%84%E6%B0%94%E5%8A%A8%E5%8A%9F%E7%8E%87%E5%92%8C%E6%83%AF%E6%80%A7%E5%8A%9F%E7%8E%87.png)

## 沿着拍打轴z_{rr}-axis的气动功率和惯性功率
![calculated results](https://github.com/xijunke/HoverEnergyConsumptionOptimizations_WGP_WKP/blob/main/WingM6_3_10parameters_fixed_C_F_2/force_moment_power_20160316/pic_png/%E6%B2%BF%E7%9D%80%E6%8B%8D%E6%89%93%E8%BD%B4z_%7Brr%7D-axis%E7%9A%84%E6%B0%94%E5%8A%A8%E5%8A%9F%E7%8E%87%E5%92%8C%E6%83%AF%E6%80%A7%E5%8A%9F%E7%8E%87.png)

## 沿着扭转轴x_{rw}-axis和拍打轴z_{rr}-axis的惯性功率
![calculated results](https://github.com/xijunke/HoverEnergyConsumptionOptimizations_WGP_WKP/blob/main/WingM6_3_10parameters_fixed_C_F_2/force_moment_power_20160316/pic_png/%E6%B2%BF%E7%9D%80%E6%89%AD%E8%BD%AC%E8%BD%B4x_%7Brw%7D-axis%E5%92%8C%E6%8B%8D%E6%89%93%E8%BD%B4z_%7Brr%7D-axis%E7%9A%84%E6%83%AF%E6%80%A7%E5%8A%9F%E7%8E%87.png)

## 针对最优翅膀形貌参数的翅扭转功率_拍打功率和总功率输出
![calculated results](https://github.com/xijunke/HoverEnergyConsumptionOptimizations_WGP_WKP/blob/main/WingM6_3_10parameters_fixed_C_F_2/force_moment_power_20160316/pic_png/%E9%92%88%E5%AF%B9%E6%9C%80%E4%BC%98%E7%BF%85%E8%86%80%E5%BD%A2%E8%B2%8C%E5%8F%82%E6%95%B0%E7%9A%84%E7%BF%85%E6%89%AD%E8%BD%AC%E5%8A%9F%E7%8E%87_%E6%8B%8D%E6%89%93%E5%8A%9F%E7%8E%87%E5%92%8C%E6%80%BB%E5%8A%9F%E7%8E%87%E8%BE%93%E5%87%BA.png)