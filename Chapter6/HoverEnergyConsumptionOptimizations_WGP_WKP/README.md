# Wing geometry and kinematic parameters (WGP and WKP) optimization of flapping wing hovering flight for minimum energy

## Here, these codes are written for the following papers:

**[1] Xijun Ke, Weiping Zhang, Xuefei Cai, and Wenyuan Chen,"*Wing geometry and kinematic parameter optimization of flapping wing hovering flight for minimum energy*", ***Aerospace Science and Technology***, 2017(64), 192-203. https://doi.org/10.1016/j.ast.2017.01.019. (IF: 5.6)**

**[2] Xijun Ke and Weiping Zhang,"*Wing geometry and kinematic parameter optimization of flapping wing hovering flight*", ***Applied Sciences***, 2016, 6,390,1-35. https://doi.org/10.3390/app6120390. (IF: 2.7)** 

## Abstract and Keywords

**Abstract**: *The optimizations of wing geometry parameters (WGP) and wing kinematic parameters (WKP) to minimize the energy consumption of flapping wing hovering flight are performed by using a revised quasi-steady aerodynamic model and hybrid genetic algorithm (hybrid-GA). The parametrization method of dynamically scaled wing with the non-dimensional conformal feature of fruit fly's wing is firstly developed for the optimization involving the WGP. And the objective function of optimization is formed on basis of the power density model with the additional penalty items of lift-to-weight ratio, boundary constraints, aspect ratio (AR) and Reynolds number (Re). The obtained optimal WGP and WKP are separately substituted into the power density model to estimate the instantaneous forces and the power output again. The lower power density, flapping frequency and larger WGP for the combined optimal WGP and WKP are obtained in comparison with the estimated values for hovering fruit fly. These results might arise from the effect of strong coupling relationship between WGP and WKP via AR and Re on minimization of power density under the condition of lift balancing weight. Moreover, the optimal flapping angle manifests the harmonic profile, and the optimal pitch angle possesses the round trapezoidal profile with certain faster time scale of pitch reversal. The conceptual model framework of combined optimization provides a useful way to design fundamental parameters of biomimic flapping wing micro aerial vehicle.*

**Keywords**: *Flapping wing micro aerial vehicle, combined optimization, and nonlinear couple.*

---------------------------------------------------------------------------------------------------------   

## Highlights:

Wing geometry and kinematic parameter optimization of flapping wing hovering flight for minimum energy

**(1) The combined optimizations of wing geometry and kinematic parameters are firstly performed to minimize the energy of flapping wing hovering flight.**

**(2) The parametrization description of dynamically scaled wing with non-dimensional conformal feature is firstly developed by including the parametrization of wing leading-edge profiles, the definition of pitch axis and mass properties.**

**(3) The revised quasi-steady aerodynamic model is developed on basis of previous aerodynamic model by additionally introducing the rotational circulation moments and aerodynamic damping moment along the chordwise axis of wing planform.**

---------------------------------------------------------------------------------------------------------   

## Wing kinematic angle definition, coordination axis convention and stroke plane of wing motion pattern for Diptera Fruit Fly

![Fig1_Left_wing_body_model](https://github.com/xijunke/HoverEnergyConsumptionOptimizations_WGP_WKP/blob/main/pic_png/Fig1_Left_wing_body_model_s1_4_12_2.png)

![Fig1_Right_wing_body_model](https://github.com/xijunke/HoverEnergyConsumptionOptimizations_WGP_WKP/blob/main/pic_png/Fig1_Right_wing_body_model_s1_4_13_2.png)

***Figure 1: Coordinate systems and definition of right-wing Euler angles relative to the stroke plane in right wing root frame (x_{rr}y_{rr}z_{rr}) [9, 10].***

## All the Figures relavant to this paper and chapter's content have been listed as following:

### 1. Here,the optimal post-processed resuts are for ***fixed aerodynamic lift and drag coefficients (C_F)***

Notes: Aerodynamic lift and drag coefficients (C_F) do not vary with wing morphology parameters

#### For the combined optimal wing shape parameters and kinematic parameters, aerodynamic force, inertial force, aerodynamic power, inertial power, wing pitch and flapping motion power, and total power output are listed here.

#### Aerodynamic forces and inertia force
```
Aerodynamic force inbcludes normal force to wing plane, which are decomposed into aerodynamic lift and thrust again. And inertial force is also given out.
```

#### Aerodynamic power and inertial power for pitch and flapping axis and total power output

#### Figures for above mentioned aerodynamic force and power, inertial force and power, and wing pitch and flapping motion power, and total power output

##### aerodynamic_forces_inertia_force_lift_and_thrust
![calculated results](https://github.com/xijunke/HoverEnergyConsumptionOptimizations_WGP_WKP/blob/main/WingM6_3_10parameters_fixed_C_F_2/force_moment_power_20160316/pic_png/aerodynamic_forces_inertia_force_lift_and_thrust.png)

##### 沿着扭转轴x-axis的气动力矩
![calculated results](https://github.com/xijunke/HoverEnergyConsumptionOptimizations_WGP_WKP/blob/main/WingM6_3_10parameters_fixed_C_F_2/force_moment_power_20160316/pic_png/%E6%B2%BF%E7%9D%80%E6%89%AD%E8%BD%AC%E8%BD%B4x-axis%E7%9A%84%E6%B0%94%E5%8A%A8%E5%8A%9B%E7%9F%A9.png)

##### 沿着拍打轴z-axis的气动力矩
![calculated results](https://github.com/xijunke/HoverEnergyConsumptionOptimizations_WGP_WKP/blob/main/WingM6_3_10parameters_fixed_C_F_2/force_moment_power_20160316/pic_png/%E6%B2%BF%E7%9D%80%E6%8B%8D%E6%89%93%E8%BD%B4z-axis%E7%9A%84%E6%B0%94%E5%8A%A8%E5%8A%9B%E7%9F%A9.png)

##### 沿着扭转轴x_{rw}-axis的气动功率和惯性功率
![calculated results](https://github.com/xijunke/HoverEnergyConsumptionOptimizations_WGP_WKP/blob/main/WingM6_3_10parameters_fixed_C_F_2/force_moment_power_20160316/pic_png/%E6%B2%BF%E7%9D%80%E6%89%AD%E8%BD%AC%E8%BD%B4x_%7Brw%7D-axis%E7%9A%84%E6%B0%94%E5%8A%A8%E5%8A%9F%E7%8E%87%E5%92%8C%E6%83%AF%E6%80%A7%E5%8A%9F%E7%8E%87.png)

##### 沿着拍打轴z_{rr}-axis的气动功率和惯性功率
![calculated results](https://github.com/xijunke/HoverEnergyConsumptionOptimizations_WGP_WKP/blob/main/WingM6_3_10parameters_fixed_C_F_2/force_moment_power_20160316/pic_png/%E6%B2%BF%E7%9D%80%E6%8B%8D%E6%89%93%E8%BD%B4z_%7Brr%7D-axis%E7%9A%84%E6%B0%94%E5%8A%A8%E5%8A%9F%E7%8E%87%E5%92%8C%E6%83%AF%E6%80%A7%E5%8A%9F%E7%8E%87.png)

##### 沿着扭转轴x_{rw}-axis和拍打轴z_{rr}-axis的惯性功率
![calculated results](https://github.com/xijunke/HoverEnergyConsumptionOptimizations_WGP_WKP/blob/main/WingM6_3_10parameters_fixed_C_F_2/force_moment_power_20160316/pic_png/%E6%B2%BF%E7%9D%80%E6%89%AD%E8%BD%AC%E8%BD%B4x_%7Brw%7D-axis%E5%92%8C%E6%8B%8D%E6%89%93%E8%BD%B4z_%7Brr%7D-axis%E7%9A%84%E6%83%AF%E6%80%A7%E5%8A%9F%E7%8E%87.png)

##### 针对最优翅膀形貌参数的翅扭转功率_拍打功率和总功率输出
![calculated results](https://github.com/xijunke/HoverEnergyConsumptionOptimizations_WGP_WKP/blob/main/WingM6_3_10parameters_fixed_C_F_2/force_moment_power_20160316/pic_png/%E9%92%88%E5%AF%B9%E6%9C%80%E4%BC%98%E7%BF%85%E8%86%80%E5%BD%A2%E8%B2%8C%E5%8F%82%E6%95%B0%E7%9A%84%E7%BF%85%E6%89%AD%E8%BD%AC%E5%8A%9F%E7%8E%87_%E6%8B%8D%E6%89%93%E5%8A%9F%E7%8E%87%E5%92%8C%E6%80%BB%E5%8A%9F%E7%8E%87%E8%BE%93%E5%87%BA.png)

### 2. Here,the optimal post-processed resuts are for ***aerodynamic lift and drag coefficients (C_F) varying with wing morphology parameters***

Notes: Aerodynamic lift and drag coefficients (C_F) vary with wing morphology parameters

#### For the combined optimal wing shape parameters and kinematic parameters, aerodynamic force, inertial force, aerodynamic power, inertial power, wing pitch and flapping motion power, and total power output are listed here.

#### Aerodynamic forces and inertia force
```
Aerodynamic force inbcludes normal force to wing plane, which are decomposed into aerodynamic lift and thrust again. And inertial force is also given out.
```

#### Aerodynamic power and inertial power for pitch and flapping axis and total power output

#### Figures for above mentioned aerodynamic force and power, inertial force and power, and wing pitch and flapping motion power, and total power output

##### aerodynamic_forces_inertia_force_lift_and_thrust
![calculated results](https://github.com/xijunke/HoverEnergyConsumptionOptimizations_WGP_WKP/blob/main/WingM6_3_10parameters_variable_C_F_2/force_moment_power_20160316/pic_png/aerodynamic_forces_inertia_force_lift_and_thrust.png)

##### 沿着扭转轴x-axis的气动力矩
![calculated results](https://github.com/xijunke/HoverEnergyConsumptionOptimizations_WGP_WKP/blob/main/WingM6_3_10parameters_variable_C_F_2/force_moment_power_20160316/pic_png/%E6%B2%BF%E7%9D%80%E6%89%AD%E8%BD%AC%E8%BD%B4x-axis%E7%9A%84%E6%B0%94%E5%8A%A8%E5%8A%9B%E7%9F%A9.png)

##### 沿着拍打轴z-axis的气动力矩
![calculated results](https://github.com/xijunke/HoverEnergyConsumptionOptimizations_WGP_WKP/blob/main/WingM6_3_10parameters_variable_C_F_2/force_moment_power_20160316/pic_png/%E6%B2%BF%E7%9D%80%E6%8B%8D%E6%89%93%E8%BD%B4z-axis%E7%9A%84%E6%B0%94%E5%8A%A8%E5%8A%9B%E7%9F%A9.png)

## 沿着扭转轴x_{rw}-axis的气动功率和惯性功率
![calculated results](https://github.com/xijunke/HoverEnergyConsumptionOptimizations_WGP_WKP/blob/main/WingM6_3_10parameters_variable_C_F_2/force_moment_power_20160316/pic_png/%E6%B2%BF%E7%9D%80%E6%89%AD%E8%BD%AC%E8%BD%B4x_%7Brw%7D-axis%E7%9A%84%E6%B0%94%E5%8A%A8%E5%8A%9F%E7%8E%87%E5%92%8C%E6%83%AF%E6%80%A7%E5%8A%9F%E7%8E%87.png)

##### 沿着拍打轴z_{rr}-axis的气动功率和惯性功率
![calculated results](https://github.com/xijunke/HoverEnergyConsumptionOptimizations_WGP_WKP/blob/main/WingM6_3_10parameters_variable_C_F_2/force_moment_power_20160316/pic_png/%E6%B2%BF%E7%9D%80%E6%8B%8D%E6%89%93%E8%BD%B4z_%7Brr%7D-axis%E7%9A%84%E6%B0%94%E5%8A%A8%E5%8A%9F%E7%8E%87%E5%92%8C%E6%83%AF%E6%80%A7%E5%8A%9F%E7%8E%87.png)

##### 沿着扭转轴x_{rw}-axis和拍打轴z_{rr}-axis的惯性功率
![calculated results](https://github.com/xijunke/HoverEnergyConsumptionOptimizations_WGP_WKP/blob/main/WingM6_3_10parameters_variable_C_F_2/force_moment_power_20160316/pic_png/%E6%B2%BF%E7%9D%80%E6%89%AD%E8%BD%AC%E8%BD%B4x_%7Brw%7D-axis%E5%92%8C%E6%8B%8D%E6%89%93%E8%BD%B4z_%7Brr%7D-axis%E7%9A%84%E6%83%AF%E6%80%A7%E5%8A%9F%E7%8E%87.png)

##### 针对最优翅膀形貌参数的翅扭转功率_拍打功率和总功率输出
![calculated results](https://github.com/xijunke/HoverEnergyConsumptionOptimizations_WGP_WKP/blob/main/WingM6_3_10parameters_variable_C_F_2/force_moment_power_20160316/pic_png/%E9%92%88%E5%AF%B9%E6%9C%80%E4%BC%98%E7%BF%85%E8%86%80%E5%BD%A2%E8%B2%8C%E5%8F%82%E6%95%B0%E7%9A%84%E7%BF%85%E6%89%AD%E8%BD%AC%E5%8A%9F%E7%8E%87_%E6%8B%8D%E6%89%93%E5%8A%9F%E7%8E%87%E5%92%8C%E6%80%BB%E5%8A%9F%E7%8E%87%E8%BE%93%E5%87%BA.png)

---------------------------------------------------------------------------------------------------------   

## The extended quasi-steady aerodynamic and inertial forces/moments model are derived in the following paper:

**Xijun Ke, Weiping Zhang, Jinhao Shi and Weidong Chen,"*The numerical solution for flapping wing hovering wingbeat dynamics*", ***Aerospace Science and Technology***, 110(2021), 106474. (IF: 5.6)**

**https://doi.org/10.1016/j.ast.2020.106474**

**https://www.sciencedirect.com/science/article/abs/pii/S1270963820311561**

### The relevant codes for this paper has opened in the following URL:

**https://github.com/xijunke/HoverWingbeatDynamics**

--------------------------------------------------------------------------------------------------------- 