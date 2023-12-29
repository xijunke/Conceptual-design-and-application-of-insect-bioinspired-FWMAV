# 说明注释区别

本文件夹的程序原始名称为: ***WingM5_4_10variable_group***

hybrid_GA_fminsearch_WingM4_4_2不同于hybrid_GA_fminsearch_WingM4_4_1

0.hybrid_GA_fminsearch_WingM4_4_2由hybrid_GA_fminsearch_WingM4_4_1进化修改而来;

1.该文件夹下为10变量GA混合优化之程序—fmincon—搜索算法;

2.含翅膀形貌学和运动学参数(人为设计谐波运动学规律)共计10个变量;

3.功率调用函数和气动力调用函数一起调用同一个函数Aero_F_M_fruitfly——程序较简洁;

4.变量约束区间改小——注意数据的上下限修改——直指频率约束f∈[0,300];

5.kenimatics_wing_and_AoA_fruitfly_sim输出(1000*9)矩阵.


![pic_forces_for_10_optimized_parameters](https://github.com/xijunke/HoverEnergyConsumptionOptimizations_WGP_WKP/blob/main/WingM5_4_10parameters/pic_forces_for_10_optimized_parameters/Forces_for_optimized_combined_wing_shape_motions_parameters.png)

![pic_power_for_10_optimal_parameters](https://github.com/xijunke/HoverEnergyConsumptionOptimizations_WGP_WKP/blob/main/WingM5_4_10parameters/pic_power_for_10_optimal_parameters/Power_for_optimized_combined_wing_shape_motions_parameters.png)


