# $Experiment 5：$$FTBS$格式

## 杨乐园  PB18010496

#### 问题描述

  1. 针对下述偏微分方程初值问题：
     $$
     \left\{
     \begin{array}{rcl}
     u_{t}+u_{x} &=& 0,\ \ -\infin<x<+\infin,t>0  \\
     u(x,0)&=&
     \left\{
     \begin{array}{ll}
     1, & 0.4\leq x \leq0.6 \\
     0, & else
     \end{array}
     \right.
     \end{array}
     \right.
     $$

构造其$FTBS$格式，用其分别计算$t=1.0,2.0,5.0$时刻的数值解以及该方程的精确解，并绘图，与精确解作比较，给出相应的评论（针对耗散性、色散性）。其中$r=\frac{\Delta t}{\Delta x}$分别取$0.2,0.8$，空间步长$\Delta x=0.05$。

#### 数值方法

​		记$v_{j}^{n}\thickapprox u(x_{j},t_{n})$，得到$FTBS$格式：$v_{j}^{n+1}=v_{j}^{n}-\frac{\Delta t}{\Delta x}(v_{j}^{n}-v_{j-1}^{n})$。

#### 数值结果

​		我们有如下数值求解结果：

从上到下依次为$T=1.0,2.0,5.0$，从左到右依次为真解图像、$r=0.2$数值解图像、$r=0.8$数值解图像。

<center class="half">    
    <img src="E:\study_materials\PDEns\05\ori1.png" width="215" height="180"/>    
    <img src="E:\study_materials\PDEns\05\lam11.png" width="215" height="180"/> 
    <img src="E:\study_materials\PDEns\05\lam21.png" width="215" height="180"/> 
</center>

   <center class="half">    
    <img src="E:\study_materials\PDEns\05\ori2.png" width="215" height="222"/>    
    <img src="E:\study_materials\PDEns\05\lam12.png" width="215" height="222"/> 
    <img src="E:\study_materials\PDEns\05\lam22.png" width="215" height="222"/> 
</center>

 <center class="half">    
    <img src="E:\study_materials\PDEns\05\ori3.png" width="215" height="222"/>    
    <img src="E:\study_materials\PDEns\05\lam13.png" width="215" height="222"/> 
    <img src="E:\study_materials\PDEns\05\lam23.png" width="215" height="222"/> 
</center>

   ​		通过观察上述数值求解结果与方程真解对比我们发现，该方程具有耗散性。并且可以看到$r=0.2$时耗散更迅速。

#### 代码

​		其中数值求解代码与绘图代码参见附件！