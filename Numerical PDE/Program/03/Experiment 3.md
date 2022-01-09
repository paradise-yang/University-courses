# $Experiment 3：$$FTCS$近似求解

## 杨乐园  PB18010496

#### 问题描述

​		针对下述偏微分方程初值问题：
$$
\left \{\begin{array}{ll}u_{t}=u_{x}, & -\infin<x<+\infin,t>0 \\u(x,0)=\sin 2\pi x , & -\infin<x<+\infin \\\text{Periodic boundary condition}, & T = 1 \end{array}\right.
$$
该方程的精确解为$u(x, t) = sin(2π(x + t)) $，对时空区域均匀剖分，其中$x_{j} = j · ∆x, j = 0, 1, 2, . . . , J$，空间步长$ ∆x =\frac{1}{J}$，令$\lambda=\frac{\Delta t}{\Delta x}$。取$\lambda=0.5,J=80$，分别取终止时间$T=0.1,0.4,0.8,1.0$。用$FTCS$方法计算其数值解，绘制出最大误差随时间变化图，并给出相应评论。

#### 数值方法

​		记$v_{j}^{n}\thickapprox u(x_{j},t_{n})$，有导数近似$u_{t} \approx \frac{u(x,t+\Delta t)-u(x,t)}{\Delta t}$，$u_{x} \approx \frac{u(x+\Delta x,t)-u(x-\Delta x,t)}{\Delta x}$，从而由偏微分方程得到相应的离散方程如下：
$$
v_{j}^{n+1}=v_{j}^{n}+\frac{\Delta t}{2\Delta x}(v_{j+1}^{n}-v_{j-1}^{n})
$$
其中定解条件为：初始条件：$v_{j}^{0}=sin2\pi x_{j}$，边界条件：$v_{j}^{n}=v_{j+J}^{n}$。

#### 数值结果

1. $t=0.1$

   ![01](E:\study_materials\PDEns\03\01.png)

2. $t=0.4$

   ![04](E:\study_materials\PDEns\03\04.png)

3. $t=0.8$

   ![08](E:\study_materials\PDEns\03\08.png)

4. $t=1.0$

   ![10](E:\study_materials\PDEns\03\10.png)

5. 误差输出以及随时间变化图

    ![error](E:\study_materials\PDEns\03\error.png)
    
    ![tt](E:\study_materials\PDEns\03\tt.png)

#### 讨论

​		通过观察误差随时间变化图可以发现：最大误差随时间增长逐渐变大。