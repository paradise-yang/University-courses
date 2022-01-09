# $Experiment 1：$有限差分法近似求解

## 杨乐园  PB18010496

#### 问题描述

​		利用有喜爱能查分法求下述偏微分方程初值问题在时刻$t = 0.3$的近似解：  
$$
\left 
\{
\begin{array}{ll}
u_{t}=u_{x}, & -\infin<x<+\infin,t>0 \\
u(x,0)=\sin 2\pi x , & -\infin<x<+\infin \\
\text{Periodic boundary condition}, & T = 1 
\end{array}
\right.
$$
该方程的精确解为$u(x, t) = sin(2π(x + t)) $，对时空区域$[0,1]\times[0,1]$均匀剖分如下：
$$
\text{时间：} t_{n} = n · ∆t,n = 0, 1, 2, . . . , N, 时间步长 ∆t = \frac{1}{N} \\
空间: x_{j} = j · ∆x, j = 0, 1, 2, . . . , J, 时间步长 ∆x =\frac{1}{J}
$$
​		问题如下：

$ 1.1: 取 ∆x = 0.02, ∆t = 0.01, $求上述偏微分方程初值问题在时刻$ t = 0.3 $的近似解，并画图比较精确解和精确解 (画图)。
$ 1.2 : 取 ∆x = 0.02, ∆t = 0.03, $求上述偏微分方程初值问题在时刻$ t = 0.3 $的近似解，并画图比较精确解和精确解 (画图)。
$ 1.3 : $对上述两种实验结果进行描述, 分析并评论。  

#### 数值方法

​		记$v_{j}^{n}\thickapprox u(x_{j},t_{n})$，有导数近似$u_{t} \approx \frac{u(x,t+\Delta t)-u(x,t)}{\Delta t}$，$u \approx \frac{u(x+\Delta x,t)-u(x,t)}{\Delta x}$，从而由偏微分方程得到相应的离散方程如下：
$$
v_{j}^{n}=v_{j}^{n}+\frac{\Delta t}{\Delta x}(v_{j+1}^{n}-v_{j}^{n})
$$
其中定解条件为：初始条件：$v_{j}^{0}=sin2\pi x_{j}$，边界条件：$v_{j}^{n}=v_{j+J}^{n}$。

#### 数值结果

1. $\Delta t=0.01$

   ![001](E:\study_materials\PDEns\01\001.png)

<center class="half">
    <img src="E:\study_materials\PDEns\01\ori.jpg" width="210" height="170"/>
    <img src="E:\study_materials\PDEns\01\pde1.jpg" width="210" height="170"/>
    <img src="E:\study_materials\PDEns\01\pde11.jpg" width="210" height="170"/>
</center>

​                    精确解                                           数值解                                        对比图

2. $\Delta t=0.03$

   ![003](E:\study_materials\PDEns\01\003.png)

<center class="half">
    <img src="E:\study_materials\PDEns\01\ori.jpg" width="210" height="170"/>
    <img src="E:\study_materials\PDEns\01\pde2.jpg" width="210" height="170"/>
    <img src="E:\study_materials\PDEns\01\pde22.jpg" width="210" height="170"/>
</center>

​                    精确解                                           数值解                                        对比图

#### 讨论

​		通过对比两结果可以发现：对于$\Delta t=0.01$时，数值求解结果在区间边界处整体偏小，而在区间中段时却整体偏大；而对于$\Delta t=0.03$时，却正好相反，即数值求解结果在区间边界处整体偏大，而在区间中段时却整体偏小。
