# $Experiment 4：$对比多种方法求解结果

## 杨乐园  PB18010496

#### 问题描述

  1. 针对下述偏微分方程初值问题：
     $$
     \left \{\begin{array}{ll}u_{t}=u_{x}, & -\infin<x<+\infin,t>0 \\u(x,0)=\sin 2\pi x , & -\infin<x<+\infin \\\text{Periodic boundary condition}, & T = 1 \end{array}\right.
     $$

该方程的精确解为$u(x, t) = sin(2π(x + t)) $，对时空区域均匀剖分，其中$x_{j} = j · ∆x, j = 0, 1, 2, . . . , J$，空间步长$ ∆x =\frac{1}{J}$，令$\lambda=\frac{\Delta t}{\Delta x}$。

$\textbf{问1.1}：$取$\lambda=0.5,J=80$，分别取终止时间$T=0.1,0.4,0.8,1.0$。分别用$FTCS$、$Lax-Friedrich$和$Lax-Wendroff$方法计算其数值解，绘制出最大误差随时间变化图，并给出相应评论。

$\textbf{问1.2：}$取$\lambda=0.5,T=1.0$，分别取终止时间$J=10,20,40,80,160$。用$Lax-Wendroff$方法计算其数值解，并与精确解画在同一张图上进行比较，并给出相应评论。

2. 针对下述偏微分方程初值问题：
   $$
   \left \{\begin{array}{ll}u_{t}=u_{x}, & -\infin<x<+\infin,t>0 \\u(x,0)=\sin 2\pi x , & -\infin<x<+\infin \\\text{Periodic boundary condition}, & T = 1 \end{array}\right.
   $$

该方程的精确解为$u(x, t) = sin(2π(x + t)) $，对时空区域均匀剖分，其中$x_{j} = j · ∆x, j = 0, 1, 2, . . . , J$，空间步长$ ∆x =\frac{1}{J}$，令$\lambda=\frac{\Delta t}{\Delta x}$。

$\textbf{问2.1}：$取$T=1.0,J=80$，分别取$\lambda=0.5,1.5$。用$CTCS$格式（$v_{j}^{1}$用$FTFS$格式）计算其数值解，并与精确解画在同一张图上进行比较，并给出相应评论。

$\textbf{问2.2：}$取$\lambda=0.5,T=1.0$，分别取$。用$CTCS$格式（$v_{j}^{1}$用$FTFS$格式）计算其数值解，并与精确解画在同一张图上进行比较，并给出相应评论。

$\textbf{问2.3：}$取$\lambda=0.5,J=80$，分别取$T=0.2,0.5$。用$FTBS$格式计算其数值解，并与精确解画在同一张图上进行比较，并给出相应评论。

#### 数值方法

​		记$v_{j}^{n}\thickapprox u(x_{j},t_{n})$，根据不同格式的导数近似以及偏微分方程得到相应的格式：

1. $FTCS$：$v_{j}^{n+1}=v_{j}^{n}+\frac{\Delta t}{2\Delta x}(v_{j+1}^{n}-v_{j-1}^{n})$

2. $Lax-Friedrich$：$v_{j}^{n+1}=(\frac{\Delta t}{2\Delta x}+\frac{1}{2})v_{j+1}^{n}+(-\frac{\Delta t}{2\Delta x}+\frac{1}{2})v_{j-1}^{n}$

3. $Lax-Wendroff$：$v_{j}^{n+1}=(\frac{\Delta t}{2\Delta x}+\frac{\Delta t^{2}}{2\Delta x^{2}})v_{j+1}^{n}+(1-\frac{\Delta t^{2}}{\Delta x^{2}})v_{j}^{n}+(-\frac{\Delta t}{2\Delta x}+\frac{\Delta t^{2}}{2\Delta x^{2}})v_{j-1}^{n}$

4. $CTCS$：$v_{j}^{n+1}=v_{j}^{n-1}+\frac{\Delta t}{\Delta x}(v_{j+1}^{n}-v_{j-1}^{n})$

5. $FTBS$：$v_{j}^{n+1}=v_{j}^{n}+\frac{\Delta t}{\Delta x}(v_{j}^{n}-v_{j-1}^{n})$

   ​	其中定解条件为：初始条件：$v_{j}^{0}=sin2\pi x_{j}$，边界条件：$v_{j}^{n}=v_{j+J}^{n}$。

#### 数值结果

1. $\textbf{问1.1}：$$\lambda=0.5,J=80$，并且$T=0.1,0.4,0.8,1.0$。

   ![](E:\study_materials\PDEns\04\11.png)

   ​		通过观察误差随终止时间的增长，我们可以看到，误差逐渐增大；三种方法中，$Lax-Wendroff$的误差结果最小，其逼近效果更好一些。

2. $\textbf{问1.2}：$$\lambda=0.5,T=1.0$，并且$J=10,20,40,80,160$，格式为$Lax-Wendroff$。

   <center class="half">    
       <img src="E:\study_materials\PDEns\04\1210.png" width="220" height="170"/>    
       <img src="E:\study_materials\PDEns\04\1220.png" width="220" height="170"/> 
       <img src="E:\study_materials\PDEns\04\1240.png" width="220" height="170"/> 
       <img src="E:\study_materials\PDEns\04\1280.png" width="220" height="170"/> 
       <img src="E:\study_materials\PDEns\04\12160.png" width="220" height="170"/> 
   </center>

   注:从上到下、从左到右分别为$J=10,20,40,80,160$的数值解与真解对比，红色为真解。

   ![12error](E:\study_materials\PDEns\04\12error.png)

   ​		通过观察不同$J$时数值解与真解对比图，以及相应模最大误差随$J$变化图像，明显可以看出，模最大误差随空间离散程度$J$的增大而猪价减小。

3. $\textbf{问2.1}：$

   $T=1.0,J=80,\lambda=0.5$，格式为$CTCS$。

   ![211](E:\study_materials\PDEns\04\211.png)

   $T=1.0,J=80,\lambda=1.5$，格式为$CTCS$。

   <center class="half">    
       <img src="E:\study_materials\PDEns\04\2121.png" width="220" height="170"/>    
       <img src="E:\study_materials\PDEns\04\2122.png" width="220" height="170"/> 
       <img src="E:\study_materials\PDEns\04\2123.png" width="220" height="170"/> 
   </center>

   输出如下：

   ![21](E:\study_materials\PDEns\04\21.png)

   ​		通过对比不同$\lambda$值时数值解与真解图像以及相应误差输出，可以看到，$\lambda=0.5$时数值逼近结果仍可以接受，误差较小；但$\lambda=1.5$时则直接不稳定，误差爆炸！

4. $\textbf{问2.2}：$$T=1.0,\lambda=0.5$，分别取$J=10,20,40,80,160$，格式为$CTCS$。

   <center class="half">    
       <img src="E:\study_materials\PDEns\04\221.png" width="220" height="170"/>    
       <img src="E:\study_materials\PDEns\04\222.png" width="220" height="170"/> 
       <img src="E:\study_materials\PDEns\04\223.png" width="220" height="170"/> 
       <img src="E:\study_materials\PDEns\04\224.png" width="220" height="170"/> 
       <img src="E:\study_materials\PDEns\04\225.png" width="220" height="170"/> 
   </center>

   误差输出结果：

   ![22](E:\study_materials\PDEns\04\22.png)

   ​		通过观察不同$J$时数值解与真解对比图，以及相应模最大误差随$J$变化图像，明显可以看出，模最大误差随空间离散程度$J$的增大而猪价减小。

5. $\textbf{问2.3}：$

   $\lambda=0.5,J=80,T=0.2$，格式为$FTBS$。

   ![231](E:\study_materials\PDEns\04\231.png)

   $\lambda=0.5,J=80,T=0.5$，格式为$FTBS$。

   <center class="half">    
       <img src="E:\study_materials\PDEns\04\2321.png" width="220" height="170"/>    
       <img src="E:\study_materials\PDEns\04\2322.png" width="220" height="170"/> 
       <img src="E:\study_materials\PDEns\04\2323.png" width="220" height="170"/> 
   </center>

误差输出结果如下：

![23](E:\study_materials\PDEns\04\23.png)

​		通过对比不同$T$值时数值解与真解图像以及相应误差输出，可以看到，$T=0.2$时数值逼近结果仍可以接受，误差较小；但$T=0.5$时则直接不稳定，误差爆炸！可见$FTBS$格式对于该方程不稳定，其不满足$CFL$条件。

#### 代码

​		其中数值求解代码与绘图代码参见附件！