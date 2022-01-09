# $Experiment 7：$初-边值近似问题

## 杨乐园  PB18010496

#### 问题描述

$\mathbf{HW\ 1.5.9\ and\ 1.5.10}$

​		针对下述初边值问题：
$$
\left \{
\begin{array}{ll}
v_t=\nu v_{xx}+F(x,t) & x\in(0,1),t>0  \\
v(x,0)=f(x) & x\in [0,1] \\
v_x(0,t)=a(t) & t\geq0 \\
v(1,t)=b(t) & t\geq0
\end{array}\right.
$$
分别使用一阶与二阶的$Neumann$边界近似方法，求解上述初-边值问题。其中
$$
\nu =0.1,\quad f(x)=x(1-x),\quad a(t)=10\sin t,\quad b(t)=4\sin 6t,\quad F(x,t)=\sin 2\pi x\sin 4\pi t
$$
并且分别在数值差分条件$M=10,\Delta t=0.05$、$M=20,\Delta t=0.01$、$M=40,\Delta t=0.002$下，求取时间$T=0.1,0.9,2.0$处数值解。

$\mathbf{HW\ 2.2.2}$

​		针对下述初边值问题：
$$
\left \{
\begin{array}{ll}
v_t=\nu v_{xx} & x\in(0,1),t>0  \\
v(x,0)=\sin 4\pi x & x\in [0,1] \\
v(0,t)=v(1,t)=0 & t\geq0
\end{array}\right.
$$
利用$FTCS$格式，其中$\nu=0.1$，分别在条件$(i)\Delta x=0.1,\Delta t=0.05$、$(ii)\Delta x=0.05,\Delta t=0.0125$、$(ii)\Delta x=0.01,\Delta t=0.0005$，求取时间$T=0.05,0.1$处的数值解。

#### 数值方法

$\mathbf{HW\ 1.5.9\ and\ 1.5.10}$

​		一阶近似格式如下：
$$
\left \{
\begin{array}{ll}
u_k^{n+1}=u_{k}^{n}+\nu\frac{\Delta t}{\Delta x^2}(u_{k+1}^{n}-2u_{k}^{n}+u_{k-1}^{n})+\Delta tF_k^n & k=1,...,M-1\\
u_{k}^{0}=f_k & k=0,1,...,M \\
u_0^n=u_1^n-\Delta x a^n & n=0,1,...,T/\Delta t\\
u_M^n=b_M & n=0,1,...,T/\Delta t
\end{array}
\right.
$$
​		二阶近似格式如下：
$$
\left \{
\begin{array}{ll}
u_k^{n+1}=u_{k}^{n}+\nu\frac{\Delta t}{\Delta x^2}(u_{k+1}^{n}-2u_{k}^{n}+u_{k-1}^{n})+\Delta tF_k^n & k=1,...,M-1\\
u_{k}^{0}=f_k & k=0,1,...,M \\
u_0^{n+1}=u_0^n-2\nu\frac{\Delta t}{\Delta x^2}(u_{1}^{n}-u_{0}^{n})+-2\nu\frac{\Delta t}{\Delta x} a^n & n=0,1,...,T/\Delta t\\
u_M^n=b_M & n=0,1,...,T/\Delta t
\end{array}
\right.
$$
$\mathbf{HW\ 2.2.2}$

​		$FTCS$格式如下：
$$
\left \{
\begin{array}{ll}
u_k^{n+1}=u_{k}^{n}+\nu\frac{\Delta t}{\Delta x^2}(u_{k+1}^{n}-2u_{k}^{n}+u_{k-1}^{n}) & k=1,...,M-1,\ \ n=1,...,T/\Delta t \\
u_{k}^{0}= \sin 4 \pi k\Delta x & k=1,...,M \\
u_{0}^{n}=u_M^n=0 & n=1,...,T/\Delta t
\end{array}
\right.
$$

#### 数值结果

$\mathbf{HW\ 1.5.9\ and\ 1.5.10}$

从上到下依次为$T=0.1,0.9,2.0$。

$\bullet M=10,\Delta t=0.05$。从左到右依次为一阶近似、二阶近似、对比、差的绝对值。

<center class="half">    
    <img src="E:\study_materials\PDEns\07\picture\111.png" width="160" height="125"/>    
    <img src="E:\study_materials\PDEns\07\picture\112.png" width="160" height="125"/> 
    <img src="E:\study_materials\PDEns\07\picture\11all.png" width="160" height="125"/> 
    <img src="E:\study_materials\PDEns\07\picture\11sub.png" width="160" height="125"/> 
</center>

<center class="half">    
    <img src="E:\study_materials\PDEns\07\picture\121.png" width="160" height="125"/>    
    <img src="E:\study_materials\PDEns\07\picture\122.png" width="160" height="125"/> 
    <img src="E:\study_materials\PDEns\07\picture\12all.png" width="160" height="125"/> 
    <img src="E:\study_materials\PDEns\07\picture\12sub.png" width="160" height="125"/> 
</center>

<center class="half">    
    <img src="E:\study_materials\PDEns\07\picture\131.png" width="160" height="125"/>    
    <img src="E:\study_materials\PDEns\07\picture\132.png" width="160" height="125"/> 
    <img src="E:\study_materials\PDEns\07\picture\13all.png" width="160" height="125"/> 
    <img src="E:\study_materials\PDEns\07\picture\13sub.png" width="160" height="125"/> 
</center>

$\bullet M=20,\Delta t=0.01$。从左到右依次为一阶近似、二阶近似、对比、差的绝对值。

<center class="half">    
    <img src="E:\study_materials\PDEns\07\picture\2111.png" width="160" height="125"/>    
    <img src="E:\study_materials\PDEns\07\picture\2112.png" width="160" height="125"/> 
    <img src="E:\study_materials\PDEns\07\picture\211all.png" width="160" height="125"/> 
    <img src="E:\study_materials\PDEns\07\picture\211sub.png" width="160" height="125"/> 
</center>

<center class="half">    
    <img src="E:\study_materials\PDEns\07\picture\2121.png" width="160" height="125"/>    
    <img src="E:\study_materials\PDEns\07\picture\2122.png" width="160" height="125"/> 
    <img src="E:\study_materials\PDEns\07\picture\212all.png" width="160" height="125"/> 
    <img src="E:\study_materials\PDEns\07\picture\212sub.png" width="160" height="125"/> 
</center>

<center class="half">    
    <img src="E:\study_materials\PDEns\07\picture\2131.png" width="160" height="125"/>    
    <img src="E:\study_materials\PDEns\07\picture\2132.png" width="160" height="125"/> 
    <img src="E:\study_materials\PDEns\07\picture\213all.png" width="160" height="125"/> 
    <img src="E:\study_materials\PDEns\07\picture\213sub.png" width="160" height="125"/> 
</center>

$\bullet M=40,\Delta t=0.02$。从左到右依次为一阶近似、二阶近似、对比、差的绝对值。

<center class="half">    
    <img src="E:\study_materials\PDEns\07\picture\2211.png" width="160" height="125"/>    
    <img src="E:\study_materials\PDEns\07\picture\2212.png" width="160" height="125"/> 
    <img src="E:\study_materials\PDEns\07\picture\221all.png" width="160" height="125"/> 
    <img src="E:\study_materials\PDEns\07\picture\221sub.png" width="160" height="125"/> 
</center>

<center class="half">    
    <img src="E:\study_materials\PDEns\07\picture\2221.png" width="160" height="125"/>    
    <img src="E:\study_materials\PDEns\07\picture\2222.png" width="160" height="125"/> 
    <img src="E:\study_materials\PDEns\07\picture\222all.png" width="160" height="125"/> 
    <img src="E:\study_materials\PDEns\07\picture\222sub.png" width="160" height="125"/> 
</center>

<center class="half">    
    <img src="E:\study_materials\PDEns\07\picture\2231.png" width="160" height="125"/>    
    <img src="E:\study_materials\PDEns\07\picture\2232.png" width="160" height="125"/> 
    <img src="E:\study_materials\PDEns\07\picture\223all.png" width="160" height="125"/> 
    <img src="E:\study_materials\PDEns\07\picture\223sub.png" width="160" height="125"/> 
</center>
$\mathbf{HW\ 2.2.2}$

$\bullet \Delta x=0.1,\Delta t=0.05$。从左到右依次$T=0.05,0.1$。其中棕色曲线为真实解。

<center class="half">    
    <img src="E:\study_materials\PDEns\07\picture\311.png" width="200" height="150"/>    
    <img src="E:\study_materials\PDEns\07\picture\312.png" width="200" height="150"/>  
</center>

$\bullet \Delta x=0.05,\Delta t=0.0125$。从左到右依次$T=0.05,0.1$。其中棕色曲线为真实解。

<center class="half">    
    <img src="E:\study_materials\PDEns\07\picture\321.png" width="200" height="150"/>    
    <img src="E:\study_materials\PDEns\07\picture\322.png" width="200" height="150"/>  
</center>

$\bullet \Delta x=0.01,\Delta t=0.0005$。从左到右依次$T=0.05,0.1$。其中棕色曲线为真实解。

<center class="half">    
    <img src="E:\study_materials\PDEns\07\picture\331.png" width="200" height="150"/>    
    <img src="E:\study_materials\PDEns\07\picture\332.png" width="200" height="150"/>  
</center>

#### 代码

​		其中数值求解代码与绘图代码参见附件！
