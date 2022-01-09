# $Experiment 2：$函数逼近

## 杨乐园  PB18010496

#### 问题描述

​		对函数$f(x)=\frac{\pi - x}{2},(0<x\leq 2\pi)$分别利用如下两个函数，取$N=10,100$进行逼近：
$$
f_{N}(x)=\sum_{w=1}^{N} \frac{\sin wx}{w} \\
\widetilde{f}_{N}(x)=\sum_{w=1}^{N} \frac{\sin \frac{w\pi}{N}}{\frac{w\pi}{N}} \frac{\sin wx}{w}
$$
通过离散变量绘制相应图像并比较，其中$\Delta x=\frac{2\pi}{m},m=20,160$。

#### 数值结果

1. 函数$f_{N}(x)=\sum_{w=1}^{N} \frac{\sin wx}{w}$在$N=10,m=20,160$结果，其中最左侧为原函数图像：

<center class="half">
    <img src="E:\study_materials\PDEns\02\ori.jpg" width="210" height="170"/>
    <img src="E:\study_materials\PDEns\02\1020.jpg" width="210" height="170"/>
    <img src="E:\study_materials\PDEns\02\10160.jpg" width="210" height="170"/>
</center>

​	 函数$f_{N}(x)=\sum_{w=1}^{N} \frac{\sin wx}{w}$在$N=100,m=20,160$结果，其中最左侧为原函数图像：

<center class="half">
    <img src="E:\study_materials\PDEns\02\ori.jpg" width="210" height="170"/>
    <img src="E:\study_materials\PDEns\02\10020.jpg" width="210" height="170"/>
    <img src="E:\study_materials\PDEns\02\100160.jpg" width="210" height="170"/>
</center>

2. 函数$\widetilde{f}_{N}(x)=\sum_{w=1}^{N} \frac{\sin \frac{w\pi}{N}}{\frac{w\pi}{N}} \frac{\sin wx}{w}$在$N=10,m=20,160$结果，其中最左侧为原函数图像：

<center class="half">
    <img src="E:\study_materials\PDEns\02\ori.jpg" width="210" height="170"/>
    <img src="E:\study_materials\PDEns\02\t1020.jpg" width="210" height="170"/>
    <img src="E:\study_materials\PDEns\02\t10160.jpg" width="210" height="170"/>
</center>

 函数$\widetilde{f}_{N}(x)=\sum_{w=1}^{N} \frac{\sin \frac{w\pi}{N}}{\frac{w\pi}{N}} \frac{\sin wx}{w}$在$N=100,m=20,160$结果，其中最左侧为原函数图像：

<center class="half">
    <img src="E:\study_materials\PDEns\02\ori.jpg" width="210" height="170"/>
    <img src="E:\study_materials\PDEns\02\t10020.jpg" width="210" height="170"/>
    <img src="E:\study_materials\PDEns\02\t100160.jpg" width="210" height="170"/>
</center>
3. 差值$f(x)-\widetilde{f}_{N}(x)$在$N=10,m=20,160$结果：

   <center class="half">
       <img src="E:\study_materials\PDEns\02\chazhi1.jpg"/>
       <img src="E:\study_materials\PDEns\02\chazhi2.jpg"/>
   </center>

差值$f(x)-\widetilde{f}_{N}(x)$在$N=100,m=20,160$结果：

<center class="half">
    <img src="E:\study_materials\PDEns\02\chazhi3.jpg"/>
    <img src="E:\study_materials\PDEns\02\chazhi4.jpg"/>
</center>


#### 讨论

​		通过对比两结果可以发现：两函数在逼近结果上，随着求和项数$N$的增大逼近结果越好，即越趋近于原函数；在$N$相同的情况下，函数$\widetilde{f}_{N}(x)=\sum_{w=1}^{N} \frac{\sin \frac{w\pi}{N}}{\frac{w\pi}{N}} \frac{\sin wx}{w}$的逼近结果要优于函数$f_{N}(x)=\sum_{w=1}^{N} \frac{\sin wx}{w}$的逼近结果。