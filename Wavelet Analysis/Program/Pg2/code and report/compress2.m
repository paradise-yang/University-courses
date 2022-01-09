clc;clear;close all;
load julia;whos;
%crit为使用小波包进行分解时所选取的熵函数类型
[thr,sorh,keepapp,crit] = ddencmp('cmp','wp',X);
%返回输入X压缩后的XD。输出参数TREED是XD的最佳小波包分解树；
%PERFL2和PERF0是恢复和压缩L2的能量百分比。PERFL2=100*(X的小波包系数范数/X的小波包系数)^2；
%函数使用由字符串CRIT定义的熵和阈值参数thr*2实现最佳分解。
[Xcomp,treed,PERF0,PERFL2] = wpdencmp(X,sorh,2,'db4',crit,thr*2,keepapp);
n=size(map,1);
colormap(pink(n));
subplot(121);image(wcodemat(X,n));title('origin image','FontSize',20);
axis off
subplot(122);
image(wcodemat(Xcomp,n));title('compress image','FontSize',20);
axis off
plot(treed)