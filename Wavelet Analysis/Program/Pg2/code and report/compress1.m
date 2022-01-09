%利用Haar小波对图像压缩
load wbarb;whos;
%对图像X用Daubechies1阶消失矩小波基函数实现2层分解
[C, S] = wavedec2(X,2,'db1');
%利用ddencmp()函数自动生成小波压缩的阈值选取方案，用小波分解
%thr为阈值，sorth为选择阈值的方式（sorth=s为软阈值，=h为硬阈值），keepapp=0/1决定是否对近似分量进行阈值处理
[thr,sorh,keepapp] = ddencmp('cmp','wv',X);
%对每层采取同一个阈值进行处理
%XC是压缩后的信号,[CXC,LXC]是XC的小波分解结构,PERF0和PERFL2是恢复和压缩L^2的范数百分比, 是用百分制表明降噪或压缩所保留的能量成分
[Xcomp, CXC, LXC, PERF0, PERFL2] = wdencmp('gbl',C,S,'db1',2,thr,sorh,keepapp);
colormap(map);
subplot(121);image(X);title('origin image','FontSize',20);
axis square
axis off
subplot(122);image(Xcomp);title('compress image','FontSize',20);
axis square
axis off