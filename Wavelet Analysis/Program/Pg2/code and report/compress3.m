clc;clear;close all;load porche
colormap(pink(255))
subplot(2,2,1); image(X);
axis square,axis off;
title('origin image','FontSize',20)
%返回压缩比CR与Bit-Per-Pixel ratio，采用EZW方法，Haar小波
[CR,BPP] = wcompress('c',X,'mask.wtc','ezw','maxloop',12,'wname','haar');
Xc = wcompress('u','mask.wtc');
subplot(2,2,2);image(Xc);
axis square,axis off;
title({['EZW - Haar'],['ratio: ' num2str(CR,'%1.2f %%'),',BPP: ' num2str(BPP,'%3.2f')]},'FontSize',20)

%采用EZW方法，bior4.4小波
[CR,BPP] = wcompress('c',X,'mask.wtc','ezw','maxloop',12,'wname','bior4.4');
Xc = wcompress('u','mask.wtc');
colormap(pink(255))
subplot(223); image(Xc);
axis square,axis off;
title({['EZW - Bio4.4'],['ratio:' num2str(CR,'%1.2f %%'),',BPP: ' num2str(BPP,'%3.2f')]},'FontSize',20)

%采用SPIHT方法，bior4.4小波
[CR,BPP] = wcompress('c',X,'mask.wtc','spiht','maxloop',12,'wname','bior4.4');
Xc = wcompress('u','mask.wtc');
colormap(pink(255))
subplot(224); image(Xc);
axis square,axis off;
title({['SPINT - Bio4.4'],['ratio: ' num2str(CR,'%1.2f %%'),',BPP: ' num2str(BPP,'%3.2f')]},'FontSize',20)