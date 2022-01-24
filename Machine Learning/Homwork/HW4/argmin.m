function [w] = argmin(z,L,lambda)
w=z;
[~,column]=size(w);
temp=lambda/L;
for i=1:column
    if z(i)<-temp
        w(i)=z(i)+temp;
    elseif z(i)>temp
        w(i)=z(i)-temp;
    else
        w(i)=0;
    end
end