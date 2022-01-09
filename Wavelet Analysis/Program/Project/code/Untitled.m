tic;
% this method is transform from Galerkin method 
%also call it as finit method
%is used for solving two point BVP which is the first and second term.
%this code was writen by HU.D.dong in February 11th 2017
%MATLAB 7.0
clear;
clc;
N=50;
h=1/N;
X=0:h:1;
f=inline('(0.5*pi^2)*sin(0.5*pi.*x)');
%以下是右端向量:
for i=2:N
    fun1=@(x) pi^2/2.*sin(pi/2.*x).*(1-(x-X(i))/h);
    fun2=@(x) pi^2/2.*sin(pi/2.*x).*((x-X(i-1))/h);
    f_phi(i-1,1)=quad(fun1,X(i),X(i+1))+quad(fun2,X(i-1),X(i));
end
funN=@(x) pi^2/2.*sin(pi/2.*x).*(x-X(N))/h;
f_phi(N)=quad(funN,X(N),X(N+1));
%以下是刚度矩阵：
A11=quad(@(x) 2/h+0.25*pi^2*h.*(1-2*x+2*x.^2),0,1);
A12=quad(@(x) -1/h+0.25*pi^2*h.*(1-x).*x,0,1);
ANN=quad(@(x) 1/h+0.25*pi^2*h*x.^2,0,1);
A=diag([A11*ones(1,N-1),ANN],0)+diag(A12*ones(1,N-1),1)+diag(A12*ones(1,N-1),-1);
Numerical_solution=A\f_phi;
Numerical_solution=[0;Numerical_solution];
%Accurate solution on above以下是精确解
%%
for i=1:length(X)
    Accurate_solution(i,1)=sin((pi*X(i))/2)/2 - cos((pi*X(i))/2)/2 + exp((pi*X(i))/2)*((exp(-(pi*X(i))/2)*cos((pi*X(i))/2))/2 + (exp(-(pi*X(i))/2)*sin((pi*X(i))/2))/2);
end 
figure(1);
grid on; 
subplot(1,2,1);
plot(X,Numerical_solution,'ro-',X,Accurate_solution,'b^:');
title('Numerical solutions vs Accurate solutions');
legend('Numerical_solution','Accurate_solution');
subplot(1,2,2);
plot(X,Numerical_solution-Accurate_solution,'b x');
legend('error_solution');
title('error');
toc;