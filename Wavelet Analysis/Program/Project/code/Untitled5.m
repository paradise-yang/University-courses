hx=0.01;
ht=0.000051;
x=zeros(1,100);
s=zeros(1,50);
d=zeros(1,50);
for i=1:100
    x(i)=(i-1)*hx;
end
u=zeros(1,100);
f=zeros(1,100);
h=dbwavf('db3');
r=zeros(50);
c=zeros(50);
b=zeros(50);
a=zeros(50);
g=[0.0249,0.0604,-0.0955,-0.3252,0.5706,-0.2352];
g1=0.0249+0.0604-0.0955-0.3252+0.5706-0.2352;
h1=0.2353+0.5706+0.3252-0.0955-0.0604+0.0249;
r0=-5900/1123;r1=11792/3369;r2=-2944/3369;r3=128/1123;r4=6/1123;
c0=(h(1)*(-r1)+h(2)*r0+h(3)*r1+h(4)*r2+h(5)*r3+h(6)*r4)*g1;
c1=(h(1)*(r1)+h(2)*r2+h(3)*r3+h(4)*r4)*g1;
c2=(h(1)*(r3)+h(2)*r4)*g1;

a0=(g(1)*(-r1)+g(2)*r0+g(3)*r1+g(4)*r2+g(5)*r3+g(6)*r4)*g1;
a1=(g(1)*(r1)+g(2)*r2+g(3)*r3+g(4)*r4)*g1;
a2=(g(1)*(r3)+g(2)*r4)*g1;

b0=(g(1)*(-r1)+g(2)*r0+g(3)*r1+g(4)*r2+g(5)*r3+g(6)*r4)*h1;
b1=(g(1)*(r1)+g(2)*r2+g(3)*r3+g(4)*r4)*h1;
b2=(g(1)*(r3)+g(2)*r4)*h1;
for i=1:50
    for j=1:50
        t=i-j;
        if t==0
            r(i,j)=-5900/1123;
            c(i,j)=-2.6061;
            a(i,j)=a0;
            b(i,j)=b0;
        end
        if t==1
            r(i,j)=11792/3369;
            c(i,j)=0.3612;
            a(i,j)=a1;
            b(i,j)=b1;
        end
        if t==2
            r(i,j)=-2944/3369;
            c(i,j)=0.0299;
            a(i,j)=a2;
            b(i,j)=b2;
        end
        if t==3
            r(i,j)=128/1123;
        end
        if t==4
            r(i,j)=6/1123;
        end
        if t==-1
            r(i,j)=-11792/3369;
            c(i,j)=-0.3612;
            a(i,j)=-a1;
            b(i,j)=-b1;
        end
        if t==-2
            r(i,j)=2944/3369;
            c(i,j)=-0.0299;
            a(i,j)=-a2;
            b(i,j)=-b2;
        end
        if t==-3
            r(i,j)=-128/1123;
        end
        if t==-4
            r(i,j)=-6/1123;
        end
    end
end

for i=1:100
    f(i)=2*(hx*(i-1)+10);
    u(i)=hx*(i-1)+10;
    if hx*(i-1)>=0.5
        f(i)=2*(hx*(i-1)-10);
        u(i)=hx*(i-1)-10;
    end
end
for i=1:50
    for j=1:6
        if j+2*i-2>100||j+2*i-2<1
            u(j+2*i-2)=0;
        end
        s(i)=s(i)+h(j)*u(j+2*i-2);
        d(i)=d(i)+g(j)*u(j+2*i-2);
    end
end

m=zeros(1,50);
n=zeros(1,50);
for k=1:500
    m=0.125*(s*r+d*c);
    n=0.125*(s*b+d*a);
    for i=1:50
        u(i)=u(i)+ht*(f(i)+m(i)-u(i));
        u(i+50)=u(i+50)+ht*(f(i+50)+n(i)-u(i+50));
    end
    for i=1:50
        s(i)=u(i);
        d(i)=u(i+50);
    end
    for j=1:100
        f(j)=(hx*(j-1)+10)*(ht*k+2);
        if hx*(j-1)>=0.5
            f(j)=(hx*(j-1)-10)*(ht*k+2);
        end
    end
end

uz=zeros(1,100);
for i=1:50
    uz(i)=((i-1)*0.01+10)*(1+k*ht);
    %uz(i)=uz(i)-u(i);
end
for i=51:100
    uz(i)=((i-1)*0.01+10)*(1+k*ht);
    %uz(i)=uz(i)-u(i);
end

plot(x,uz,x,u(1:100),'r')
%plot(x,uz(1:100))
xlabel('t=0.005时的算法')
