% bad1d.m
for N=1:50
    [Ah,Bh,Ch,Dh,zh,wh] = SEBhat(N); % Standard GLL points
    [Au,Bu,Cu,Du,zu,wu] = unihat(N); % Uniform points
    No = N+2; % Overintgration
    [Ao,Bo,Co,Do,zo,wo] = SEBhat(No); 
    J = interp_mat(zo,zu);
    Bu = J'*Bo*J; Au = Du'*Bu*Du; % Full mass and stiffness matrices

    n=size(Ah,1);
    Ah=Ah(1:n-1,1:n-1); % Dirichlet at x=1; Neumann at x=-1;
    Au=Au(1:n-1,1:n-1);
    NN(N)=N; 
    cu(N)=cond(Au); 
    cg(N)=cond(Ah);
end

nmax=25; c=zeros(nmax,1); N=zeros(nmax,1);
for n=1:25
    A=zeros(n,n);
    for i=1:n
        for j=1:n
            A(i,j) = (i*j)/(i+j-1); % Stiffness matrix for phi_i = x^i on (0,1]
        end
    end
    N(n) = n; c(n) = cond(A);
end
semilogy(NN,cu,'r+',NN,cu,'r-',NN,cg,'bo',NN,cg,'b-',N,c,'kx',N,c,'k-')
axis([0 50 1 1.e15]);
title('Condition of 1D Laplacian, u(0)=u''(1)=0');
xlabel('Polynomial Order N'); ylabel('Condition Number');
print -deps ’bad1d.ps’

[Ah,Bh,Ch,Dh,z,w]=semhat(N);

for N=1:50
    [Ah,Bh,Ch,Dh,z,rho] = semhat(N);
    a=0;b=1; x=a+0.5*(b-a)*(z+1); Lx=(b-a);
    p=x; % This is p(x)
    Ab = (2./Lx)*Dh'*diag(rho.*p)*Dh;
    Bb = 0.5*Lx*Bh;
    P = eye(N+1); P=P(:,1:N); % Prolongation matrix
    R = P';
    
    A=R*Ab*P;
    
    f=x;
    rhs = R*Bb*f;
    
    u = P*(A\rhs); % Extend numerical solution by 0 at Dirichlet boundaries.
    ue = 0.25*(1-x.*x);
    
    plot(x,u,'r.-',x,ue,'g.-')
    pause(1);
    
    eN(N)=max(abs(u-ue));
    nN(N)=N;
end