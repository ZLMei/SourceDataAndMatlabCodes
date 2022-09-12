clear all;clc;close all
redblue=[0:0.024:0.96,0.96:0.001:1; 0:0.024:0.96,0.96:-0.024:0;1:-0.001:0.96,0.96:-0.024:0 ]';

%%   Data Generation and Identification

stp = 1050;

N = 10;    

A = [0,0,0,0,1,0,0,1,0,1;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,1,1,1,0,0;
    0,0,0,0,0,0,0,1,0,0;
    1,0,0,0,0,1,0,0,1,0;
    0,0,1,0,1,0,0,1,0,1;
    0,1,1,0,0,0,0,1,1,0;
    1,0,1,1,0,1,1,0,1,0;
    0,0,0,0,1,0,1,1,0,0;
    1,0,0,0,0,1,0,0,0,0];

D = diag(sum(A,2));
L = D -A;

I =@(x) 0.5*(1-cos(x));

x = zeros(stp,N);

rl = 400;
range = 2;   %  range of initial values

beta = 1.45*pi;
delta = 0.525;

xp = (rand(1,N) - 0.5) * range; 


for i = 1:stp
    x(i,:) = xp;
    xp = beta * I(xp) + delta+ I(xp)*L' ;
end

kc = 180;                 %  numbers of RBF observables
[~,center] = kmeans(x,kc,'maxIter',1e4);   %  generater centers for observables using k-means clustering

clear x

x = zeros(stp,N);


tst = rand(stp,N);
tnn = zeros(stp,N);

for kk = 1:max(size(tst))
    tnn(kk,:) = beta * I(tst(kk,:)) + delta ;
end
C = tst'*pinv(Psi(tst,center));           % obtain  C  s.t. x = C * Psi(x)

An = C * Psi(tnn,center) *pinv(Psi(tst,center));

[k,~]=size(Psi(tnn,center));    % total number of observables

clear tst tnc tnn


rho = 5;       %  step size
lambda = 0.5;      %  l1 gain
rep = 1;


xr = zeros(rl,N);
yr = zeros(rl,N);


xp = (rand(1,N) - 0.5) * range;


a1 = rand(k,k);
z = rand(k,k);
w = rand(k,k);

a2 = zeros(N,k);

Ar = zeros(k*N, stp);
Br = zeros(k*N, stp);


for cc=1:stp
    
    x(cc,:) = xp;
    
    %%%
    xr(2:rl,:) = xr(1:rl-1,:);
    xr(1,:) = xp;
    
    %%%
    
    
    xp = beta * I(xp) + delta + I(xp)*L';
    
    yr(2:rl,:) = yr(1:rl-1,:);
    yr(1,:) = xp;
    
    
    
    X1 = Psi(xr,center);
    Y1 = Psi(yr,center);

    
    
    a1 = -0.5*(X1*X1'+rho/2*eye(k))^-1 *(-2*X1*Y1' + w-rho*( z));
    z = sth(a1+1/rho*w,lambda/rho);
    w = w+ rho*(a1-z);

    
    
    
    
    ar = C*a1';
    
    
    a2(:,1) =sum(ar(:,1))/N *ones(N,1);
    sar = sum(ar,1);
    for jj = 1:(k-1)/N
        
        a2(:,1+N*(jj-1)+(1:N)) = diag(sar(:,1+N*(jj-1)+(1:N) ));
        
    end
    
    
    
    Ar(:,cc) = reshape(ar,[],1);
    
    
    
    if mod(cc,150)==0
        disp([num2str(cc/stp*100),'%']);
    end
end



%%  Plots and Verification


plotrange = 4;
stepsize = 0.1;
bias = 2;
px = (-plotrange:stepsize:plotrange ) + bias;
py = (-plotrange:stepsize:plotrange ) + bias;



ztrue1 = zeros(length(px),length(py));
zid1   = zeros(length(px),length(py));



figure(1)
for s = 2:N
    for i = 1:length(px)
        for j = 1:length(py)
            zid1(i,j) = real([1,zeros(1,N-1)] * (ar-a2) * Psi([px(i),zeros(1,s-2),py(j),zeros(1,N-s)],center));
            ztrue1(i,j) =D(1,1)*I(px(i)) - I(py(j));
        end
    end
    
    
    subplot(2,round(N/2),s)
    surf(px,py,zid1)
    shading flat
    colormap(redblue);
    ylabel('$$x_{1}$$','interpreter','latex')
    xlabel(['$$x_{',num2str(s),'}$$'],'interpreter','latex')
    view([0 0 1])
end



s = 5;
N = 10;

plotrange = 6;
stepsize = 0.1;
bias = 2;
px = (-plotrange:stepsize:plotrange ) + bias;
py = (-plotrange:stepsize:plotrange ) + bias;



ztrue1 = zeros(length(px),length(py));
ztrue2 = zeros(length(px),length(py));
zid1   = zeros(length(px),length(py));
zid2   = zeros(length(px),length(py));


for i = 1:length(px)
    for j = 1:length(py)
        zid1(i,j) = real([1,zeros(1,N-1)] * (ar-a2) * Psi([px(i),zeros(1,s-2),py(j),zeros(1,N-s)],center));
        ztrue1(i,j) =D(1,1)*I(px(i)) - I(py(j));
    end
end



figure(2)

surf(px,py,zid1-ztrue1)
shading flat
colormap(redblue);
caxis([-1 1])
ylabel('$$x_{1}$$','interpreter','latex')
xlabel(['$$x_{',num2str(s),'}$$'],'interpreter','latex')
axis([min(px),max(px),min(py),max(py),-inf inf])
view([0 0 1])
colorbar()
hold on
scatter3(x(:,1),x(:,5),10+x(:,1)*0,'kx')


%%


function y=Psi(x,~)
x = x';
[N,m] = size(x);

T = [0.1*pi:0.1*pi:3*pi];
y = zeros(length(T)*2,m);
for i=1:length(T)
    y((i-1)*2*N+(1:N),:) = sqevn(x,T(i));
    y((i-1)*2*N+N+(1:N),:) = sqodd(x,T(i));
end
y = [ones(1,m);x;y];
end

function y = sth(x,s)
y=(abs(x)>s).*(x-sign(x)*s);
end

function y = sqevn(x,T)
y = sqodd(x-T/4,T);
end

function y = sqodd(x,T)   
plus  = (mod(x,T) < T/2) + 0;
minus = (mod(x,T) >= T/2) + 0;
y = plus-minus;
end


