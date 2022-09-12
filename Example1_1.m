clear all;clc;close all
redblue=[0:0.024:0.96,0.96:0.001:1; 0:0.024:0.96,0.96:-0.024:0;1:-0.001:0.96,0.96:-0.024:0 ]';
load('Source_data_for_example_1-1.mat');
% loaded data:
% A  : adjacency matrix
% D  : Degree matrix
% At : Theoretically obtained 'true' A1 matrix
% x  : trajectories of 10 oscillators   [step * node index];

stp = 250;
ts = 0.01;

N = 10;

I =@(x) 0.5*(1-cos(x));



m = 35;

k = 31;


rho = 3;       
lambda = 1e-4;      %  l1 gain

xr = rand(m,N);   % initilize a data storage to store data for m steps
C = xr'*pinv(Psi(xr));  % plift: calculate the Psi(z) observables

xr = zeros(m,N);  % re-initialize the data storages
yr = zeros(m,N);

an = zeros(N,k);

[q,~]=size(Psi(xr));



A0 = rand(q,q);  % initialize the variables for ADMMD
Z = rand(q,q);
W = rand(q,q);

Ar = zeros(q*N, stp);  %  to record the time evolution of entries in A1
Br = zeros(q*N, stp);  %  to record the time evolution of entries in A1p

Nr = 0.1;   % noise

for i=2:stp
    xr = update_data(xr, x(i-1,:));   % update the data of m steps stroed in the storage
    yr = update_data(yr, x(i,  :));
    
    X = Psi(xr+(rand(size(xr))-0.5)*Nr);   % calculate the data matrices using the observable set Psi(z)
    Y = Psi(yr+(rand(size(xr))-0.5)*Nr);
    
    A0 = -0.5*(X*X'+rho/2*eye(q))^-1 *(-2*X*Y' + W-rho*( Z));  % performing ADMM
    Z = sth(A0+1/rho*W,lambda/rho);
    W = W+ rho*(A0-Z);
    
    A1 = C*A0';       % the proposed method
    A1p = C*(Y*pinv(X));  % the pseudo-inverse-based method
    
    an(:,1) =sum(A1(:,1))/N *ones(N,1);
    sar = sum(A1,1);
    for jj = 1:(k-1)/N
        
        an(:,1+N*(jj-1)+(1:N)) = diag(sar(:,1+N*(jj-1)+(1:N) ));
        
    end
    
    bn(:,1) =sum(A1p(:,1))/N *ones(N,1);
    sarp = sum(A1p,1);
    for jj = 1:(k-1)/N
        
        bn(:,1+N*(jj-1)+(1:N)) = diag(sarp(:,1+N*(jj-1)+(1:N) ));
        
    end
    
    
    
    
    Ar(:,i) = reshape(A1,[],1);
    Br(:,i) = reshape(A1p,[],1);

    
end

A2 = getA2(A1);


%%%%

plotrange = 6;
stepsize = 0.1;
bias = 0;
px = (-plotrange:stepsize:plotrange ) + bias;
py = (-plotrange:stepsize:plotrange ) + bias;

ztrue1 = zeros(length(px),length(py));
ztrue2 = zeros(length(px),length(py));
zid1   = zeros(length(px),length(py));
zid2   = zeros(length(px),length(py));





figure(1)
for s = 2:N
    for i = 1:length(px)
        for j = 1:length(py)
            zid1(i,j) = real([1,zeros(1,N-1)] * (A1-an) * Psi([px(i),zeros(1,s-2),py(j),zeros(1,N-s)]));
            
            zid2(i,j) = [1,zeros(1,N-1)] *( A1p-bn) * Psi([px(i),zeros(1,s-2),py(j),zeros(1,N-s)]);
            ztrue1(i,j) =D(1,1)*I(px(i)) - I(py(j));
        end
    end

    subplot(2,round(N/2),s)
    surf(px,py,zid1)
    shading flat
    colormap(redblue);
    %             caxis([-1 1])
    ylabel('$$x_{11}$$','interpreter','latex')
    if s <=9
        xlabel(['$$x_{',num2str(s),'1}$$'],'interpreter','latex')
    else
        xlabel(['$$x_{(',num2str(s),')1}$$'],'interpreter','latex')
    end
    
    view([0 0 1])
    
end





figure(2)

s = 5;
for i = 1:length(px)
    for j = 1:length(py)
        zid1(i,j) = real([1,zeros(1,N-1)] * (A1-an) * Psi([px(i),zeros(1,s-2),py(j),zeros(1,N-s)]));
        
        zid2(i,j) = [1,zeros(1,N-1)] *( A1p-bn) * Psi([px(i),zeros(1,s-2),py(j),zeros(1,N-s)]);
        ztrue1(i,j) =D(1,1)*I(px(i)) - I(py(j));
    end
end

subplot(1,4,1)
surf(px,py,zid1)
shading flat
colormap(redblue);
caxis([-4 4])
ylabel({'$$(a)$$','$$x_{1}$$'},'interpreter','latex')
xlabel(['$$x_{',num2str(s),'}$$'],'interpreter','latex')

view([0 0 1])

subplot(1,4,2)
surf(px,py,zid2)
shading flat
colormap(redblue);
caxis([-4 4])
ylabel({'$$(b)$$','$$x_{1}$$'},'interpreter','latex')
xlabel(['$$x_{',num2str(s),'}$$'],'interpreter','latex')
view([0 0 1])

subplot(1,4,3)
surf(px,py,zid1-ztrue1)
shading flat
colormap(redblue);
caxis([-4 4])
ylabel({'$$(c)$$','$$x_{1}$$'},'interpreter','latex')
xlabel(['$$x_{',num2str(s),'}$$'],'interpreter','latex')

view([0 0 1])
subplot(1,4,4)
surf(px,py,zid2-ztrue1)
shading flat
colormap(redblue);
caxis([-4 4])
ylabel({'$$(d)$$','$$x_{1}$$'},'interpreter','latex')
xlabel(['$$x_{',num2str(s),'}$$'],'interpreter','latex')
view([0 0 1])


disp([ norm(zid1-ztrue1,'f')/norm(ztrue1,'f'),norm(zid2-ztrue1,'f')/norm(ztrue1,'f')])














figure(3)
subplot(2,1,1)
plot(Ar')
ylabel('$$(a)$$','interpreter','latex')
axis([-inf inf -6 5])
subplot(2,1,2)
plot(Br')
ylabel('$$(b)$$','interpreter','latex')
xlabel('step $$[k]$$','interpreter','latex')
axis([-inf inf -6 5])


figure(4)
subplot(2,4,[1,2])
imagesc(A1)
colormap(redblue);colorbar();caxis([-5 5])
ylabel('$$(a)$$','interpreter','latex')

subplot(2,4,[3,4])
imagesc(A2)
colormap(redblue);colorbar();caxis([-5 5])
ylabel('$$(b)$$','interpreter','latex')

subplot(2,4,[6,7])
imagesc((A1-A2))
colormap(redblue);colorbar();caxis([-5 5])
ylabel('$$(c)$$','interpreter','latex')


figure(5)
phi=0:2*pi/N:2*pi/N*(N-1);
nx=cos(phi);
ny=sin(phi);
hold on
scatter(nx,ny);
for i=1:N
    text(1.2*nx(i),1.2*ny(i),num2str(i));
    for j=1:N
        if A(i,j)
            plot([nx(i),nx(j)],[ny(i),ny(j)],'b')  ;
        end
    end
end
axis off


function y=Psi(x)
x = x';
[N,m] = size(x);
y = [ones(1,m);x/2/pi;cos(x);sin(x)];
end


function y = sth(x,s)
y=(abs(x)>s).*(x-sign(x)*s);
end


function y = update_data(X,x)
X(2:end,:,:)=X(1:end-1,:,:);
X(1,:,:)=x;
y = X;
end


function y = getA2(A)
A(:,12:21) = diag(sum(A(:,12:21),1));
A(:,1) = ones(10,1)*sum(A(:,1))/10;
y = A;
end
