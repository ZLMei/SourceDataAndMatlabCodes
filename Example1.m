clear all;clc;close all
redblue=[0:0.024:0.96,0.96:0.001:1; 0:0.024:0.96,0.96:-0.024:0;1:-0.001:0.96,0.96:-0.024:0 ]';
load('Source_data_for_example_1.mat');

stp = 1000;
ts = 0.01;

N = 10;

I =@(x) 0.5*(1-cos(x));



rl = 35;


rho = 3;       %  step size
lambda = 1e-4;      %  l1 gain

xr = rand(rl,N);   % initilize a data storage to store data from a m-step-long horizon
C = xr'*pinv(plift(xr));  % plift: calculate the Psi(z) observables

xr = zeros(rl,N);  % re-initialize the data storages
yr = zeros(rl,N);

[q,~]=size(plift(xr));

er = zeros(stp,2);

A0 = rand(q,q);  % initialize the variables for ADMMD
Z = rand(q,q);
W = rand(q,q);

Ar = zeros(q*N, stp);  %  to record the time evolution of entries in A1
Br = zeros(q*N, stp);  %  to record the time evolution of entries in A1p

Nr = 0.1;   % noise

for i=2:stp
    xr = update_data(xr, x(i-1,:));   % update the data of m steps stroed in the storage
    yr = update_data(yr, x(i,  :));
    
    X = plift(xr+(rand(size(xr))-0.5)*Nr);   % calculate the data matrices using the observable set Psi(z)
    Y = plift(yr+(rand(size(xr))-0.5)*Nr);

        A0 = -0.5*(X*X'+rho/2*eye(q))^-1 *(-2*X*Y' + W-rho*( Z));  % performing ADMM
        Z = sth(A0+1/rho*W,lambda/rho);
        W = W+ rho*(A0-Z);

    A1 = C*A0';       % the proposed method
    A1p = C*(Y*pinv(X));  % the pseudo-inverse-based method

    Ar(:,i) = reshape(A1,[],1);
    Br(:,i) = reshape(A1p,[],1);
    
    er(i,:) = [norm(A1-At,'f')/norm(At,'f')*100,norm(A1p-At,'f')/norm(At,'f')*100];  % calculate and record the error indexes \mathcal{I}(x)
    
end

A2 = getA2(A1);

figure(2)
subplot(2,1,1)
plot(Ar')
ylabel('$$(a)$$','interpreter','latex')
axis([-inf inf -6 5])
subplot(2,1,2)
plot(Br')
ylabel('$$(b)$$','interpreter','latex')
xlabel('step $$[k]$$','interpreter','latex')
axis([-inf inf -6 5])

figure(3)
se1=semilogy(er);
hold on
se2=semilogy(5+0*(1:1:stp),'k--');
% axis([-inf inf -5 150])
ylabel('error(\%) $$\mathcal{I}[k]$$', 'interpreter','latex' )
legend([se1(1),se1(2),se2],'proposed method','pseudo-inverse-based method','5% error')
grid on
xlabel('step $$[k]$$','interpreter','latex')
axis([-inf inf 1e0 1e2])

figure(1)
subplot(2,4,[1,2])
imagesc(A1)
colormap(redblue);colorbar();caxis([-5 5])
ylabel('$$(a)$$','interpreter','latex')

subplot(2,4,[3,4])
imagesc(A2)
colormap(redblue);colorbar();caxis([-5 5])
ylabel('$$(b)$$','interpreter','latex')

subplot(2,4,[5,6])
imagesc((A1-A2))
colormap(redblue);colorbar();caxis([-5 5])
ylabel('$$(c)$$','interpreter','latex')

subplot(2,4,7)

phi=0:2*pi/N:2*pi/N*(N-1);
nx=cos(phi); 
ny=sin(phi);
hold on
scatter(nx,ny);
for i=1:N
    text(1.2*nx(i),1.2*ny(i),num2str(i));
    for j=1:N
        if A(i,j)
%             quiver(nx(j),ny(j),nx(i)-nx(j),ny(i)-ny(j),0,'k')           %  black : large index to small index
            plot([nx(i),nx(j)],[ny(i),ny(j)],'b')  ;                    % undirected
        end
    end
end
tx=text(-1.2,0,{'$$(d)$$'},'interpreter','latex','Rotation',90);
% set(tx','Rotation','90')
axis off
% ylabel('$$(d)$$','interpreter','latex')
% ylabel on

subplot(2,4,8)
imagesc(L)
colormap(redblue);colorbar();caxis([-5 5])
ylabel('$$(e)$$','interpreter','latex')



% figure(3)
% hold on
% scatter(reshape(At,[],1),reshape(ar,[],1),'rx')
% plot(min(min(At))-0.5:max(max(At))+0.5,min(min(At))-0.5:max(max(At))+0.5,'k')


function y=plift(x)
x = x';
[N,m] = size(x);

% y = [ones(1,m);x;cos(x);sin(x);x.^2;x.^3;exp(x)];
y = [ones(1,m);x/2/pi;cos(x);sin(x)];

% X=zeros(N^2,m);
% for i =1:m
%     X(:,i) = reshape(cos(x(:,i))*sin(x(:,i))',[],1);
% end
% y = [ones(1,m);x;X];
% y = [ones(1,m);x];
end



function A = graph(N,s)
A = rand(N,N);
A = (A+A')/2 > s;
A = A.*(ones(N,N)-eye(N));
end

function y = sth(x,s)
y=(abs(x)>s).*(x-sign(x)*s);
end



function y = f(x,A,w)
N = length(x);

u=zeros(1,N);
for i =1:N
    for j=1:N
        u(1,i) = u(1,i) + A(i,j)*sin(x(j)-x(i));
        %         u(1,i) = u(1,i) + A(i,j)*(x(j)-x(i));
    end
end
y = w + u;
end

function y = update_data(X,x)  
X(2:end,:,:)=X(1:end-1,:,:);
X(1,:,:)=x;
y = X;
end

function y = rk(x,A,w,ts)
k1 = f(x,A,w);
k2 = f(x+0.5*ts*k1,A,w);
k3 = f(x+0.5*ts*k2,A,w);
k4 = f(x+ts*k3,A,w);

y = x + 1/6*ts*(k1+2*k2+2*k3+k4);


end

function y = getA2(A)
A(:,12:21) = diag(sum(A(:,12:21),1));
A(:,1) = ones(10,1)*sum(A(:,1))/10;
y = A;
end
