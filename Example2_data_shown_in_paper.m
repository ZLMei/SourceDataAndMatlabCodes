clear all;clc;close all
redblue=[0:0.024:0.96,0.96:0.001:1; 0:0.024:0.96,0.96:-0.024:0;1:-0.001:0.96,0.96:-0.024:0 ]';

load('Source_data_for_example_2.mat')
% loaded data:
%  Arecord :   the adjacency matrices
%  x       :   the trajectories of the nodes


time = 120;
ts = 0.01;
tc = round((30:30:150)/ts);
stp = time/ts;

N = 10;

spars = 0.6;


A = Arecord(:,:,1);
D = diag(sum(A,2));
L = D - A;

sumD=sum(sum(D));



m =450;    % recording length (amount of data)
range = 20;   % initial values

rho =1e-3;       %  step size
lambda = 1e-7;      %  l1 gain


xr = rand(m,N,3);
X = Psi(xr(:,:,1),xr(:,:,2),xr(:,:,3));
[k,~]=size(X);
C = [xr(:,:,1)';xr(:,:,2)';xr(:,:,3)']*pinv(X);
clear xr X

a2 = zeros(3*N,k);

xr = zeros(m,N,3);
yr = zeros(m,N,3);

con = @(A) (A-eye(length(A)))/ts;


xpc = xint;
xp0 = xpc;

rankX = zeros(stp,1);
rankY = zeros(stp,1);




a1 = zeros(k,k);
z = zeros(k,k);
w = zeros(k,k);



Ar = zeros( stp, k*N*3);
Br = zeros( stp, k*N*3);
er = zeros(stp,1);


cter = 2;
sigma = 0.5;

Noise = 0.000;
for i =1:stp
    xpc = x(i,:,:);
    xr = update_data(xr,xpc);
    
    new_x = rk(xpc(:,:,1),xpc(:,:,2),xpc(:,:,3),sigma*A,ts);
    
    xpc(1,:,1) = new_x(1,:);
    xpc(1,:,2) = new_x(2,:);
    xpc(1,:,3) = new_x(3,:);
    
    
    yr = update_data(yr, xpc);
    
    
    
    X = Psi(xr(:,:,1),xr(:,:,2),xr(:,:,3));
    Y = Psi(yr(:,:,1),yr(:,:,2),yr(:,:,3));
    
    rankX(i) = sum(eig(X*X')>1e-2);
    
    
    a1 = -0.5*(X*X'+rho/2*eye(k))^-1 *(-2*X*Y' + w-rho*z);
    z = sth(a1+1/rho*w,lambda/rho);
    w = w+ rho*(a1-z);
    
    
    ar = C*(z)';
    br = C*(Y*pinv(X));
    Ar(i,:) = reshape(ar-a2,1,[]);
    
    
    sar1 = sum(ar(1:N,:),1);
    sar2 = sum(ar(N+1:2*N,:),1);
    sar3 = sum(ar(2*N+1:3*N,:),1);
    sar = [sar1;sar2;sar3];
    for ss = 1:3
        
        a2(N*(ss-1) + (1:N),1) =sum(ar(N*(ss-1) + (1:N),1))/N *ones(N,1);
        for jj = 1:(k-1)/N
            
            
            a2(N*(ss-1) + (1:N),1+N*(jj-1)+(1:N)) = diag(sar(ss,1+N*(jj-1)+(1:N) ));
        end
        
    end
    
    
    if max(i==tc)
        
        A = Arecord(:,:,cter);
        cter =cter+1;
        
        
        disp([num2str(i/stp*100),'%'])
        
    end
    
    
end


figure(1)

subplot(3,1,1)
plot(ts:ts:time, x(:,:,2));
ylabel('$$x_{i,2}$$','interpreter','latex')


subplot(3,1,2)
plot(ts:ts:time,rankX)
ylabel('$$card\{eig(XX^*)>10^{-2}\}$$','interpreter','latex')


subplot(3,1,3)
hold on
c13 = Ar(:,401);
c31 = Ar(:,343);

gr=plot(ts:ts:time,Ar','g');
rd=plot(ts:ts:time,c13,'--r','linewidth',2);
bl=plot(ts:ts:time,c31,'--b','linewidth',2);
ylabel('entries in $$A_1-A_2$$','interpreter','latex')
axis([-inf inf -0.08 0.08])
text(2,-0.05,'$$(c)$$','interpreter','latex')

xlabel('time [s]')


figure(2)

subplot(2,4,[1,2])
imagesc(C*a1')
ylabel('$$(a)$$','interpreter','latex')
colormap(redblue);colorbar();caxis([-0.03 0.03])


subplot(2,4,[3,4])
imagesc(a2)
ylabel('$$(b)$$','interpreter','latex')
colormap(redblue);colorbar();caxis([-0.03 0.03])

subplot(2,4,[6,7])
imagesc(C*a1'-a2)
ylabel('$$(c)$$','interpreter','latex')
colormap(redblue);colorbar();caxis([-0.03 0.03])






function y=Psi(x1,x2,x3)
[m,~]=size(x1);
x = [x1';x2';x3'];
y = [ones(1,m);x;x1'.*x2';x1'.*x3';x2'.*x3'];

end


function y = sth(x,s)
y=(abs(x)>s).*(x-sign(x)*s);
end


function dd=f(xp, yp, zp, A)
D = diag(sum(A,2));
L = D - A;
u = -yp*L';

dx = -10*xp+ 10*yp;                     %% Lorenz
dy = -xp.*zp + 28*xp - yp + u;
dz = xp.*yp - 8/3*zp;

dd = [dx;dy;dz];
end


function y = rk(x1, x2, x3, A, ts)
k1 = f(x1, x2, x3, A);
% k2 = f(x1 + 0.5*ts*k1(1,:), x2 + 0.5*ts*k1(2,:), x3 + 0.5*ts*k1(3,:), A);
% k3 = f(x1 + 0.5*ts*k2(1,:), x2 + 0.5*ts*k2(2,:), x3 + 0.5*ts*k2(3,:), A);
% k4 = f(x1 + ts*k3(1,:), x2 + ts*k3(2,:), x3 + ts*k3(3,:), A);

% y = [x1;x2;x3] + 1/6*ts*(k1+2*k2+2*k3+k4);
y = [x1;x2;x3] + ts*(k1);
end

function y = update_data(X,x)
X(2:end,:,:)=X(1:end-1,:,:);
X(1,:,:)=x;

y = X;

end
