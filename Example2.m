clear all;clc;close all
redblue=[0:0.024:0.96,0.96:0.001:1; 0:0.024:0.96,0.96:-0.024:0;1:-0.001:0.96,0.96:-0.024:0 ]'; % colorbar

load('Source_data_for_example_2.mat');
% loaded data:
% A  : adjacency matrix
% D  : Degree matrix
% a1t: Theoretically obtained 'true' A1 matrix
% x  : trajectories of 10 Lorenz oscillators   [step * node index * state component];

time = 80;
ts = 0.01;   % length of a time step
tc = 40;     % time of topology change in the network

stp = time/ts;

N = 10;    % number of nodes

m =350;    % recording length (amount of data)

rho =1e-3;       %  parameters for ADMM
lambda = 1e-7;


xr = rand(m,N,3);    % initilize a data storage to store data from a m-step-long horizon

X = plift(xr(:,:,1),xr(:,:,2),xr(:,:,3));   % plift: the Psi(z) observable set
[q,~]=size(X);
C = [xr(:,:,1)';xr(:,:,2)';xr(:,:,3)']*pinv(X);

clear xr X

xr = zeros(m,N,3);     % re-initialize the storages
yr = zeros(m,N,3);


A0 = zeros(q,q);   % initialize variables of ADMM
Z = zeros(q,q);
W = zeros(q,q);

Ar = zeros( stp, q*N*3);    %  to record the time evolution of entries in A1-A2
er = zeros(stp,1);

for i =2:stp
    xr = update_data(xr, x(i-1,:,:));   % update the data of m steps stroed in the storage
    yr = update_data(yr, x(i,  :,:));
   
    X = plift(xr(:,:,1),xr(:,:,2),xr(:,:,3));   %  calculate the data matrices using the observable set Psi(z)
    Y = plift(yr(:,:,1),yr(:,:,2),yr(:,:,3));
    
    A0 = -0.5*(X*X'+rho/2*eye(q))^-1 *(-2*X*Y' + W-rho*Z);  % update variables via ADMM
    Z = sth(A0+1/rho*W,lambda/rho);
    W = W+ rho*(A0-Z);
    
    A1 = C*(A0)';
    A2 = getA2(A1);
    
    Ar(i,:) = reshape(A1-A2,1,[]);     
    
    er(i,:) = norm(A1-a1t,'f')/norm(a1t,'f')*100;  % calculate and record the error index \mathcal{I}(x)

    if i==tc/ts
        A(1,3)=0;      % cut off the connection between node 1 and node 3 
        A(3,1)=0;
        D = diag(sum(A,2));
        L = D - A;
        
        a1t(11,14) = 0;
        a1t(13,12) = 0;
        a1t(11,12) = 0.98;
        a1t(13,14)=0.97;
    end
end

figure(1)
subplot(2,1,1)
hold on
c13 = Ar(:,401);
c31 = Ar(:,343);

gr=plot(ts:ts:time,Ar','g');
rd=plot(ts:ts:time,c13,'--r','linewidth',2);
bl=plot(ts:ts:time,c31,'--b','linewidth',2);
legend([rd,bl,gr(1)],{'$$c_{13}$$','$$c_{31}$$','other entries'},'interpreter','latex')
% xlabel('time[s], $$h=0.01$$s','interpreter','latex')
ylabel('entries in $$A_1-A_2$$','interpreter','latex')
axis([-inf inf -0.08 0.08])
text(2,-0.05,'$$(a)$$','interpreter','latex')

subplot(2,1,2)
semilogy(ts:ts:time,er)
hold on
semilogy(ts:ts:time,er*0+5,'--k')
legend({'error index $$\mathcal{I}(A_1)$$','$$5\%$$ error'},'interpreter','latex')
ylabel('error(\%) $$\mathcal{I}[k]$$','interpreter','latex')
xlabel('time[s], $$h=0.01$$s','interpreter','latex')
yticks([1e-2 1e-1 1e-0 1e1 1e2 1e3])
text(2,0.1,'$$(b)$$','interpreter','latex')
axis([-inf inf 1e-3 1e3])


function y=plift(x1,x2,x3)   % calculate the observables
[m,~]=size(x1);
x = [x1';x2';x3'];
y = [ones(1,m);x;x1'.*x2';x1'.*x3';x2'.*x3'];
end

function y = sth(x,s)   % soft-thresholding
y=(abs(x)>s).*(x-sign(x)*s);
end

function y = update_data(X,x)  
X(2:end,:,:)=X(1:end-1,:,:);
X(1,:,:)=x;
y = X;
end

function y = getA2(A)
A(11:20,12:21) = diag(sum(A(11:20,12:21),1));
y = A;
end
