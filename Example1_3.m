clear all;clc;close all
redblue=[0:0.024:0.96,0.96:0.001:1; 0:0.024:0.96,0.96:-0.024:0;1:-0.001:0.96,0.96:-0.024:0 ]';

stp = 1000;

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



range = 2;
beta = 1.45*pi;
delta = 0.525;




tst = rand(stp,N);
tnn = zeros(stp,N);

for kk = 1:max(size(tst))
    tnn(kk,:) = beta * I(tst(kk,:)) + delta ;
end
C = tst'*pinv(plift(tst));

An = C * plift(tnn) *pinv(plift(tst));

[k,~]=size(plift(tnn));

clear tst tnc tnn

rho = 4;       %  step size
lambda = 0.5;      %  l1 gain
rep = 1;

%%%%%%%%%

dataamount = 100:50:800;
Rep = 10;

ERR = zeros(length(dataamount),Rep);

cter = 1;

for DA = 1:length(dataamount)
    for rrep = 1:Rep

        rl = dataamount(DA);
        stp = rl*3;
        
        
        xr = zeros(rl,N);
        yr = zeros(rl,N);
        
        x = zeros(stp,N);
        xp = (rand(1,N) - 0.5) * range;
        
        er = zeros(stp,2);
        
        a1 = rand(k,k);
        
        z = rand(k,k);
        w = rand(k,k);
        
        a2 = zeros(N,k);
        
        for cc=1:stp
            
            x(cc,:) = xp;
            
            xr(2:rl,:) = xr(1:rl-1,:);
            xr(1,:) = xp;
            
            xp = beta * I(xp) + delta + I(xp)*L';
            
            yr(2:rl,:) = yr(1:rl-1,:);
            yr(1,:) = xp;
            
            
            X1 = plift(xr);
            Y1 = plift(yr);
            
            if cc>rl
                
                a1 = -0.5*(X1*X1'+rho/2*eye(k))^-1 *(-2*X1*Y1' + w-rho*( z));
                z = sth(a1+1/rho*w,lambda/rho);
                w = w+ rho*(a1-z);
                
                ar = C*a1';
                
                a2(:,1) =sum(ar(:,1))/N *ones(N,1);
                sar = sum(ar,1);
                for jj = 1:(k-1)/N
                    
                    a2(:,1+N*(jj-1)+(1:N)) = diag(sar(:,1+N*(jj-1)+(1:N) ));
                    
                end
                
            end
            
        end
        
        
        plotrange = 4;
        stepsize = 0.1;
        bias = 2;
        px = (-plotrange:stepsize:plotrange ) + bias;
        py = (-plotrange:stepsize:plotrange ) + bias;
        
        
        %
        ztrue1 = zeros(length(px),length(py));
        err   = zeros(length(px),length(py));
        
        s = 5;
        
        for i = 1:length(px)
            for j = 1:length(py)
                err(i,j) =[1,zeros(1,N-1)] * (ar-a2) * plift([px(i),zeros(1,s-2),py(j),zeros(1,N-s)]) - (D(1,1)*I(px(i)) - I(py(j)));
                ztrue1(i,j) =D(1,1)*I(px(i)) - I(py(j));
                
            end
        end
        
        
        
        ERR(DA,rrep) = norm(err,'f')/norm(ztrue1,'f');
        
        disp([num2str(cter/Rep/length(dataamount)*100),'%'])
        cter = cter+1;
        
    end
end


ERRav = sum(ERR,2)/Rep;
plot(dataamount/k,ERRav*100);
xlabel('data amount as multiples of $$q$$','interpreter','latex')
ylabel('error $$\mathcal{I}(A_1)$$','interpreter','latex')
grid on




function y=plift(x)
x = x';
[N,m] = size(x);

%
T = [0.1*pi:0.1*pi:3*pi];
y = zeros(length(T)*2,m);
for i=1:length(T)
    y((i-1)*2*N+(1:N),:) = sqevn(x,T(i));
    y((i-1)*2*N+N+(1:N),:) = sqodd(x,T(i));
    %     y((i-1)*2*N+N+(1:N),:) = cos(x);
end
% y = [ones(1,m);x;x.^2;sign(x).*abs(x).^0.5;y];

y = [ones(1,m);x;y];
end





function A = graph(N,s)
A = rand(N,N);
A = (A+A')/2 > s;
A = A.*(ones(N,N)-eye(N));
end

function y = sth(x,s)
y=(abs(x)>s).*(x-sign(x)*s);
end



function y = rk(x,A,w,ts)
k1 = f(x,A,w);
k2 = f(x+0.5*ts*k1,A,w);
k3 = f(x+0.5*ts*k2,A,w);
k4 = f(x+ts*k3,A,w);

y = x + 1/6*ts*(k1+2*k2+2*k3+k4);


end



function y = sqevn(x,T)

y = sqodd(x-T/4,T);

end


function y = sqodd(x,T)

plus  = (mod(x,T) < T/2) + 0;
minus = (mod(x,T) >= T/2) + 0;

y = plus-minus;

end


