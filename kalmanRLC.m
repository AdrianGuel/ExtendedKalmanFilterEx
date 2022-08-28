%kalman filter example
%Aplication in a RLC system , we measure voltage on C and estimate current on L
%Author: Adrian-Josue Guel-Cortez 2022

%Solution to the model
R=10; L=1; C=0.1;
T=1e-3;
t=0:T:10;
A=[-1/(C*R),1/C;-1/L,0];
B=[0;1/L];
Vi=0*sin(0.3*t);
y=zeros(2,length(t));
y(:,1)=[1;0.1];
for j=2:length(t)
    y(:,j)=y(:,j-1)+T*(model(y(:,j-1),Vi(j-1),R,L,C)+normrnd(0,1e-5,[2,1]));
end

%Kalman filter
n=2;
P_previous=[1e-5, 0;0 , 1e-3];
Q=((1e-5)^2)*eye(n,n);
R=(1e-2)^2; H=[1,0]; 
alpha=1e-4;
epsilon=1e-3;
x_previous=[1;0.1];
x=zeros(n,length(t));
for k=2:length(t)
    z_k=y(1,k)+normrnd(0,sqrt(R));
    x_model=A*x_previous+B*Vi(k);
    P_kminus=A*P_previous*A'+Q;
    K_k=P_kminus*H'*((H*P_kminus*H'+R)^(-1));
    x_k=x_model+K_k*(z_k-H*x_model);
    P_k=(eye(n,n)-K_k*H)*P_kminus;
    %Q=(1-alpha)*Q+epsilon*K_k*(z_k-H*x_model)*(z_k-H*x_model)'*K_k';
    x(:,k)=x_k;
    x_previous=x_k;
    P_previous=P_k;
end

figure
set(gcf,'color','w');
subplot(2,1,1)
plot(t,y(1,:),'b')
hold on
plot(t,x(1,:),'r')
xlabel('t')
ylabel('V_c')

subplot(2,1,2)
plot(t,y(2,:))
hold on
plot(t,x(2,:),'r')
xlabel('t')
ylabel('i_L')
legend('true value','estimation')

function aux=model(x,Vi,R,L,C)
aux=[-x(1)/(C*R)+x(2)/C;
    -x(1)/L+Vi/L];
end
