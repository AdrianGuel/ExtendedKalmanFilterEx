%kalman filter example
%Aplication in a RLC system 
%Author: Adrian-Josue Guel-Cortez 2021

%Solution to the model
R=0.5; L=0.7; C=1;
T=1e-3;
t=0:T:15;
A=[-1/(C*R),1/C;-1/L,0];
B=[0;1/L];
Vi=sin(t);
y=zeros(2,length(t));
y(:,1)=[1;0.1];
for j=2:length(t)
    y(:,j)=y(:,j-1)+T*model(y(:,j-1),Vi(j-1),R,L,C);
end

%Kalman filter
n=2;
P_previous=1e-3*eye(n,n);
Q=eye(n,n);
R=1e-3; H=[1,0]; 

x_previous=[1;0.1];
x=zeros(n,length(t));
for k=2:length(t)
    x_model=A*x_previous+B*Vi(k);
    P_kminus=A*P_previous*A'+Q;
    K_k=P_kminus*H'*((H*P_kminus*H'+R)^(-1));
    x_k=x_model+K_k*(y(1,k)-H*x_model);
    P_k=(eye(n,n)-K_k*H)*P_kminus;
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

function aux=model(x,Vi,R,L,C)
aux=[-x(1)/(C*R)+x(2)/C;
    -x(1)/L+Vi/L];
end
