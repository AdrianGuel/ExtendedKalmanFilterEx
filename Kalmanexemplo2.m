%Kalman filter example
% Aplication in a simple mass-spring-damper system
%Author: Adrian-Josue Guel-Cortez 2022
%In this case, we know m but we do not know k and b.

%Solution to the model
t=0:1e-3:15;
F=zeros(1,length(t));
%F=0.1*sin(t);
k=0.5; b=0.7; m=1;
h=1e-3;
y=zeros(2,length(t));
y(:,1)=[1;0.1];
for j=2:length(t)
    y(:,j)=y(:,j-1)+h*model(y(:,j-1),F(j-1),m,b,k);
end
plot(t,y(2,:))

%Kalman filter
n=2;
A=[0,1;-k/m,-b/m];
Ad=eye(n,n)-A;
P_previous=50000*eye(n,n);
Q=zeros(n,n);
Q(1,1)=50; Q(2,2)=50;
R=0.01; H=[1,0]; V=1;
W=eye(n,n);

x_previous=[0.3,0.5];
x=zeros(l*length(t),n);
for i=2:l*length(t)
    x_model=f(x_previous,Ad,F(i-1),[0 0 0 0 0]);
    P_kminus=Ad*P_previous*Ad'+Q;
    K_k=P_kminus*H'*((H*P_kminus*H'+R)^(-1));
    x_k=x_model+K_k*(y(1,i)-H*x_model);
    P_k=(eye(n,n)-K_k*H)*P_kminus;
    x(i,:)=x_k';
    x_previous=x_k;
    P_previous=P_k;
end


figure
set(gcf,'color','w');
subplot(2,1,1)
plot(t,y(1,:),'b')
hold on
plot(t,x(:,1),'--r','LineWidth',2)
subplot(2,1,2)
plot(t,y(2,:),'b')
hold on
plot(t,x(:,2),'--r','LineWidth',2)


function aux=f(x,Ad,F,w)
T=1e-3;
aux=Ad*x+T*B*F+w';
end


function aux=model(x,F,m,b,k)
aux=[x(2);
    -(k/m)*x(1)-(b/m)*x(2)+F+rand()];
end
