%Extended kalman filter example
% Aplication in a simple mass-spring-damper system
%Author: Adrian-Josue Guel-Cortez 2020
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

%Kalman filter
n=4;
P_previous=50000*eye(n,n);
Q=zeros(n,n);
Q(1,1)=5000; Q(2,2)=5000;
R=0.01; H=[1,0,0,0]; V=1;
W=eye(n,n);

l=1;
x_previous=[1,0.1,0.2,0.4,0.2];
x=zeros(l*length(t),n);
k=2;
for i=2:l*length(t)
    x_model=f(x_previous,F(i-1),[0 0 0 0 0],m);
    P_kminus=A(x_model,m)*P_previous*A(x_model,m)'+W*Q*W';
    K_k=P_kminus*H'*((H*P_kminus*H'+V*R*V')^(-1));
    x_k=x_model+K_k*(y(1,i)+0.01*randn(1,1)-x_model(1));
    P_k=(eye(n,n)-K_k*H)*P_kminus;
    x(i,:)=x_k';
    x_previous=x_k;
    P_previous=P_k;
%         k=k+1;
%     if k==length(t)
%        k=2;
%        P_previous(1:3,1:3)=50*eye(3,3);
%     end
end

figure
set(gcf,'color','w');
subplot(3,2,[1 2])
plot(t,y(1,:),'b')
hold on
plot(t,x(:,1),'--r','LineWidth',2)
subplot(3,2,3)
plot(t,y(2,:),'b')
hold on
plot(t,x(:,2),'--r','LineWidth',2)
subplot(3,2,4)
plot(t,x(:,3))
subplot(3,2,5)
plot(t,x(:,4))
k_estimated=x(end,3);
b_estimated=x(end,4);

function aux=A(x,m)
T=1e-3;
aux=[1 T 0 0;
    -T*x(3)/m 1-T*x(4)/m -T*x(1)/m -T*x(2)/m;
    0 0 1 0;
    0 0 0 1];
end

function aux=f(x,F,w,m)
T=1e-3;
aux=[x(1)+T*x(2)+T*w(1);
    x(2)-T*x(3)*x(1)/m-T*x(4)*x(2)/m+T*F/m+T*w(2);
    x(3)+T*w(3);
    x(4)+T*w(4)];
end

function aux=model(x,F,m,b,k)
aux=[x(2);
    -(k/m)*x(1)-(b/m)*x(2)+F];
end
