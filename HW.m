clc;clear;close all;
T=1;%雷达扫描周期
N=100/T;%总的采样次数
X=zeros(4,N);%目标真实位置、速度
X(:,1)=[-100,2,200,20];%目标初始位置、速度
S(:,1)=[-100,2,200,20];
Z=zeros(2,N);%传感器对位置的观测
Z(:,1)=[X(1,1),X(3,1)];%观测初始化

delta_w=1e-2; %如果增大这个参数，目标真实轨迹就是曲线了
Q=delta_w*diag([0.5,1,0.5,1]);%过程噪声均值
R=eye(2);%观测噪声均值

phi=[1,T,0,0;0,1,0,0;0,0,1,T;0,0,0,1];%状态转移矩阵
H=[1,0,0,0;0,0,1,0];%观测矩阵
Xn=zeros(4,0);
for i=2:N
    S(:,i)=phi*S(:,i-1);%目标理论轨迹
    X(:,i)=phi*X(:,i-1)+sqrtm(Q)*randn(4,1);%目标真实轨迹
    Z(:,i)=H*X(:,i)+sqrtm(R)*randn(2,1);%对目标的观测
end

% Kalman 滤波
Xkf=zeros(4,N);
Xkf(:,1)=X(:,1);%卡尔曼滤波状态初始化
M(1,:)=Xkf(:,1);
P0=100e-2*eye(4);%协方差阵初始化


for i=2:N
    Xn=phi*Xkf(:,i-1);%预测
    M(i,:)=Xn;

    P1=phi*P0*phi'+Q;%预测误差协方差
    K=P1*H'*inv(H*P1*H'+R)%增益
    Xkf(:,i)=Xn+K*(Z(:,i)-H*Xn);%状态更新
    P0=(eye(4)-K*H)*P1;             %滤波误差协方差更新
end
% 误差分析
for i=1:N

    Err_Observation(i)=RMSx(Z(:,i),X(:,i));%滤波前的误差
    Err_KalmanFilter(i)=RMSx(X(:,i),Xkf(:,i));%滤波后的误差
end


figure
hold on;box on;
plot(S(1,:),S(3,:),'g','LineWidth',1);%理想轨迹
plot(X(1,:),X(3,:),'b','LineWidth',1);%真实轨迹
plot(Z(1,:),Z(2,:),'r','LineWidth',1);%观测轨迹
plot(Xkf(1,:),Xkf(3,:),'c','LineWidth',1);%卡尔曼滤波轨迹
plot(M(1,:),M(3,:),'k','LineWidth',1);%一步预测轨迹
legend('理论轨迹','实际运动轨迹','观测轨迹','滤波后轨迹','一步预测轨迹');
xlabel('横坐标 X/m');
ylabel('纵坐标 Y/m');
 
figure
hold on;box on;
plot(Err_Observation);
plot(Err_KalmanFilter);
legend('滤波前误差','滤波后误差');
xlabel('观测时间/s');
ylabel('误差值');

% 计算欧氏距离子函数
function dist=RMS(X1,X2)
    if length(X2)<=2
        dist=sqrt((X1(1)-X2(1))^2+(X1(3)-X2(2))^2);
    else
        dist=sqrt((X1(1)-X2(1))^2+(X1(3)-X2(3))^2);
    end
end
function dist=RMSx(X1,X2)
        dist=(X1(1)-X2(1))^2;
end
