# SOC-Estimation-of-Li-ion-battery-using-Kalman-filter
1) This paper was written as a part of coursework for UTA EE 5322 under Dr Frank Lewis
2) Refer the Conference paper for reference. 
3) Access the Matlab file for validation. The following code was used for validation. 

data = xlsread('0deg1.xls');
k = 7000;
A1 = eye(2);
A2 = eye(2);
H1 = [1,1];
H2 = [1,1];
x01 = [4.15;0];
x02 = [0;0];
xhat1=[4.26;0];
xhat2=[0.99;0];
P1 = 0.000000001*eye(2);
P2 = 0.07*eye(2);
Q1 = 0.0000001*eye(2);
Q2 = 0.05*eye(2);
G1 = 0.0000008*eye(2);
G2 = 0.0015*eye(2);
w= rand(2,7000);
vk1 = rand(1,2);
vk2 = rand(1,2);
r1 = [0.5];
r2 = [0.8];
for j = 1:k
x1(:,j+1)=A1*x01+G1*w(:,j);
x2(:,j+1)=A2*x02+G2*w(:,j);
x01=x1(:,j+1);
x02=x2(:,j+1);
z1(:,j+1) = H1*x01 + vk1;
z2(:,j+1) = H2*x02 + vk2;
p1mi= A1*P1*A1'+ G1*Q1*G1';
p2mi= A2*P2*A2'+ G2*Q2*G2';
xhm1(:,j+1) = A1*x01;
xhm2(:,j+1) = A2*x02;
K1 = p1mi*H1' * inv(H1*p1mi*H1'+r1);
K2 = p2mi*H2' * inv(H2*p1mi*H2'+r2);
t1(:,j+1) = K1;
t2(:,j+1) = K2;
P1 = (eye(2)-K1*H1)*p1mi;
P2 = (eye(2)-K2*H2)*p2mi;
xhat1(:,j+1) = x01 + K1 .* (z1(:,j+1) - H1*xhm1(:,j+1)) ;
xhat2(:,j+1) = x02 + K2 .* (z2(:,j+1) - H2*xhm2(:,j+1)) ;
 end
plot(xhat1(1,:))
hold on
plot(data(1:7000,8), 'g')
legend('estimated mV','original mV')
xlabel('time stamp')
ylabel('Voltage in mV')
title('States and estimates of xk vs xhatk for Voltage')
figure
plot(xhat2(1,:))
hold on
plot(data(60:7060,7), 'r')
legend('estimated mA','original mA')
xlabel('time stamp')
ylabel('Current in mA')
title('States and estimates of xk vs xhatk for current')
