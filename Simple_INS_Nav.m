clear all;
close all;
obsdata=dlmread('IMU.txt');
obsdata(:,5:6)=-obsdata(:,5:6);
obsdata(:,2:3)=-obsdata(:,2:3);
[no_epochs,no_columns] = size(obsdata);
% Parameters
deg_to_rad = pi/180;
omega_ie = 7.292115E-5;  % Earth rotation rate (rad/s)

a = 6378137.0000;	% Earth radius in meters
b = 6356752.3142;	% Earth semiminor in meters	
f = (a-b)/a;% Eccentricity 

%初始化
% yaw = 93.87999*deg_to_rad;
% roll= 0.69465*deg_to_rad;
% pitch= 1.75973*deg_to_rad;
yaw = 95.71560*deg_to_rad;
roll= 0.66252*deg_to_rad;
pitch= 1.74943*deg_to_rad;
t_yaw=zeros(3,3);
t_yaw(1,1)=cos(yaw); 
t_yaw(2,2)=t_yaw(1,1);
t_yaw(1,2)=sin(yaw);
t_yaw(2,1)=-t_yaw(1,2);
t_yaw(3,3)=1;
t_pitch=zeros(3,3);
t_pitch(1,1)=1;
t_pitch(2,2)=cos(pitch);
t_pitch(3,3)=t_pitch(2,2);
t_pitch(2,3)=sin(pitch);
t_pitch(3,2)=-t_pitch(2,3);
t_roll=zeros(3,3);
t_roll(1,1)=cos(roll);
t_roll(3,3)=t_roll(1,1);
t_roll(2,2)=1;
t_roll(1,3)=-sin(roll);
t_roll(3,1)=-t_roll(1,3);






%初始坐标
preRbeXYZ=[-2279083.2976     5008645.5754     3214152.3198]';
% initLat=30.4568478102*deg_to_rad;
% initLng=114.4669524500*deg_to_rad;
[initLat,initLng,h]=ecef2llh(preRbeXYZ');



C_N_e=zeros(3,3);
C_N_e(1, 1) =- sin(initLng);C_N_e(1, 2) = -sin(initLat)*cos(initLng);C_N_e(1, 3) = cos(initLat)*cos(initLng);
C_N_e(2, 1) =  cos(initLng);C_N_e(2, 2) = -sin(initLat)*sin(initLng);C_N_e(2, 3) = cos(initLat)*sin(initLng);
                            C_N_e(3,2) = cos(initLat);               C_N_e(3, 3) = sin(initLat);

% C_N_e(1, 1) =- sin(initLng);C_N_e(1, 2) = -sin(initLat)*cos(initLng);C_N_e(1, 3) = cos(initLat)*cos(initLng);
% C_N_e(2, 1) =  cos(initLng);C_N_e(2, 2) = -sin(initLat)*sin(initLng);C_N_e(2, 3) = cos(initLat)*sin(initLng);
%                             C_N_e(3,2) = cos(initLat);               C_N_e(3, 3) = sin(initLat);

C_B_N=t_yaw'*t_pitch'*t_roll';

prCbe=C_N_e*C_B_N;
%初始速度
preVbe=[0,0,0]';
%初始时刻
preWs=obsdata(1,1);
vbe=zeros(no_epochs,3);
vbe(1,:)=preVbe';
ckvbe(1,:)=preVbe';
% geoPoses=zeros(no_epochs,2);
% geoPoses(1,:)=[30.4568478102*deg_to_rad, 114.4669524500*deg_to_rad];
rbeXYZs=zeros(no_epochs,3);
rbeXYZs(1,:)=preRbeXYZ';
cbes=zeros(no_epochs,3);

cbes(1,1)=pitch;
cbes(1,2)=roll;
cbes(1,3)=yaw;
[L_b,lambda_b,h_b,v_eb_n,C_b_n] = (ECEF_to_NED(preRbeXYZ,preVbe,prCbe));
cbes(1,:)=CTM_to_Euler([0,1,0;1,0,0;0,0,1]*C_b_n)';
ZVPT_THREDHOLD=0.5;
ZVPT_WINDOW=20;
C_B_E_S=zeros(no_epochs,9);
C_B_E_S(1,:)=reshape(prCbe,1,9);
idz=ZVPT_WINDOW;
%迭代
for epoch = 2:no_epochs
    ws = obsdata(epoch,1);
    gyroX = obsdata(epoch,2);
    gyroY = obsdata(epoch,3);
    gyroZ = obsdata(epoch,4);
    acceX = obsdata(epoch,5);
    acceY = obsdata(epoch,6);
    acceZ = obsdata(epoch,7);
    %时间间隔
    tor_i=ws-preWs;
    
    %比力
    f_ib_b=[acceX,acceY,acceZ]';
    %旋转角
    omega_ib_b=[gyroX,gyroY,gyroZ]';
    %状态更新
    [r_eb_e,v_eb_e,C_b_e] = Nav_equations_ECEF(tor_i,preRbeXYZ,...
        preVbe,prCbe,f_ib_b,omega_ib_b);
     ckvbe(epoch,:)=v_eb_e';
     rbeXYZs(epoch,:)=r_eb_e';
     C_B_E_S(epoch,:)=reshape(C_b_e,1,9);
     if epoch>idz
         vel_m=mean(ckvbe(epoch-idz:epoch,:));
         if(norm(vel_m)<=ZVPT_THREDHOLD)
             rbeXYZs(epoch,:)=mean(rbeXYZs(epoch-idz:epoch,:))';
             r_eb_e= rbeXYZs(epoch,:)';
             C_B_E_S(epoch,:)=mean(C_B_E_S(epoch-idz:epoch,:))';
             C_b_e=reshape(C_B_E_S(epoch,:),3,3);
%              ckvbe(epoch,:)=vel_m';
%              v_eb_e=ckvbe(epoch,:)';
         end
     end
    
    
    preWs=ws;
    preRbeXYZ=r_eb_e;
    preVbe=v_eb_e;
    prCbe=C_b_e;
    
    
    
    %更新
    [L_b,lambda_b,h_b,v_eb_n,C_b_n] = ECEF_to_NED(r_eb_e,v_eb_e,C_b_e);



   
    vbe(epoch,:)=([0,1,0;1,0,0;0,0,-1]*v_eb_n)';
   cbes(epoch,:)=CTM_to_Euler([0,1,0;1,0,0;0,0,1]*C_b_n)';
end


% figure
% plot(geoPoses(:,2)/deg_to_rad,geoPoses(:,1)/deg_to_rad);
pkdata=dlmread('check.txt');
figure
plot(vbe(:,1),'r');
hold on
plot(pkdata(:,12),'g');
figure
plot(vbe(:,2),'r');
hold on
plot(pkdata(:,13),'g');
figure
plot(vbe(:,3),'r');
hold on
plot(pkdata(:,14),'g');
figure
plot(cbes(:,2)/deg_to_rad,'r');
hold on
plot(pkdata(:,15),'g');
figure
plot(cbes(:,1)/deg_to_rad,'r');
hold on
plot(pkdata(:,16),'g');
figure
plot(cbes(:,3)/deg_to_rad,'r');
hold on
plot(pkdata(:,17),'g');
figure
plot3(rbeXYZs(:,1),rbeXYZs(:,2),rbeXYZs(:,3),'r')
hold on
plot3(pkdata(:,9),pkdata(:,10),pkdata(:,11),'g')
% hold on
% cdata=dlmread('vel.txt');
% plot3(cdata(:,1),cdata(:,2),cdata(:,3),'b+')
figure
plot(rbeXYZs(:,1),'r')
hold on
plot(pkdata(:,9),'g')
figure
plot(rbeXYZs(:,2),'r')
hold on
plot(pkdata(:,10),'g')

figure
plot(rbeXYZs(:,3),'r')
hold on
plot(pkdata(:,11),'g')


