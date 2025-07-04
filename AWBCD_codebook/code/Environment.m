

global C M N dy dz

%% __________________________________________________________________________________________________________

%房间坐标
room_x_w=20;%房间沿x轴长度
room_y_w=20;%房间沿y轴长度
room_z_w=5;%房间沿z轴长度
room_x = 15;  
room_y = 15;  
room_z = room_z_w/2;
x_room= [room_x - room_x_w/2, room_x + room_x_w/2, room_x + room_x_w/2, room_x - room_x_w/2, room_x - room_x_w/2, room_x + room_x_w/2, room_x + room_x_w/2, room_x - room_x_w/2];  
y_room = [room_y - room_y_w/2, room_y - room_y_w/2, room_y + room_y_w/2, room_y + room_y_w/2, room_y - room_y_w/2, room_y - room_y_w/2, room_y + room_y_w/2, room_y + room_y_w/2];  
z_room = [room_z - room_z_w/2, room_z - room_z_w/2, room_z - room_z_w/2, room_z - room_z_w/2,room_z + room_z_w/2, room_z + room_z_w/2, room_z + room_z_w/2, room_z + room_z_w/2];

%绘图
% room_plot


%% __________________________________________________________________________________________________________
%频率
f=30*10^9;
lambda=3*10^8/f;

%基站坐标
x_BS=0;
y_BS=0;
z_BS=10;
M=1;
space=lambda/2;
AN_x=x_BS+zeros(M,1);
AN_y = y_BS - (M-1)/2*space : space : y_BS + (M-1)/2*space;
AN_z = z_BS;
[AN_y, AN_z] = meshgrid(AN_y, AN_z);


%旋转
y_theta =60;  
Rz = [cosd(y_theta),-sind(y_theta),0;sind(y_theta),cosd(y_theta),0;0,0,1];   
coords_AN = [AN_x(:)-x_BS, AN_y(:)-y_BS, AN_z(:)-z_BS];    
rotated_coords_AN = (Rz * coords_AN')';
AN_x_rotated = rotated_coords_AN(:,1)+x_BS;  
AN_y_rotated = rotated_coords_AN(:,2)+y_BS;  
AN_z_rotated = rotated_coords_AN(:,3)+z_BS; 

%法向量    
% AN = [AN_x_rotated(1) - AN_x_rotated(2); AN_y_rotated(1) - AN_y_rotated(2); AN_z_rotated(1)- AN_z_rotated(2)];  
% AM = [AN_x_rotated(3) - AN_x_rotated(1); AN_y_rotated(3) - AN_y_rotated(1); AN_z_rotated(3)- AN_z_rotated(1)];  
% AN_v = cross(AN, AM)/norm(cross(AN, AM));

%功率频率
P_tx_all_dbm=44;
P_tx=db2pow(P_tx_all_dbm)*10^-3/M*ones(M,1);
G_tx_db=0*ones(M,1);
G_tx=db2pow(G_tx_db);

%画图
% transmitting_antenna_plot


%% __________________________________________________________________________________________________________
%RIS坐标
J=3;
dy=lambda/2;
dz=lambda/2;
N_y=nn;
N_z=nn;
RIS_y_w=N_y*dy;
RIS_z_w=N_z*dz;
RIS1_x=room_x - room_x_w/2;
RIS1_y=15;
RIS1_z=room_z_w*3.5/5;
RIS2_x=room_x - room_x_w/2;
RIS2_y=10;
RIS2_z=room_z_w*3.5/5;
RIS3_x=room_x - room_x_w/2;
RIS3_y=20;
RIS3_z=room_z_w*3.5/5;
for m=1:N_z
    for n=1:N_y
        x_cell1(m,n)=RIS1_x;
        y_cell1(m,n)=RIS1_y-RIS_y_w/2+(n-1/2)*dy;%沿y轴电磁元坐标
        z_cell1(m,n)=RIS1_z+RIS_z_w/2-(m-1/2)*dz;%沿z轴电磁元坐标
    end
end
for m=1:N_z
    for n=1:N_y
        x_cell2(m,n)=RIS2_x;
        y_cell2(m,n)=RIS2_y-RIS_y_w/2+(n-1/2)*dy;%沿y轴电磁元坐标
        z_cell2(m,n)=RIS2_z+RIS_z_w/2-(m-1/2)*dz;%沿z轴电磁元坐标
    end
end
for m=1:N_z
    for n=1:N_y
        x_cell3(m,n)=RIS3_x;
        y_cell3(m,n)=RIS3_y-RIS_y_w/2+(n-1/2)*dy;%沿y轴电磁元坐标
        z_cell3(m,n)=RIS3_z+RIS_z_w/2-(m-1/2)*dz;%沿z轴电磁元坐标
    end
end
x_cell=[x_cell1 x_cell2 x_cell3];
y_cell=[y_cell1 y_cell2 y_cell3];
z_cell=[z_cell1 z_cell2 z_cell3];
x_cell_vector=x_cell(:);
y_cell_vector=y_cell(:);
z_cell_vector=z_cell(:);
x_RIS=[RIS1_x RIS2_x RIS3_x];
y_RIS=[RIS1_y RIS2_y RIS3_y];
z_RIS=[RIS1_z RIS2_z RIS3_z];

% x_cell=[x_cell1];
% y_cell=[y_cell1];
% z_cell=[z_cell1];
% x_cell_vector=x_cell(:);
% y_cell_vector=y_cell(:);
% z_cell_vector=z_cell(:);
% x_RIS=[RIS1_x];
% y_RIS=[RIS1_y];
% z_RIS=[RIS1_z];

%法向量
RIS1_v=[1 0 0];


%RIS电磁性质
G=8;%电磁元增益
D=sqrt(RIS_y_w^2+RIS_z_w^2);%阵列最大几何尺寸
L=2*D^2/lambda;%弗劳恩霍夫距离
A_T=0.9;
G_rx_db=0;%接收天线增益
G_rx=db2pow(G_rx_db);

%初始化
phi=2*pi*rand(1,J*N_y*N_z);
%phi=0*rand(1,J*N_y*N_z);
N=J*N_y*N_z;
Gamma_initial=sqrt(A_T).*exp(1j*phi);

%绘图
% RIS_plot
if bit~=0
    Qphi=2.*pi.*[0:2^bit]./2^bit;
end
%% __________________________________________________________________________________________________________ 
%目标区域
gap=0.1;
aim_wide=3;
aim_wx=aim_wide;
aim_wy=aim_wide;

% %中心坐标

% % aim_z=RIS1_z;
aim_x=10;
aim_y=10;
aim_z=0.5;


% % % 圆环参数
% r1 = 8; 
% r2 = 10;
% 
% theta1=deg2rad(-5);
% theta2=deg2rad(5);
% theta_cen=(theta1+theta2)/2;
% r_cen=(r1+r2)/2;
% aim_x=r_cen*cos(theta_cen)+RIS1_x;
% aim_y=r_cen*sin(theta_cen)+RIS1_y;
% theta = linspace(theta1,theta2, 60);%60
% radii = linspace(r1, r2, 30);%30
% % 构建圆环点坐标
% [Theta, Radii] = meshgrid(theta, radii);
% X = Radii .* cos(Theta) + RIS1_x; % x坐标
% Y = Radii .* sin(Theta) + RIS1_y; % y坐标
% x_aim=X(:);
% y_aim=Y(:);
% z_aim=aim_z*ones(length(x_aim),1);




x_coords=15-2+gap/2:gap:15+2-gap/2;
y_coords=15-3+gap/2:gap:15-2-gap/2;
z_coords=aim_z;
[X_grid, Y_grid, Z_grid] = meshgrid(x_coords,y_coords,z_coords);  
x_aim=X_grid(:);
y_aim=Y_grid(:);
z_aim=Z_grid(:);


x_coords1=15-2+gap/2:gap:15-1-gap/2;
y_coords1=15-2+gap/2:gap:15+2-gap/2;
z_coords1=aim_z;
[X_grid1, Y_grid1, Z_grid1] = meshgrid(x_coords1,y_coords1,z_coords1);  
x_aim1=X_grid1(:);
y_aim1=Y_grid1(:);
z_aim1=Z_grid1(:);
x_aim=[x_aim1;x_aim];
y_aim=[y_aim1;y_aim];
z_aim=[z_aim1;z_aim];

A_r=lambda^2/(4*pi);
G_r=A_r*4*pi/lambda^2;
weight=ones(1,size(x_aim,1));
noise_dbm=-105;
noise=db2pow(noise_dbm)*10^-3;
C = length(x_aim);

%绘图
%aim_plot

%%__________________________________________________________________________________________________________
%计算G
r_BS_cell = sqrt(sum(( reshape([AN_x_rotated, AN_y_rotated,AN_z_rotated], [M, 1, 3])-reshape([x_cell_vector, y_cell_vector, z_cell_vector], [1,length(x_cell_vector), 3]) ).^2, 3));
theta_cell_BS = compute_angles(x_cell_vector, y_cell_vector, z_cell_vector,AN_x_rotated, AN_y_rotated, AN_z_rotated, RIS1_v);
F_tx=ones(M,length(x_cell_vector));
F_r=cosd(theta_cell_BS).^3; 

%计算H
r_cell_aim = sqrt(sum(( reshape([x_cell_vector, y_cell_vector, z_cell_vector], [length(x_cell_vector), 1, 3]) - reshape([x_aim, y_aim,z_aim], [1, length(x_aim), 3])).^2, 3));
theta_cell_aim = compute_angles(x_aim,y_aim,z_aim, x_cell_vector, y_cell_vector, z_cell_vector,RIS1_v);
F_t=cosd(theta_cell_aim).^3;
F_rx=ones(length(x_cell_vector),length(x_aim));

%计算A F
F_tx_exp = repmat(F_tx, [1, 1, C]);
F_r_exp = repmat(F_r, [1, 1, C]);
F_t_exp = permute(repmat(F_t, [1, 1, M]), [3, 1, 2]);
F_rx_exp = permute(repmat(F_rx, [1, 1, M]), [3, 1, 2]);
F = F_tx_exp .* F_r_exp .* F_t_exp .* F_rx_exp;
clear F_tx_exp F_r_exp F_t_exp F_rx_exp
P_tx_exp = repmat(P_tx, [1, N, C]);
G_tx_exp = repmat(G_tx, [1, N, C]);
r_BS_cell_exp = repmat(r_BS_cell, [1, 1, C]);
r_cell_aim_exp = permute(repmat(r_cell_aim, [1, 1, M]), [3, 1, 2]);
A = sqrt(F .* P_tx_exp .* G_tx_exp) ./ (r_BS_cell_exp .* r_cell_aim_exp);


% %%__________________________________________________________________________________________________________
% %场景
% scenario_plot
% toc

