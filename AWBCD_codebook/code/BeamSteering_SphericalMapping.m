Beta=2;
Gamma=Gamma_initial;

v_x = aim_x - RIS1_x;
v_y = aim_y - RIS1_y;
v_z = aim_z - RIS1_z;
R=sqrt(v_x^2+v_y^2+v_z^2);
theta0  = atan2(v_y, v_x);
phi0 = acos(-abs(v_z)/R);

%%--------------------------------------------------------

X_r2 = r2 .* cos(theta_cen) + RIS1_x; % x坐标
Y_r2 = r2 .* sin(theta_cen) + RIS1_y; % y坐标
Z_r2 = aim_z;

X_r1 = r1 .* cos(theta_cen) + RIS1_x; % x坐标
Y_r1 = r1 .* sin(theta_cen) + RIS1_y; % y坐标
Z_r1 = aim_z;

v2=[X_r2-RIS1_x Y_r2-RIS1_y Z_r2-RIS1_z];
v1=[X_r1-RIS1_x Y_r1-RIS1_y Z_r1-RIS1_z];

dot_product12 = dot(v1, v2);
norm_v1 = norm(v1);
norm_v2 = norm(v2);
phi_rad12 = acos(dot_product12 / (norm_v1 * norm_v2));
phi_rad12 = mod(phi_rad12, pi / 2);
phi_rad12_f=phi_rad12/Beta;
phi_rad_12=linspace(phi0-phi_rad12_f/2,phi0+phi_rad12_f/2,N_z);

%%--------------------------------------------------------

X_r3 = r_cen .* cos(theta1) + RIS1_x; % x坐标
Y_r3 = r_cen .* sin(theta1) + RIS1_y; % y坐标
Z_r3 = aim_z;

X_r4 = r_cen .* cos(theta2) + RIS1_x; % x坐标
Y_r4 = r_cen .* sin(theta2) + RIS1_y; % y坐标
Z_r4 = aim_z;

v3=[X_r4-RIS1_x Y_r3-RIS1_y Z_r3-RIS1_z];
v4=[X_r4-RIS1_x Y_r4-RIS1_y Z_r4-RIS1_z];

dot_product34 = dot(v3, v4);
norm_v3 = norm(v3);
norm_v4 = norm(v4);
theta_rad34 = acos(dot_product34 / (norm_v3 * norm_v4));
theta_rad34 = mod(theta_rad34, pi / 2);
theta_rad34_f=theta_rad34/Beta;
theta_rad34=linspace(theta0-theta_rad34_f/2,theta0+theta_rad34_f/2,N_y);

[theta_grid,phi_grid]=meshgrid(theta_rad34,phi_rad_12);

% 计算直角坐标
x = R * cos(theta_grid) .* sin(phi_grid) + RIS1_x;
y = R * sin(theta_grid) .* sin(phi_grid) + RIS1_y;
z = R * cos(phi_grid) + RIS1_z;

x0 = R * cos(theta0) .* sin(phi0) + RIS1_x;
y0 = R * sin(theta0) .* sin(phi0) + RIS1_y;
z0 = R * cos(phi0) + RIS1_z;

% 绘制图形
% figure;
% scatter3(x, y, z); % 绘制球面
% hold on;
% scatter3(x0, y0, z0,'black'); % 绘制球面
% % 添加球心
% scatter3(RIS1_x, RIS1_y, RIS1_z, 100, 'filled', 'MarkerFaceColor', 'r'); % 球心标记
% 
% % 设置图形参数
% axis equal;
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% title('球面上的点分布');
% grid on;

%%---------------------------------------------------------------------------
d=sqrt((x_cell-x).^2+(y_cell-y).^2+(z_cell-z).^2)+reshape(r_BS_cell,[N_z,N_y]);


%bit
% teco=wrapTo2Pi(-d*2*pi/lambda);%最大
% for ci=1:size(teco,1)
%     for cj=1:size(teco,2)
%         index1=find(abs(teco(ci,cj)-Qphi)<=2*pi/2^qbit/2);
%         ph(ci,cj)=Qphi(index1);
%     end
% end
% Gamma=A_T.*exp(1j*ph);

%continue
Gamma=A_T.*exp(-1j*d*2*pi/lambda);

Gamma_qiu=Gamma(:)';
[P_rx,P_rx_mean,P_rx_dbm,P_rx_mean_dbm]=rx_power(r_BS_cell_exp,r_cell_aim_exp,G,G_rx,Gamma_qiu,lambda,F,P_tx_exp,G_tx_exp);
R=log2(1+P_rx/noise);
%%--------------------
% Save_com(:,count)=R;


% figure
% [n, m] = size(Gamma);
% for s=[0]
%     hold on
%     noisy_theta = s*randn(n, m);
%     Gamma_noisy = Gamma .*exp(1j*noisy_theta);
% 
%     %计算功率
%     [P_rx_finall,P_rx_finall_mean,P_rx_finall_mean_dbm]=rx_power(r_BS_cell_exp,r_cell_aim_exp,G,G_rx,Gamma_noisy,lambda,F,P_tx_exp,G_tx_exp);
%     mean(log2(1+P_rx_finall/noise))
%     [f_finall, x_finall] = ksdensity(log2(1+P_rx_finall/noise),'Width',0.1);
%     cdf_values_f_finall = cumsum(f_finall) / sum(f_finall);
%     [f_initial, x_initial] = ksdensity(log2(1+P_rx_initial/noise),'Width',0.1);
%     cdf_values_f_initial = cumsum(f_initial) / sum(f_initial);
%     plot(x_finall,cdf_values_f_finall,LineWidth=1.25);
%     plot(x_initial,cdf_values_f_initial,'-k',LineWidth=1.25);
%     hold off
%     ylim([0 1])
%     grid on
% end
% 
% 
% 
% figure
% [n, m] = size(Gamma);
% for s=[0]
%     hold on
%     noisy_theta = s*randn(n, m);
%     Gamma_noisy = Gamma .*exp(1j*noisy_theta);
% 
%     %计算功率
%     [P_rx_finall,P_rx_finall_mean,P_rx_finall_mean_dbm]=rx_power(r_BS_cell_exp,r_cell_aim_exp,G,G_rx,Gamma_noisy,lambda,F,P_tx_exp,G_tx_exp);
%     mean(pow2db(10^3*P_rx_finall))
%     [f_finall, x_finall] = ksdensity(pow2db(10^3*P_rx_finall),'Width',0.1);
%     cdf_values_f_finall = cumsum(f_finall) / sum(f_finall);
%     [f_initial, x_initial] = ksdensity(pow2db(10^3*P_rx_initial),'Width',0.1);
%     cdf_values_f_initial = cumsum(f_initial) / sum(f_initial);
%     plot(x_finall,cdf_values_f_finall,LineWidth=1.25);
%     plot(x_initial,cdf_values_f_initial,'-k',LineWidth=1.25);
%     hold off
%     ylim([0 1])
%     grid on
% end