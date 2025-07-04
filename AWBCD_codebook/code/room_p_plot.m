
tic;
gap_room=3*gap;

% x_r=room_x-room_x_w/2+gap_room/2:gap_room:room_x+room_x_w/2-gap_room/2;
% y_r=room_y-room_y_w/2+gap_room/2:gap_room:room_y+room_y·_w/2-gap_room/2;
x_left=aim_x-aim_wx/2+gap/2:-gap_room:room_x-room_x_w/2;
x_right=aim_x-aim_wx/2+gap/2+gap_room:gap_room:room_x+room_x_w/2;
x_r=[x_left(end:-1:1) x_right];
y_left=aim_y-aim_wy/2+gap/2:-gap_room:room_y-room_y_w/2;
y_right=aim_y-aim_wy/2+gap/2+gap_room:gap_room:room_y+room_y_w/2;
y_r=[y_left(end:-1:1) y_right];
z_r=aim_z;
[X_grid_r, Y_grid_r, Z_grid_r] = meshgrid(x_r,y_r,z_r);  
x_room=X_grid_r(:);
y_room=Y_grid_r(:);
z_room=Z_grid_r(:);
C_room=length(x_room);

%计算
r_cell_room = sqrt(sum(( reshape([x_cell_vector, y_cell_vector, z_cell_vector], [length(x_cell_vector), 1, 3]) - reshape([x_room, y_room, z_room], [1, length(x_room), 3])).^2, 3));
theta_cell_room = compute_angles(x_room, y_room, z_room, x_cell_vector, y_cell_vector, z_cell_vector,RIS1_v);
F_t_room=cosd(theta_cell_room).^3;
F_rx_room=ones(length(x_cell_vector),length(x_room));

P_tx_exp = repmat(P_tx, [1, N, C_room]);
G_tx_exp = repmat(G_tx, [1, N, C_room]);
r_BS_cell_exp = repmat(r_BS_cell, [1, 1, C_room]);
r_cell_room_exp = permute(repmat(r_cell_room, [1, 1, M]), [3, 1, 2]);
F_tx_room_exp = repmat(F_tx, [1, 1, C_room]);
F_r_room_exp = repmat(F_r, [1, 1, C_room]);
F_te=F_tx_room_exp .* F_r_room_exp ;
clear F_tx_room_exp F_r_room_exp
F_t_room_exp=permute(repmat(F_t_room, [1, 1, M]), [3, 1, 2]);
F_rx_room_exp = permute(repmat(F_rx_room, [1, 1, M]), [3, 1, 2]);
F_room = F_te.* F_t_room_exp .* F_rx_room_exp;
clear F_t_room_exp F_rx_room_exp
toc

% clearvars -except Gamma_iteration
%%随机----------------------------------------------------------------------------------------------
% 获取 Gamma 的尺寸
Gamma=Gamma_iteration_weight;

%%----------------------------------------------------------------------------------------------
%计算最初功率
% [P_rx_room_initial,P_rx_room_initial_mean,P_rx_room_initial_mean_dbm]=rx_power(r_BS_cell_exp,r_cell_room_exp,G,G_rx,Gamma_initial,lambda,F_room,P_tx_exp,G_tx_exp);
% room_initial_dbm=reshape(pow2db(10^3*P_rx_room_initial),[length(x_r),length(y_r)]);

%计算最后功率
[P_rx_room_finall,P_rx_room_finall_mean,P_rx_room_finall_dbm,P_rx_room_finall_mean_dbm]=rx_power(r_BS_cell_exp,r_cell_room_exp,G,G_rx,Gamma,lambda,F_room,P_tx_exp,G_tx_exp);
room_finall_dbm=reshape(pow2db(10^3*P_rx_room_finall),[length(x_r),length(y_r)]);
room_finall_bps=reshape(log2(1+P_rx_room_finall/noise),[length(x_r),length(y_r)]);
toc

% figure
% pixel=[1000 1000];
% affine=pixel(1)/room_x_w;
% initial_data = imresize(flipud(room_initial_dbm), pixel, 'lanczos3');
% imagesc(initial_data);
% xticks(linspace(1,pixel(1),11));
% yticks(linspace(1,pixel(2),11));
% xticklabels(room_x - room_x_w/2:room_x + room_x_w/2);
% yticklabels(room_y + room_x_w/2:-1:room_x - room_x_w/2);
% % colormap('jet');
% colorbar
% grid off
% caxis([-80 -50]);
% rectangle('Position', [aim_x-aim_wx/2-(room_x-room_x_w/2), room_y_w-(aim_y+aim_wy/2-(room_y-room_y_w/2)), aim_wx, aim_wy].*affine, 'EdgeColor', 'k', 'LineWidth', 2);
% hold off;

% figure
% pixel=5*size(room_finall_dbm);
% affine=pixel(1)/room_x_w;
% finall_data = imresize(flipud(room_finall_dbm), pixel, 'lanczos3');
% imagesc(finall_data);
% xticks(linspace(1,pixel(1),21));
% yticks(linspace(1,pixel(2),21));
% xticklabels(room_x - room_x_w/2:room_x + room_x_w/2);
% yticklabels(room_y + room_x_w/2:-1:room_x - room_x_w/2);
% colorbar
% grid off
% caxis([-55 -33]);
% % rectangle('Position', [aim_x-aim_wx/2-(room_x-room_x_w/2), room_y_w-(aim_y+aim_wy/2-(room_y-room_y_w/2)), aim_wx, aim_wy].*affine, 'EdgeColor', 'k', 'LineWidth', 2);
% hold off;

size(room_finall_bps)
pixel=[1000 1000];
affine=pixel(1)/room_x_w;
finall_data = imresize(room_finall_bps, pixel, 'lanczos3');

%%
% finall_data=finall_data*0.98;

figure
hold on
imagesc(finall_data);
xticks(linspace(1,pixel(1),21));
yticks(linspace(1,pixel(2),21));
xticklabels(room_x - room_x_w/2:room_x + room_x_w/2);
yticklabels(room_y + room_x_w/2:-1:room_x - room_x_w/2);
axis([0, 1000, 0, 1000]);
colorbar
grid off
colormap('jet');
my_handle=colorbar;
my_handle.Label.String = 'Spectrum Effectiveness (bps/Hz)';
my_handle.Label.FontSize = 16;
my_handle.FontSize = 16;
my_handle.FontName = 'Times New Roman';
my_handle.Label.FontName ='Times New Roman';
caxis([0 20]);
xticks([1 250 500 750 1000]);
xticklabels({'0', '5', '10', '15', '20'});
yticks([1 250 500 750 1000]);
yticklabels({'0', '5', '10', '15', '20'});
xlabel('x-axis(m)')
ylabel('y-axis(m)')
ax = gca;
ax.XLabel.FontSize = 16;
ax.YLabel.FontSize = 16;
ax.XLabel.FontName = 'Times New Roman';
ax.YLabel.FontName = 'Times New Roman';
% Theta=linspace(theta1, theta2, 50);%30
% X1 = r1 .* cos(Theta) + RIS1_x; % x坐标
% Y1 = r1 .* sin(Theta) + RIS1_y; % y坐标
% x_h1=X1(:);
% y_h1=Y1(:);
% plot((x_h1-5)./20.*1000,(y_h1-5)./20.*1000, 'k-', 'LineWidth', 2);
% X2 = r2 .* cos(Theta) + RIS1_x; % x坐标
% Y2 = r2 .* sin(Theta) + RIS1_y; % y坐标
% x_h2=X2(:);
% y_h2=Y2(:);
% plot((x_h2-5)./20.*1000,(y_h2-5)./20.*1000, 'k-', 'LineWidth', 2);
% radii = linspace(r1, r2, 30);%30
% X3 = radii .* cos(theta1) + RIS1_x; % x坐标
% Y3 = radii .* sin(theta1) + RIS1_y; % y坐标
% x_h3=X3(:);
% y_h3=Y3(:);
% plot((x_h3-5)./20.*1000,(y_h3-5)./20.*1000, 'k-', 'LineWidth', 2);
% radii = linspace(r1, r2, 60);%30
% X4 = radii  .* cos(theta2) + RIS1_x; % x坐标
% Y4 = radii  .* sin(theta2) + RIS1_y; % y坐标
% x_h4=X4(:);
% y_h4=Y4(:);
% plot((x_h4-5)./20.*1000,(y_h4-5)./20.*1000, 'k-', 'LineWidth', 2);
% hold off
% set(gca,'FontName','Times New Rome','FontSize',12);

% ax = gca;
% ax.XAxis.TickLabel.FontSize = 2;
% ax.XAxis.TickLabel.FontName = 'Times New Roman';
% ax.YAxis.TickLabel.FontSize = 22;
% ax.YAxis.TickLabel.FontName = 'Times New Roman';
%rectangle('Position', [aim_x-aim_wx/2-(room_x-room_x_w/2), room_y_w-(aim_y+aim_wy/2-(room_y-room_y_w/2)), aim_wx, aim_wy].*affine, 'EdgeColor', 'k', 'LineWidth', 2);
%rectangle('Position', [3.4, 3.4, 3, 3].*affine, 'EdgeColor', 'k', 'LineWidth', 2);
%rectangle('Position', [13.5, 13.5, 3, 3].*affine, 'EdgeColor', 'k', 'LineWidth', 2);

% line([7.5 7.5].* affine, [10, 11].* affine,'color','k','LineWidth', 2);
% line([7.5 12.5].* affine, [11, 11].* affine,'color','k','LineWidth', 2);
% line([12.5 12.5].* affine, [10, 11].* affine,'color','k','LineWidth', 2);
% line([7.5 9.5].* affine, [10, 10].* affine,'color','k','LineWidth', 2);
% line([10.5 12.5].* affine, [10, 10].* affine,'color','k','LineWidth', 2);
% line([9.5 10.5].* affine, [7, 7].* affine,'color','k','LineWidth', 2);
% line([9.5 9.5].* affine, [7, 10].* affine,'color','k','LineWidth', 2);
% line([10.5 10.5].* affine, [7, 10].* affine,'color','k','LineWidth', 2);

line([8 8].* affine, [7, 12].* affine,'color','k','LineWidth', 2);
line([9 9].* affine, [8 12].* affine,'color','k','LineWidth', 2);
line([8 12].* affine, [7, 7].* affine,'color','k','LineWidth', 2);
line([9 12].* affine, [8, 8].* affine,'color','k','LineWidth', 2);
line([12 12].* affine, [7, 8].* affine,'color','k','LineWidth', 2);
line([8 9].* affine, [12, 12].* affine,'color','k','LineWidth', 2);