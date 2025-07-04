Gamma=Gamma_initial;
% 计算目标方向单位向量
aim_vector = [aim_x - RIS1_x, aim_y - RIS1_y, aim_z - RIS1_z];
norm_aim_vector = sqrt(sum(aim_vector.^2));
u = aim_vector / norm_aim_vector;

% 计算每个单元到阵列中心的向量并计算相位
phi_dft = zeros(N_z, N_y);
for i = 1:N_z
    for j = 1:N_y
        r_ij = [x_cell(i,j) - RIS1_x, y_cell(i,j) - RIS1_y, z_cell(i,j) - RIS1_z];
        % 计算相位差
        phi_dft(i, j) = -2 * pi * dot(u, r_ij) / lambda;
    end
end


%con
Gamma=A_T.*exp(-1j*reshape(r_BS_cell,[N_z,N_y])*2*pi/lambda).*exp(-1j*phi_dft);

%bit
% teco=wrapTo2Pi(-reshape(r_BS_cell,[N_z,N_y])*2*pi/lambda+-phi_dft);%最大
% for ci=1:size(teco,1)
%     for cj=1:size(teco,2)
%         index1=find(abs(teco(ci,cj)-Qphi)<=2*pi/2^qbit/2);
%         ph(ci,cj)=Qphi(index1);
%     end
% end
% Gamma=A_T.*exp(1j*ph);


Gamma_dft=Gamma(:)';
[P_rx,P_rx_mean,P_rx_dbm,P_rx_mean_dbm]=rx_power(r_BS_cell_exp,r_cell_aim_exp,G,G_rx,Gamma_dft,lambda,F,P_tx_exp,G_tx_exp);
R=log2(1+P_rx/noise);

% Save_com(:,count)=R;