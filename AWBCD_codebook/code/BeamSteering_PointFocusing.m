Gamma=Gamma_initial;



d=sqrt((x_cell-aim_x).^2+(y_cell-aim_y).^2+(z_cell-aim_z).^2)+reshape(r_BS_cell,[N_z,N_y]);
Gamma=A_T.*exp(-1j*d*2*pi/lambda);
Gamma_dian=Gamma(:)';
[P_rx,P_rx_mean,P_rx_dbm,P_rx_mean_dbm]=rx_power(r_BS_cell_exp,r_cell_aim_exp,G,G_rx,Gamma_dian,lambda,F,P_tx_exp,G_tx_exp);
R=log2(1+P_rx/noise);
% Save_com(:,count)=R;