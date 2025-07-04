figure
[n, m] = size(Gamma_iteration);
for s=[0]
    hold on

% 

    [f_initial, x_initial] = ksdensity(log2(1+P_rx_initial/noise),'Width',0.1);
    cdf_values_f_initial = cumsum(f_initial) / sum(f_initial);
    plot(x_initial,cdf_values_f_initial,'-k',LineWidth=1.5);
%     %计算功率
% 
    [P_rx_finall,P_rx_finall_mean,P_rx_finall_mean_dbm]=rx_power(r_BS_cell_exp,r_cell_aim_exp,G,G_rx,Gamma_iteration,lambda,F,P_tx_exp,G_tx_exp);

    mean_iteration_bps=mean(log2(1+P_rx_finall/noise))
    [f_finall, x_finall] = ksdensity(log2(1+P_rx_finall/noise),'Width',0.1);
    cdf_values_f_finall = cumsum(f_finall) / sum(f_finall);
    plot(x_finall,cdf_values_f_finall,'-r',LineWidth=1.5);
% 
    [P_rx_finall,P_rx_finall_mean,P_rx_finall_mean_dbm]=rx_power(r_BS_cell_exp,r_cell_aim_exp,G,G_rx,Gamma_qiu,lambda,F,P_tx_exp,G_tx_exp);

    mean_qiu_bps=mean(log2(1+P_rx_finall/noise))
    [f_finall, x_finall] = ksdensity(log2(1+P_rx_finall/noise),'Width',0.1);
    cdf_values_f_finall = cumsum(f_finall) / sum(f_finall);
    plot(x_finall,cdf_values_f_finall,'b',LineWidth=1.5);
    % 
    [P_rx_finall,P_rx_finall_mean,P_rx_finall_mean_dbm]=rx_power(r_BS_cell_exp,r_cell_aim_exp,G,G_rx,Gamma_dian,lambda,F,P_tx_exp,G_tx_exp);

    mean_dian_bps=mean(log2(1+P_rx_finall/noise))
    [f_finall, x_finall] = ksdensity(log2(1+P_rx_finall/noise),'Width',0.1);
    cdf_values_f_finall = cumsum(f_finall) / sum(f_finall);
    plot(x_finall,cdf_values_f_finall,'g',LineWidth=1.5);
    % 
    [P_rx_finall,P_rx_finall_mean,P_rx_finall_mean_dbm]=rx_power(r_BS_cell_exp,r_cell_aim_exp,G,G_rx,Gamma_dft,lambda,F,P_tx_exp,G_tx_exp);
TE=log2(1+P_rx_finall/noise)
    mean_dian_bps=mean(log2(1+P_rx_finall/noise))
    [f_finall, x_finall] = ksdensity(log2(1+P_rx_finall/noise),'Width',0.1);
    cdf_values_f_finall = cumsum(f_finall) / sum(f_finall);
    plot(x_finall,cdf_values_f_finall,'cyan',LineWidth=1.5);

    % [P_rx_finall,P_rx_finall_mean,P_rx_finall_mean_dbm]=rx_power(r_BS_cell_exp,r_cell_aim_exp,G,G_rx,Gamma_iteration_weight,lambda,F,P_tx_exp,G_tx_exp);
    % mean_iteration_bps=mean(log2(1+P_rx_finall/noise))
    % [f_finall, x_finall] = ksdensity(log2(1+P_rx_finall/noise),'Width',0.1);
    % cdf_values_f_finall = cumsum(f_finall) / sum(f_finall);
    % plot(x_finall,cdf_values_f_finall,'-m',LineWidth=1.5);

    hold off

    ylim([0 1])
    grid on
end


figure

hold on

[f_finall, x_finall] = ksdensity(R(:,1),'Width',0.1);
cdf_values_f_finall = cumsum(f_finall) / sum(f_finall);
plot(x_finall,cdf_values_f_finall,LineWidth=1.5);

[f_finall, x_finall] = ksdensity(R(:,2),'Width',0.1);
cdf_values_f_finall = cumsum(f_finall) / sum(f_finall);
plot(x_finall,cdf_values_f_finall,LineWidth=1.5);
% 
[f_finall, x_finall] = ksdensity(R(:,3),'Width',0.1);
cdf_values_f_finall = cumsum(f_finall) / sum(f_finall);
plot(x_finall,cdf_values_f_finall,LineWidth=1.5);
% 
[f_finall, x_finall] = ksdensity(R(:,4),'Width',0.1);
cdf_values_f_finall = cumsum(f_finall) / sum(f_finall);
plot(x_finall,cdf_values_f_finall,LineWidth=1.5);

[f_finall, x_finall] = ksdensity(R(:,5),'Width',0.1);
cdf_values_f_finall = cumsum(f_finall) / sum(f_finall);
plot(x_finall,cdf_values_f_finall,LineWidth=1.5);

[f_finall, x_finall] = ksdensity(R(:,6),'Width',0.1);
cdf_values_f_finall = cumsum(f_finall) / sum(f_finall);
plot(x_finall,cdf_values_f_finall,LineWidth=1.5);

hold off









% figure
% for s=[0]
%     hold on
% 
%     [f_initial, x_initial] = ksdensity(pow2db(10^3*P_rx_initial),'Width',0.2);
%     cdf_values_f_initial = cumsum(f_initial) / sum(f_initial);
%     plot(x_initial,cdf_values_f_initial,'-k',LineWidth=1.5);
%     % 计算功率
%     [P_rx_finall,P_rx_finall_mean,P_rx_finall_dbm,P_rx_finall_mean_dbm]=rx_power(r_BS_cell_exp,r_cell_aim_exp,G,G_rx,Gamma_iteration,lambda,F,P_tx_exp,G_tx_exp);
%     mean_iteration_dbm=mean(pow2db(10^3*P_rx_finall))
%     [f_finall, x_finall] = ksdensity(pow2db(10^3*P_rx_finall),'Width',0.1);
%     cdf_values_f_finall = cumsum(f_finall) / sum(f_finall);
%     plot(x_finall,cdf_values_f_finall,'r',LineWidth=1.5);
% 
%     % 
%     [P_rx_finall,P_rx_finall_mean,P_rx_finall_mean_dbm]=rx_power(r_BS_cell_exp,r_cell_aim_exp,G,G_rx,Gamma_qiu,lambda,F,P_tx_exp,G_tx_exp);
%     mean_qiu_dbm=mean(pow2db(10^3*P_rx_finall))
%     [f_finall, x_finall] = ksdensity(pow2db(10^3*P_rx_finall),'Width',0.1);
%     cdf_values_f_finall = cumsum(f_finall) / sum(f_finall);
%     plot(x_finall,cdf_values_f_finall,'b',LineWidth=1.5);
%     % 
%     [P_rx_finall,P_rx_finall_mean,P_rx_finall_mean_dbm]=rx_power(r_BS_cell_exp,r_cell_aim_exp,G,G_rx,Gamma_dian,lambda,F,P_tx_exp,G_tx_exp);
%     mean_dian_dbm=mean(pow2db(10^3*P_rx_finall))
%     [f_finall, x_finall] = ksdensity(pow2db(10^3*P_rx_finall),'Width',0.2);
%     cdf_values_f_finall = cumsum(f_finall) / sum(f_finall);
%     plot(x_finall,cdf_values_f_finall,'g',LineWidth=1.5);
% 
%     [P_rx_finall,P_rx_finall_mean,P_rx_finall_mean_dbm]=rx_power(r_BS_cell_exp,r_cell_aim_exp,G,G_rx,Gamma_dft,lambda,F,P_tx_exp,G_tx_exp);
%     mean_dian_dbm=mean(pow2db(10^3*P_rx_finall))
%     [f_finall, x_finall] = ksdensity(pow2db(10^3*P_rx_finall),'Width',0.2);
%     cdf_values_f_finall = cumsum(f_finall) / sum(f_finall);
%     plot(x_finall,cdf_values_f_finall,'cyan',LineWidth=1.5);
% 
%     % [P_rx_finall,P_rx_finall_mean,P_rx_finall_dbm,P_rx_finall_mean_dbm]=rx_power(r_BS_cell_exp,r_cell_aim_exp,G,G_rx,Gamma_iteration_weight,lambda,F,P_tx_exp,G_tx_exp);
%     % mean_iteration_w_dbm=mean(pow2db(10^3*P_rx_finall))
%     % [f_finall, x_finall] = ksdensity(pow2db(10^3*P_rx_finall),'Width',0.2);
%     % cdf_values_f_finall = cumsum(f_finall) / sum(f_finall);
%     % plot(x_finall,cdf_values_f_finall,'m',LineWidth=1.5);
% 
%     hold off
%     % xlim([-50 -30])
%     % legend('random','iteration','sphere mapping','point focusing','DFT codebook','iteration_weight')
%     xlabel('equivalent power (dbm)')
%     ylabel('probability')
%     ylim([0 1])
%     grid on
% end


