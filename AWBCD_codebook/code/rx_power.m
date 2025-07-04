function [P_rx,P_rx_mean,P_rx_dbm,P_rx_mean_dbm]=rx_power(r_BS_cell_exp,r_cell_aim_exp,G,G_rx,Gamma,lambda,F,P_tx_exp,G_tx_exp)

global C M N dy dz


te = sqrt(P_tx_exp .* G_tx_exp .* F) .* Gamma ./ (r_BS_cell_exp .* r_cell_aim_exp).* exp(-1j * 2 * pi ./ lambda .* (r_BS_cell_exp + r_cell_aim_exp));
P_rx = G*G_rx*dy*dz*lambda^2/(64*pi^3)*abs(sum(sum(te,1),2)).^2;
P_rx=permute(P_rx,[3,2,1]);
P_rx_mean=mean(P_rx);
P_rx_dbm=pow2db(P_rx*10^3);
P_rx_mean_dbm=mean(P_rx_dbm);