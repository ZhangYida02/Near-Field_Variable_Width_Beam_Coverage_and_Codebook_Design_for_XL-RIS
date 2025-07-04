r_thr1=17.5; 
r_thr2=17; 
r_thr3=16.5;
r_thr4=16; 
%%——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
SD2_save=[];
R_save=[];
pcount1_save=[];
pcount2_save=[];
pcount3_save=[];
pcount4_save=[];
%% __________________________________________________________________________________________________________
%初始平均接收功率
[P_rx_initial,P_rx_initial_mean,P_rx_initial_dbm,P_rx_initial_mean_dbm]=rx_power(r_BS_cell_exp,r_cell_aim_exp,G,G_rx,Gamma_initial,lambda,F,P_tx_exp,G_tx_exp);


R=log2(1+P_rx_initial/noise);

%存
logical_vector = R < r_thr1;
pcount1_save = [pcount1_save sum(logical_vector)/length(R)];
logical_vector = R < r_thr2;
pcount2_save = [pcount2_save sum(logical_vector)/length(R)];
logical_vector = R < r_thr3;
pcount3_save = [pcount3_save sum(logical_vector)/length(R)];
logical_vector = R < r_thr4;
pcount4_save = [pcount4_save sum(logical_vector)/length(R)]

R_save=[R_save mean(R)]

%% __________________________________________________________________________________________________________
% 
%计算all_PQab
phi_b_exp=repmat(phi, [M, 1, C]);
exp_component = exp(1j * (-2 * pi / lambda * (r_BS_cell_exp + r_cell_aim_exp) + phi_b_exp));
te_matrix = A .* exp_component; 
te_sum = sum(te_matrix, 1); 
all_PQab= permute(sum(te_sum, 2),[1,3,2]);
toc
Weight=[];
%% __________________________________________________________________________________________________________
%优化
bar= waitbar(0, '请等待...');
max_u=50;
P=zeros(N,C);
Q=zeros(N,C);
U=zeros(N,C);
V=zeros(N,C);
beta=zeros(N,C);
alpha=zeros(N,C);
weight=100*ones(1,C);
Gamma=Gamma_initial;
gama=0.99;

for u1=1:max_u
    for k=1:N

        if mod(k,100)==0
            waitbar((k+(u1-1)*N)/(max_u*N), bar, [num2str(100*(k+(u1-1)*N)/(max_u*N)),'%']);
        end

        %求P beta
        exp_component= exp(1j * (-2 * pi / lambda * (r_BS_cell(:, k) + r_cell_aim(k, :))));
        A_slice= squeeze(A(:, k, :));
        Su=sum(A_slice);
        all_Pb=Su*exp_component;

        %求Q alpha
        all_Qa=all_PQab-all_Pb*exp(1j*phi(k));
        P(k,:)=abs(all_Pb);
        beta(k,:)=wrapTo2Pi(angle(all_Pb));
        Q(k,:)=abs(all_Qa);
        alpha(k,:)=wrapTo2Pi(angle(all_Qa));

        %U V
        U(k,:)=P(k,:).^2+Q(k,:).^2;
        V(k,:)=2.*P(k,:).*Q(k,:);

        %计算Phi_b Phi_a       
        Phi_a = -sum(weight.*V(k,:) ./ U(k,:) .* sin(beta(k, :) - alpha(k, :)));
        Phi_b = sum(weight.*V(k,:) ./ U(k,:) .* cos(beta(k, :) - alpha(k, :)));

        % Phi_a = -sum(weight.*Q(k,:) ./ P(k,:) .* sin(beta(k, :) - alpha(k, :)));
        % Phi_b = sum(weight.*Q(k,:) ./ P(k,:) .* cos(beta(k, :) - alpha(k, :)));

        %调相增强
        target='enhance';
        if Phi_a>=0

            if bit~=0
                te=wrapTo2Pi(pi/2-atan(Phi_b/Phi_a));%最大
                index1=find(abs(te-Qphi)<=2*pi/2^bit/2);
                phi(k)=Qphi(index1);
            else
                phi(k)=wrapTo2Pi(pi/2-atan(Phi_b/Phi_a));%最大
            end

        else
            if bit~=0
                te=wrapTo2Pi(-pi/2-atan(Phi_b/Phi_a));%最大
                index1=find(abs(te-Qphi)<=2*pi/2^bit/2);
                phi(k)=Qphi(index1);
            else
                phi(k)=wrapTo2Pi(-pi/2-atan(Phi_b/Phi_a));%最大
            end

        end
        %
        %调相减弱
        % target='recede';
        % if Phi_a>=0
        %    phi(k)=3*pi/2-atan(Phi_b/Phi_a);%最小
        % else
        %    phi(k)=pi/2-atan(Phi_b/Phi_a);%最小
        % end

        %更新Gamma
        Gamma(k)=A_T*exp(1j*phi(k));

        %更新all_PQab
        all_PQab=all_Qa+all_Pb*exp(1j*phi(k));

        % 计算过程功率


    end
    if mod(u1,1)==0
        [P_rx,P_rx_mean,P_rx_dbm,P_rx_mean_dbm]=rx_power(r_BS_cell_exp,r_cell_aim_exp,G,G_rx,Gamma,lambda,F,P_tx_exp,G_tx_exp);
        R=log2(1+P_rx/noise);

        change=double((mean(R)-R) > 0);
        weight=weight+change';
        Weight=[Weight;weight];

        logical_vector = R < r_thr1;
        pcount1_save = [pcount1_save sum(logical_vector)/length(R)];
        logical_vector = R < r_thr2;
        pcount2_save = [pcount2_save sum(logical_vector)/length(R)];
        logical_vector = R < r_thr3;
        pcount3_save = [pcount3_save sum(logical_vector)/length(R)];
        logical_vector = R < r_thr4;
        pcount4_save = [pcount4_save sum(logical_vector)/length(R)];

        R_save=[R_save mean(R)];
        SD2_save=[SD2_save sqrt(var(R))];


    end
    R_save=[R_save mean(log2(1+P_rx/noise))];
end

close(bar);
clear all_PQab all_Pb all_Qa
Gamma_iteration_weight=Gamma;

plot(1:length(R_save),R_save,'-*')