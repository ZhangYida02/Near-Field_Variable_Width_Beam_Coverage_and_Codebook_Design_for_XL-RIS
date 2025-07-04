clear
clc

e_save=[];
p_save=[];
count_bit=0;

tic
for bit=[0]
    count_bit=count_bit+1
    count_nn=0;
    for nn=[20]

        
        count_nn=count_nn+1

        clearvars -except nn e_save p_save bit count_bit count_nn
        close all
        rng(1)

        Environment

         BeamSteering_Iteration
        % Beamsteering_WeightUpdateIteration
        % BeamSteering_SphericalMapping
        % BeamSteering_PointFocusing
        % BeamSteering_DFT
        % BeamSteering_Compare

        % e_save(count_bit,count_nn)=meanR1_save(end);
        % p_save(count_bit,count_nn)=pcount4_save(end);
        % 
        % e_save
        % p_save
    end
end
toc
% figure
% hold on
% plot([0:size(e_save,2)-1],p_save,'-*',Color='b')
% plot([0:length(pcount2_save)-1],pcount2_save,'-*',Color='c')
% plot([0:length(pcount3_save)-1],pcount3_save,'-*',Color='y')
% plot([0:length(pcount4_save)-1],pcount4_save,'-*',Color='k')
% hold off



