% **********************************************
% This script launches BOATS in Monte Carlo mode
% **********************************************
clear all

addpath('mc/');
parameters

% Set parameter distributions for undetermined distributions
%boats.param.main.MC_param_Pel=set_MC('Pelagic');
%boats.param.main.MC_param_Dem=set_MC('Demersal');
%param_loops=boats.param.main.nrun;
load('best_param_164.mat')
param_loops = 5;

% Loop ensemble of simulations -> spin-up + transient
for indr = 1:param_loops
    ens_param = [];
    % Parameter selection
%    if (boats.param.ecology.demersal)&&(boats.param.ecology.pelagic)
%        ens_param=load_MC(boats.param.main.MC_param_Pel,'Pelagic' , 1,ens_param);
%        ens_param=load_MC(boats.param.main.MC_param_Dem,'Demersal', 1,ens_param);
%    else
%        ens_param=load_MC(boats.param.main.MC_param_Pel,'Pelagic',[]);
%    end
    ens_param = best_param(indr);
    % run
    run_boats('frc/Ecological_depth_dist_Zeuph.mat','frc/Economical_LME_q5.mat',ens_param,'PP','nh','final')
    run_boats('frc/Ecological_depth_dist_Zeuph.mat','frc/Economical_LME_q5.mat',ens_param,'restart','h','snap10year')
end 
