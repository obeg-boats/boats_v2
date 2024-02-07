% **********************************************
% This script launches BOATS in Monte Carlo mode
% **********************************************
clear all

addpath('mc/');

% % Set parameter distributions for undetermined distributions %%%%%%%%%%%%%
% % FOR MONTE CARLO ENSEMBLE
% parameters
% boats.param.main.MC_param_Pel=set_MC('Pelagic');
% boats.param.main.MC_param_Dem=set_MC('Demersal');
% param_loops=boats.param.main.nrun;
% FOR 5 BEST ENSEMBLES
load('best_param_5v2.mat')
param_loops = 5;

% Loop ensemble of simulations -> spin-up + transient
for indr = 1:param_loops
    ens_param = [];
    % Parameter selection %%%%%%%%%%%%%%%%%%%%%%%
%     % FOR MONTE CARLO ENSEMBLE
%     if (boats.param.ecology.demersal)&&(boats.param.ecology.pelagic)
%         ens_param=load_MC(boats.param.main.MC_param_Pel,'Pelagic' , 1,ens_param);
%         ens_param=load_MC(boats.param.main.MC_param_Dem,'Demersal', 1,ens_param);
%     else
%         ens_param=load_MC(boats.param.main.MC_param_Pel,'Pelagic',[]);
%     end
    % FOR 5 BEST ENSEMBLES
    ens_param = best_param(indr);
    
    % run %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    run_boats('frc/Ecological.mat','frc/Economical.mat',ens_param,'PP','nh','final')
    run_boats('frc/Ecological.mat','frc/Economical.mat',ens_param,'restart','hd','annual')
end 
