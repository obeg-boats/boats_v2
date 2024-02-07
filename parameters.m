%**************************************************************************
% BOATS PARAMETERS
% Default parameters and tunable variables, paths, of the BOATS model
%**************************************************************************
 

%**************************************************************************
% MAIN TUNABLE PARAMETERS
%**************************************************************************
% Paths *******************************************
 boats.param.path.wrkdir = [''];
 boats.param.path.outdir = ['output/'];
% Names and Switches ******************************
 boats.param.main.sim_type     = 'nh';                                     % No Harvest 'nh' or Harvest simulations 'h'
 boats.param.main.sim_init     = 'PP';                                     % Initialisation from 'PP' or else 'restart'
 boats.param.main.sim_name     = 'Boats_V2';                               % Simulation name
 boats.param.main.model_version= 'V2';                                     % Model version
 boats.param.main.save_restart = 1;                                        % Save restart: yes=1 ; no=0
 boats.param.main.save_output  = 1;                                        % Save output: yes=1 ; no=0
% Simulation features *****************************
 boats.param.main.run_length   = 300;                                      % Simulation length in years 
 boats.param.main.dtt          = 30;                                       % days per timestep
 boats.param.main.nforcing     = 12;                                       % number of forcing to loop
 boats.param.main.param_ens    = 1;                                        % Use parameters ensembles: yes=1 ; no=0
 boats.param.main.dataset_ens  = 'ensemble_parameters.mat';                % if param_ens=1 name of ensemble parameters
% Monte Carlo ensemble ***************************
 boats.param.main.nrun         = 1000;                                       % Number of Replicates for ensemble of Monte Carlo simulations
%**************************************************************************
% END MAIN TUNABLE PARAMETERS
%**************************************************************************
 

%**************************************************************************
% CONVERSION FACTORS & EPSILON
%**************************************************************************
 boats.param.conversion.sperd        = 3600*24;                            % seconds per day
 boats.param.conversion.spery        = boats.param.conversion.sperd*360;   % seconds per year
 boats.param.conversion.sperfrc      = boats.param.conversion.sperd...
                                       *boats.param.main.dtt...
                                       *boats.param.main.nforcing;         % seconds per forcing
 boats.param.conversion.gC_2_wetB    = 10;                                 % grams of wet fish biomass per gram of fish carbon
 boats.param.conversion.mmolC_2_wetB = (12*boats.param.conversion.gC_2_wetB)/1000; % grams of wet fish biomass per mmol of fish carbon
 boats.param.conversion.C_2_K        = 273.15;                             % deg C to Kelvin
 boats.param.conversion.epsln        = 1e-50;                              % small epsilon
 boats.param.main.dtts               = 30*24*3600;                         % seconds per timestep

 
%**************************************************************************
% PARAMS RELATED TO THE ENVIRONMENT
%**************************************************************************
% Temperature *************************************
 boats.param.environment.E_activation_A = [0.45 0.45];                            % Activation energy of metabolism (growth A) (eV) (Savage et al., 2004)
 boats.param.environment.E_activation_m = [0.45 0.45];                            % Activation energy of metabolism (mortality) (eV) (Savage et al., 2004)
 boats.param.environment.k_Boltzmann    = 8.617e-05;                       % Boltzmann Constant (eV K-1)
 boats.param.environment.temp_ref_A     = 10 + 273.15;                     % Reference temperature (K) (Andersen and Beyer, 2013, p. 18)
% Primary production ******************************
 boats.param.environment.kappa_eppley = [0.06 0];                             % Eppley constant (degC-1)
 boats.param.environment.Prod_star    = [0.35 0];                              % Pivotal primary production (m mol C m-3 d-1)
 boats.param.environment.mc_phy_l     = 5.6234132519e-06;                  % mass of typical large phytoplankton (g)
 boats.param.environment.mc_phy_s     = 5.6234132519e-15;                  % mass of typical small phytoplankton (g)
 boats.param.environment.cap_npp      = 10000;                             % limit on npp (m mol C m-2 d-)
 boats.param.environment.mc_benthic   = 1/5000*0.02;                  % mass of typical benthic organism (b) 
 
%**************************************************************************
% PARAMS RELATED TO THE ECOLOGICAL MODULE 
%**************************************************************************
 boats.param.ecology.pelagic     = 1;                                      % Pelagic community (0 or 1)
 boats.param.ecology.demersal    = 1;					   % Demersal community (0 or 1)
% Spectrum ****************************************
 boats.param.ecology.te          = [0.125 0.125];                                  % trophic efficiency
 boats.param.ecology.ppmr        = [5000 5000];                                   % predator to prey mass ratio
 boats.param.ecology.tro_sca     = log10(boats.param.ecology.te)./log10(boats.param.ecology.ppmr); % trophic scaling
 boats.param.ecology.b_allo      = [0.66 0.66];                                   % allometric scaling
 boats.param.ecology.zeta1       = [0.57 0.57];                                   % constant mortality scaling
 boats.param.ecology.h_allo      = [0.5 0.5];                                    % mass scaling of mortality
 boats.param.ecology.eff_a       = 0.8;                                    % efficiency of activity (Andersen and Beyer, 2013, p. 4)
 boats.param.ecology.A00         = [4.46 4.46];                                   % allometric growth rate (Andersen and Beyer, 2013, p. 4)
 boats.param.ecology.fmass_0 	 = 10;                                     % initial mass class (g)
 boats.param.ecology.fmass_e 	 = 1e5;                                    % final mass class (g)
% Reproduction ************************************
 boats.param.ecology.m_egg       = 5.2e-4;                                 % egg mass (g)
 boats.param.ecology.frac_fem    = 0.5;                                    % fraction of individuals that allocate energy to reproduction (females)
 boats.param.ecology.egg_surv    = [0.01 0.01];                                   % egg survival
 boats.param.ecology.rep_slope   = 5;                                      % slope parameter of sigmoidal allocation to reproduction function
 boats.param.ecology.rep_pos     = 1;                                      % position parameter of sigmoidal allocation to reproduction function as fraction of malpha
% Size class **************************************
 boats.param.ecology.nfmass 	 = 50;                                     % number of fish mass classes
 boats.param.ecology.minf        = [0.01*(30/0.95)^3 0.01*(90/0.95)^3 1e5];% asymptotic mass
 boats.param.ecology.eta_alpha   = 0.25;						           % mass at maturity as fraction of asymptotic mass 
 boats.param.ecology.malpha      = boats.param.ecology.eta_alpha*boats.param.ecology.minf; % maturity mass
 boats.param.ecology.nfish       = length(boats.param.ecology.minf);       % number of fish groups
% Iron Limitation *********************************
 boats.param.ecology.kfe         = 5;                                      % coefficient determining sensitivity of fish biomass to iron (No3) concentration
                                                                           % set high (1e6) to have no limitation!

%**************************************************************************
% PARAMS RELATED TO THE ECONOMICAL MODULE 
%**************************************************************************
 boats.param.economy.landedvalue_global = 8.4233e+10;                      % SAUP 1990-2006 average ($)
 boats.param.economy.yield_global       = 7.9963e+13;                      % SAUP 1990-2006 average (g)
 boats.param.economy.price_global       = boats.param.economy.landedvalue_global/boats.param.economy.yield_global;       % Global price ($ g-1)
 boats.param.economy.cost_global        = boats.param.economy.landedvalue_global;    % Global total cost ($) Assume C = R (matches Lam et al., 2011)
 boats.param.economy.effort_global      = 14.6229e9;                       % Global effort (W) (Watson et al. (2012) 1990-2006 average)
 boats.param.economy.cost_effort_0      = boats.param.economy.cost_global/(boats.param.economy.effort_global*boats.param.conversion.spery); % Cost per unit effort ($ W-1)
 boats.param.economy.k_e                = 1e-6;                            % Fleet dynamic parameter (W $-1 s-1) 1 W of effort per dollar of net revenue per second
 boats.param.economy.sel_pos_1          = 1;                               % Selectivity position shift 1
 boats.param.economy.sel_pos_2          = 1;%0.5;                             % Selectivity position shift 2
 boats.param.economy.sel_pos_3          = 1;%0.25;                            % Selectivity position shift 3
 boats.param.economy.sel_pos_scale      = [1 1];                               % Selectivity position scale
 boats.param.economy.sel_slope          = [17.8 17.8];                              % Selectivity slope
 boats.param.economy.harvest_start      = 0;                               % Year of starting harvest [y]
 boats.param.economy.qcatch0            = 1e-5 * 1.05^(-100);              % Base catchability
 boats.param.economy.price_0            = boats.param.economy.price_global;% Base price (constant)
 boats.param.economy.rc_q_ini           = 0;                               % Discount rate of catchability to determine initial value
 boats.param.economy.q_discount_y       = 0;                               % Number of years to discount catchability
 boats.param.economy.costpar		= 1; 				   % 0 or 1 for depth dependent costs
 boats.param.economy.catchpar           = 1;                               % 0 or 1 for depth dependent catchabilities

%**************************************************************************************************************
% END OF SCRIPT


 
