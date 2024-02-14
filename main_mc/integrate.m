%**************************************************************************************************************
% FUNCTION integrate.m
%**************************************************************************************************************

%-----------------------------------------------------------------------------------------
% Authors
%-----------------------------------------------------------------------------------------

% David A. Carozza  (david.carozza@gmail.com)
% Daniele Bianchi   (dbianchi@atmos.ucla.edu)
% Eric D. Galbraith (eric.d.galbraith@gmail.com)

%-----------------------------------------------------------------------------------------
% Introduction
%-----------------------------------------------------------------------------------------

% Bioeconomic Open-Access Trophic Size-Spectrum (BOATS) model based on the
% McKendrick-von Foerster model for a mass-spectrum of fish biomass
% with an open access economic framework to calculate effort and harvest.
% Forced with monthly primary production and temperature data (from model or observations)

% Primary production in units of mmolC m-2 s-1
% Core unit of biomass is mmolC for fish and harvest plots in model
% Convert these to grams of wet biomass (g wetB) using mmolC_2_wetB
% Time in seconds
% Time step is 30 days (1 month)

% 3 fish groups defined by different asymptotic masses
% 1 effort group and 1 selectivity function per fish group
% 2-dimensional version of the code

%-----------------------------------------------------------------------------------------
% Other comments
%-----------------------------------------------------------------------------------------

% We set a lower limit on effort by adding a small constant (epsln = 1e-15)
% to each application of effort when calculating harvest, cost, and effort change.
% This guarantees that as effort approaches zero, the effort change equation approaches 
% the analytical simplification. This prevents dividing by zero.
% We cannot use the analytical simplification of the change in effort directly because we
% limit harvest and we need harvest directly to calculate revenue.

% dbianchi 04/16/2016 : output routine changed to allow for general and flexible 
% output "modes". Each mode specifies the variables to be saved, the processing of variables
% (e.g. averages, integrals) and the interval time used for averaging variables (e.g. annual,
% decadal, final timestep, user defined etc.) Additionally, the main calculations were 
% optimized for fast array calculations using "bsxfun", and tested with matlab's built-in 
% profiler.

%**************************************************************************************************************
% MAIN CODE
%**************************************************************************************************************
function boats = integrate(boats)

 % start timer
 tic

 %---------------------------------
 % Aliases for variables category for readability
 MAIN=boats.param.main;
 CONV=boats.param.conversion;
 ENVI=boats.param.environment;
 ECOL=boats.param.ecology;
 ECON=boats.param.economy;
 FORC=boats.forcing;
 STRU=boats.structure;
 INIT=boats.initial; 
 
 %--------------------------------------------------------------------------------------------------------------
 % Defines Ecologic simulations
 %--------------------------------------------------------------------------------------------------------------
 
 %-------------------------------------------------------------------------------------
 % Time integration parameters
 time                   = [MAIN.dtts:MAIN.dtts:MAIN.run_length*CONV.spery]';		   % seconds
 ntime                  = length(time);

 %-----------------------------------------------------------------------------------------
 % Update ecological parameters
 A0                     = ECOL.A00/CONV.spery;                             % growth rate per second (Andersen and Beyer, 2013, p. 4)

 %-----------------------------------------------------------------------------------------
 % Initialize biological arrays
 if (ECOL.pelagic)&&(ECOL.demersal)
     en_input_P             = nan(FORC.nvec,2*ECOL.nfish,ECOL.nfmass);
     en_input_vb            = nan(FORC.nvec,2*ECOL.nfish,ECOL.nfmass);
     en_input               = nan(FORC.nvec,2*ECOL.nfish,ECOL.nfmass);
     gamma                  = nan(FORC.nvec,2*ECOL.nfish,ECOL.nfmass);
     flux_in                = nan(FORC.nvec,2*ECOL.nfish,ECOL.nfmass);
     flux_out               = nan(FORC.nvec,2*ECOL.nfish,ECOL.nfmass);
     flux_fish_growth       = nan(FORC.nvec,2*ECOL.nfish,ECOL.nfmass);
     flux_in_rep            = nan(FORC.nvec,2*ECOL.nfish);
     flux_in_P              = nan(FORC.nvec,2*ECOL.nfish);
     flux_in_num_eggs       = nan(FORC.nvec,2*ECOL.nfish);
     ena_regime             = nan(FORC.nvec,2*ECOL.nfish,ECOL.nfmass); 
     mortality              = nan(FORC.nvec,2*ECOL.nfish,ECOL.nfmass);
 else
     en_input_P             = nan(FORC.nvec,ECOL.nfish,ECOL.nfmass);
     en_input_vb            = nan(FORC.nvec,ECOL.nfish,ECOL.nfmass);
     en_input               = nan(FORC.nvec,ECOL.nfish,ECOL.nfmass);
     gamma                  = nan(FORC.nvec,ECOL.nfish,ECOL.nfmass);
     flux_in                = nan(FORC.nvec,ECOL.nfish,ECOL.nfmass);
     flux_out               = nan(FORC.nvec,ECOL.nfish,ECOL.nfmass);
     flux_fish_growth       = nan(FORC.nvec,ECOL.nfish,ECOL.nfmass);
     flux_in_rep            = nan(FORC.nvec,ECOL.nfish);
     flux_in_P              = nan(FORC.nvec,ECOL.nfish);
     flux_in_num_eggs       = nan(FORC.nvec,ECOL.nfish);
     ena_regime             = nan(FORC.nvec,ECOL.nfish,ECOL.nfmass); 
     mortality              = nan(FORC.nvec,ECOL.nfish,ECOL.nfmass);
 end

 %-----------------------------------------------------------------------------------------
 % Biomass initial condition
 dfish                  = boats.initial.dfish;
 % Partition of primary production between groups
 part_PP_b              = 1/ECOL.nfish;                                    % partition of primary production at boundary condition (recruitment)
 part_PP_g              = 1/ECOL.nfish;                                    % partition of primary production of growth 
 
 
 %--------------------------------------------------------------------------------------------------------------
 % Defines fileds for Economics simulations
 %--------------------------------------------------------------------------------------------------------------
 if (strcmp(MAIN.sim_type,'hd')||strcmp(MAIN.sim_type,'hf'))

   %-----------------------------------------------------------------------------------------
   % Initialize economics arrays
   if (ECOL.pelagic)&&(ECOL.demersal)
       dharvest             = nan(FORC.nvec,2*ECOL.nfish,ECOL.nfmass);
       dfish_temp           = nan(FORC.nvec,2*ECOL.nfish,ECOL.nfmass);
       effort               = nan(FORC.nvec,2*ECOL.nfish);
       effort_change        = nan(FORC.nvec,2*ECOL.nfish);
       cost                 = nan(FORC.nvec,2*ECOL.nfish);
       revenue              = nan(FORC.nvec,2*ECOL.nfish);
   else
       dharvest             = nan(FORC.nvec,ECOL.nfish,ECOL.nfmass);
       dfish_temp           = nan(FORC.nvec,ECOL.nfish,ECOL.nfmass);
       effort               = nan(FORC.nvec,ECOL.nfish);
       effort_change        = nan(FORC.nvec,ECOL.nfish);
       cost                 = nan(FORC.nvec,ECOL.nfish);
       revenue              = nan(FORC.nvec,ECOL.nfish);
   end
  
   %-----------------------------------------------------------------------------------------
   % Initialize some additional output array
   catchability_used 	 = nan(1,ntime);
   price_used            = nan(1,ntime);
   cost_effort_used 	 = nan(1,ntime);

   %-----------------------------------------------------------------------------------------
   % Effort initial condition
   effort(:,:)      =  boats.initial.effort;
   
 end


 %--------------------------------------------------------------------------------------------------------------
 % Prepare output subroutine
 %--------------------------------------------------------------------------------------------------------------
 outmode = boats.output;
 outmode.noutm = length(outmode.modes);
 % Prepare averaging time bounds
 % This is for cases where time bounds are not specified
 for indm=1:outmode.noutm
    tmode = outmode.modes{indm};
    outmode.(tmode).nvar = length(outmode.(tmode).var_name);
    %----------------------------------------------------------------------
    % Defines the time intervals for averaging
    % here anly the defaults are included
    switch outmode.modes{indm}
    case 'all'
       % NOTE: For output every timestep, no checks/averaging required
       odt = (30*CONV.sperd) /MAIN.dtts; 
       oendt = time(end)/MAIN.dtts;
       outmode.(tmode).it_bounds = [ [1:odt:oendt]' [odt:odt:oendt]'];
    case 'annual'
       % BOATS assumption is 12 months of 30 days 
       odt = (30*12*CONV.sperd) /MAIN.dtts; 
       oendt = time(end)/MAIN.dtts;
       outmode.(tmode).it_bounds = [ [1:odt:oendt]' [odt:odt:oendt]'];
    case 'decadal'
       odt = (10*30*12*sperd) /dtts; 
       oendt = time(end)/dtts;
       outmode.(tmode).it_bounds = [ [1:odt:oendt]' [odt:odt:oendt]'];
    case 'snap5year'
       % one year every 10 years
       odt  = (5*30*12*sperd) /dtts; 
       odt1 = (1*30*12*sperd) /dtts; 
       oendt = time(end)/dtts;
       vec2 = [odt:odt:oendt]';
       vec1 = vec2-odt1+1;
       outmode.(tmode).it_bounds = [vec1 vec2];
    case 'snap10year'
       % one year every 10 years
       odt  = (10*30*12*CONV.sperd) /MAIN.dtts; 
       odt1 = (1*30*12*CONV.sperd) /MAIN.dtts; 
       oendt = time(end)/MAIN.dtts;
       vec2 = [odt:odt:oendt]';
       vec1 = vec2-odt1+1;
       outmode.(tmode).it_bounds = [vec1 vec2];
    case 'final'
       % Final year only
       odt = (30*12*CONV.sperd) /MAIN.dtts; 
       oendt = time(end)/MAIN.dtts;
       outmode.(tmode).it_bounds = [oendt-odt+1 oendt];
    otherwise
       % Uses user-defined time bounds and convert from years
       outmode.(tmode).it_bounds = outmode.(tmode).t_bounds*spery/dtts;
       outmode.(tmode).it_bounds(outmode.(tmode).it_bounds==0) = 1;
    end
    %----------------------------------------------------------------------
    % correct for case when only 1 time integral is given, but it's valid
    outmode.(tmode).t_bounds = time(outmode.(tmode).it_bounds);
    if (size(outmode.(tmode).t_bounds,1)>1)&(size(outmode.(tmode).t_bounds,2)==1); 
       outmode.(tmode).t_bounds = outmode.(tmode).t_bounds';
    end
    outmode.(tmode).time = mean(outmode.(tmode).t_bounds,2);
    outmode.(tmode).year = outmode.(tmode).time/CONV.spery;
    outmode.(tmode).ndt = diff(outmode.(tmode).it_bounds,[],2)+1;
    outmode.(tmode).ntime = length(outmode.(tmode).time);
    % Does a check of variables - if variables have not been defined, removes them from output 
    ivbad = [];   
    for indv=1:outmode.(tmode).nvar
       if ~exist(outmode.(tmode).var_name{indv})
          ivbad = [ivbad indv];
       end
    end
    outmode.(tmode).var_name(ivbad) = [];
    outmode.(tmode).var_type(ivbad) = [];
    outmode.(tmode).var_proc(ivbad) = [];
    outmode.(tmode).var_outn(ivbad) = [];
    outmode.(tmode).nvar = length(outmode.(tmode).var_name);
    if isfield(outmode.(tmode),'var_derv')
       outmode.(tmode).var_derv(ivbad) = [];
    end
 end
 % Removes all invalid output modes:
 % (1) no variables to save
 % (2) no timesteps to save
 imbad = [];
 for indm=1:outmode.noutm
    nvar = outmode.(outmode.modes{indm}).nvar;
    itsz = size(outmode.(outmode.modes{indm}).it_bounds);
    if nvar==0|itsz(2)==1|any(itsz==0)
       imbad = [imbad indm];
    end
 end
 outmode = rmfield(outmode,outmode.modes(imbad));
 outmode.modes(imbad) = [];
 outmode.noutm = length(outmode.modes);

 
 %--------------------------------------------------------------------------------------------------------------
 % Initialize output arrays
 %--------------------------------------------------------------------------------------------------------------
 % This require hardwiring the dimensions of output arrays based on the
 % variable original size/type (4D,3D, etc) and processing (2di, si, gi)
 disp(['Initializing output...']);
 
 if (ECOL.pelagic)&&(ECOL.demersal)
     nfish = 2*ECOL.nfish;
 else
     nfish = ECOL.nfish;
 end
 
 for indm=1:outmode.noutm
    ontime = outmode.(outmode.modes{indm}).ntime;
    for indv=1:outmode.(outmode.modes{indm}).nvar
      % Defines output size according to variable type
       switch outmode.(outmode.modes{indm}).var_type{indv}
       case '4D'
          % Further defines output size according to variable processing
          switch outmode.(outmode.modes{indm}).var_proc{indv}
          case 'none'
          % Saves the entire variable without any processing
             outmode.(outmode.modes{indm}).(outmode.(outmode.modes{indm}).var_outn{indv}) = ...
                                           zeros(ontime,FORC.nlat,FORC.nlon,nfish,ECOL.nfmass);
          case 'si'
          % Saves the variable after integral over size dimensions
             outmode.(outmode.modes{indm}).(outmode.(outmode.modes{indm}).var_outn{indv}) = ...
                                           zeros(ontime,FORC.nlat,FORC.nlon,nfish);
          case 'si100'
          % Saves the variable after integral over size dimensions
             outmode.(outmode.modes{indm}).(outmode.(outmode.modes{indm}).var_outn{indv}) = ...
                                           zeros(ontime,FORC.nlat,FORC.nlon,nfish);
          case 'si1000'
          % Saves the variable after integral over size dimensions
             outmode.(outmode.modes{indm}).(outmode.(outmode.modes{indm}).var_outn{indv}) = ...
                                           zeros(ontime,FORC.nlat,FORC.nlon,nfish);
          case 'si10000'
          % Saves the variable after integral over size dimensions
             outmode.(outmode.modes{indm}).(outmode.(outmode.modes{indm}).var_outn{indv}) = ...
                                           zeros(ontime,FORC.nlat,FORC.nlon,nfish);
          case 'si100000'
          % Saves the variable after integral over size dimensions
             outmode.(outmode.modes{indm}).(outmode.(outmode.modes{indm}).var_outn{indv}) = ...
                                           zeros(ontime,FORC.nlat,FORC.nlon,nfish);
          case '2di_si'
          % Saves the variable after integral over size dimensions and over 2D domain
             outmode.(outmode.modes{indm}).(outmode.(outmode.modes{indm}).var_outn{indv}) = ...
                                           zeros(ontime,nfish);
          case 'LMEi_si'
          % Saves the variable after integral over size dimensions and over 2D domain only across LME
             outmode.(outmode.modes{indm}).(outmode.(outmode.modes{indm}).var_outn{indv}) = ...
                                           zeros(ontime,nfish);
          case '2di_si_gi'
          % Saves the variable summed over the fish groups after integral over size dimensions and over 2D
             outmode.(outmode.modes{indm}).(outmode.(outmode.modes{indm}).var_outn{indv}) = ...
                                           zeros(ontime,1);
          otherwise
             error(['processing case ' outmode.(outmode.modes{indm}).var_proc{indv} ' not specified']);
          end
       case '3D'
          % Further defines output size according to variable processing
          switch outmode.(outmode.modes{indm}).var_proc{indv}
          case 'none'
          % Saves the entire variable without any processing
             outmode.(outmode.modes{indm}).(outmode.(outmode.modes{indm}).var_outn{indv}) = ...
                                           zeros(ontime,FORC.nlat,FORC.nlon,nfish);
          case '2di'
          % Saves the variable after integral over 2D domain
             outmode.(outmode.modes{indm}).(outmode.(outmode.modes{indm}).var_outn{indv}) = ...
                                           zeros(ontime,nfish);
          case 'LMEi'
          % Saves the variable after integral over 2D domain only across LME
             outmode.(outmode.modes{indm}).(outmode.(outmode.modes{indm}).var_outn{indv}) = ...
                                           zeros(ontime,nfish);
          case '2di_gi'
          % Saves the variable summed over groups after integral over 2D domain 
             outmode.(outmode.modes{indm}).(outmode.(outmode.modes{indm}).var_outn{indv}) = ...
                                           zeros(ontime,1);
          otherwise
             error(['processing case ' outmode.(outmode.modes{indm}).var_proc{indv} ' not specified']);
          end
       case 'DER'
          % Initialization generally not be required for derived variables
       otherwise
          error(['processing case ' outmode.(outmode.modes{indm}).var_type{indv} ' not specified']);
       end
    end
 end

 
%--------------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------------
% MAIN LOOP
%--------------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------------
for indt = 1:ntime

  % local month
  local_month = ceil((mod(time(indt),CONV.sperfrc)/CONV.sperfrc)*MAIN.nforcing);
  if local_month==0;local_month=MAIN.nforcing;end
  disp(['indt : ' num2str(indt) ' / ' num2str(ntime) ' - local month : ' num2str(local_month)]);

  %--------------------------------------------------------------------------------------------------
  % Physical and ecological forcings (NPP / PFB, in mmolC m-2 s-1, NPP_ED in mmolC m-2 d-1, T in C/K)
  %--------------------------------------------------------------------------------------------------
  npp         = squeeze(FORC.npp_vec(:,local_month));
  npp_ed      = squeeze(FORC.npp_ed_vec(:,local_month));
  pfb         = squeeze(FORC.pfb_vec(:,local_month));
  temp_phyto  = squeeze(FORC.temperature_pel_vec(:,local_month));
  temp_fish_pel_A = squeeze(FORC.temperature_pel_K_vec(:,local_month));
  temp_fish_pel_m = temp_fish_pel_A;  
  temp_fish_dem_A = squeeze(FORC.temperature_dem_K_vec(:,local_month));
  temp_fish_dem_m = temp_fish_dem_A;

  % For heterogenous catch
  if ECON.catchpar
     catch_profile = FORC.catcha_profile;
%     catch_profile = squeeze(FORC.catcha_profile(:,local_month)); 
  end

  %-----------------------------------------------------------------------------------------------------------
  % Large fraction of phytoplankton and representative phytoplankton mass (Dunne)
  %-----------------------------------------------------------------------------------------------------------
  s_over_p = ( -1.0 + ( 1.0 + 4.0 .* npp_ed ./ (exp(ENVI.kappa_eppley(1).*temp_phyto) .* ...
    ENVI.Prod_star(1)) ).^0.5) .* 0.5;
  frac_lg_du = s_over_p ./ (1.0 + s_over_p);                               % large fraction of PP as in Dunne et al. (2005)
  mphyto = (ENVI.mc_phy_l.^frac_lg_du) .* (ENVI.mc_phy_s.^(1.0 - frac_lg_du));
  
  %-----------------------------------------------------------------------------------------------------------
  % Representative benthis mass
  %-----------------------------------------------------------------------------------------------------------
  mbentho = ENVI.mc_benthic(1);

  %-----------------------------------------------------------------------------------------------------------
  % Growth rate
  %-----------------------------------------------------------------------------------------------------------
  % growth rate = production distribution * mass / biomass distribution
  % multiply by growth partition function (part_PP_g)
  %-------------------------------------------------------------------------------------
% % Original code - changed to optimize calculation
% en_input_P = repmat(npp./mphyto,[1 1 nfish nfmass]) .* ...
%              ((fmass_4d ./ repmat(mphyto,[1 1 nfish nfmass])).^(tro_sca-1)) .* ...
%              fmass_4d ./ squeeze(dfish + epsln) * part_PP_g;
  % Optimized code by using "bsxfun" instead of repmat

  if (ECOL.pelagic)&&(ECOL.demersal)
%       en_input_P(:,1:3,:) = repmat(npp./mphyto,[1 ECOL.nfish ECOL.nfmass]) .* ...
%                             STRU.fmass_4d_vec(:,1:3,:).^(repmat(INIT.tro_sca(:,1),[1 ECOL.nfish ECOL.nfmass])-1)    ./ repmat(mphyto.^(INIT.tro_sca(:,1)-1) ,[1 ECOL.nfish ECOL.nfmass]) .* ...
%                             STRU.fmass_4d_vec(:,1:3,:) ./ squeeze(dfish(:,1:3,:) + CONV.epsln) * part_PP_g;

      en_input_P(:,1:3,:) = bsxfun(@times,npp./mphyto, ...
                   (bsxfun(@rdivide,STRU.fmass_4d_vec(:,1:3,:).^(repmat(INIT.tro_sca(:,1),[1 ECOL.nfish ECOL.nfmass])-1),repmat(mphyto.^(INIT.tro_sca(:,1)-1),[1 ECOL.nfish ECOL.nfmass])) .* ... 
                   STRU.fmass_4d_vec(:,1:3,:) ./ squeeze(dfish(:,1:3,:) + CONV.epsln) * part_PP_g));
      en_input_P(:,4:6,:) = bsxfun(@times,pfb./mbentho, ...
                   (bsxfun(@rdivide,STRU.fmass_4d_vec(:,4:6,:).^(repmat(INIT.tro_sca(:,2),[1 ECOL.nfish ECOL.nfmass])-1),repmat(mbentho.^(INIT.tro_sca(1,2)-1),[size(pfb,1) ECOL.nfish ECOL.nfmass])) .* ... 
                   STRU.fmass_4d_vec(:,4:6,:) ./ squeeze(dfish(:,4:6,:) + CONV.epsln) * part_PP_g));

%       en_input_P(:,1:3,:) = bsxfun(@times,npp./mphyto, ...
%                    (STRU.fmass_4d_vec(:,1:3,:)./repmat(mphyto,[1 ECOL.nfish ECOL.nfmass])).^(repmat(INIT.tro_sca(:,1),[1 ECOL.nfish ECOL.nfmass])-1) .* ... 
%                    STRU.fmass_4d_vec(:,1:3,:) ./ squeeze(dfish(:,1:3,:) + CONV.epsln) * part_PP_g);
%       en_input_P(:,4:6,:) = bsxfun(@times,pfb./mbentho, ...
%                    (STRU.fmass_4d_vec(:,4:6,:)./repmat(mbentho,[size(pfb,1) ECOL.nfish ECOL.nfmass])).^(INIT.tro_sca(1,2)-1) .* ... 
%                    STRU.fmass_4d_vec(:,4:6,:) ./ squeeze(dfish(:,4:6,:) + CONV.epsln) * part_PP_g);
  else
      en_input_P = bsxfun(@times,npp./mphyto, ...
                   (bsxfun(@rdivide,STRU.fmass_4d_vec.^(repmat(INIT.tro_sca(:,1),[1 ECOL.nfish ECOL.nfmass])-1),repmat(mphyto.^(INIT.tro_sca(:,1)-1),[1 ECOL.nfish ECOL.nfmass])) .* ... 
                   STRU.fmass_4d_vec ./ squeeze(dfish + CONV.epsln) * part_PP_g));
  end
   
  %---------------------------------------------------------------------------------------
  % Based on allometric scaling (von Bertalanffy)
  % calculate temperature dependencies, then growth rate, the activity loss constant (ka)  
  if (ECOL.pelagic)&&(ECOL.demersal)
      temp_dep_A_P  = repmat(exp( (-ENVI.E_activation_A(1)/ENVI.k_Boltzmann) .* (1./temp_fish_pel_A - 1/ENVI.temp_ref_A)),[1 ECOL.nfish ECOL.nfmass]);
      temp_dep_m_P  = repmat(exp( (-ENVI.E_activation_m(1)/ENVI.k_Boltzmann) .* (1./temp_fish_pel_m - 1/ENVI.temp_ref_A)),[1 ECOL.nfish ECOL.nfmass]);
      A_P 	      = A0(1) .* temp_dep_A_P;
      ka_P 	      = A_P * ECOL.eff_a .* STRU.minf_4d_p_bm1_P_vec;                      %  A * eff_a .* minf_4d.^(b_allo-1) (s-1)
      en_input_vb(:,1:3,:) = A_P .* STRU.fmass_4d_p_b_P_vec - ka_P .* STRU.fmass_4d_vec(:,1:3,:);              % A .* fmass_4d.^b_allo - ka .* fmass_4d;
      temp_dep_A_D  = repmat(exp( (-ENVI.E_activation_A(2)/ENVI.k_Boltzmann) .* (1./temp_fish_dem_A - 1/ENVI.temp_ref_A)),[1 ECOL.nfish ECOL.nfmass]);
      temp_dep_m_D  = repmat(exp( (-ENVI.E_activation_m(2)/ENVI.k_Boltzmann) .* (1./temp_fish_dem_m - 1/ENVI.temp_ref_A)),[1 ECOL.nfish ECOL.nfmass]);
      A_D 	      = A0(2) .* temp_dep_A_D;
      ka_D 	      = A_D * ECOL.eff_a .* STRU.minf_4d_p_bm1_D_vec;                      %  A * eff_a .* minf_4d.^(b_allo-1) (s-1)
      en_input_vb(:,4:6,:) = A_D .* STRU.fmass_4d_p_b_D_vec - ka_D .* STRU.fmass_4d_vec(:,4:6,:);              % A .* fmass_4d.^b_allo - ka .* fmass_4d;
      
      temp_dep_m = cat(2,temp_dep_m_P,temp_dep_m_D);
  else
      temp_dep_A_P  = repmat(exp( (-ENVI.E_activation_A(1)/ENVI.k_Boltzmann) .* (1./temp_fish_pel_A - 1/ENVI.temp_ref_A)),[1 ECOL.nfish ECOL.nfmass]);
      temp_dep_m_P  = repmat(exp( (-ENVI.E_activation_m(1)/ENVI.k_Boltzmann) .* (1./temp_fish_pel_m - 1/ENVI.temp_ref_A)),[1 ECOL.nfish ECOL.nfmass]);
      A_P 	      = A0(1) .* temp_dep_A_P;
      ka_P 	      = A_P * ECOL.eff_a .* STRU.minf_4d_p_bm1_P_vec;                      %  A * eff_a .* minf_4d.^(b_allo-1) (s-1)
      en_input_vb = A_P .* STRU.fmass_4d_p_b_P_vec - ka_P .* STRU.fmass_4d_vec;              % A .* fmass_4d.^b_allo - ka .* fmass_4d;
      
      temp_dep_m = temp_dep_m_P;
  end
  
  %---------------------------------------------------------------------------------------
  % Input energy (energy available to growth and reproduction)
  % minimum of en_input_P and en_input_vb
  en_input                     = min(en_input_P,en_input_vb);
  en_input(STRU.mask_notexist_4d_vec)   = NaN;

  %---------------------------------------------------------------------------------------
  % Somatic growth rate (growth rate of body tissues)
  gamma                        = (1 - STRU.rep_alloc_frac_vec).*en_input;
    
  %---------------------------------------------------------------------------------------
  % Flux out of a mass class
  flux_out   = gamma .* squeeze(dfish+CONV.epsln) ./ STRU.delfm_4d_vec;
  
  %---------------------------------------------------------------------------------------
  % Biomass growth
  % increase in biomass in growth to next mass class
  % arises in conversion from abundance-based MFVF to biomass-based
  flux_fish_growth   = gamma .* squeeze(dfish+CONV.epsln) ./ STRU.fmass_4d_vec;
  
  %---------------------------------------------------------------------------------------
  % boundary condition (flux in to first mass class)
  %---------------------------------------------------------------------------------------    
  % Boundary condition based on primary production
  % multiply by boundary condition partition function (part_PP_b)
  if (ECOL.pelagic)&&(ECOL.demersal)
      flux_in_P(:,1:3) = part_PP_b * (repmat(npp./mphyto,[1 ECOL.nfish])) .* (STRU.fmass_bc./repmat(mphyto,[1 ECOL.nfish])).^(repmat(INIT.tro_sca(:,1),[1 ECOL.nfish])-1) * STRU.fmass_bc / STRU.delfm_4d_vec(1,1,1);
      flux_in_P(:,4:6) = part_PP_b * (repmat(pfb./mbentho,[1 ECOL.nfish])) .* (STRU.fmass_bc./repmat(mbentho,[size(pfb,1) ECOL.nfish])).^(repmat(INIT.tro_sca(:,2),[1 ECOL.nfish])-1) * STRU.fmass_bc / STRU.delfm_4d_vec(1,4,1);
  else
      flux_in_P = part_PP_b * (repmat(npp./mphyto,[1 ECOL.nfish])) .* (STRU.fmass_bc./repmat(mphyto,[1 ECOL.nfish])).^(repmat(INIT.tro_sca(:),[1 ECOL.nfish])-1) * STRU.fmass_bc / STRU.delfm_4d_vec(1,1,1);
  end
  %---------------------------------------------------------------------------------------
  % Flux in of number of eggs produced
  if (ECOL.pelagic)&&(ECOL.demersal)
      flux_in_num_eggs(:,1:3) = (ECOL.frac_fem/ECOL.m_egg) .* ...
        nansum( STRU.rep_alloc_frac_vec(:,1:3,:) .* en_input(:,1:3,:) .* squeeze(dfish(:,1:3,:)) .* STRU.delfm_4d_over_fmass_4d_vec(:,1:3,:),3) / STRU.delfm_4d_vec(1,1,1);
      flux_in_num_eggs(:,4:6) = (ECOL.frac_fem/ECOL.m_egg) .* ...
        nansum( STRU.rep_alloc_frac_vec(:,4:6,:) .* en_input(:,4:6,:) .* squeeze(dfish(:,4:6,:)) .* STRU.delfm_4d_over_fmass_4d_vec(:,4:6,:),3) / STRU.delfm_4d_vec(1,4,1);
  else
      flux_in_num_eggs = (ECOL.frac_fem/ECOL.m_egg) .* ...
        nansum( STRU.rep_alloc_frac_vec .* en_input .* squeeze(dfish) .* STRU.delfm_4d_over_fmass_4d_vec,3) / STRU.delfm_4d_vec(1,1,1);
  end
  %---------------------------------------------------------------------------------------
  % Boundary condition based on recruitment (production and survival of eggs)
  % If the flux of number of eggs is less than 0.001 eggs m-2 y-1
  % then the flux due to the production and survival of eggs is zero
  mask_eggs_low = (flux_in_num_eggs < 0.001/CONV.spery);
  if (ECOL.pelagic)&&(ECOL.demersal)
      flux_in_rep(:,1:3) = (STRU.fmass_bc*ECOL.egg_surv(1)) .* flux_in_num_eggs(:,1:3);
      flux_in_rep(:,4:6) = (STRU.fmass_bc*ECOL.egg_surv(2)) .* flux_in_num_eggs(:,4:6);
  else
      flux_in_rep = (STRU.fmass_bc*ECOL.egg_surv(1)) .* flux_in_num_eggs;
  end
  flux_in_rep(mask_eggs_low) = 0;
 
  %---------------------------------------------------------------------------------------
  % Boundary condition (Beverton-Holt form)
  flux_in(:,:,1) = flux_in_P .* ((flux_in_rep + CONV.epsln) ./ (flux_in_P + flux_in_rep + CONV.epsln));
    
  %---------------------------------------------------------------------------------------
  % Flux in of other mass classes
  flux_in(:,:,2:end) = gamma(:,:,1:end-1) .* squeeze(dfish(:,:,1:end-1)) ./ STRU.delfm_2end_4d_vec;
  
  %---------------------------------------------------------------------------------------
  % Input energy (available energy to growth and reproduction regime
  % 1 if en_input_P less - then en_input_P determines input energy
  ena_regime = (en_input_P < en_input_vb);

  %---------------------------------------------------------------------------------------
  % Mortality
  %---------------------------------------------------------------------------------------
  % Mass-specific mortality rate 
  % Calculate associated growth rate with mortality temperature dependence temp_dep_m
  % calculate mortality rate mortality0
  % mortality00 is the exp(zeta_1) term from the model description
  if (ECOL.pelagic)&&(ECOL.demersal)
      A           = A0(1)*temp_dep_m(:,1:3,:);
      mortality0(:,1:3,:)  = (exp(ECOL.zeta1(1))/3)*A;
      A           = A0(2)*temp_dep_m(:,4:6,:);
      mortality0(:,4:6,:)  = (exp(ECOL.zeta1(2))/3)*A;
  else
      A           = A0(1)*temp_dep_m;
      mortality0  = (exp(ECOL.zeta1(1))/3)*A;
  end
  
  %---------------------------------------------------------------------------------------  
  % Mortality rate 
  % Charnov et al. (2013)
  % minf_4d_p_hplusbm1 is minf_4d.^(h_allo + b_allo - 1)
  mortality = mortality0 .* STRU.minf_4d_p_hplusbm1_vec .* STRU.fmass_4d_p_mh_vec .* squeeze(dfish);
    
  %---------------------------------------------------------------------------------------
  % Integrate fish and effort
  %---------------------------------------------------------------------------------------
  if strcmp(MAIN.sim_type,'nh') % No economic harvesting
  
    %-------------------------------------------------------------------------------------
    % Integrate dfish
    %-------------------------------------------------------------------------------------
    dfish  = squeeze(dfish) + ( flux_in - flux_out + flux_fish_growth - mortality ) * MAIN.dtts;
    mask_dfish_neg = (squeeze(dfish) < 0);
    dfish(mask_dfish_neg)  = 0;
     
  elseif strcmp(MAIN.sim_type,'hd')  % Economic harvesting
     
    %-------------------------------------------------------------------------------------
    % Update fish to calculate harvest
    %-------------------------------------------------------------------------------------
    dfish_temp  = squeeze(dfish) + ( flux_in - flux_out + flux_fish_growth - mortality ) * MAIN.dtts;
    mask_dfish_temp_neg = (squeeze(dfish_temp) < 0);
    dfish_temp(mask_dfish_temp_neg)  = 0;
   
    %-------------------------------------------------------------------------------------
    % Catchability 
    %-------------------------------------------------------------------------------------   
    if (time(indt) < ECON.harvest_start*CONV.spery)                        
        % Catchability zero before start year (adjusted for dharvest calculation)
        if (ECOL.pelagic)&&(ECOL.demersal)
    	    qcatch = 0 * ones(1,2*ECOL.nfish);
	else
            qcatch = 0 * ones(1,ECOL.nfish);
	end
    else  
        % Set catchability forcing scenario (adjusted for dharvest calculation)     
        catchability_used(1,indt) = FORC.catchability(indt);
        if (ECOL.pelagic)&&(ECOL.demersal)
            qcatch = 4*FORC.catchability(indt) * ones(1,2*ECOL.nfish);
        else
            qcatch = FORC.catchability(indt) * ones(1,ECOL.nfish);
        end
    end

    % Spatially variable catchability
    if ECON.catchpar  
	qcatch_star = repmat(catch_profile,[1,3]).*repmat(qcatch(1,1:3),[size(catch_profile,1),1]);
	qcatch_star = cat(2,qcatch_star,0.9*repmat(qcatch(1,1:3),[size(catch_profile,1),1]));
    end 

    %-------------------------------------------------------------------------------------
    % dharvest [nlat,nlon,nfish,nfmass]
    %-------------------------------------------------------------------------------------
    % qcatch * effort * selectivity * dfish_temp
    % Set upper limit (min) so that no more than the fish present can be caught (dfish_temp/dtts)
    %-------------------------------------------------------------------------------------
%   % Original code - changed to optimize calculation
%   dharvest = min(squeeze(dfish_temp)/dtts, permute(repmat(qcatch(:),[1 nlat nlon nfmass]),[2 3 1 4]) .* ...
%              repmat(squeeze(effort+epsln),[1 1 1 nfmass]) .* selectivity_4d .* squeeze(dfish_temp));
    % Optimized code by using "bsxfun" instead of repmat
    if ECON.catchpar
        dharvest = min(squeeze(dfish_temp)/MAIN.dtts, ...
               bsxfun(@times,bsxfun(@times,qcatch_star,squeeze(effort+CONV.epsln)), ...
               STRU.selectivity_4d_vec .* squeeze(dfish_temp)));
    else
        dharvest = min(squeeze(dfish_temp)/MAIN.dtts, ...
               bsxfun(@times,bsxfun(@times,permute(qcatch(:),[2 1]),squeeze(effort+CONV.epsln)), ...
               STRU.selectivity_4d_vec .* squeeze(dfish_temp)));
   % dharvest = min(squeeze(dfish_temp)/MAIN.dtts, ...
   %            bsxfun(@times,bsxfun(@times,qcatch,squeeze(effort+CONV.epsln)), ...
   %            STRU.selectivity_4d_vec .* squeeze(dfish_temp)));
    end

    mask_dharvest_neg = (squeeze(dharvest) < 0);
    dharvest(mask_dharvest_neg)  = 0;
    
    %----------------------------------------------------------------------------------
    % Price forcing
    % input price ($ g-1), multiply by mmolC_2_wetB to get ($ mmolC-1)
    %----------------------------------------------------------------------------------
    % Set price forcing scenario (adjusted for dharvest calculation)
    price_used(1,indt) = FORC.price(indt);
    if (ECOL.pelagic)&&(ECOL.demersal)
        price       = FORC.price(indt) * ones(2*ECOL.nfish,ECOL.nfmass);
    else
        price       = FORC.price(indt) * ones(ECOL.nfish,ECOL.nfmass);
    end
    
    %-------------------------------------------------------------------------------------
    % Cost per unit effort forcing
    % input cost per unit effort ($ W-1 s-1)
    %-------------------------------------------------------------------------------------
    % Set cost forcing scenario (adjusted for dharvest calculation)
    cost_effort_used(1,indt) = FORC.cost(indt);
    if (ECOL.pelagic)&&(ECOL.demersal)
	if ECON.costpar
            cost_effort = [repmat(FORC.cost(indt) * boats.forcing.cost_profile2,[1,3]),repmat(FORC.cost(indt) * (boats.forcing.cost_profile1+boats.forcing.cost_profile2)/2,[1,3])];
            %cost_effort = FORC.cost(indt) * ones(length(boats.forcing.cost_profile),ECOL.nfish);
            %cost_effort = [cost_effort,repmat(FORC.cost(indt) * boats.forcing.depth_profile,[1,3])];
	else
	    cost_effort = FORC.cost(indt) * ones(1,2*ECOL.nfish);
	end
    else
        cost_effort = FORC.cost(indt) * ones(1,ECOL.nfish);
    end

    %-------------------------------------------------------------------------------------
    % revenue [nlat,nlon,nfish]
    %-------------------------------------------------------------------------------------
    % sum over mass (price * dharvest * delfm)
    %-------------------------------------------------------------------------------------
%   % Original code - changed to optimize calculation
%   revenue = nansum( permute(repmat(price,[1 1 nlat nlon]),[3 4 1 2]) .* squeeze(dharvest) .* delfm_4d, 4);
    % Optimized code by using "bsxfun" instead of repmat
    revenue = nansum(bsxfun(@times,permute(price,[3 1 2]),squeeze(dharvest) .* STRU.delfm_4d_vec),3);
    
    %-------------------------------------------------------------------------------------
    % cost [nlat,nlon,nfish]
    %-------------------------------------------------------------------------------------
    % cost_effort * effort
    %-------------------------------------------------------------------------------------
%   % Original code - changed to optimize calculation
%   cost =  permute(repmat(cost_effort(:),[1 nlat nlon]),[2 3 1]) .* squeeze(effort + epsln);
    % Optimized code by using "bsxfun" instead of repmat
%    cost =  bsxfun(@times,permute(cost_effort(:),[2 1]),squeeze(effort + CONV.epsln));
    cost =  bsxfun(@times,cost_effort,squeeze(effort + CONV.epsln));

    %-------------------------------------------------------------------------------------
    % effort_change [nlat,nlon,nfish]
    %-------------------------------------------------------------------------------------
    effort_change = ECON.k_e * (revenue - cost) ./ (effort + CONV.epsln);
    
    %-------------------------------------------------------------------------------------
    % integrate dfish [nlat,nlon,nfish,nfmass]
    %-------------------------------------------------------------------------------------
    dfish  = squeeze(dfish_temp) - squeeze(dharvest) * MAIN.dtts;
    mask_dfish_neg = (squeeze(dfish) < 0);
    dfish(mask_dfish_neg)  = 0;

    %-------------------------------------------------------------------------------------
    % integrate effort [nlat,nlon,nfish]
    %-------------------------------------------------------------------------------------
    effort = squeeze(effort) + effort_change * MAIN.dtts;
    mask_effort_neg = (squeeze(effort) < 0);
    effort(mask_effort_neg)  = 0;

  elseif strcmp(MAIN.sim_type,'hf')  % Economic harvesting

    %-------------------------------------------------------------------------------------
    % Update fish to calculate harvest
    %-------------------------------------------------------------------------------------
    dfish_temp  = squeeze(dfish) + ( flux_in - flux_out + flux_fish_growth - mortality ) * MAIN.dtts;
    mask_dfish_temp_neg = (squeeze(dfish_temp) < 0);
    dfish_temp(mask_dfish_temp_neg)  = 0;

    %-------------------------------------------------------------------------------------
    % Catchability
    %-------------------------------------------------------------------------------------
    if (time(indt) < ECON.harvest_start*CONV.spery)
        % Catchability zero before start year (adjusted for dharvest calculation)
        if (ECOL.pelagic)&&(ECOL.demersal)
            qcatch = 0 * ones(1,2*ECOL.nfish);
        else
            qcatch = 0 * ones(1,ECOL.nfish);
        end
    else
        % Set catchability forcing scenario (adjusted for dharvest calculation)
        catchability_used(1,indt) = FORC.catchability(indt);
        if (ECOL.pelagic)&&(ECOL.demersal)
            qcatch = 4*FORC.catchability(indt) * ones(1,2*ECOL.nfish);
        else
            qcatch = FORC.catchability(indt) * ones(1,ECOL.nfish);
        end
    end

    %-------------------------------------------------------------------------------------
    % dharvest [nlat,nlon,nfish,nfmass]
    %-------------------------------------------------------------------------------------
    % qcatch * effort * selectivity * dfish_temp
    % Set upper limit (min) so that no more than the fish present can be caught (dfish_temp/dtts)
    %-------------------------------------------------------------------------------------
%   % Original code - changed to optimize calculation
%   dharvest = min(squeeze(dfish_temp)/dtts, permute(repmat(qcatch(:),[1 nlat nlon nfmass]),[2 3 1 4]) .* ...
%              repmat(squeeze(effort+epsln),[1 1 1 nfmass]) .* selectivity_4d .* squeeze(dfish_temp));
    % Optimized code by using "bsxfun" instead of repmat
    dharvest = min(squeeze(dfish_temp)/MAIN.dtts, ...
               bsxfun(@times,bsxfun(@times,permute(qcatch(:),[2 1]),squeeze(effort+CONV.epsln)), ...
               STRU.selectivity_4d_vec .* squeeze(dfish_temp)));

    mask_dharvest_neg = (squeeze(dharvest) < 0);
    dharvest(mask_dharvest_neg)  = 0;

    %-------------------------------------------------------------------------------------
    dfish  = squeeze(dfish_temp) - squeeze(dharvest) * MAIN.dtts;
    mask_dfish_neg = (squeeze(dfish) < 0);
    dfish(mask_dfish_neg)  = 0;

  end % end integrate
  
  
%-----------------------------------------------------------------------------------------  
% Output : calculate and save output variables
%-----------------------------------------------------------------------------------------
% Fills in output variables inside previously defined output "modes" structures. 
% Handles the time averaging dynamically by dividing the temporal output by the number 
% of timesteps and summing. Sadly, this part requires a few "eval statements"
 for indm=1:outmode.noutm
    % Decides whether the current timestep should be saved out 
    % finds the index corresponding to the current time averaging interval - if it exists
    prodt = prod(outmode.(outmode.modes{indm}).it_bounds - indt,2);
    iit = find(prodt<=0,1,'first');
    if ~isempty(iit);
       % Processes the current output timestep
       for indv=1:outmode.(outmode.modes{indm}).nvar
          tvar_outn = outmode.(outmode.modes{indm}).var_outn; 
          % Defines output size according to variable type
          switch outmode.(outmode.modes{indm}).var_type{indv}
          case '4D'
             tmpvar = eval(outmode.(outmode.modes{indm}).var_name{indv});
             switch outmode.(outmode.modes{indm}).var_proc{indv}
             case 'none'
             % Saves the entire variable without any processing
                % Convert vectors to maps
                for ifish = 1:nfish
                    for ifmass = 1:ECOL.nfmass
                        [var_map(:,:,ifish,ifmass)]    = function_vec_2_map(tmpvar(:,ifish,ifmass),FORC.indlat,FORC.indlon,FORC.mask(:,:,1));
                    end
                end
                outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit,:,:,:,:)  = ...
                squeeze(outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit,:,:,:,:))  + ...
                                              var_map/outmode.(outmode.modes{indm}).ndt(iit) + ...
                                              STRU.mask_land_g_s_nan;
             case 'si'
             % Saves the variable after integral over size dimensions
                tmpvar1 = squeeze(nansum( tmpvar .* STRU.delfm_4d_vec,3)) * CONV.mmolC_2_wetB + STRU.mask_land_g_nan_vec;
                for ifish = 1:nfish
                    [var_map(:,:,ifish)]    = function_vec_2_map(tmpvar1(:,ifish),FORC.indlat,FORC.indlon,FORC.mask(:,:,1));
                end
                outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit,:,:,:)  = ...
                squeeze(outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit,:,:,:))  + ...
                                              var_map/outmode.(outmode.modes{indm}).ndt(iit);
             case 'si100'
             % Saves the variable after integral over size dimensions
                temp_delfm_4d = STRU.delfm_4d_vec;
                temp_delfm_4d(:,:,1:12) = temp_delfm_4d(:,:,1:12)*NaN;
                temp_delfm_4d(:,:,13) = temp_delfm_4d(:,:,13)*0.52;
                tmpvar1 = squeeze(nansum( tmpvar .* temp_delfm_4d,3)) * CONV.mmolC_2_wetB + STRU.mask_land_g_nan_vec;
                for ifish = 1:nfish
                    [var_map(:,:,ifish)]    = function_vec_2_map(tmpvar1(:,ifish),FORC.indlat,FORC.indlon,FORC.mask(:,:,1));
                end
                outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit,:,:,:)  = ...
                squeeze(outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit,:,:,:))  + ...
                                              var_map/outmode.(outmode.modes{indm}).ndt(iit);
             case 'si1000'
             % Saves the variable after integral over size dimensions
                temp_delfm_4d = STRU.delfm_4d_vec;
                temp_delfm_4d(:,:,1:25) = temp_delfm_4d(:,:,1:25)*NaN;
                tmpvar1 = squeeze(nansum( tmpvar .* temp_delfm_4d,3)) * CONV.mmolC_2_wetB + STRU.mask_land_g_nan_vec;
                for ifish = 1:nfish
                    [var_map(:,:,ifish)]    = function_vec_2_map(tmpvar1(:,ifish),FORC.indlat,FORC.indlon,FORC.mask(:,:,1));
                end
                outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit,:,:,:)  = ...
                squeeze(outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit,:,:,:))  + ...
                                              var_map/outmode.(outmode.modes{indm}).ndt(iit);
             case 'si10000'
             % Saves the variable after integral over size dimensions
                temp_delfm_4d = STRU.delfm_4d_vec;
                temp_delfm_4d(:,:,1:37) = temp_delfm_4d(:,:,1:37)*NaN; 
                temp_delfm_4d(:,:,38) = temp_delfm_4d(:,:,38)*0.52;
                tmpvar1 = squeeze(nansum( tmpvar .* temp_delfm_4d,3)) * CONV.mmolC_2_wetB + STRU.mask_land_g_nan_vec;
                for ifish = 1:nfish
                    [var_map(:,:,ifish)]    = function_vec_2_map(tmpvar1(:,ifish),FORC.indlat,FORC.indlon,FORC.mask(:,:,1));
                end
                outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit,:,:,:)  = ...
                squeeze(outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit,:,:,:))  + ...
                                              var_map/outmode.(outmode.modes{indm}).ndt(iit);
             case 'si100000'
             % Saves the variable after integral over size dimensions
                temp_delfm_4d = STRU.delfm_4d_vec;

                tmpvar1 = squeeze(nansum( tmpvar .* temp_delfm_4d,3)) * CONV.mmolC_2_wetB + STRU.mask_land_g_nan_vec;
                for ifish = 1:nfish
                    [var_map(:,:,ifish)]    = function_vec_2_map(tmpvar1(:,ifish),FORC.indlat,FORC.indlon,FORC.mask(:,:,1));
                end
                outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit,:,:,:)  = ...
                squeeze(outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit,:,:,:))  + ...
                                              var_map/outmode.(outmode.modes{indm}).ndt(iit);
             case 'sis'
             % Saves the variable after integral over size dimensions
                tmpvar1 = squeeze(nansum( tmpvar .* STRU.selectivity_4d_vec,3)) * CONV.mmolC_2_wetB + STRU.mask_land_g_nan_vec;
                for ifish = 1:nfish
                    [var_map(:,:,ifish)]    = function_vec_2_map(tmpvar1(:,ifish),FORC.indlon,FORC.mask(:,:,1));
                end
                outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit,:,:,:)  = ...
                   squeeze(outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit,:,:,:))  + ...
                   var_map/outmode.(outmode.modes{indm}).ndt(iit);
             case '2di_si'
             % Saves the variable after integral over size dimensions and over 2D domain
                tmpvar1 = squeeze(nansum(tmpvar .* STRU.delfm_4d_vec,3)) * CONV.mmolC_2_wetB;
                tmpvar2 = nan([1 size(tmpvar1,2)]);
                for indgg=1:size(tmpvar1,2)
                   tmpvar2(1,indgg) = integrate_2d(squeeze(tmpvar1(:,indgg)),FORC.surf_vec);
                end
                outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit,:)  = ...
                squeeze(outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit,:))  + ...
                                              tmpvar2/outmode.(outmode.modes{indm}).ndt(iit);
             case 'LMEi_si'
             % Saves the variable after integral over size dimensions and over 2D domain only across LME
                tmpvar1 = squeeze(nansum(tmpvar .* STRU.delfm_4d,4)) * CONV.mmolC_2_wetB;
                tmpvar2 = nan([1 3]);
                for indgg=1:3
                   tmpvar2(1,indgg) = integrate_2d(mask_LME_nan+squeeze(tmpvar1(:,:,indgg)),boats.surf);
                end
                outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit,:)  = ...
                squeeze(outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit,:))  + ...
                                              tmpvar2/outmode.(outmode.modes{indm}).ndt(iit);
             case '2di_si_gi'
             % Saves the variable summed over the fish groups after integral over size dimensions and over 2D
                tmpvar1 = squeeze(nansum(tmpvar .* delfm_4d,4)) * mmolC_2_wetB;
                tmpvar2 = nan([1 3]);
                for indgg=1:3
                   tmpvar2(1,indgg) = integrate_2d(squeeze(tmpvar1(:,:,indgg)),boats.surf);
                end
                tmpvar3 = nansum(tmpvar2,2);
                outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit)  = ...
                squeeze(outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit))  + ...
                                              tmpvar3/outmode.(outmode.modes{indm}).ndt(iit);
             otherwise
                error(['processing case ' outmode.(outmode.modes{indm}).var_proc{indv} ' not specified']);
             end
          case '3D'
             tmpvar = eval(outmode.(outmode.modes{indm}).var_name{indv});
             % Further defines output size according to variable processing
             switch outmode.(outmode.modes{indm}).var_proc{indv}
             case 'none'
             % Saves the entire variable without any processing
                for ifish = 1:nfish
                    [var_map(:,:,ifish)]    = function_vec_2_map(tmpvar(:,ifish),FORC.indlat,FORC.indlon,FORC.mask(:,:,1));
                end
                outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit,:,:,:)  = ...
                squeeze(outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit,:,:,:))  + ...
                                              var_map/outmode.(outmode.modes{indm}).ndt(iit) + ...
                                              STRU.mask_land_g_nan;
             case '2di'
             % Saves the variable after integral over 2D domain
                tmpvar1 = nan([1 size(tmpvar1,2)]);
                for indgg=1:size(tmpvar1,2)
                   tmpvar1(1,indgg) = integrate_2d(squeeze(tmpvar(:,indgg)),FORC.surf_vec);
                end
                outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit,:)  = ...
                squeeze(outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit,:))  + ...
                                              tmpvar1/outmode.(outmode.modes{indm}).ndt(iit);
             case 'LMEi'
             % Saves the variable after integral over 2D domain only across LME
                tmpvar1 = nan([1 3]);
                for indgg=1:3
                   tmpvar1(1,indgg) = integrate_2d(mask_LME_nan+squeeze(tmpvar(:,:,indgg)),boats.surf);
                end
                outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit,:)  = ...
                squeeze(outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit,:))  + ...
                                              tmpvar1/outmode.(outmode.modes{indm}).ndt(iit);
             case '2di_gi'
             % Saves the variable summed over groups after integral over 2D domain
                tmpvar1 = nan([1 3]);
                for indgg=1:3
                   tmpvar1(1,indgg) = integrate_2d(squeeze(tmpvar(:,:,indgg)),boats.surf);
                end
                tmpvar2 = nansum(tmpvar1,2);
                outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit,:)  = ...
                squeeze(outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit,:))  + ...
                                              tmpvar2/outmode.(outmode.modes{indm}).ndt(iit);
             otherwise
                error(['processing case ' outmode.(outmode.modes{indm}).var_proc{indv} ' not specified']);
             end
          case 'DER'
             % Calculations for derived variables may not be required
          otherwise
             error(['processing case ' outmode.(outmode.modes{indm}).var_type{indv} ' not specified']);
          end
       end
    end
 end
end % for indt = 1:ntime

%--------------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------------
% END OF MAIN LOOP
%--------------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------------


%--------------------------------------------------------------------------------------------------------------
% Save last timestep of dfish and effort for restart
%--------------------------------------------------------------------------------------------------------------
 if MAIN.save_restart==1
   boats.dfish      = dfish;
   if (strcmp(MAIN.sim_type,'hd')||strcmp(MAIN.sim_type,'hf'))
     boats.effort   = effort;
   end
 end

%-----------------------------------------------------------------------------------------
% Output : additional derived variables needed for output - 
% these are calculated outside the time-stepping loop
%-----------------------------------------------------------------------------------------
 for indm=1:outmode.noutm
    % Derived variables for output
    % Generally avoid doing calculation and integrals within 
    % the main  time loop if it can be avoided
    ider = find(strcmp(outmode.(outmode.modes{indm}).var_type,'DER'));
    for indv=1:length(ider)
       vproc = outmode.(outmode.modes{indm}).var_proc{ider(indv)};
       switch vproc
       case 'gi';
          % sums over groups the deriving variable
          % for 'gi' option, assumes group dimension is the last
          derv = outmode.(outmode.modes{indm}).var_derv{ider(indv)};
          tempvar = outmode.(outmode.modes{indm}).(derv);
          idim = length(size(tempvar));
          % uses a nanmask to preserve nans - here adpts minimum nan number
          tmpnan = squeeze(all(isnan(tempvar),idim));
          tempvar1 = nansum(tempvar,idim);
          tempvar1(tmpnan) = nan;
          outn = outmode.(outmode.modes{indm}).var_outn{ider(indv)}; 
          outmode.(outmode.modes{indm}).(outn) = tempvar1;
       otherwise
          error(['processing mode ' vproc ' not valid']);
       end
    end
 end
      
%-----------------------------------------------------------------------------------------
% Adds in output structure
%-----------------------------------------------------------------------------------------
 boats.output = outmode;

%-----------------------------------------------------------------------------------------
% Save forcing arrays
%-----------------------------------------------------------------------------------------
 if (strcmp(MAIN.sim_type,'hd')||strcmp(MAIN.sim_type,'hf'))
   boats.forcing_used.catchability = catchability_used;
   boats.forcing_used.price        = price_used;
   boats.forcing_used.cost_effort  = cost_effort_used;
 end

%-----------------------------------------------------------------------------------------
% Save time taken to run this script
%-----------------------------------------------------------------------------------------
  runtime_integrate       = toc;
  boats.runtime_integrate = runtime_integrate;
 
 end
 
%**************************************************************************************************************
% END OF FUNCTION


%**************************************************************************************************************
% SUBFUNCTION
%**************************************************************************************************************
 function out = integrate_2d(tvar,area);
    tmp = area .* tvar;
    out = nansum(tmp(:));
 end
 
%**************************************************************************************************************
% END OF SCRIPT
