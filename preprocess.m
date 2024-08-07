%**************************************************************************
% BOATS PRE-PROCESS
% Definition of the forcing for the BOATS model ('run_boats.m'). 
% This preprocessing step allows to load and convert different forcing
% sources (Observations/ROMS/...) in order to convert them into suitable
% inputs for BOATS.
%**************************************************************************
% BOATS inputs generated :
% Ecological input 'Ecological.mat' (for no harvest & harvest simulations):
%   - mask [boolean] (mask of land = 1 and ocean = 0 cells in the domain)
%   - lon/lat [degrees] (coordinates of the computation nodes)
%   - surface [m^2] (surface of the computation cells for integration)
%   - npp [mmolC m^-2 s^-1] (net primary production)
%   - npp_ed [mmolC m^-3 d^-1] (net primary production averaged on euphotic zone depth)
%   - pfb [mmolC m^-2 s^-1] (particule flux at bottom)
%   - temperature pelagic / demersal [deg C] (water temperature)
% Economic input 'Economic.mat' (for harvest simulations):
%   - price [$ g^-1] (fish price history)
%   - cost [$ W^-1] (exploitation cost history)
%   - catchability [m^2 W^-1 s^-1] (catchability history)
%   - effEtar [s-1] (effective effort target)
%   - societenf [??] (societal enforcement)
%**************************************************************************
addpath('preprocess')
clear all


%**************************************************************************
% DEFINE FORCING CHARACTERISTICS
%**************************************************************************
% Preprocess options *******************************
plot_input = 1;                                     % (yes 1 or no 0)
create_ecology = 0;                                 % (yes 1 or no 0)
create_economy = 0;                                 % (yes 1 or no 0)
create_regulation = 1;                              % (yes 1 or no 0)

% General forcing paths and characteristics ********
nlat = 180;                                         % nlat = m nodes along latitude
nlon = 360;                                         % nlon = n nodes along longitude
ntime = 12;                                         % ntime = t distinct time steps per forcing
ngroup = 3;                                         % ngroup = g fish groups per forcing (facultativ)
nensemble = 1;                                      % nensemble = e ensembles per forcing (facultativ)
% mask
mask_path = 'frc/mask_notlme_high_nan.mat';                     % Path of forcing dataset where is the mask
mask_var = 'mask_notlme_high_nan';                 % Name of mask variable in mask_path
%mask_path = 'frc/geo_time.mat';                     % Path of forcing dataset where is the mask
%mask_var = 'geo_time.mask_land_2d';                 % Name of mask variable in mask_path
wet=1;                                              % Value for ocean cells in mask_var
dry=NaN;                                            % Value for continent cells in mask_var
%wet=0;                                              % Value for ocean cells in mask_var
%dry=1;                                            % Value for continent cells in mask_var
mask_dim = [nlat nlon 1 1 1];                       % Dimension of the mask forcing generated
% coordinates lon/lat
lon_path = 'frc/geo_time.mat';                      % Path of forcing dataset where is the longitude
lon_var = 'geo_time.lon';                           % Name of longitude variable in lon_var
lat_path = 'frc/geo_time.mat';                      % Path of forcing dataset where is the latitude
lat_var = 'geo_time.lat';                           % Name of latitude variable in lat_var
lon_dim = [nlat nlon 1 1 1];                        % Dimension of the longitude array generated
lat_dim = [nlat nlon 1 1 1];                        % Dimension of the latitude array generated
% surface
surf_path = 'frc/geo_time.mat';                     % Path of forcing dataset where is the surface
surf_var = 'geo_time.surf';                         % Name of surface variable in surf_path
surf_unit = '[m^2]';                                % surf_var unit ([km^2] or [m^2] or ??)
surf_dim = [nlat nlon 1 1 1];                       % Dimension of the surface array generated

% Ecological forcing paths and characteristics *****
% depth
depth_path = 'frc/data_monthly_orig.mat';           % Path of forcing dataset where is the depth
depth_var = 'data_monthly.depth';                   % Name of depth variable in depth_path
depth_type = 'value';                               % Type of depth definition (map or value)
depth_unit = '[m]';                                 % depth_var unit ([m] or [km] or ??)
ed=75;                                              % User defined depth of euphotic zone
depth_dim = [nlat nlon 1 1 1];                      % Dimension of the depth array generated
% npp
npp_path = 'frc/data_monthly_orig.mat';             % Path of forcing dataset where is the primary production
npp_var = 'data_monthly.npp';                       % Name of primary porduction variable in npp_path
npp_unit = '[mmolC m^-2 d^-1]';                     % npp_var unit ([mgC m^-2 d^-1] or ??)
npp_dim = [nlat nlon ntime 1 1];                    % Dimension of the npp array generated
% pfb
b 	= -0.8;					    % Exponent for Martin's curve attenuation (typical value = -0.8)
% Iron
no3min_path = 'frc/woa05_min_no3.mat';              % Path of forcing dataset where is min no3
no3min_var = 'min_no3';                             % Name of min no3 variable in npp_path
no3min_unit = '[??]';                               % no3min_var unit (??)
no3min_dim = [nlat nlon 1 1 1];                     % Dimension of the no3 min array generated
% Topo, high and/or low resolution
depth_path_high    = 'frc/ETOPO_depth.mat';         % Path of forcing dataset where is the depth at high resolution 
depth_var_high 	   = 'ETOPO';                       % Name of depth variable in depth_path_high
lon_depth_var_high = 'LON';			           % Longitude of depth variable in depth_path_high
lat_depth_var_high = 'LAT';			    % Latitude of depth variable in depth_path_high
depth_dim_high = [1800 3600 1 1 1];	            % Dimension of the high resolution depth array generated
depth_path_low    = 'frc/data_monthly_orig.mat';         % Path of forcing dataset where is the depth at high resolution 
depth_var_low 	   = 'data_monthly.etopo';                       % Name of depth variable in depth_path_high
depth_dim_low = [180 360 1 1 1];	            % Dimension of the high resolution depth array generated
% temperature pelagic
temp_path_pel = 'frc/data_monthly_orig.mat';        % Path of forcing dataset where is the pelagic temperature
temp_var_pel = 'data_monthly.temp75';               % Name of temperature variable in temp_path         
temp_unit_pel = '[degC]';                           % temp_var unit ([degC])
temp_dim_pel = [nlat nlon ntime 1 1];               % Dimension of the temp array generated
% temperature demersal
temp_path_dem = 'frc/temp_bot.mat';                     % Path of forcing dataset where is the demersal temperature
temp_var_dem = 'temp_bot_c';                          % Name of temperature variable in temp_path
temp_unit_dem = '[degC]';                           % temp_var unit ([degC])
temp_dim_dem = [nlat nlon ntime 1 1];               % Dimension of the temp array generated

% Economic forcing paths and characteristics *****
% price
price_path = '/Users/jguiet/OneDrive - University of California/BOATS/input/economy/price_forcing.mat'; % Path of forcing dataset where is the price
price_var = 'price_forcing.price_average_pd_1850_2100';                                 % Name of price variable in price_path or user defined (udef)
price_type = 'cst';                                 % Type of user defined price if price_var = 'udef' (constant cst or ??)
price_ref1 = 1.264079636832035e-04;                 % First parameter for user defined price forcing
price_dim=[1 1 3000 1 1];                           % Dimension of price forcing (out)
price_unit = '[$ g^-1]';                            % price_var unit ([$ g^-1] or ??)
% cost
cost_path = '/Users/jguiet/OneDrive - University of California/BOATS/input/economy/cost_effort_forcing.mat'; % Path of forcing dataset where is the cost
cost_var = 'cost_effort_forcing.cost_effort_average_pd_1850_2100';                      % Name of cost variable in cost_path or user defined (udef)
cost_type = 'cst';                                  % Type of user defined cost if cost_var = 'udef' (constant cst or ??)
cost_ref1 = 1.8520e-07;                             % First parameter for user defined cost forcing
cost_dim=[1 1 3000 1 1];                            % Dimension of user defined cost forcing
cost_unit = '[$ W^-1]';                             % cost_var unit ([$ W^-1] or ??)
% catchability
catch_path = '/Users/jguiet/OneDrive - University of California/BOATS/input/economy/catchability_forcing.mat'; % Path of forcing dataset where is the catchability
catch_var = 'catchability_forcing.catch_p05_250y';  % Name of catchability variable in catch_path or user defined (udef)
catch_type = 'rate';                                % Type of user defined cost if cost_var = 'udef' (constant cst or rate or ??)
catch_ref1 = 7.6045e-08;                            % First parameter for user defined catchability forcing
catch_ref2 = 0.05;                                  % Second parameter for user defined catchability forcing
catch_dim=[1 1 3000 1 1];                           % Dimension of user defined catchability forcing
catch_unit = '[m^2 W^-1 s^-1]';                     % catch_var unit ([m^2 W^-1 s^-1] or ??)

% Regulation forcing paths and characteristics *****
% MSY/MEY effective effort target
effEtarg_path = '/Users/jguiet/OneDrive - University of California/BOATS/input/regulation/E_eff_MSY_true_14028B_.mat'; % Path of forcing dataset where is the effort target
effEtarg_var = 'E_eff_MSY_true';
effEtarg_dim = [nlat nlon 1 2*ngroup nensemble];      % Dimension of the effort target array generated
effEtarg_unit = '[s^-1]';                           % efftarg_var unit ([s^-1])
% Societal enforcement
societenf_path = 'frc/Se_mask_CPI.mat'; % Path of forcing dataset where is the societal enforcement
societenf_var = 'udef';                             % Name of catchability variable in socenf_path ('Se_mask') or user defined (udef)
societenf_type = 'cst';                             % Type of user defined socenf if socenf_var = 'udef' (constant cst or rate or ??)
societenf_ref1 = 1;                                 % First parameter for user defined socenf forcing (e.g. if Se constant = 1, put this here!
societenf_dim = [nlat nlon 1 1 1];
societenf_unit = '[]';
%**************************************************************************
% END DEFINE FORCING CHARACTERISTICS
%**************************************************************************


%**************************************************************************
% LOAD and CONVERT
%**************************************************************************
disp('Lets go !!!')

% Load and convert ecological forcing *************
if create_ecology
    disp('Create ecological forcing')
    % Create mask *********************************
    mask=get_var(mask_path,mask_var,mask_dim);
    mask=arrange_mask(mask,wet,dry);
    % Create coordinates lon/lat ******************
    lon=get_var(lon_path,lon_var,lon_dim);
    lat=get_var(lat_path,lat_var,lat_dim);
    %ATT JG [lon lat]=arrange_coord(lon,lat);
    % Create surface ******************************
    surface=get_var(surf_path,surf_var,surf_dim);
    surface=arrange_surf(surface,surf_unit,lon,lat);
    % Create npp *********************************
    npp=get_var(npp_path,npp_var,npp_dim);
    depth=get_var(depth_path,depth_var,depth_dim);
    [npp npp_ed]=arrange_npp(npp,npp_unit,depth,depth_unit,npp_dim,ed,depth_type);
    % Create temperature *************************
    temperature_pel=get_var(temp_path_pel,temp_var_pel,temp_dim_pel);
    temperature_dem=get_var(temp_path_dem,temp_var_dem,temp_dim_dem);
    depth_low=get_var(depth_path_low,depth_var_low,depth_dim_low);
    temperature_dem=correctdem(temperature_pel,temperature_dem,depth_low,ed);
    % Create pfb *********************************
    depth_high=get_var(depth_path_high,depth_var_high,depth_dim_high);
    tmp = depth_high;
    depth_high(:,1801:3600)=tmp(:,1:1800);
    depth_high(:,1:1800)=tmp(:,1801:3600);
    lon_high=get_var(depth_path_high,lon_depth_var_high,depth_dim_high);
    lon_high =lon_high + 180;
    lat_high=get_var(depth_path_high,lat_depth_var_high,depth_dim_high);
    martin_att  = martin_attenuation(depth_high,ed,b);
    martin_att_bin = bin_var(martin_att,lon_high,lat_high,lon,lat);
    pep     = particle_export(npp,temperature_pel,ed);
    pfb		= pep.*repmat(martin_att_bin,[1,1,size(npp,3)]);
    % Create te_no3 *********************************
    no3min=get_var(no3min_path,no3min_var,no3min_dim);
    % Plot forcings ******************************
    if plot_input
        % surface
        plot_domain2D(surface,lon,lat,mask,'surface [m^2]',1)
        % npp
        plot_domain2D(npp,lon,lat,mask,'npp [mmolC m^{-2} s^{-1}]',2)
        plot_domain2D(npp_ed,lon,lat,mask,'npp_{ed} [mmolC m^{-3} d^{-1}]',4)
        % temperature
        plot_domain2D(temperature_pel,lon,lat,mask,'temperature pelagic [deg C]',3)
        plot_domain2D(temperature_dem,lon,lat,mask,'temperature demersal [deg C]',3)
        % pfb
        plot_domain2D(pfb,lon,lat,mask,'pfb [mmolC m^{-2} s^{-1}]',2)
    end
    % Save forcing
    Ecological.mask=mask;
    Ecological.lon=lon;
    Ecological.lat=lat;
    Ecological.surface=surface;
    Ecological.npp=npp;
    Ecological.npp_ed=npp_ed;
    Ecological.pfb=pfb;
    Ecological.no3min=no3min;
    Ecological.temperature_pel=temperature_pel;
    Ecological.temperature_dem=temperature_dem;
    Ecological.depth=depth_low;
    save('Ecological.mat','Ecological','-v7.3')
end

% Load and convert economic forcing *************  
if create_economy
    disp('Create economic forcing')
    % Create price ********************************
    price=get_var(price_path,price_var,price_dim);
%     price=udef_var(price_type, price_dim, price_ref1);
    % ATT JG arrange
    % Create price ********************************
    cost=get_var(cost_path,cost_var,cost_dim);
%     cost=udef_var(cost_type, cost_dim, cost_ref1);
    % ATT JG arrange
    % Create price ********************************
    catchability=get_var(catch_path,catch_var,catch_dim);
%     catchability=udef_var(catch_type, catch_dim, catch_ref1, catch_ref2);
    % ATT JG arrange 
    % Plot forcings ******************************
    if plot_input
        % price
        plot_domain1D(1:length(price),price,'time','price [$ g^-1]',5)
        % cost
        plot_domain1D(1:length(cost),cost,'time','cost [$ W^-1]',6)
        % catchability
        plot_domain1D(1:length(catchability),catchability,'time','catchability [m^2 W^-1 s^-1]',7)
    end
    % Save forcing
    Economic.price=price;
    Economic.cost=cost;
    Economic.catchability=catchability;
    save('Economic.mat','Economic','-v7.3')
end

% Load and convert regulation forcing ************  
if create_regulation
    disp('Create regulation forcing')
    % Create effort target ***********************
    effEtarg=get_var(effEtarg_path,effEtarg_var,effEtarg_dim);
    % Create social enforcement ******************
    societenf=get_var(societenf_path,societenf_var,societenf_dim);         % SOCIETAL ENFORCEMENT AS A MAP
    societenf=udef_var(societenf_type, societenf_dim, societenf_ref1);     % SOCIETAL ENFORCEMENT AS A CONSTANT    
    % Plot forcings ******************************
    if plot_input
        plot_domain1D(1:length(societenf),societenf,'time','societenf [m^2 W^-1 s^-1]',8)
    end
    % Save forcing
    if exist('frc/Economic.mat')==2
        load('frc/Economic.mat')
    end
    Economic.effEtarg=effEtarg;
    Economic.societenf=societenf;
    save('frc/Economic.mat','Economic','-v7.3')
end
disp('Done !')



