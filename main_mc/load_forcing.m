%**************************************************************************************************************
% FUNCTION load_forcing.m
% Load forcings preprocessed by "preprocess.m" for the simulation :
% Ecological.mat
% Economical.mat
%**************************************************************************************************************
function forcing = load_forcing(boats,forcing_ecological,forcing_economical)

%---------------------------------
% Forcing ecological:
if exist(forcing_ecological,'file')
    load(forcing_ecological);
else
    disp('Hum, double-check the path for ecological forcing:')
    disp(forcing_ecological)
end
forcing.mask=repmat(Ecological.mask,[1 1 size(Ecological.npp,3)]);
forcing.nlat=size(forcing.mask,1);
forcing.nlon=size(forcing.mask,2);
forcing.npp=Ecological.npp;
forcing.npp(find(forcing.mask==1))=NaN;
forcing.npp_ed=Ecological.npp_ed;
forcing.npp_ed(find(forcing.mask==1))=NaN;
forcing.pfb=Ecological.pfb;
forcing.pfb(find(forcing.mask==1))=NaN;
forcing.no3min=Ecological.no3min;
forcing.temperature_pel=Ecological.temperature_pel;
forcing.temperature_pel_K=Ecological.temperature_pel+boats.param.conversion.C_2_K;
forcing.temperature_dem=Ecological.temperature_dem;
forcing.temperature_dem_K=Ecological.temperature_dem+boats.param.conversion.C_2_K;
forcing.depth=Ecological.depth;
forcing.surf=Ecological.surface;
forcing.zeuph=Ecological.zeuph;

%--------------------------------- 
% Forcing economical
if (strcmp(boats.param.main.sim_type,'hd')||strcmp(boats.param.main.sim_type,'hf'))
    if exist(forcing_economical,'file')
        load(forcing_economical);
    else
        disp('Hum, double-check the path for economical forcing:')
        disp(forcing_economical)
    end
    load(forcing_economical)
    forcing.price=Economical.price;
    forcing.cost=Economical.cost;
    forcing.catchability=Economical.catchability;
end

 
%---------------------------------
  % Convert maps to vectors
  for itime = 1:size(forcing.npp,3)
      [forcing.npp_vec(:,itime) forcing.indlat forcing.indlon]           = function_map_2_vec(squeeze(forcing.npp(:,:,itime)),squeeze(forcing.mask(:,:,1)));
      [forcing.npp_ed_vec(:,itime) forcing.indlat forcing.indlon]        = function_map_2_vec(squeeze(forcing.npp_ed(:,:,itime)),squeeze(forcing.mask(:,:,1)));
      [forcing.pfb_vec(:,itime) forcing.indlat forcing.indlon]           = function_map_2_vec(squeeze(forcing.pfb(:,:,itime)),squeeze(forcing.mask(:,:,1)));
      [forcing.temperature_pel_vec(:,itime) forcing.indlat forcing.indlon]   = function_map_2_vec(squeeze(forcing.temperature_pel(:,:,itime)),squeeze(forcing.mask(:,:,1)));
      [forcing.temperature_pel_K_vec(:,itime) forcing.indlat forcing.indlon] = function_map_2_vec(squeeze(forcing.temperature_pel_K(:,:,itime)),squeeze(forcing.mask(:,:,1)));
      [forcing.temperature_dem_vec(:,itime) forcing.indlat forcing.indlon]   = function_map_2_vec(squeeze(forcing.temperature_dem(:,:,itime)),squeeze(forcing.mask(:,:,1)));
      [forcing.temperature_dem_K_vec(:,itime) forcing.indlat forcing.indlon] = function_map_2_vec(squeeze(forcing.temperature_dem_K(:,:,itime)),squeeze(forcing.mask(:,:,1)));
      [forcing.zeuph_vec(:,itime) forcing.indlat forcing.indlon]           = function_map_2_vec(squeeze(forcing.zeuph(:,:,itime)),squeeze(forcing.mask(:,:,1)));
  end % itime
  [forcing.no3min_vec forcing.indlat forcing.indlon]              = function_map_2_vec(forcing.no3min,squeeze(forcing.mask(:,:,1)));
  [forcing.surf_vec forcing.indlat forcing.indlon]                = function_map_2_vec(forcing.surf,squeeze(forcing.mask(:,:,1)));
  [forcing.depth_vec forcing.indlat forcing.indlon]               = function_map_2_vec(forcing.depth,squeeze(forcing.mask(:,:,1)));
  forcing.nvec=size(forcing.surf_vec,1);

if (boats.param.economy.depthdep)
  % Define a depth depedent profile to exploit demersal catch
  forcing.depth_profile = forcing.depth_vec*NaN;  

  % LINEAR
  ind_cst = find(-forcing.depth_vec<=100);
  ind_var = find(-forcing.depth_vec >100);
  forcing.depth_profile(ind_cst) = 1;
  forcing.depth_profile(ind_var) = 1 + 1./5.85/100 * (-forcing.depth_vec(ind_var)-100); 

  % EXP
  %forcing.depth_profile = exp((-forcing.depth_vec-100)./966);
end

if (boats.param.economy.zeudep)
  forcing.zeu_profile = forcing.zeuph_vec*NaN;
  %1/zeu
  zmax  = 170;
  zmean = 57.0527;
  qmin  = 0.1;
  forcing.zeu_profile = qmin + (1-qmin) * (1./forcing.zeuph_vec-1/zmax)./(1/zmean-1/zmax);
  % zeu
  %zmax  = 170;
  %zmean = 65;
  %qmin  = 0.1;
  %forcing.zeu_profile = repmat(qmin + (1-qmin) * (zmax - forcing.zeuph_vec)./(zmax-zmean),[1,3]); 
end

%**************************************************************************************************************
% END FUNCTION

