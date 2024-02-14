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
forcing.depth(find(forcing.depth>-1)) = -1;
%forcing.zeuph=Ecological.zeuph;
forcing.dist=Ecological.dist;
forcing.surf=Ecological.surface;

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
%      [forcing.zeuph_vec(:,itime) forcing.indlat forcing.indlon]           = function_map_2_vec(squeeze(forcing.zeuph(:,:,itime)),squeeze(forcing.mask(:,:,1)));
  end % itime
  [forcing.no3min_vec forcing.indlat forcing.indlon]              = function_map_2_vec(forcing.no3min,squeeze(forcing.mask(:,:,1)));
  [forcing.surf_vec forcing.indlat forcing.indlon]                = function_map_2_vec(forcing.surf,squeeze(forcing.mask(:,:,1)));
  [forcing.depth_vec forcing.indlat forcing.indlon]               = function_map_2_vec(forcing.depth,squeeze(forcing.mask(:,:,1)));
  [forcing.dist_vec forcing.indlat forcing.indlon]               = function_map_2_vec(forcing.dist,squeeze(forcing.mask(:,:,1)));
  forcing.nvec=size(forcing.surf_vec,1);

if (boats.param.economy.costpar)
  % Define a depth depedent profile to exploit demersal catch
  forcing.cost_profile1 = forcing.depth_vec*NaN;
  forcing.cost_profile2 = forcing.depth_vec*NaN;

  ind_cst = find(-forcing.depth_vec<=200);
  ind_var = find(-forcing.depth_vec >200);
  forcing.cost_profile1(ind_cst) = 1;
  forcing.cost_profile1(ind_var) = 1 + 0.5/5.85/200 * (-forcing.depth_vec(ind_var)-200);

  ind_cst = find(forcing.dist_vec<=370);
  ind_var = find(forcing.dist_vec >370);
  forcing.cost_profile2(ind_cst) = 1;
  forcing.cost_profile2(ind_var) = 1 + 0.5*5.85/5.85/370 * (forcing.dist_vec(ind_var)-370);

  % LINEAR
  ind_cst = find(-forcing.depth_vec<=200);
  ind_var = find(-forcing.depth_vec >200);
  forcing.cost_profile(ind_cst) = 1;
  forcing.cost_profile(ind_var) = 1 + 0.5./5.85/200 * (-forcing.depth_vec(ind_var)-200); 
end

if (boats.param.economy.catchpar)
  % define a depth dependent catchability
  forcing.catcha_profile = forcing.depth_vec*NaN;
  
  %log10(zbot)
  zmax = 3.76;
  zmean = 3.38;
  qmin  = 0.8;
  forcing.catcha_profile = qmin + (1-qmin) * (zmax - log10(-forcing.depth_vec))./(zmax-zmean);

%  %(zbot)
%  zmax = 7473;
%  zmean = 3810;
%  qmin  = 0.8;
%  forcing.catcha_profile = qmin + (1-qmin) * (zmax + forcing.depth_vec)./(zmax-zmean);

%  %(1/zeu)
%  zmax  = 170;
%  zmean = 57.0527;
%  qmin  = 0.1;
%  forcing.catcha_profile = qmin + (1-qmin) * (1./forcing.zeuph_vec-1/zmax)./(1/zmean-1/zmax);

%  %(zeu)
%  zmax  = 170;
%  zmean = 65;
%  qmin  = 0.1;
%  forcing.catcha_profile = qmin + (1-qmin) * (zmax - forcing.zeuph_vec)./(zmax-zmean);

  forcing.catcha_profile(find(isinf(forcing.catcha_profile))) = 5*qmin;
  forcing.catcha_profile(find(forcing.catcha_profile>5*qmin)) = 5*qmin;
end

%**************************************************************************************************************
% END FUNCTION

