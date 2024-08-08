 function plot_boats_diagnostics(boats,varargin)
%**************************************************************************
A.coast = 0;		% Adds coastline
A.print = 0;		% Print figure
A.mode = 1;		% 1. Only timeseries
			% 2. Timeseries and maps 
A.outmode = 'annual';	% Print figure
A.period = [];
A.ref = [];
%--------------------------------------------------------------------------
addpath /Users/danielebianchi/AOS1/matlabpathfiles/
A = parse_pv_pairs(A, varargin);
%**************************************************************************

 % Clips data to have integer number of years
 out = boats.output.(A.outmode);

 lon = boats.output.lon;
 lat = boats.output.lat;
 year = ceil(out.year);

 fish_g_out = out.fish_g_out;
 fish_t_out = out.fish_t_out;
 fish_gi_g  = out.fish_gi_g;
 fish_gi_t  = out.fish_gi_t;

 % Samples by specific period
 if ~isempty(A.period)
    if length(A.period)==1
       indy1 = findin(A.period(1),year);
       indy2 = length(year);
    else
       indy1 = findin(A.period(1),year);
       indy2 = findin(A.period(2),year);
    end
 else
    indy1 = 1;
    indy2 = length(year);
 end
 year = year(indy1:indy2);
 fish_g_out = fish_g_out(indy1:indy2,:,:,:); 
 fish_t_out = fish_t_out(indy1:indy2,:,:);
 fish_gi_g = fish_gi_g(indy1:indy2,:);
 fish_gi_t = fish_gi_t(indy1:indy2);

 % For plotting, changes reference year
 if ~isempty(A.ref)
    iyref = findin(A.ref(1),year);
    year = year - year(iyref) + A.ref(2); 
 end

 %--------------------------------------------
 % Global maps
 %--------------------------------------------
 fish_t_out_m = log10(squeeze(nanmean(fish_t_out,1)));
 fish_t_out_s = log10(squeeze(nanstd(fish_t_out,1)));
 
 fish_g_out_m = log10(squeeze(nanmean(fish_g_out,1)));
 fish_g_out_s = log10(squeeze(nanstd(fish_g_out,1)));
 
 %--------------------------------------------
 % Makes plots
 %--------------------------------------------
 figure
 switch A.mode
 case 1
    tcl = tiledlayout(2,1);
    fsty = 'nor2';
 otherwise
    tcl = tiledlayout(2,3);
    fsty = 'nor7g';
 end
 ftz = 15;
 %----
 nexttile
 plot(year,fish_gi_t*1e-15,'-','linewidth',4,'color',[0.3 0.3 0.3]);
 title('Total','fontsize',ftz);
 ylabel('Pg','fontsize',ftz);
 xlabel('year','fontsize',ftz);
 xlim([year(1) year(end)]);
 box on;
 %----
 if A.mode==2
    nexttile
    pcolor(lon,lat,fish_t_out_m);
    shading flat
    if A.coast==1;plot_coast;end
    cb = colorbar;
    cb.Ticks = [-2:1:3];
    cb.TickLabels = {'0.01','0.1','1','10','100','1000'}; 
    caxis([-2 3]);
    cb.Label.String = 't km^-^2';
    title('Total','fontsize',ftz);
    %----
    nexttile
    pcolor(lon,lat,fish_g_out_m(:,:,1));
    shading flat
    if A.coast==1;plot_coast;end
    cb = colorbar;
    cb.Ticks = [-2:1:2];
    cb.TickLabels = {'0.01','0.1','1','10','100'}; 
    caxis([-2 2]);
    cb.Label.String = 't km^-^2';
    title('Small','fontsize',ftz);
 end
 %----
 nexttile
 hold on
 plot(year,fish_gi_g(:,1)*1e-15,'-','linewidth',4,'color',[0.8 0.4 0.4]);
 plot(year,fish_gi_g(:,2)*1e-15,'-','linewidth',4,'color',[0.4 0.8 0.4]);
 plot(year,fish_gi_g(:,3)*1e-15,'-','linewidth',4,'color',[0.4 0.4 0.8]);
 legend({'s','m','l'},'Location','NorthEast','fontsize',ftz);
 title('Group','fontsize',ftz);
 ylabel('Pg','fontsize',ftz);
 xlabel('year','fontsize',ftz);
 xlim([year(1) year(end)]);
 box on;
 %----
 if A.mode==2
    nexttile
    pcolor(lon,lat,fish_g_out_m(:,:,2));
    shading flat
    if A.coast==1;plot_coast;end
    cb = colorbar;
    cb.Ticks = [-2:1:2];
    cb.TickLabels = {'0.01','0.1','1','10','100'}; 
    caxis([-2 2]);
    cb.Label.String = 't km^-^2';
    title('Medium','fontsize',ftz);
    %----
    nexttile
    pcolor(lon,lat,fish_g_out_m(:,:,3));
    shading flat
    if A.coast==1;plot_coast;end
    cb = colorbar;
    cb.Ticks = [-2:1:2];
    cb.TickLabels = {'0.01','0.1','1','10','100'}; 
    caxis([-2 2]);
    cb.Label.String = 't km^-^2';
    title('Large','fontsize',ftz);
 end

 % Common title
 sname = boats.param.main.sim_name;
 title(tcl,sname,'fontsize',round(ftz*1.2),'interpreter','none');

 if A.print==1
    fname = ['fig_' sname];
    mprint_fig('name',fname,'for','jpeg','sty',fsty);
 end
