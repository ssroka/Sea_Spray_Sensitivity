clear;close all;clc
%{

Author: Sydney Sroka

This script calculates the enthalpy flux and spray CK for a specified set
of environmental conditions
              > relative humidity RH [%]
              > air-sea temperature differnce DT [K]
              > 10-m wind speed U10 [m/s]
Usage:
The user can specify the input parameters in the USER INPUT section below
and the resultant spray-induced enthalpy flux and spra-CK will be
calculated.

%}

%% USER INPUT

% fluid parameters
% rho_w = 1020;    % saltwater density [kg m^-3]
% s0    = 0.034;   % salinity in [kg/kg]
% cp    = 4000;    % specific heat capacity of saltwater [J kg^-1 K^-1]
% Lv    = 2434054; % latent heat of vaporization [J /kg]
% R     = 8.31447; % universal gas constant [J mol^-1 K^-1]
Nayar_flag = true;

% parameters for time-of-flight Newton-Raphson calculation
maxEr_uf = 4e-4; % maximum error in terminal velocity (uf) in [m/s]
maxIt    = 1000; % maximum number of iterations to execute to calc uf

maxEr_s  = 1e-6; % maximum error in salinity (s) in [kg/kg]

SSGF_str_cell = {'troit','OS','Zhao'};

options = odeset('reltol',1e-5,'abstol',1e-5);

MATLAB_colors = [...
    0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250];

nSSGF = length(SSGF_str_cell);

% parameters corresponding to SGF from Troitskaya et al. 2018a
eq_N = 9; % eq_N \in [5,7,9] - which equation to use for the number of bags
% 5 = derived equation
% 7 = laboratory conditions
% 9 = field conditions
Omega = 2.5; % wave-age parameter (usually between 2.5- 3.5)

% environmental conditions

% _vec quantities will be looped over holding the other parameters
% constant at their _default values
RH_vec  = [88 92 96];    % [%]
DT_vec  = [1.0 2.0 3.0];  % [K]
SST_vec  = [27 28 29];  % [K]
n_t_vec = 100; % time vector to interpolate the velocity of drop on to
t_vec = logspace(-8,4,n_t_vec);


% approximation flag
% true => use the fitted equation
use_approx_flag = false;
% false => use the full microphysical model formula
%          WARNING if set to false and using results in SST_C_29
%                  r0 MUST be \in [50:50:1200]
%                  DT MUST be \in [1:1:3]
%                  RH MUST be \in [80:5:95]
%          WARNING if set to false and using results in SST_C_27
%                  r0 MUST be \in [100:50:1200]
%                  DT MUST be \in [0.5:0.5:3.5]
%                  RH MUST be \in [88:2:96]


tof_vec = 1.0;%[0.1:0.1:2.0];

ud0 = 0; % initial condition
nu_a = 1.48e-5; % kinematic viscosity of air
rho_a = 1.2;
rho_s = 1022;

t0 = 20; % initial guess for re-entry time
subplotLetter = 'abcd';
%% ---------------- Begin -------------------------------------------------

action_str = 'add';
paths_for_calc_CK;

helpingAnonFxns;

nDT = length(DT_vec);
nRH = length(RH_vec);
nSST = length(SST_vec);
ntof = length(tof_vec);

r0_vec_sort = [31 2 13 24 25 26 27 28 29 30 32 33 34 35 36 37 38 39 40 1 3 4 5 6 7 8 9 10 11 12 14 15 16 17 18 19 20 21 22 23];
r0_vec_ALL = [50:50:2000];

r0_vec_sort = [1 2 3];
r0_vec_ALL = [400:100:600];

r0_vec = [400 500 600]; % initial drop radii [micron]

% set defaults

RH_def = RH_vec(2);
DT_def = DT_vec(2);
SST_def = SST_vec(2);
r0_def = r0_vec(2);

nr0 = length(r0_vec);

% parameters that are fixed given the directory
% arbitrarily select the first values in r0_vec, DT_vec, and RH_vec
TH = sprintf('/Users/ssroka/Documents/MATLAB/mpm/sandbox/results/res_fm_Eng/SST_27_28_29_10000/SST_27/timehistory_%dK_%d_RH',DT_vec(1)*10,RH_vec(1));
load(TH,'timehistory'); % will instantiate a variable called timehistory

s0 = timehistory(1).ic.s0;        % [kg/kg]
p0 = timehistory(1).ic.p0;        % [Pa]


Teq_mat = zeros(nr0,nRH,nDT,nSST);

t_interp = [logspace(-2,log10(5000),50) logspace(log10(6e3),log10(70e3),20)];

c = [45 149 165; 240 97 35; 95 42 129]./256;
RH_count = 1;
DT_count = 1;
SST_count = 1;
r0_count = 1;

U10 = 54;

for tof_ind = 1:ntof
    tof_fac = tof_vec(tof_ind);
    tic
    for SST_ind = 1:nSST
        SST_C = SST_vec(SST_ind);
        % microphysical model data is stored
        TH_src = sprintf('/Users/ssroka/Documents/MATLAB/mpm/sandbox/results/res_fm_Eng/SST_27_28_29_10000/SST_%d',SST_C);
        
        for DT_ind = 1:nDT
            DT = DT_vec(DT_ind);
            fprintf('DT = %4.1f  -----------------------\n',DT)
            
            % parameters that are fixed given DT
            % arbitrarily select the first values in r0_vec and RH_vec
            TH = sprintf('%s/timehistory_%dK_%d_RH',TH_src,DT*10,RH_vec(1));
            load(TH,'timehistory'); % will instantiate a variable called timehistory
            
            T_a = timehistory(1).ic.T_a;      % [deg C]
            
            for RH_ind = 1:nRH
                RH = RH_vec(RH_ind);
                fprintf('RH = %4.1f\n\n',RH)
                
                % parameters that are fixed given RH
                TH = sprintf('%s/timehistory_%dK_%d_RH',TH_src,DT*10,RH);
                load(TH,'timehistory'); % will instantiate a variable called timehistory
                
                
                
                for r0_ind = 1:nr0
                    
                    r0_file = r0_vec_sort(round(r0_vec(r0_ind))==r0_vec_ALL);
                    
                    r0_m = r0_vec(r0_ind)*1e-6;              % convert to [m]
                    
                    m_s = timehistory(r0_file).ic.m_s;      % [kg]
                    t_full = timehistory(r0_file).time_vec;
                    r_full = timehistory(r0_file).r_t;
                    T_full = timehistory(r0_file).T_s_t;
                    
                    
                    tauf(r0_ind) = compute_tauf(U10,SST_C,r0_m,m_s,s0,p0,T_a,maxEr_uf,maxIt); % meters
                    
                    
                    if (tof_ind == 1) && (RH_ind == 1) && (DT_ind == 1) && (SST_ind == 1)
                    end
                    if (tof_ind == 1)
                        [Teq_mat(r0_ind,RH_ind,DT_ind,SST_ind),indmin] = min(Celsius2Kelvin(T_full));
                        %                                 [~,ind_tauT] = min(abs(((T_full(1:indmin)-(Teq_mat(r0_ind,RH_ind,DT_ind,SST_ind)-273.15))./(T_full(1)+273.15-Teq_mat(r0_ind,RH_ind,DT_ind,SST_ind)))-exp(-1)));
                        %                                 tauT_mat(r0_ind,RH_ind,DT_ind,SST_ind) = t_full(ind_tauT);
                    end
                    
                    %                             tauT = tauT_mat(r0_ind,RH_ind,DT_ind,SST_ind);
                    T_eq = Teq_mat(r0_ind,RH_ind,DT_ind,SST_ind);
                    
                    T_interp = interp1(t_full,T_full,t_interp);
                    r_interp = interp1(t_full,r_full*1e6,t_interp); % convert to microns
                    
                    T_K = Celsius2Kelvin(T_full);
                    
                    % check defaults
                    RH_def_flag = (RH_vec(RH_ind) == RH_def);
                    DT_def_flag = (DT_vec(DT_ind) == DT_def);
                    SST_def_flag = (SST_vec(SST_ind) == SST_def);
                    r0_def_flag = (r0_vec(r0_ind) == r0_def);
                    
                    if DT_def_flag && SST_def_flag && r0_def_flag
                        subplot(2,2,1) % RH
                        yyaxis left
                        h(RH_count,1,1)= loglog(t_interp,T_interp,':',...
                            'displayname',sprintf('$T$, $RH = $ %d \\%%',...
                            RH_vec(RH_ind)),'linewidth',3,'color',c(RH_count,:));
                        hold on
                        yyaxis right
                        h(RH_count,2,1)= loglog(t_interp,r_interp,'-',...
                            'displayname',sprintf('$r,$ $RH = $ %d \\%%',...
                            RH_vec(RH_ind)),'linewidth',3,'color',c(RH_count,:));
                        hold on
                        RH_count = RH_count + 1;
                        loglog([1 1]*tauf(r0_ind),[200 610],'--','color',[1 1 1]*0.1,'linewidth',3)
                        
                    end
                    if RH_def_flag && SST_def_flag && r0_def_flag
                        subplot(2,2,2) % DT
                        yyaxis left
                        h(DT_count,1,2)= loglog(t_interp,T_interp,':',...
                            'displayname',sprintf('$T,$ $\\Delta T = $%d $^{o}$C',...
                            DT_vec(DT_ind)),'linewidth',3,'color',c(DT_count,:));
                        hold on
                        yyaxis right
                        h(DT_count,2,2)= loglog(t_interp,r_interp,'-',...
                            'displayname',sprintf('$r,$ $\\Delta T = $%d $^{o}$C',...
                            DT_vec(DT_ind)),'linewidth',3,'color',c(DT_count,:));
                        hold on
                        DT_count = DT_count + 1;
                        loglog([1 1]*tauf(r0_ind),[200 610],'--','color',[1 1 1]*0.1,'linewidth',3)
                        
                    end
                    if RH_def_flag && DT_def_flag && r0_def_flag
                        subplot(2,2,3) % SST
                        yyaxis left
                        h(SST_count,1,3)= loglog(t_interp,T_interp,':',...
                            'displayname',sprintf('$T,$ $T_s = $%d $^{o}$C',...
                            SST_vec(SST_ind)),'linewidth',3,'color',c(SST_count,:));
                        hold on
                        yyaxis right
                        h(SST_count,2,3)= loglog(t_interp,r_interp,'-',...
                            'displayname',sprintf('$r,$ $T_s = $%d $^{o}$C',...
                            SST_vec(SST_ind)),'linewidth',3,'color',c(SST_count,:));
                        hold on
                        SST_count = SST_count + 1;
                        loglog([1 1]*tauf(r0_ind),[200 610],'--','color',[1 1 1]*0.1,'linewidth',3)
                        
                    end
                    if RH_def_flag && DT_def_flag && SST_def_flag
                        subplot(2,2,4) % r0
                        yyaxis left
                        h(r0_count,1,4) = loglog(t_interp,T_interp,':',...
                            'displayname',sprintf('$T,$ $r_0 = $%d $$\\mu$$m',...
                            r0_vec(r0_ind)),'linewidth',3,'color',c(r0_count,:));
                        hold on
                        yyaxis right
                        h(r0_count,2,4) = loglog(t_interp,r_interp,'-',...
                            'displayname',sprintf('$r,$ $r_0 = $%d $$\\mu$$m',...
                            r0_vec(r0_ind)),'linewidth',3,'color',c(r0_count,:));
                        hold on
                        loglog([1 1]*tauf(r0_ind),[200 610],'--','color',c(r0_count,:),'linewidth',3)
                        r0_count = r0_count +1;
                    end
                    
                    
                    clear t_full r_full T_full
                    
                end
                
                
                
            end
        end
    end
    toc
    fprintf('completed tof %f\n',tof_ind)
end
for i = 1:4
    subplot(2,2,i)
    yyaxis left
    set(gca,'fontsize',20,'YColor','k','ytick',25:29,'ylim',[24 29],...
        'xtick',sort([10.^[0 2 5] tauf(2)]),'xticklabel',{'$10^0$','$\tau_f$','$10^2$','$10^5$'},...
        'TickLabelInterpreter','latex')
    ylh = ylabel('[$^o$C]');
    set(ylh,'rotation',0,'interpreter','latex',...
        'HorizontalAlignment','left','VerticalAlignment','bottom',...
        'units','normalized','position',[-0.1 1.1 0])
    yyaxis right
    set(gca,'fontsize',20,'YColor','k','ytick',300:100:600,'ylim',[200 610])
    ylh = ylabel('[$\mu$m]');
    set(ylh,'rotation',0,'interpreter','latex',...
        'HorizontalAlignment','left','VerticalAlignment','bottom',...
        'units','normalized','position',[1 1.1 0])
    xlabel('$t$ [s]','interpreter','latex')
    lh(i) = legend([h(1,1,i) h(2,1,i) h(3,1,i) h(1,2,i) h(2,2,i) h(3,2,i)],...
        'location','eastoutside','interpreter','latex');
    th = text(0.1,0.1,sprintf('(%s)',subplotLetter(i)));
    set(th,'units','normalized','position',[0.01 0.07],'fontsize',20)
    drawnow
end
set(gcf,'position',[132          31        1259         767],'color','w')

addpath('/Users/ssroka/Documents/MATLAB/util/')
H = gcf;
drawnow
rearrange_figure(H,lh,'2x2_4_legends')

update_figure_paper_size()
print(sprintf('imgs/compare_drops'),'-dpdf')



function [Cd] = get_Cd(Re)

w = log10(Re);

if Re < 0.01
    
    Cd = 3/16+24/Re;
    
elseif Re <= 20
    
    Cd = 24/Re*(1+0.1315*Re^(0.82-0.05*w));
    
elseif Re<= 260
    
    Cd = 24/Re*(1+0.1935*Re^(0.6305));
    
elseif Re <= 1500
    
    Cd = 10.^(1.6435-1.1242*w+0.1558*w.^2);
    
elseif Re <= 1.2e4
    
    Cd = 10.^(-2.4571+2.5558*w-0.9295*w.^2+0.1049.*w.^3);
    
else
    
    Cd = 10.^(-1.9181+0.637*w-0.0636*w.^2);
    
end

end
