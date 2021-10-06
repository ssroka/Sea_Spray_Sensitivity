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
Troit_waveage_num = 3;
Zhao_waveage_num = 3;
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
if Troit_waveage_num == 1
    Omega = 2.5; % wave-age parameter (usually between 2.5- 3.5)
elseif Troit_waveage_num == 2
    Omega = 3.5; % wave-age parameter (usually between 2.5- 3.5)
elseif Troit_waveage_num == 3
    Omega = 2.5; % used in paper
end
% environmental conditions

% _vec quantities will be looped over holding the other parameters
% constant at their _default values
RH_vec  = [80:2:98];    % [%]
DT_vec  = [0.5:0.5:3.5];  % [K]
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
cw = 4007.2; % [J/kg K] specific heat of sea water at 0.034 s, 1e5 Pa, 29 deg C

t0 = 20; % initial guess for re-entry time
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


for SSGF_ind = 1:nSSGF
    SSGF_str = SSGF_str_cell{SSGF_ind};
    % % parameters common to all environmental conditions
    % SSGF_str = 'troit';% 'troit' 'OS' 'Zhao'
    addpath('/Users/ssroka/MIT/Research/EmanuelGroup/thesis/review/');
    if Zhao_waveage_num == 1
        beta = 1.2; % Zhao only
    elseif Zhao_waveage_num == 2
        beta = 0.2;
    elseif Zhao_waveage_num == 3
        beta = 0.4;
    end
    
    if strcmp(SSGF_str,'Zhao')
        r0_vec = [50:50:500];%[30 500]  % initial drop radii [micron]
        U10_vec = [10:10:80];    % [m/s]
        SSGF_str_title = ['Zhao et al. (2006)' newline '$$r_0 = [50-500]\mu$$m'];
    elseif strcmp(SSGF_str,'troit')
        r0_vec = [50:50:2000];%[50:50:2000];  % initial drop radii [micron]
        U10_vec = [10:10:80];    % [m/s]
        SSGF_str_title = ['Troitskaya et al. (2018)a' newline '$$r_0 = [50-2000]\mu$$m'];
    elseif strcmp(SSGF_str,'OS')
        r0_vec = [100:50:1000];%[86 1036];  % initial drop radii [micron]
        U10_vec = [36 40.5 45 49.5 54];%[10:10:80];    % [m/s]
        SSGF_str_title = ['Ortiz-Suslow et al. (2016)' newline '$$r_0 = [100-1000]\mu$$m'];
    end
    
    nr0 = length(r0_vec);
    nU10 = length(U10_vec);
    
    % parameters that are fixed given the directory
    % arbitrarily select the first values in r0_vec, DT_vec, and RH_vec
    %     TH = sprintf('/Users/ssroka/Documents/MATLAB/mpm/sandbox/results/res_fm_Eng/SST_27_28_29/SST_27/timehistory_%dK_%d_RH',DT_vec(1)*10,RH_vec(1));
    TH = sprintf('/Users/ssroka/Documents/MATLAB/mpm/sandbox/results/res_fm_Eng/SST_27_28_29_fullRH/SST_27/timehistory_%dK_%d_RH',DT_vec(1)*10,RH_vec(1));
    load(TH,'timehistory'); % will instantiate a variable called timehistory
    
    s0 = timehistory(1).ic.s0;        % [kg/kg]
    S0 = timehistory(1).ic.S0;        % [ppt g/kg]
    p0 = timehistory(1).ic.p0;        % [Pa]
    
    ZERO = zeros(nDT,nRH,nU10,length(SST_vec),ntof);
    tof_mat = zeros(nr0,nU10); % the formula appears to depend on temp, but that's just for rho_a, nu_a calcs
    tau_Hsu_mat = zeros(nr0,nU10); % the formula appears to depend on temp, but that's just for rho_a, nu_a calcs
    %     tauT_mat= zeros(nr0,nRH,nDT,nSST);
    Teq_mat = zeros(nr0,nRH,nDT,nSST);
    
    ZERO_frac = zeros(nr0,nU10,nRH,nDT,nSST,ntof);
    
    frac_tauf_T = ZERO_frac;
    frac_taub_T = ZERO_frac;
    frac_tauHsu_T = ZERO_frac;
    frac_tauf_u = ZERO_frac;
    frac_taub_u = ZERO_frac;
    frac_tauHsu_u = ZERO_frac;
    frac_tauf_r = ZERO_frac;
    frac_taub_r = ZERO_frac;
    frac_tauHsu_r = ZERO_frac;
    Qk          = ZERO_frac;
    
    ZERO_ratio = zeros(nU10,nRH,nDT,nSST,ntof);
    
    RT  = ZERO_ratio;
    RbT = ZERO_ratio;
    RHsuT = ZERO_ratio;
    Ru  = ZERO_ratio;
    Rbu = ZERO_ratio;
    RHsuu = ZERO_ratio;
    Rr  = ZERO_ratio;
    Rbr = ZERO_ratio;
    RHsur = ZERO_ratio;
    
    ud_mat = zeros(n_t_vec,nr0,nU10);
    
    dFdr_vec_mat = zeros(nr0,nU10);
    t_b_mat = zeros(nU10,1);
    
    
    for tof_ind = 1:ntof
        tof_fac = tof_vec(tof_ind);
        tic
        for SST_ind = 1:length(SST_vec)
            SST_C = SST_vec(SST_ind);
            % microphysical model data is stored
            %             TH_src = sprintf('/Users/ssroka/Documents/MATLAB/mpm/sandbox/results/res_fm_Eng/SST_27_28_29/SST_%d',SST_C);
            TH_src = sprintf('/Users/ssroka/Documents/MATLAB/mpm/sandbox/results/res_fm_Eng/SST_27_28_29_fullRH/SST_%d',SST_C);
            T_s = timehistory(1).ic.T_s;
            for DT_ind = 1:nDT
                DT = DT_vec(DT_ind);
                fprintf('DT = %4.1f  -----------------------\n',DT)
                
                % parameters that are fixed given DT
                % arbitrarily select the first values in r0_vec and RH_vec
                TH = sprintf('%s/timehistory_%dK_%d_RH',TH_src,DT*10,RH_vec(1));
                load(TH,'timehistory'); % will instantiate a variable called timehistory
                
                T_a = timehistory(1).ic.T_a;      % [deg C]
                T_a_K = Celsius2Kelvin(T_a);
                for RH_ind = 1:nRH
                    RH = RH_vec(RH_ind);
                    fprintf('RH = %4.1f\n\n',RH)
                    
                    % parameters that are fixed given RH
                    TH = sprintf('%s/timehistory_%dK_%d_RH',TH_src,DT*10,RH);
                    load(TH,'timehistory'); % will instantiate a variable called timehistory
                    
                    for U10_ind = nU10:-1:1
                        U10 = U10_vec(U10_ind);
                        if (RH_ind == 1) && (DT_ind == 1) && (SST_ind == 1) && (tof_ind==1)
                            switch SSGF_str
                                case 'troit'
                                    dFdr_vec_mat(:,U10_ind) = calc_Troit_SGF(U10,eq_N,Omega,r0_vec*1e-6); % meters
                                case 'OS'
                                    dFdr_vec_mat(:,U10_ind) = calc_OS_SGF(U10,r0_vec); % meters
                                case 'Zhao'
                                    wp = 
                                    dFdr_vec_mat(:,U10_ind) = Zhao2006(r0_vec,U10,beta); % meters
                            end
                            t_b_mat(U10_ind) = fzero(@(t) U10*t*sin(pi/2)-0.5*9.81*t.^2,t0);
                            %                         V_flux = dFdr_vec_mat(U10_ind,:).*4./3.*pi.*((r0_vec).^3);
                            
                        end
                        t_b = t_b_mat(U10_ind);
                        dFdr_vec = dFdr_vec_mat(:,U10_ind);
                        
                        for r0_ind = 1:nr0
                            
                            r0_file = r0_vec_sort(round(r0_vec(r0_ind))==r0_vec_ALL);
                            
                            r0_m = r0_vec(r0_ind)*1e-6;              % convert to [m]
                            
                            m_s = timehistory(r0_file).ic.m_s;      % [kg]
                            t_full = timehistory(r0_file).time_vec;
                            r_full = timehistory(r0_file).r_t;
                            T_full = timehistory(r0_file).T_s_t;
                            
                            
                            if (tof_ind == 1) && (RH_ind == 1) && (DT_ind == 1) && (SST_ind == 1)
                                [tof_mat(r0_ind,U10_ind),u_f] = compute_tauf(U10,SST_C,r0_m,m_s,s0,p0,T_a,maxEr_uf,maxIt); % meters
                                Hs = 0.0087*(beta^1.86)*(U10^2)/9.8;% see Hsu et. al (2017)
                                tau_Hsu_mat(r0_ind,U10_ind) = 2*Hs/u_f; % meters
                                Re =@(ud) (U10-ud).*(2.*r0_m)./nu_a;
                                %                                 [~,uode] = ode45(@(t,ud) 3/8.*get_Cd(Re(ud)).*rho_a(T_a,p0)./rho_s(T_s,[],[],s0,p0).*(U10 - ud).^2./r0_m,t_vec,ud0,options) ;
                                [~,uode] = ode45(@(t,ud) 3/8.*get_Cd(Re(ud)).*rho_a./rho_s.*(U10 - ud).^2./r0_m,t_vec,ud0,options) ;
                                ud_mat(:,r0_ind,U10_ind) = uode;
                                tau_ac_mat(r0_ind,U10_ind) = interp1(uode,t_vec,U10*(1-exp(-1)));
                            end
                            if (tof_ind == 1) && (U10_ind == nU10)
                                [Teq_mat(r0_ind,RH_ind,DT_ind,SST_ind),indmin] = min(Celsius2Kelvin(T_full));
                                [~,ind_tauT] = min(abs(((T_full(1:indmin)-(Teq_mat(r0_ind,RH_ind,DT_ind,SST_ind)-273.15))./(T_full(1)+273.15-Teq_mat(r0_ind,RH_ind,DT_ind,SST_ind)))-exp(-1)));
                                tauT_mat(r0_ind,RH_ind,DT_ind,SST_ind) = t_full(ind_tauT);
                            end
                            t_Hsu = tau_Hsu_mat(r0_ind,U10_ind);

                            %                             tauT = tauT_mat(r0_ind,RH_ind,DT_ind,SST_ind);
                            T_eq = Teq_mat(r0_ind,RH_ind,DT_ind,SST_ind);
                            tof = tof_fac.*tof_mat(r0_ind,U10_ind);
                            inds = t_full<tof;
                            ufinal_tauf = interp1(t_vec,ud_mat(:,r0_ind,U10_ind),tof);
                            ufinal_taub = interp1(t_vec,ud_mat(:,r0_ind,U10_ind),t_b);
                            ufinal_tauHsu = interp1(t_vec,ud_mat(:,r0_ind,U10_ind),t_Hsu);
                            
                            Tfinal_tauf = interp1(t_full,Celsius2Kelvin(T_full),tof);
                            Tfinal_taub = interp1(t_full,Celsius2Kelvin(T_full),t_b);
                            Tfinal_tauHsu = interp1(t_full,Celsius2Kelvin(T_full),t_Hsu);
                            
                            rfinal_tauf = interp1(t_full,r_full,tof);
                            rfinal_taub = interp1(t_full,r_full,t_b);
                            rfinal_tauHsu = interp1(t_full,r_full,t_Hsu);
                            
                            t = [t_full(inds); tof];
                            
                            T_K = Celsius2Kelvin(T_full);
                            
                            %                             Qk(r0_ind,U10_ind,RH_ind,DT_ind,SST_ind,tof_ind) = rho_s*cw*4/3*pi*r_full(1).^3*(T_K(1)-Tfinal_tauf);
                            m_0 = rho_s*4/3*pi*r_full(1).^3;
                            m_f = rho_s*4/3*pi*rfinal_tauf.^3;
                            Qk(r0_ind,U10_ind,RH_ind,DT_ind,SST_ind,tof_ind) = cw*(m_0*T_K(1)-m_f*Tfinal_tauf-(m_0-m_f)*T_a_K);
                            
                            frac_tauf_T(r0_ind,U10_ind,RH_ind,DT_ind,SST_ind,tof_ind) = abs(Tfinal_tauf-T_K(1))/abs(T_eq - T_K(1));
                            frac_taub_T(r0_ind,U10_ind,RH_ind,DT_ind,SST_ind,tof_ind) = abs(Tfinal_taub-T_K(1))/abs(T_eq - T_K(1));
                            frac_tauHsu_T(r0_ind,U10_ind,RH_ind,DT_ind,SST_ind,tof_ind) = abs(Tfinal_tauHsu-T_K(1))/abs(T_eq - T_K(1));
                            
                            frac_tauf_u(r0_ind,U10_ind,RH_ind,DT_ind,SST_ind,tof_ind) = abs(ufinal_tauf-ud0)/abs(U10 - ud0);
                            frac_taub_u(r0_ind,U10_ind,RH_ind,DT_ind,SST_ind,tof_ind) = abs(ufinal_taub-ud0)/abs(U10 - ud0);
                            frac_tauHsu_u(r0_ind,U10_ind,RH_ind,DT_ind,SST_ind,tof_ind) = abs(ufinal_tauHsu-ud0)/abs(U10 - ud0);
                            
                            frac_tauf_r(r0_ind,U10_ind,RH_ind,DT_ind,SST_ind,tof_ind) = (rfinal_tauf/r0_m).^3;
                            frac_taub_r(r0_ind,U10_ind,RH_ind,DT_ind,SST_ind,tof_ind) = (rfinal_taub/r0_m).^3;
                            frac_tauHsu_r(r0_ind,U10_ind,RH_ind,DT_ind,SST_ind,tof_ind) = (rfinal_tauHsu/r0_m).^3;
                            
                            clear t_full r_full T_full
                            
                        end
                        % N = numerator
                        % D = denominator
                        % b = ballistic
                        % u = momentum
                        
                        N = trapz(r0_vec,squeeze(frac_tauf_T(:,U10_ind,RH_ind,DT_ind,SST_ind,tof_ind)).*dFdr_vec.*(r0_vec'*1e-6).^3);
                        D = trapz(r0_vec,dFdr_vec.*(r0_vec'*1e-6).^3);
                        
                        RT(U10_ind,RH_ind,DT_ind,SST_ind,tof_ind) = N/D;
                        Nb = trapz(r0_vec,frac_taub_T(:,U10_ind,RH_ind,DT_ind,SST_ind,tof_ind).*dFdr_vec.*(r0_vec'*1e-6).^3);
                        RbT(U10_ind,RH_ind,DT_ind,SST_ind,tof_ind) = Nb/D;
                        NHsu = trapz(r0_vec,frac_tauHsu_T(:,U10_ind,RH_ind,DT_ind,SST_ind,tof_ind).*dFdr_vec.*(r0_vec'*1e-6).^3);
                        RHsuT(U10_ind,RH_ind,DT_ind,SST_ind,tof_ind) = NHsu/D;
                        
                        Nu = trapz(r0_vec,frac_tauf_u(:,U10_ind,RH_ind,DT_ind,SST_ind,tof_ind).*dFdr_vec.*(r0_vec'*1e-6).^3);
                        Ru(U10_ind,RH_ind,DT_ind,SST_ind,tof_ind) = Nu/D;
                        Nbu = trapz(r0_vec,frac_taub_u(:,U10_ind,RH_ind,DT_ind,SST_ind,tof_ind).*dFdr_vec.*(r0_vec'*1e-6).^3);
                        Rbu(U10_ind,RH_ind,DT_ind,SST_ind,tof_ind) = Nbu/D;
                        NHsuu = trapz(r0_vec,frac_tauHsu_u(:,U10_ind,RH_ind,DT_ind,SST_ind,tof_ind).*dFdr_vec.*(r0_vec'*1e-6).^3);
                        RHsuu(U10_ind,RH_ind,DT_ind,SST_ind,tof_ind) = NHsuu/D;
                        
                        Nr = trapz(r0_vec,frac_tauf_r(:,U10_ind,RH_ind,DT_ind,SST_ind,tof_ind).*dFdr_vec.*(r0_vec'*1e-6).^3);
                        Rr(U10_ind,RH_ind,DT_ind,SST_ind,tof_ind) = Nr/D;
                        Nbr = trapz(r0_vec,frac_taub_r(:,U10_ind,RH_ind,DT_ind,SST_ind,tof_ind).*dFdr_vec.*(r0_vec'*1e-6).^3);
                        Rbr(U10_ind,RH_ind,DT_ind,SST_ind,tof_ind) = Nbr/D;
                        NHsur = trapz(r0_vec,frac_tauHsu_r(:,U10_ind,RH_ind,DT_ind,SST_ind,tof_ind).*dFdr_vec.*(r0_vec'*1e-6).^3);
                        RHsur(U10_ind,RH_ind,DT_ind,SST_ind,tof_ind) = NHsur/D;
                    end
                end
            end
        end
        toc
        fprintf('completed tof %f\n',tof_ind)
    end
    fprintf('saved ratios for %s\n',SSGF_str)
    save(sprintf('ratios_%s_%d_%d',SSGF_str,Zhao_waveage_num,Troit_waveage_num),...
        'RT','RbT','Ru','Rbu','Rr','Rbr','RHsuT','RHsuu','RHsur','Qk',...
        'frac_tauf_T','frac_taub_T','frac_tauHsu_T',...
        'frac_tauf_u','frac_taub_u','frac_tauHsu_u'...
        ,'frac_tauf_r','frac_taub_r','frac_tauHsu_r',...
        'r0_vec','U10_vec','RH_vec','DT_vec','SST_vec','tof_vec','tau_Hsu_mat',...
        'dFdr_vec_mat','tof_mat','t_b_mat','ud_mat','tau_ac_mat');
end

%{
    subplot(2,3,1)
    plot(U10_vec,RT,'-o','linewidth',2,'displayname',SSGF_str_title)
    hold on
    subplot(2,3,4)
    plot(U10_vec,RbT,'-o','linewidth',2,'displayname',SSGF_str_title)
    hold on
    subplot(2,3,2)
    plot(U10_vec,Ru,'-o','linewidth',2,'displayname',SSGF_str_title)
    hold on
    subplot(2,3,5)
    plot(U10_vec,Rbu,'-o','linewidth',2,'displayname',SSGF_str_title)
    hold on
    subplot(2,3,3)
    plot(U10_vec,Rr,'-o','linewidth',2,'displayname',SSGF_str_title)
    hold on
    subplot(2,3,6)
    plot(U10_vec,Rbr,'-o','linewidth',2,'displayname',SSGF_str_title)
    hold on




%}


%%
%{
for i = [1 4]
    subplot(2,3,i)
    if i == 1
        t_str = 'tau_f';
    else
        t_str = 'tau_b';
    end
    title(sprintf('$$\\frac{\\int \\frac{T(\\%s)-T_s}{T_w-T_s} \\frac{dF}{dr_0} r_0^3 dr_0 }{ \\int \\frac{dF}{dr_0} r_0^3 dr_0}$$',t_str),'interpreter','latex','rotation',0)
end

clc

for i = [2 5]
    subplot(2,3,i)
    if i == 3
        t_str = 'tau_f';
    else
        t_str = 'tau_b';
    end
    title(sprintf('$$\\frac{\\int \\frac{u(\\%s)-u_0}{U_{10}-u_0} \\frac{dF}{dr_0} r_0^3 dr_0 }{ \\int \\frac{dF}{dr_0} r_0^3 dr_0}$$',t_str),'interpreter','latex','rotation',0)
    
end

for i = [3 6]
    subplot(2,3,i)
    if i == 3
        t_str = 'tau_f';
    else
        t_str = 'tau_b';
    end
    title(sprintf('$$\\frac{\\int \\left(\\frac{r(\\%s)}{r_0}\\right)^3 \\frac{dF}{dr_0} r_0^3 dr_0 }{ \\int \\frac{dF}{dr_0} r_0^3 dr_0}$$',t_str),'interpreter','latex','rotation',0)
    
end
set(gcf,'position',[1  35   1404   770],'color','w')

drawnow
% add a bit space to the figure
fig = gcf;
for i = 1:6
    set(gca,'fontsize',20)
    xlabel('$U_{10}$ [m/s]','interpreter','latex')
    h(i) = subplot(2,3,i);
    drawnow
end
% add legend
subplot(2,3,5)
lh = legend('show','fontsize',15,'interpreter','latex','location','southeast');

% addpath('/Users/ssroka/Documents/MATLAB/util/')
% rearrange_figure(h,lh,'2x3_1_legend')

% update_figure_paper_size()
% print(sprintf('imgs/mass_flux_Tu'),'-dpdf')
%}
%%
%{
figure(2)
subplot(2,3,1)
imagesc(U10_vec,r0_vec,pct_mass_flux)
xlabel('U10 [m/s]')
ylabel('r0')
title('% mass flux')
colorbar
subplot(2,3,2)
imagesc(U10_vec,r0_vec,frac_tauf)
xlabel('U10 [m/s]')
ylabel('r0')
title('(Tf-T0)/(Tw-T0)')
colorbar
subplot(2,3,3)
imagesc(U10_vec,r0_vec,(frac_tauf>=0.9).*pct_mass_flux)
colorbar
xlabel('U10 [m/s]')
ylabel('r0')
title('(Tf-T0)/(Tw-T0)>0.9  X  % mass flux')
subplot(2,3,4)
plot(U10_vec,sum(pct_mass_flux),'-o','linewidth',2)
colorbar
xlabel('U10 [m/s]')
title('\Sigma % mass flux')
subplot(2,3,5)
axis off
title(SSGF_str_title)
subplot(2,3,6)
imagesc(U10_vec,r0_vec,frac_tauf)
colorbar
plot(U10_vec,sum((frac_tauf>=0.9).*pct_mass_flux),'-o','linewidth',2)
xlabel('U10 [m/s]')
title('\Sigma % mass flux')

for i  = [1:6]
    subplot(2,3,i)
    colormap(othercolor('RdYlGn10'))
    set(gca,'fontsize',25)
end
set(gcf,'position',[70         107        1259         698])
savefig([SSGF_str '_mass_flux'])
%%
action_str = 'remove';
paths_for_calc_CK;

%}

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
