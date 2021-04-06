clear;close all;clc


SSGF_str_cell = {'troit'};
nSSGF = length(SSGF_str_cell);


R_titles = ...
    {'$\frac{\int \frac{T(\tau_f)-T_s}{T_w-T_s} \frac{dF}{dr_0} r_0^3 dr_0 }{ \int \frac{dF}{dr_0} r_0^3 dr_0}$',...
    '$\frac{\int \frac{u(\tau_f)-u_0}{U_{10}-u_0} \frac{dF}{dr_0} r_0^3 dr_0 }{ \int \frac{dF}{dr_0} r_0^3 dr_0}$',...
    '$\frac{\int \left(\frac{r(\tau_f)}{r_0}\right)^3 \frac{dF}{dr_0} r_0^3 dr_0 }{ \int \frac{dF}{dr_0} r_0^3 dr_0}$',...
    '$\frac{\int \frac{T(\tau_b)-T_s}{T_w-T_s} \frac{dF}{dr_0} r_0^3 dr_0 }{ \int \frac{dF}{dr_0} r_0^3 dr_0}$',...
    '$\frac{\int \frac{u(\tau_b)-u_0}{U_{10}-u_0} \frac{dF}{dr_0} r_0^3 dr_0 }{ \int \frac{dF}{dr_0} r_0^3 dr_0}$',...
    '$\frac{\int \left(\frac{r(\tau_b)}{r_0}\right)^3 \frac{dF}{dr_0} r_0^3 dr_0 }{ \int \frac{dF}{dr_0} r_0^3 dr_0}$',...
    };

MATLAB_colors = [...
    0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250];

% colors = distinguishable_colors(8);

colors =[...
    1.0000         0         0
    1.0000    0.1034    0.7241
    1.0000    0.8276         0
    0    1.0000         0
    0    0.3448         0
    0         0    0.1724
    0         0    1.0000
    0.5172    0.5172    1.0000];

for SSGF_ind = 1:nSSGF
    SSGF_str = SSGF_str_cell{SSGF_ind};
    
    if strcmp(SSGF_str,'Zhao')
        SSGF_str_title = ['Zhao et al. (2006)' newline '$$r_0 = [50-500]\mu$$m'];
    elseif strcmp(SSGF_str,'troit')
        SSGF_str_title = ['Troitskaya et al. (2018)a' newline '$$r_0 = [50-2000]\mu$$m'];
    elseif strcmp(SSGF_str,'OS')
        SSGF_str_title = ['Ortiz-Suslow et al. (2016)' newline '$$r_0 = [100-1000]\mu$$m'];
    end
    
    load(sprintf('ratios_%s',SSGF_str));
    
    nr0 = length(r0_vec);
    nU10 = length(U10_vec);
    nDT = length(DT_vec);
    nRH = length(RH_vec);
    nSST = length(SST_vec);
    ntof = length(tof_vec);
end

figure(1)
count = 1;
for i = [1 2 10 20 40]
    semilogy(U10_vec,tof_mat(i,:),'-','color',colors(count,:),'linewidth',3,...
        'displayname',sprintf('$r_0$ = %d $$\\mu$$m',r0_vec(i)))
    hold on
    count = count + 1;
end
% plot(U10_vec,t_b_mat,'k-','linewidth',3,'displayname','Projectile, All Drop Sizes')

legend('location','southeast','interpreter','latex')
xlabel('$U_{10}$ [m/s]','interpreter','latex')
title('$\tau_f$ [s]','interpreter','latex')
set(gca,'fontsize',20);
set(gcf,'position',[148   401   638   382],'color','w')
update_figure_paper_size()
print(sprintf('imgs/res_time'),'-dpdf')


