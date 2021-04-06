clear;close all;clc


SSGF_str_cell = {'troit'};
nSSGF = length(SSGF_str_cell);


R_titles = ...
    {'$\frac{T(\tau_f)-T_s}{T_w-T_s}$',...
    '$\frac{u(\tau_f)-u_0}{U_{10}-u_0}$',...
    '$\left(\frac{r(\tau_f)}{r_0}\right)^3 $',...
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
    
    tof_ind = 1;
    for SST_ind = 2
        for DT_ind = 4
            for RH_ind = 5
                for U10_ind = 1:nU10
                    frac_tauf_T_vec = squeeze(frac_tauf_T(:,U10_ind,RH_ind,DT_ind,SST_ind,tof_ind));
                    %     frac_taub_T(r0_ind,U10_ind,RH_ind,DT_ind,SST_ind,tof_ind);
                    %
                    frac_tauf_u_vec = squeeze(frac_tauf_u(:,U10_ind,RH_ind,DT_ind,SST_ind,tof_ind));
                    %     frac_taub_u(r0_ind,U10_ind,RH_ind,DT_ind,SST_ind,tof_ind);
                    %
                    frac_tauf_r_vec = squeeze(frac_tauf_r(:,U10_ind,RH_ind,DT_ind,SST_ind,tof_ind));
                    %     frac_taub_r(r0_ind,U10_ind,RH_ind,DT_ind,SST_ind,tof_ind);
                    %     plot(r0_vec,[frac_tauf_T_vec frac_tauf_u_vec frac_tauf_r_vec])
                    subplot(1,3,1)
                    semilogx(r0_vec,[frac_tauf_T_vec],'color',colors(U10_ind,:),'linewidth',3,'displayname',sprintf('$U_{10} =$ %d m/s',U10_vec(U10_ind)))
                    hold on
                    subplot(1,3,2)
                    semilogx(r0_vec,[frac_tauf_u_vec],'color',colors(U10_ind,:),'linewidth',3,'displayname',sprintf('$U_{10} =$ %d m/s',U10_vec(U10_ind)))
                    hold on
                    subplot(1,3,3)
                    semilogx(r0_vec,[frac_tauf_r_vec],'color',colors(U10_ind,:),'linewidth',3,'displayname',sprintf('$U_{10} =$ %d m/s',U10_vec(U10_ind)))
                    hold on
                end
            end
        end
    end
end
%%
subplot(1,3,1)
title(R_titles{1},'interpreter','latex','fontsize',20)
subplot(1,3,2)
title(R_titles{2},'interpreter','latex','fontsize',20)
subplot(1,3,3)
title(R_titles{3},'interpreter','latex','fontsize',20)

for i = 1:3
    subplot(1,3,i)
    xlabel('$r_{0}$ [$\mu$m]','interpreter','latex')
    set(gca,'fontsize',20,'xtick',[100 500 1000 2000])
    %     drawnow
end
% look like contour labels
% th1 = text(180,0.5,sprintf('$U_{10}$ = %d',U10_vec(1)));
% th2 = text(1200,0.5,sprintf('$U_{10}$ = %d',U10_vec(nU10)));
% set(th1,'fontsize',15,'interpreter','latex')
% set(th2,'fontsize',15,'interpreter','latex')
subplot(1,3,3)
lh = legend('-dynamiclegend');
set(lh,'NumColumns',1,'fontsize',20,'interpreter','latex','location','southeast')

set(gcf,'position',[440    36   884   762],'color','w')

% addpath('/Users/ssroka/Documents/MATLAB/util/')
% rearrange_figure(h,lh,'2x3_1_legend')

% update_figure_paper_size()
% print(sprintf('imgs/cmp_drops'),'-dpdf')

% compare residence times
%%
%{
figure(2)
for i = [1:5 10 25 40]
        plot(U10_vec,tof_mat(i,:),'--','linewidth',3,'displayname',sprintf('Andreas 1992, r0 = %d $$\\mu$$ m',r0_vec(i)))
hold on
end
plot(U10_vec,t_b_mat,'k-','linewidth',3,'displayname','Projectile, All Drop Sizes')

legend('location','best','interpreter','latex')
xlabel('$U_{10}$ [m/s]','interpreter','latex')
ylabel('residence time [s]','interpreter','latex')
set(gca,'fontsize',20);
set(gcf,'position',[1   233   982   572],'color','w')
update_figure_paper_size()
print(sprintf('imgs/res_time_3'),'-dpdf')
%}

