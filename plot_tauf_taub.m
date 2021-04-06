clear;close all;clc


SSGF_str_cell = {'troit','OS','Zhao'};
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

colors =[...
    1.0000   0   0
    0      0.5   0
    0        0   0];

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
    
    % [ratio](U10_ind,RH_ind,DT_ind,SST_ind,tof_ind)
    mu = zeros(1,nU10);
    sig = zeros(1,nU10);
    
    R_cell = {RT,Ru,Rr,RbT,Rbu,Rbr};
    
    for R_ind = 1:length(R_cell)
        R = R_cell{R_ind};
        
        for U10_ind = 1:nU10
            % get ratio entries for this U10
            R_vec = reshape(R(U10_ind,:,:,:,:),nRH*nDT*nSST*ntof,1);
            mu(U10_ind) = mean(R_vec);
            sig(U10_ind) = std(R_vec);
        end
        subplot(2,3,R_ind)
        if R_ind == 5
        h(SSGF_ind) = plot(U10_vec,mu,'-o',...
            'color',colors(SSGF_ind,:),...
            'displayname',SSGF_str_title,...
            'linewidth',3);
        else
            plot(U10_vec,mu,'-o',...
            'color',colors(SSGF_ind,:),...
            'displayname',SSGF_str_title,...
            'linewidth',3);
        end
        hold on
        ebh = errorbar(U10_vec,mu,sig);
        set(ebh,'linewidth',2,'color',colors(SSGF_ind,:))
        title(R_titles{R_ind},'interpreter','latex')
    end
    
end


for i = 1:6
    subplot(2,3,i)
    set(gca,'fontsize',20)
    xlabel('$U_{10}$ [m/s]','interpreter','latex')
    drawnow
end
% add legend
subplot(2,3,4)
lh = legend([h(1) h(2) h(3)]);
set(lh,'fontsize',15,'interpreter','latex','location','southeast');
set(gcf,'position',[1  35   1404   770],'color','w')

% addpath('/Users/ssroka/Documents/MATLAB/util/')
% rearrange_figure(h,lh,'2x3_1_legend')

update_figure_paper_size()
print(sprintf('imgs/mass_flux_Tu'),'-dpdf')

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

