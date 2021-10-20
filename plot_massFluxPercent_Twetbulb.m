clear;close all;clc

Troit_waveage_num = 3;
Zhao_waveage_num = 3;

SSGF_str_cell = {'troit','OS','Zhao'};
nSSGF = length(SSGF_str_cell);

% R_titles = ...
%     {'$\frac{\int \frac{T(\tau_f)-T_s}{T_w-T_s} \frac{dF}{dr_0} r_0^3 dr_0 }{ \int \frac{dF}{dr_0} r_0^3 dr_0}$',...
%      '$\frac{\int \frac{u(\tau_f)-u_0}{U_{10}-u_0} \frac{dF}{dr_0} r_0^3 dr_0 }{ \int \frac{dF}{dr_0} r_0^3 dr_0}$',...
%      '$\frac{\int \left(\frac{r(\tau_f)}{r_0}\right)^3 \frac{dF}{dr_0} r_0^3 dr_0 }{ \int \frac{dF}{dr_0} r_0^3 dr_0}$',...
%      '$\frac{\int \frac{T(\tau_b)-T_s}{T_w-T_s} \frac{dF}{dr_0} r_0^3 dr_0 }{ \int \frac{dF}{dr_0} r_0^3 dr_0}$',...
%      '$\frac{\int \frac{u(\tau_b)-u_0}{U_{10}-u_0} \frac{dF}{dr_0} r_0^3 dr_0 }{ \int \frac{dF}{dr_0} r_0^3 dr_0}$',...
%      '$\frac{\int \left(\frac{r(\tau_b)}{r_0}\right)^3 \frac{dF}{dr_0} r_0^3 dr_0 }{ \int \frac{dF}{dr_0} r_0^3 dr_0}$',...
%     };

R_titles = ...
    {'$\gamma_1 = \frac{\int \frac{T_s-T(\tau_f)}{T_s-T_w} \frac{dF}{dr_0} r_0^3 dr_0 }{ \int \frac{dF}{dr_0} r_0^3 dr_0}$',...
    '$\gamma_2 = \frac{\int \left(\frac{r(\tau_f)}{r_0}\right)^3 \frac{dF}{dr_0} r_0^3 dr_0 }{ \int \frac{dF}{dr_0} r_0^3 dr_0}$',...
    '$\gamma_3 = \frac{\int \frac{u(\tau_f)}{U_{10}} \frac{dF}{dr_0} r_0^3 dr_0 }{ \int \frac{dF}{dr_0} r_0^3 dr_0}$',...
    '$\gamma_1 = \frac{\int \frac{T_s-T(\tau_b)}{T_s-T_w} \frac{dF}{dr_0} r_0^3 dr_0 }{ \int \frac{dF}{dr_0} r_0^3 dr_0}$',...
    '$\gamma_2 = \frac{\int \left(\frac{r(\tau_b)}{r_0}\right)^3 \frac{dF}{dr_0} r_0^3 dr_0 }{ \int \frac{dF}{dr_0} r_0^3 dr_0}$',...
    '$\gamma_3 = \frac{\int \frac{u(\tau_b)}{U_{10}} \frac{dF}{dr_0} r_0^3 dr_0 }{ \int \frac{dF}{dr_0} r_0^3 dr_0}$',...
    '$\gamma_1 = \frac{\int \frac{T_s-T(\tau_{Hsu})}{T_s-T_w} \frac{dF}{dr_0} r_0^3 dr_0 }{ \int \frac{dF}{dr_0} r_0^3 dr_0}$',...
    '$\gamma_2 = \frac{\int \left(\frac{r(\tau_{Hsu})}{r_0}\right)^3 \frac{dF}{dr_0} r_0^3 dr_0 }{ \int \frac{dF}{dr_0} r_0^3 dr_0}$',...
    '$\gamma_3 = \frac{\int \frac{u(\tau_{Hsu})}{U_{10}} \frac{dF}{dr_0} r_0^3 dr_0 }{ \int \frac{dF}{dr_0} r_0^3 dr_0}$',...
    };

MATLAB_colors = [...
    0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250];

colorsF4 =[...
    184 51 106
    247 176 91
    133 203 51
    16 126 125
    49 133 252
    77 255 243
    ]./256;

colors =[...
    1.0000   0   0
    0      0.5   0
    0        0   0];

frac_var = 'Tur';

for SSGF_ind = 1:nSSGF
    SSGF_str = SSGF_str_cell{SSGF_ind};
    
    if strcmp(SSGF_str,'Zhao')
        SSGF_str_title = ['Zhao et al. (2006)' newline '$$r_0 = [50-500]\mu$$m'];
    elseif strcmp(SSGF_str,'troit')
        SSGF_str_title = ['Troitskaya et al. (2018a)' newline '$$r_0 = [50-2000]\mu$$m'];
    elseif strcmp(SSGF_str,'OS')
        SSGF_str_title = ['Ortiz-Suslow et al. (2016)' newline '$$r_0 = [100-1000]\mu$$m'];
    end
    
    load(sprintf('ratios_%s_%d_%d',SSGF_str,Zhao_waveage_num,Troit_waveage_num));
    
    nr0 = length(r0_vec);
    nU10 = length(U10_vec);
    nDT = length(DT_vec);
    nRH = length(RH_vec);
    nSST = length(SST_vec);
    ntof = length(tof_vec);
    
    % [ratio](U10_ind,RH_ind,DT_ind,SST_ind,tof_ind)
    mu = zeros(1,nU10);
    sig = zeros(1,nU10);
    mu_Qk = zeros(nr0,nU10);
    sig_Qk = zeros(nr0,nU10);
    
    R_cell = {RT,Rr,Ru,RbT,Rbr,Rbu,RHsuT,RHsur,RHsuu};
    count_plot = 1; % because R_ind is in the 'wrong' order for how I want to plot these
    for R_ind = [1:9]
        
        R = R_cell{R_ind};
        
        for U10_ind = 1:nU10
            % get ratio entries for this U10
            R_vec = reshape(R(U10_ind,:,:,:,:),nRH*nDT*nSST*ntof,1);
            mu(U10_ind) = mean(R_vec);
            sig(U10_ind) = std(R_vec);
        end
        % select figure(1), figure(2), or figure(3)
        if R_ind<=3
            figure(1);
        elseif R_ind>=7
            figure(3);
            if count_plot>3; count_plot = 1; end % Hack the counter
        else
            figure(2);
            if count_plot>3; count_plot = 1; end % Hack the counter
        end
        subplot(3,1,count_plot)
        if ismember(R_ind,[3 6 9])
            line_h(SSGF_ind,R_ind) = plot(U10_vec,mu,'-o',...
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
        ylh(R_ind) = ylabel(R_titles{R_ind},'interpreter','latex','fontsize',30,'rotation',0);
        fprintf('\n%s %s\navg sig = %f\n',frac_var(count_plot),SSGF_str_title,mean(sig))
        count_plot = count_plot + 1;
        xlabel('$U_{10}$ [m/s]','interpreter','latex')
        set(gcf,'position',[36  102 1257  688],'color','w')
        set(gca,'fontsize',30)
    end
    
    if  strcmp(SSGF_str,'troit')
        figure(4)
        count = 1;
        for r0_ind = [1 2 5 10 20 40]
            for U10_ind = 1:nU10
                % get ratio entries for this U10
                Qk_vec = reshape(Qk(r0_ind,U10_ind,:,:,:,:),nRH*nDT*nSST*ntof,1);
                mu_Qk(r0_ind,U10_ind) = mean(Qk_vec);
                sig_Qk(r0_ind,U10_ind) = std(Qk_vec);
            end
            semilogy(U10_vec,mu_Qk(r0_ind,:),'-o',...
                'displayname',sprintf('$r_0 = %d \mu $m',r0_vec(r0_ind)),...
                'linewidth',3,'color',colorsF4(count,:));
            hold on
            ebh = errorbar(U10_vec,mu_Qk(r0_ind,:),sig_Qk(r0_ind,:));
            set(ebh,'linewidth',2,'color',colorsF4(count,:))
            th = text(U10_vec(end-1)+2,mu_Qk(r0_ind,end-1),sprintf('$r_0$ = %d $$\\mu$$m',r0_vec(r0_ind)));
            set(th,'interpreter','latex','BackgroundColor', [1 1 1],'fontsize',15)
            count = count + 1;
        end
        title('$Q$ [J]','interpreter','latex')
    end
    
    xlabel('$U_{10}$ [m/s]','interpreter','latex')
    set(gcf,'position',[36  102 1257  688],'color','w')
    set(gca,'fontsize',30)
end

plot_letter = 'abc';
for j = 1:3
    figure(j)
    for i = 1:3
        subplot(3,1,i)
        set(gca,'fontsize',18,'ytick',[0.1 0.5 0.8 1],'yticklabel',{'10%','50%','80%','100%'})
        plot(U10_vec,0.5*ones(nU10,1),'--','linewidth',2,'color',[1 1 1]*0.5)
        plot(U10_vec,0.8*ones(nU10,1),':','linewidth',2,'color',[1 1 1]*0.5)
        set(ylh(i+3*(j-1)),'fontsize',30)
        if i == 3
            xlabel('$U_{10}$ [m/s]','interpreter','latex')
        else
            xlabel('')
            set(gca,'xtick',U10_vec,'XTickLabel',{})
        end
        %     txth = text(0.1,0.1,sprintf('(%s)',plot_letter(i)));
        %     set(txth,'units','normalized','position',[0.95 0.05],'fontsize',20,'interpreter','latex')
        drawnow
        set(gcf,'position',[58 1 1001 802],'color','w')
    end
    % add legend
    % figure(2)
    text_array_subplot = {'a)','b)','c)'};
    t_sp_loc = [0.01,0.2;0.01,0.1;0.05,0.1];
    for i = 1:3
        h(i) = subplot(3,1,i);
        t_sp(i) = text(t_sp_loc(i,1),t_sp_loc(i,2),text_array_subplot{i},'fontsize',30,'units','normalized','BackgroundColor', [1 1 1]);
    end
    lh = legend([line_h(3,3+3*(j-1)) line_h(2,3+3*(j-1)) line_h(1,3+3*(j-1))],'location','southeast','interpreter','latex');
    % set(lh,'fontsize',20,'interpreter','latex','location','east');
    
    % addpath('/Users/ssroka/Documents/MATLAB/util/')
    rearrange_figure(h,lh,'3x1_1_legend_3_yl',t_sp,ylh([1:3]+(3*(j-1))))
end
v = 'Tur';
figure(1)
% update_figure_paper_size()
% print(sprintf('imgs/mass_flux_%s_%d_%d',v,Zhao_waveage_num,Troit_waveage_num),'-dpdf')
% savefig(sprintf('imgs/mass_flux_%s_%d_%d',v,Zhao_waveage_num,Troit_waveage_num))

% figure(4)
% update_figure_paper_size()
% print(sprintf('imgs/Qk_%d_%d',Zhao_waveage_num,Troit_waveage_num),'-dpdf')

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