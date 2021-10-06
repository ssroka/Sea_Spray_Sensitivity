clear;close all;clc


SSGF_str_cell = {'troit'};
nSSGF = length(SSGF_str_cell);

TH_src = '/Users/ssroka/Documents/MATLAB/mpm/sandbox/results/res_fm_Eng/';

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
    184 51 106
    247 176 91
    133 203 51
    16 126 125
    49 133 252
    77 255 243
    ]./256;

r0_vec_sort = [31 2 13 24 25 26 27 28 29 30 32 33 34 35 36 37 38 39 40 1 3 4 5 6 7 8 9 10 11 12 14 15 16 17 18 19 20 21 22 23];
r0_vec_ALL = [50:50:2000];



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






% find tau_T
tau_T = zeros(nr0,1);
for i = 1:nr0
    f = ls(sprintf('%s/RH_92_DT_2_SST_28_long_run/r_0_%d_*.mat',TH_src,r0_vec_ALL(i)));
    load(f(1:end-1),'T_s_t','time_vec');% rm newline, also a tau_T is in here, don't load it
    [Teq,Teq_ind] = min(T_s_t);
    %{
    loglog([1e-2 1e4],[1 1]*T_s_t(Teq_ind))
    hold on
    loglog(time_vec(1:end-1),T_s_t(1:end-1))
    
    %}
    T_tau_T = (T_s_t(1)-Teq)*exp(-1)+Teq;
    tau_T(i) = interp1(T_s_t(1:Teq_ind),time_vec(1:Teq_ind),T_tau_T);
end


% find tau_r
tau_r = zeros(nr0,1);
for i = 1:nr0
    f = ls(sprintf('%s/RH_92_DT_2_SST_28_long_run/r_0_%d_*.mat',TH_src,r0_vec_ALL(i)));
    load(f(1:end-1),'r_t','time_vec');% rm newline, also a tau_T is in here, don't load it
    req = r_t(end);
    %{
    
    loglog([1e-2 1e4],[1 1]*r_t(req_ind))
    hold on
    loglog(time_vec(1:end-1),r_t(1:end-1))
    yyaxis right
    loglog(time_vec(1:end-1),dr)
    
    plot((log(r_t(2:end))-log(r_t(1:end-1)))./(log(time_vec(2:end))-log(time_vec(1:end-1))))

    
    %}
    r_tau_r = (r_t(1)-req)*exp(-1)+req;
    r_tau_r_ind = find(r_tau_r>r_t,1,'first');
    tau_r(i) = interp1(r_t(1:r_tau_r_ind),time_vec(1:r_tau_r_ind),r_tau_r);
    
    
end

figure(1)
count = 2;
grey = [1 1 1]*0.75;
for i = [1 2 5 10 20 40]
    for j = 1:6
        subplot(2,3,j)
        semilogy(U10_vec,tof_mat(i,:),'-','color',grey,'linewidth',3,...
            'displayname',sprintf('$\\tau_{f}, r_0$ = %d $$\\mu$$m',r0_vec(i)))
        hold on
        semilogy(U10_vec,tau_ac_mat(i,:),'--','color',grey,'linewidth',3,...
            'displayname',sprintf('$\\tau_{ac}, r_0$ = %d $$\\mu$$m',r0_vec(i)))
    end
    count = count + 1;
end
count = 1;

for i = [1 2 5 10 20 40]
    figure(2)
    semilogy(U10_vec,tof_mat(i,:),'-','color',colors(count,:),'linewidth',3,...
        'displayname',sprintf('$\\tau_{f}, r_0$ = %d $$\\mu$$m',r0_vec(i)))
    hold on
    semilogy(U10_vec,tau_ac_mat(i,:),'--','color',colors(count,:),'linewidth',3,...
        'displayname',sprintf('$\\tau_{ac}, r_0$ = %d $$\\mu$$m',r0_vec(i)))
    th = text(U10_vec(end-1)-2,tof_mat(i,end-1),sprintf('$r_0$ = %d $$\\mu$$m',r0_vec(i)));
    set(th,'interpreter','latex','BackgroundColor', [1 1 1],'fontsize',15)
    th = text(U10_vec(end-1)-2,tau_ac_mat(i,end-1),sprintf('$r_0$ = %d $$\\mu$$m',r0_vec(i)));
    set(th,'interpreter','latex','BackgroundColor', [1 1 1],'fontsize',15)
    
    figure(1)
    subplot(2,3,count)
    semilogy(U10_vec,tof_mat(i,:),'-','color',colors(count,:),'linewidth',3,...
        'displayname',sprintf('$\\tau_{f}, r_0$ = %d $$\\mu$$m',r0_vec(i)))
    hold on
    semilogy(U10_vec,tau_ac_mat(i,:),'--','color',colors(count,:),'linewidth',3,...
        'displayname',sprintf('$\\tau_{ac}, r_0$ = %d $$\\mu$$m',r0_vec(i)))
    semilogy(U10_vec,tau_T(i)*ones(nU10,1),'-o','markersize',10,'color',colors(count,:),'linewidth',2,...
        'displayname',sprintf('$\\tau_{T}, r_0$ = %d $$\\mu$$m',r0_vec(i)))
    semilogy(U10_vec,tau_r(i)*ones(nU10,1),'-^','markersize',10,'color',colors(count,:),'linewidth',2,...
        'displayname',sprintf('$\\tau_{r}, r_0$ = %d $$\\mu$$m',r0_vec(i)))
    
    figure(3)
    semilogy(U10_vec,abs(tof_mat(i,:)-tau_T(i)),'^-','color',colors(count,:),'linewidth',2)
    hold on
    th = text(U10_vec(end-1)-2,tof_mat(i,end-1)-tau_T(i),sprintf('$r_0$ = %d $$\\mu$$m',r0_vec(i)));
    set(th,'interpreter','latex','BackgroundColor', [1 1 1],'fontsize',15)
    
    semilogy(U10_vec,abs(tof_mat(i,:)-tau_r(i)),'o-','color',colors(count,:),'linewidth',2)
    th = text(U10_vec(end-1)-2,tof_mat(i,end-1)-tau_r(i),sprintf('$r_0$ = %d $$\\mu$$m',r0_vec(i)));
    set(th,'interpreter','latex','BackgroundColor', [1 1 1],'fontsize',15)
    
    count = count + 1;
end


figure(2)
h1 = semilogy(U10_vec,NaN*ones(nU10,1),'k-','displayname','$\tau_{f}$','linewidth',3);
h2 = semilogy(U10_vec,NaN*ones(nU10,1),'k--','displayname','$\tau_{ac}$','linewidth',3);
legend([h1 h2],'location','northwest','interpreter','latex')
xlabel('$U_{10}$ [m/s]','interpreter','latex')
ylabel('[s]','interpreter','latex','rotation',0)
set(gca,'fontsize',20,'ytick',10.^[-2:2:4],'xtick',20:20:80,'xlim',[10 80]);
set(gcf,'position',[148   243   794   540],'color','w')
update_figure_paper_size()
print(sprintf('imgs/res_time'),'-dpdf')
plot_letter = 'abcdef';
figure(1)
count = 1;
for i =  [1 2 5 10 20 40]
    if count < 4
        leg_loc = 'northwest';
    else
        leg_loc = 'southwest';
    end
    subplot(2,3,count)
    h1 = semilogy(U10_vec,NaN*ones(nU10,1),'^-','color',colors(count,:),'displayname','$\tau_r$','linewidth',2);
    h2 = semilogy(U10_vec,NaN*ones(nU10,1),'-','color',colors(count,:),'displayname','$\tau_{f}$','linewidth',2);
    h3 = semilogy(U10_vec,NaN*ones(nU10,1),'o-','color',colors(count,:),'displayname','$\tau_{T}$','linewidth',2);
    h4 = semilogy(U10_vec,NaN*ones(nU10,1),'--','color',colors(count,:),'displayname','$\tau_{ac}$','linewidth',2);
    lh = legend([h1 h2 h3 h4],'location',leg_loc,'interpreter','latex','numcolumns',2);
    xlabel('$U_{10}$ [m/s]','interpreter','latex')
    ylabel('[s]','interpreter','latex','rotation',0)
    title(sprintf('$r_0$ = %d $\\mu$m',r0_vec(i)),'interpreter','latex','rotation',0)
    set(gca,'fontsize',20,'ylim',[1e-3 2e5],'ytick',10.^[-2:2:4],'xtick',20:20:80,'xlim',[10 80]);
    th = text(0.9,0.05,sprintf('(%s)',plot_letter(count)),'units','normalized');
    set(th,'interpreter','latex','BackgroundColor', [1 1 1],'fontsize',15)
    
    count = count + 1;
end
set(gcf,'position',[36  102 1257  688],'color','w')
update_figure_paper_size()
print(sprintf('imgs/res_time_each_r0'),'-dpdf')


figure(3)
h1 = semilogy(U10_vec,NaN*ones(nU10,1),'k^-','displayname','$\tau_{f}-\tau_T$','linewidth',3);
h2 = semilogy(U10_vec,NaN*ones(nU10,1),'ko-','displayname','$\tau_{f}-\tau_r$','linewidth',3);
legend([h1 h2],'location','northwest','interpreter','latex')
xlabel('$U_{10}$ [m/s]','interpreter','latex')
ylabel('[s]','interpreter','latex','rotation',0)
set(gca,'fontsize',20,'ytick',10.^[-2:2:4],'xtick',20:20:80,'xlim',[10 80]);
set(gcf,'position',[148   243   794   540],'color','w')
update_figure_paper_size()
print(sprintf('imgs/res_time_minus_tauTr'),'-dpdf')


