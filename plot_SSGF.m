clear;close all;clc


SSGF_str_cell = {'troit','OS','Zhao'};
nSSGF = length(SSGF_str_cell);


colors =[...
    1.0000   0   0
    0      0.5   0
    0        0   0];
% parameters corresponding to SGF from Troitskaya et al. 2018a
eq_N = 9; % eq_N \in [5,7,9] - which equation to use for the number of bags
% 5 = derived equation
% 7 = laboratory conditions
% 9 = field conditions
Omega = 2.5; % wave-age parameter (usually between 2.5- 3.5)

beta = 1.2; % Zhao only

U10 = 54;

action_str = 'add';
paths_for_calc_CK;

LS = {'-','--','-.'};

for SSGF_ind = 1:nSSGF
    SSGF_str = SSGF_str_cell{SSGF_ind};
    load(sprintf('ratios_%s',SSGF_str),'r0_vec');
    switch SSGF_str
        case 'troit'
            dFdr_vec = calc_Troit_SGF(U10,eq_N,Omega,r0_vec*1e-6); % meters
            SSGF_str_title = ['Troitskaya et al. (2018)a' newline '$$r_0 = [50-2000]\mu$$m'];
        case 'OS'
            dFdr_vec = calc_OS_SGF(U10,r0_vec); % meters
            SSGF_str_title = ['Ortiz-Suslow et al. (2016)' newline '$$r_0 = [100-1000]\mu$$m'];
        case 'Zhao'
            addpath('/Users/ssroka/MIT/Research/EmanuelGroup/thesis/review/');
            dFdr_vec = Zhao2006(r0_vec,U10,beta); % meters
            SSGF_str_title = ['Zhao et al. (2006)' newline '$$r_0 = [50-500]\mu$$m'];
    end
    nr0 = length(r0_vec);
    
    % normalized mass
    total_mass = trapz(r0_vec,dFdr_vec.*4./3.*pi.*r0_vec.^3);
    dr0 = 50;
    dFdr_vec_norm = (dFdr_vec*dr0.*[0.5 ones(1,nr0-2) 0.5]).*(4./3.*pi.*r0_vec.^3)./total_mass;
    
    
    subplot(1,2,1)
    h(SSGF_ind,1) = loglog(r0_vec,dFdr_vec,LS{SSGF_ind},...
        'color','k',...
        'displayname',SSGF_str_title,...
        'linewidth',3);
    hold on
    
    subplot(1,2,2)
    h(SSGF_ind,2) = loglog(r0_vec,dFdr_vec_norm,LS{SSGF_ind},...
        'color','k',...
        'displayname',sprintf('%s',SSGF_str_title),...
        'linewidth',3);
    hold on
end

for i = 1:2
    subplot(1,2,i)
    set(gca,'fontsize',20)
    xlabel('$r_{0}$ [$\mu$m]','interpreter','latex')
    set(gca,'fontsize',20,'xtick',[100 500 1000 2000])
    drawnow
    % add legend
    if i == 1
        lh = legend([h(1,i) h(2,i) h(3,i)]);
        set(lh,'fontsize',15,'interpreter','latex','location','southwest');
        title('$\frac{dF}{dr_0}$ [m$^{-2}$s$^{-1}\mu$m$^{-1}$]','interpreter','latex')
    else
        title('$\frac{F\frac{4}{3}\pi r_0^3 }{ \int \frac{dF}{dr_0} \frac{4}{3}\pi r_0^3 dr_0}  $ [m$^{-2}$s$^{-1}$]','interpreter','latex')
    end
    set(gcf,'position',[35         400        1365         394],'color','w')
    
end
% addpath('/Users/ssroka/Documents/MATLAB/util/')
% rearrange_figure(h,lh,'2x3_1_legend')

update_figure_paper_size()
print(sprintf('imgs/SSGF'),'-dpdf')

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

