clear;close all;clc

% SSGF_str = 'troit';% 'troit' 'OS' 'Zhao'
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

beta = 1./Omega; % Zhao shows 0.2 and 1.2

U10_fig1 = 54;

action_str = 'add';
paths_for_calc_CK;

LS = {'-','--','-.'};

a = zeros(nSSGF,1);
b = zeros(nSSGF,1);
R2 = zeros(nSSGF,1);

for SSGF_ind = 1:nSSGF
    SSGF_str = SSGF_str_cell{SSGF_ind};
    load(sprintf('ratios_%s',SSGF_str),'r0_vec');
    
    
    
    % % parameters common to all environmental conditions
    
    if strcmp(SSGF_str,'Zhao')
        addpath('/Users/ssroka/MIT/Research/EmanuelGroup/thesis/review/');
        beta = 0.4; % Zhao only
        r0_vec = [50:50:500];%[30 500]  % initial drop radii [micron]
        U10_vec = [10:10:80];    % [m/s]
        SSGF_str_title = ['Zhao et al. (2006)' newline '$$r_0 = [50-500]\mu$$m'];
        U10_inds = true(1,length(U10_vec));
        
    elseif strcmp(SSGF_str,'troit')
        r0_vec = [50:50:2000];%[50:50:2000];  % initial drop radii [micron]
        U10_vec = [10:10:80];    % [m/s]
        SSGF_str_title = ['Troitskaya et al. (2018a)' newline '$$r_0 = [50-2000]\mu$$m'];
        U10_inds = U10_vec>20;
        
    elseif strcmp(SSGF_str,'OS')
        r0_vec = [100:50:1000];%[86 1036];  % initial drop radii [micron]
        U10_vec = [36 40.5 45 49.5 54];%[10:10:80];    % [m/s]
        SSGF_str_title = ['Ortiz-Suslow et al. (2016)' newline '$$r_0 = [100-1000]\mu$$m'];
        U10_inds = true(1,length(U10_vec));
        
    end
    nU10 = length(U10_vec);
    total_mass = zeros(nU10+1,1);
    for U10_ind = 1:(nU10+1)
        
        if U10_ind>nU10
            U10 = U10_fig1;
        else
            U10 = U10_vec(U10_ind);
        end
        
        switch SSGF_str
            case 'troit'
                dFdr_vec = calc_Troit_SGF(U10,eq_N,Omega,r0_vec*1e-6); % meters
                %             SSGF_str_title = ['Troitskaya et al. (2018)a' newline '$$r_0 = [50-2000]\mu$$m'];
            case 'OS'
                dFdr_vec = calc_OS_SGF(U10,r0_vec); % meters
                %             SSGF_str_title = ['Ortiz-Suslow et al. (2016)' newline '$$r_0 = [100-1000]\mu$$m'];
            case 'Zhao'
                addpath('/Users/ssroka/MIT/Research/EmanuelGroup/thesis/review/');
                dFdr_vec = Zhao2006(r0_vec,U10,beta); % meters
                %             SSGF_str_title = ['Zhao et al. (2006)' newline '$$r_0 = [50-500]\mu$$m'];
        end
        nr0 = length(r0_vec);
        
        % normalized mass
        total_mass(U10_ind) = trapz(r0_vec,dFdr_vec.*4./3.*pi.*r0_vec.^3);
        dr0 = 50;
        dFdr_vec_norm = (dFdr_vec*dr0.*[0.5 ones(1,nr0-2) 0.5]).*(4./3.*pi.*r0_vec.^3)./total_mass(U10_ind);
        
        if U10_ind>nU10
            figure(1)
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
    end
    X = [ones(sum(U10_inds),1) log(U10_vec(U10_inds)')];
    Y = log(total_mass(U10_inds));
    [B,~,~,~,STATS] = regress(Y,X);
    a(SSGF_ind) = exp(B(1));
    b(SSGF_ind) = B(2);
    R2(SSGF_ind) = STATS(1);
    
    figure(2)
    h2(SSGF_ind) = loglog(U10_vec,total_mass(1:nU10)*(1e-6).^3,...
        'color','k','linestyle',LS{SSGF_ind},...
        'displayname',SSGF_str_title,...
        'linewidth',3);
    hold on
end

figure(2)
% semilogy(U10_vec,U10_vec.^12.*1e-20,'b','linewidth',3)
h2(SSGF_ind+1)=semilogy(U10_vec,U10_vec.^1.5.*1e-5,'b','linewidth',3,'displayname','$U_{10}^{1.5}$');
h2(SSGF_ind+2)=semilogy(U10_vec,U10_vec.^9.*1e-20,'r','linewidth',3,'displayname','$U_{10}^9$');
h2(SSGF_ind+3)=semilogy(U10_vec,U10_vec.^5.*1e-11,'g','linewidth',3,'displayname','$U_{10}^5$');
lh = legend([h2(nSSGF:-1:1) h2(nSSGF+3:-1:nSSGF+1)]);
set(lh,'fontsize',15,'interpreter','latex','location','southeast','numcolumn',2);
title('$\int \frac{4}{3}\pi r_0^3 \frac{dF}{dr_0} dr_0$ [m$^{1}$s$^{-1}$]','interpreter','latex')
set(gcf,'color','w')
xlabel('$U_{10}$ [m/s]','interpreter','latex')
set(gca,'fontsize',20,'xtick',10:10:80)
update_figure_paper_size()
print(sprintf('imgs/SSGF_totMass'),'-dpdf')
% TotalMass = a (U10^b)
% log(TotalMass) = log(a)+blog(U10)
a
b
R2

plot_letter = 'ab';
txt_xy = [0.025 0.05];
figure(1)
for i = 1:2
    subplot(1,2,i)
    set(gca,'fontsize',20)
    xlabel('$r_{0}$ [$\mu$m]','interpreter','latex')
    set(gca,'fontsize',20,'xtick',[100 500 1000 2000])
    drawnow
    % add legend
    if i == 1
        th = text(txt_xy(1),txt_xy(2),sprintf('(%s)',plot_letter(i)),'units','normalized');
        set(th,'interpreter','latex','BackgroundColor', [1 1 1],'fontsize',20)
        
        lh = legend([h(3,i) h(2,i) h(1,i)]);
        set(lh,'fontsize',15,'interpreter','latex','location','northeast');
        title('$\frac{dF}{dr_0}$ [m$^{-2}$s$^{-1}\mu$m$^{-1}$]','interpreter','latex')
    else
        th = text(txt_xy(1),txt_xy(2),sprintf('(%s)',plot_letter(i)),'units','normalized');
        set(th,'interpreter','latex','BackgroundColor', [1 1 1],'fontsize',20)
        title('$\frac{F\frac{4}{3}\pi r_0^3 }{ \int \frac{dF}{dr_0} \frac{4}{3}\pi r_0^3 dr_0}  $','interpreter','latex')
    end
    set(gcf,'position',[35         344        1268         450],'color','w')
    
end
% addpath('/Users/ssroka/Documents/MATLAB/util/')
% rearrange_figure(h,lh,'2x3_1_legend')
subplot(1,2,1)
set(gca,'ylim',[1e-2 1e7],'ytick',10.^[-2:2:6])
update_figure_paper_size()
print(sprintf('imgs/SSGF'),'-dpdf')
figure(2)

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

