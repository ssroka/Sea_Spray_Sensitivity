% script to evaluate the SGF from Troitskaya 2018 JPO
clear
close all
clc

rho_w = 1020;       % saltwater density [kg m^-3]
N_rim = 8.3;        %
H10   = 10;         % height [m]
g     = 9.81;       % gravity[m s^-2]
sig   = 0.07;       % N/m surface tension
alpha = 0.0057;     % mixing length??
kappa = 0.4;        % von karman constant
L     = 20;         % [m] Lhuissier and Villermaux 2012
Q0    = 9.27e2;     % [m^-4 s]
U0    = 2;          % [m/s]
M0    = 2.58e-4;    % [m^-2 s^-1]
M1    = 6.93e5;     % [none]
nu    = 1.5e-5;     % air [m^2/s]
Omega_u0 = 12.4;    % [m^1/2 s^-3/2]
gam   = 0.5;

z0        =@(u_star) alpha.*u_star.^2./g;
U10       =@(u_star) (u_star+0.14)/0.051;
R2_u_star =@(u_star) 9.6/u_star*1e-3;
C_D       =@(u_star) (U10(u_star)<=40)*(0.057-0.48/U10(u_star)).^2 ...
                     + (U10(u_star)>40)*(2.37/U10(u_star)-0.012).^2;
% lab
omega_p   =@(u_star) Omega_u0.*u_star.^(-gam);
% field
ReB       =@(u_star,Omega) (0.051*U10(u_star)-0.14)^2*U10(u_star)./g./nu./Omega;

%
% eq 5 [m^-2 s^-1]
N5 =@(u_star) Q0*u_star^2*exp(-U0^2/u_star^2);
% % eq 7 - for real experiemnes with fetch
N7 =@(u_star,Omega)M0.*(ReB(u_star,Omega)).^(3./2).*exp(-M1./((ReB(u_star,Omega)).^(3./2)));
% % eq 9 - for field conditions

%%%%%%%%%%%%%%%%%%%%
% TYPO IN PAPER!!!!!!!!!!!!! correct eq 9 here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N9 =@(u_star,Omega)Q0*(u_star.^2.).*(12.4*(u_star.^-gam*U10(u_star))./(g*Omega)).^(1.5).*...
    exp(-U0.^2./u_star.^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% appendix B 
th =@(u_star) 0.001*(R2_u_star(u_star))^(4/3)*L^(-1/3);

% appendix C 
TH =@(u_star) 0.0021*(R2_u_star(u_star));

% eq 24
dfdr =@(N,r,th,TH,u_star)...
    N.*(3.3e-9./L.*(rho_w.*g.*L.^2./sig).^(1.18).*(r./th).^7.3.*...
    exp(-5.2.*sinh(3./7.*log(r./th)))+...
    1.5e-4.*N_rim./TH.*(r./TH).^4.5.*...
    exp(-3.94.*sinh(0.5.*log(r./TH))));
    
T1 =@(r,th,TH) 3.3e-9./L.*(rho_w.*g.*L.^2./sig).^(1.18).*(r./th).^7.3.*...
     exp(-5.2.*sinh(3./7.*log(r./th)));
T2 =@(r,th,TH) 1.5e-4.*N_rim./TH.*(r./TH).^4.5.*...
    exp(-3.94.*sinh(0.5.*log(r./TH)));


% fig 12 eq 24
figure(12)

for u_star=1:.1:2
r = logspace(log10(6),log10(2000))*1e-6;
term1 = T1(r,th(u_star),TH(u_star));
term2 = T2(r,th(u_star),TH(u_star));
SGF = dfdr(N5(u_star),r,th(u_star),TH(u_star),u_star);

subplot(1,3,1)
loglog(r*1e6,SGF*1e-6,'k')
hold on
xlabel('r_0')
ylabel('dF/dr')
set(gca,'fontsize',20)
if u_star == 1.5
subplot(1,3,2)
loglog(r*1e6,N5(u_star)*term1*1e-6,'-.')
hold on
loglog(r*1e6,N5(u_star)*term2*1e-6,'-.')
xlabel('r_0')
ylabel('dF/dr')
set(gca,'fontsize',20)
end
subplot(1,3,3)
loglog(r*1e6,SGF*1e-6.*4./3.*pi.*(r*1e6).^3,'k')
hold on
xlabel('r_0')
ylabel('dF/dr')
set(gca,'fontsize',20)
end
 set(gcf,'color','w')
 
 
 % fig 13
 figure(13)
 u_star_vec = 1.6:0.2:2.4;
 for i = 1:length(u_star_vec)
     u_star = u_star_vec(i);
     SGF = dfdr(N5(u_star),r,th(u_star),TH(u_star),u_star);
     if ismember(u_star,[2 2.4])
         LS = '-';
     else
         LS = '--';
     end
     loglog(r*1e6,SGF*1e-6,'linewidth',4,'color',[1 1 1]*0.5,'linestyle',LS)
     hold on
 end
 set(gca,'xlim',[20 2000],'ylim',[1e-3 1e4])
 
 
 
 

 
 
 
 % fig 14 eq 24
figure(14)
% u_star = 1.65;
% u_star = 1.39;
U10_desired = [30 35];
for i_U10 = 1:2
u_star = fzero(@(u_star) U10(u_star)-U10_desired(i_U10),1.5);
subplot(1,2,i_U10)
for Omega     =[2.5 3.5]
    switch Omega
        case 2.5
            c = 'k';
        case 3.5
            c = [1 1 1]*0.5;
    end
    for eq_N = [7 9]
        switch eq_N
            case 7
                N = N7(u_star,Omega);
                LS = '-';
            case 9
                N = N9(u_star,Omega);
                LS = '--';

        end
        r = logspace(log10(30),log10(2000))*1e-6;
        term1 = T1(r,th(u_star),TH(u_star));
        term2 = T2(r,th(u_star),TH(u_star));
        SGF = dfdr(N,r,th(u_star),TH(u_star),u_star);
        
        loglog(r*1e6,SGF*1e-6,'color',c,'linewidth',3,'linestyle',LS)
        hold on
        set(gca,'ylim',[1e0 1e5],'xlim',[30 3000])
    end
end

 set(gcf,'color','w')

set(gca,'fontsize',20)

xlabel('$$r_0$$','interpreter','latex')
ylabel('$$\frac{dF}{dr} \left[m^{-2}s^{-1} \mu m^{-1}\right]$$','interpreter','latex')
end
title('SSGF','interpreter','latex')
% 75 4414

