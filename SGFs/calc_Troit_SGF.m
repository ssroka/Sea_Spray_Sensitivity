function [SGF,u_star,ReB] = calc_Troit_SGF(U10_desired,eq_N,Omega,r)
%{
 Author: Sydney Sroka
 Implementation of equation 24 in Troitskaya et al. (2018)
 Usage:
     [SGF] = calc_Troit_SGF(U10_desired,eq_N,Omega,r)
     
    INPUT:
       r     = an Nx1 vector of drop radii in [m]
       Omega = wave age parameter usually between [ 2.5 - 3.5 ]
        eq_N = number of bags observed in
                   laboratory conditions (7)
                   field conditions (9)
 U10_desired = a scalar wind speed in [m/s] which identifies the SGF

    OUTPUT:
        SGF = the function in units of
              the number of drops * m^-2 s^-1 micron^-1
%}

% constants
rho_w = 1020;       % saltwater density [kg m^-3]
N_rim = 8.3;        % constant
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

% functions of friction velocity u_star
z0        =@(u_star) alpha.*u_star.^2./g;
U10       =@(u_star) (u_star+0.14)/0.051;
R2_u_star =@(u_star) 9.6/u_star*1e-3; % convert from mm to m
C_D       =@(u_star) (U10(u_star)<=40)*(0.057-0.48/U10(u_star)).^2 ...
    + (U10(u_star)>40)*(2.37/U10(u_star)-0.012).^2;

if eq_N == 7
    % lab
    omega_p   =@(u_star) Omega_u0.*u_star.^(-gam);
elseif eq_N == 9
    % field
    omega_p   =@(U10) g*Omega/U10;
end

% windsea Reynolds number
ReB       =@(u_star,Omega) (0.051*U10(u_star)-0.14)^2*U10(u_star)./g./nu./Omega;


% equations for the number of bags <N> in untis of [m^-2 s^-1]

% eq 5 - for experiments, without ReB
N5 =@(u_star) Q0*u_star^2*exp(-U0^2/u_star^2);

% eq 7 - for real experiments with fetch
N7 =@(u_star,Omega)M0.*(ReB(u_star,Omega)).^(3./2).*exp(-M1./((ReB(u_star,Omega)).^(3./2)));

% % eq 9 - for field conditions
%%%%%%%%%%%%%%%%%%%%
% TYPO IN PAPER!!!!!!!!!!!!! correct eq 9 here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N9 =@(u_star,Omega)Q0*(u_star.^2.).*(12.4*(u_star.^-gam*U10(u_star))./(g*Omega)).^(1.5).*...
    exp(-U0.^2./u_star.^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% appendix B units of [m]
th =@(u_star) 0.001*(R2_u_star(u_star))^(4/3)*L^(-1/3);

% appendix C units of [m]
TH =@(u_star) 0.0021*(R2_u_star(u_star));

% eq 24 dfdr in units of [number of drops m^-2 s^-1 micron^-1]
dfdr =@(N,r,th,TH,u_star)...
    N.*(3.3e-9./L.*(rho_w.*g.*L.^2./sig).^(1.18).*(r./th).^7.3.*...
    exp(-5.2.*sinh(3./7.*log(r./th)))+...
    1.5e-4.*N_rim./TH.*(r./TH).^4.5.*...
    exp(-3.94.*sinh(0.5.*log(r./TH))));

% u_star is in m/s
u_star = fzero(@(u_star) U10(u_star)-U10_desired,1.5);

switch eq_N
    case 5
        N = N5(u_star);
    case 7 % laboratory
        N = N7(u_star,Omega);
    case 9 % field
        N = N9(u_star,Omega);
end

SGF = dfdr(N,r,th(u_star),TH(u_star),u_star)*1e-6;


% 1533











end