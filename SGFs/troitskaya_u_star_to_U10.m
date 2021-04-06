% Troitskaya 2018 a (JPO)
ccc
% neutral stability before eq 1

kappa = 0.4; % von karman
H10 = 10; % m
alpha = 0.0057;
g = 9.81;% m/s/s

z0 = @(u_star)alpha*u_star^2/g;

U10 = @(u_star) u_star/kappa*log(H10/z0(u_star))


u_star1 = 1.1; % where bag breakup becomes the dominant mechanism in experiments
U10(u_star1)

u_star2 = 1.51; % max wind speed tested
U10(u_star2)

u_star = 3.25; 
U10(u_star)

%% Andreas 2011
ccc

k = 0.4; % von karman
H10 = 10; % m
alpha = 0.0185;
nu = 1.48e-5; % kniematic viscosity of air
g = 9.81;% m/s/s

u_star = linspace(0.1,3);

z0 =@(u_star) 0.135.*nu./u_star+alpha.*u_star.^2./g;

U10 = @(u_star) u_star./k.*log(H10./z0(u_star));

CD10 = (k./log(H10./z0(u_star))).^2;


plot(U10(u_star),CD10)
