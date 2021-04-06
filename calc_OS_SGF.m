function [SSGF] = calc_OS_SGF(U10,r)
% r in microns
U10_fm_paper = [36 40.5 45 49.5 54];

if ~(ismember(U10,U10_fm_paper))
    error(' U10 must be one of 36 40.5 45 49.5 54 for Ortiz-Suslow')
elseif any(r<86)
    error('radii must be > 86 microns for OrtizSuslow 2016')
end
    
Table1_data = [...
36    -2.1984 7.1125 -4.0758
40.5  -2.1556 7.1772 -3.9759
45    -1.9576 6.4917 -3.1024
49.5  -2.1235 7.3997 -3.9673
54    -2.1009 7.3857 -3.8341];

R = log10(r);
p = Table1_data(U10==Table1_data(:,1),2:4);
S=@(R) p(1).*R.^2+p(2).*R+p(3); % S is the log-scaled source fxn
SSGF = 10.^(S(R)); % convert to regular source function

    















