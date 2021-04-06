paths_2_add = {
    '~/Documents/MATLAB/mpm/sandbox/util';
    '~/Documents/MATLAB/util';
    '~/Documents/MATLAB/mpm/sandbox/mm_vars/';
    '~/Documents/MATLAB/mpm/sandbox/Q_s/';
    '~/Documents/MATLAB/mpm/sandbox/NayarFxns/';
    '~/Documents/MATLAB/util/distinguishable_colors';
    '~/Documents/MATLAB/util/othercolor';
    '~/Documents/MATLAB/mpm_figures/';
    '/Users/ssroka/MIT/Research/EmanuelGroup/thesis/chapter2/SGFs';
    '/Users/ssroka/Documents/MATLAB/mpm/sandbox/dropletAcceleration';
    '/Users/ssroka/Documents/MATLAB/mpm/sandbox/approximationEndpoints';
    '/Users/ssroka/Documents/MATLAB/mpm/sandbox/mm_integration'
    };
switch action_str
    case 'add'
        for ipath = 1:length(paths_2_add)
            addpath(paths_2_add{ipath})
        end
    case 'remove'
        for ipath = 1:length(paths_2_add)
            rmpath(paths_2_add{ipath})
        end
    otherwise
        error('Please select a valid path action, \n''add'' or ''remove''')
end


