classdef mpm_sim
  properties
    name
    r_t
    T_s_t
    time_vec
    ic
  end
  methods
    function obj = mpm_sim(name,r_t,T_s_t,time_vec,ic)
       obj.name = name;
       obj.r_t = r_t;
       obj.T_s_t = T_s_t;
       obj.time_vec = time_vec;
       obj.ic = ic;
    end
  end
end
