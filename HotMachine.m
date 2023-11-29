function [T_act_hot,T_out,p_out] = HotMachine(p_in, p_target, T_in, gamma)
%COMPRESSOR
turbine_eff = 1;
compressor_eff = 1;

% take an increment of gas and compress
p_out = p_target; % This could be updated later
r = (p_out/p_in);
T_out = T_in * r^((gamma-1)/gamma);%use pressures for efficiencies

%conditional loop for compression and expansion
if p_target < p_in
    T_act_hot = T_in - turbine_eff*(T_in - T_out);
    
else
    T_act_hot = (T_out - T_in)/compressor_eff + T_in;
    
end


end