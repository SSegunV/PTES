function [T_out,p_out] = HX1(p_in,T_in,T_in_f)
%HEX
eps = 1;
% take an increment of gas and cool/heat
if T_in < T_in_f
% heat
T_out = eps*(T_in_f - T_in) + T_in;
% cool
else
T_out = T_in - eps*(T_in - T_in_f);
end
p_out = p_in;