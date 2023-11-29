% Initial model of PTES system
clc; clear all;
% Parameters
rad_hot_tank = 0.5; %m
height_hot_tank = 4; %m
rad_cold_tank = 0.5; %m
height_cold_tank = 4; %m

Vol_cold_tank = rad_cold_tank^2*pi*height_cold_tank; %m3
Vol_hot_tank = rad_hot_tank^2*pi*height_hot_tank; %m3; 
R_gas = 287.05; % J/kgK
density_TES = 5000; % kg/m3
HeatCap_TES = 1; % kJ/kgK
void_fraction = 0.5;

% initial conditions
p_cold_tank = 201325;
p_hot_tank = 501325;

T_cold_tank = 273.15+100;
T_hot_tank = 273.15+20;

mass_cold_tank = p_cold_tank*Vol_cold_tank*void_fraction/(R_gas*T_cold_tank);
mass_hot_tank = p_hot_tank*Vol_hot_tank*void_fraction/(R_gas*T_hot_tank);

% we will have six state points.
T_gas = ones(6,1)*290;
p_gas = ones(6,1)*201325; % 1 bar gauge

% consider an increment of air leaving cold store
dm = 0.002; % kg of air
tau = 1; % s

% Easiest to discretise the thermals stores
n_slices = 400;
L_slice = (height_hot_tank/n_slices);
fprintf('The length of each slice is %.2f m \n', L_slice)
gas_vol_slice = void_fraction*rad_cold_tank^2*pi*height_cold_tank/n_slices;
heat_cap_slice = ((1-void_fraction)*rad_cold_tank^2*pi*height_cold_tank/n_slices)*density_TES*HeatCap_TES;
fprintf('The heat capacity per slice is %.2f kJ/K \n', heat_cap_slice)
Hot_TES = zeros(n_slices, 1) + T_hot_tank;
Cold_TES = zeros(n_slices, 1) + T_cold_tank;
A_slice = 7.4416;
h = 40;

cp_gas = 1; %kJ/kgK
cv_gas = 0.707; %kJ/kgK
gamma = cp_gas/cv_gas;

mg_slice = (p_hot_tank*gas_vol_slice)/(R_gas*T_hot_tank);
ms_slice = ((1-void_fraction)*rad_cold_tank^2*pi*height_hot_tank/n_slices)*density_TES;
%%% CHARGING

W_net_charge = 0;

mass_circulated = 0;


%%%%%%%%%%%%%%%%%%%%%%
% In order to stop the charging, specify a threshold temperature 
T_threshold = 320; %Kelvin
% This is broadly how the system would work in practice, controlled by a
% thermocouple. So when the temperature at the far end of the hot store
% reaches the threshold temperature the charging process stops.
%%%%%%%%%%%%%%%%%%%%%%

while Hot_TES(380)<T_threshold
    % stop the loop once the desired temperature is achieved

    mass_cold_tank = mass_cold_tank - dm;

    T_gas(1) = T_cold_tank;
    p_gas(1) = p_cold_tank;
    % new cold tank pressure and temp calculated at the end

    %%%%%%%%%%%%%% gas enters heat exchanger 1 %%%%%%%%%%%%%%
    T_inlet = T_gas(1);
    p_inlet = p_gas(1);
    [T_outlet, p_outlet] = Hex1(p_inlet, T_inlet);
    % set gas to the HEX outlet
    T_gas(2) = T_outlet;
    p_gas(2) = p_outlet;

    %%%%%%%%%%%%%% gas enters the hot machine %%%%%%%%%%%%%%
    T_inlet = T_gas(2);
    p_inlet = p_gas(2);
    [T_actual, T_outlet, p_outlet] = HotMachine(p_inlet, p_hot_tank, T_inlet, gamma);
    Work_hot = dm*cp_gas*(T_actual-T_inlet);
    % set gas to the compressor outlet
    T_gas(3) = T_actual;
    p_gas(3) = p_outlet;

    %%%%%%%%%%%%%% gas enters the hot store %%%%%%%%%%%%%%
    mass_hot_tank = mass_hot_tank + dm;
    
    % check vol of dm is less than vol slice
    vol_dm = dm*R_gas*T_gas(3)/p_gas(3);
    if vol_dm > gas_vol_slice
        fprintf('Warning - gas volume greater than slice volume \n')
        fprintf('Mass increment too large for store discretisation \n')
    end
    
    T_in_slice = T_gas(3);
    
    % gas exchanges heat with hot TES
    for j = 1:n_slices
        % update gas temp
        T_gas_equ = (1000*dm*cp_gas*T_in_slice + mg_slice*cp_gas*1000*Hot_TES(j))/((mg_slice+dm)*cp_gas*1000);
        % update solid temp
        T_solid_equ = (ms_slice*HeatCap_TES*1000*Hot_TES(j) + h*A_slice*(T_gas_equ - Hot_TES(j))*tau)/(ms_slice*HeatCap_TES*1000);
        % update gas temp
        T_gas_equ = ((mg_slice+dm)*cp_gas*1000*T_gas_equ - h*A_slice*(T_gas_equ - Hot_TES(j))*tau)/((mg_slice+dm)*cp_gas*1000);      
        % update temp in next slice
        Hot_TES(j) = T_solid_equ;
        T_in_slice = T_gas_equ;
        %update dm
        dm = dm + ((p_hot_tank*gas_vol_slice)/R_gas)*(1/T_hot_tank - 1/T_gas_equ);
    end
    
    % dm removed from hot tank
    mass_hot_tank = mass_hot_tank - dm;
    T_gas(4) = T_in_slice;
    p_gas(4) = p_gas(3);
    
    %%%%%%%%%%%%%% gas enters heat exchanger 2 %%%%%%%%%%%%%%
    T_inlet = T_gas(4);
    p_inlet = p_gas(4);
    [T_outlet, p_outlet] = Hex2(p_inlet, T_inlet);
    % set gas to the HEX outlet
    T_gas(5) = T_outlet;
    p_gas(5) = p_outlet;

    %%%%%%%%%%%%%% gas enters Cold Machine %%%%%%%%%%%%%%
    T_inlet = T_gas(5);
    p_inlet = p_hot_tank;
    [T_actual, T_outlet, p_outlet] = ColdMachine(p_inlet, p_cold_tank, T_inlet, gamma);
    Work_cold = dm*cp_gas*(T_actual-T_inlet);
    % set gas to the compressor outlet
    T_gas(6) = T_actual;
    p_gas(6) = p_outlet; 
    
    %%%%%%%%%%%%%% gas enters the cold store %%%%%%%%%%%%%%
    mass_cold_tank = mass_cold_tank + dm;
    
    T_in_slice = T_gas(6);
    
    % gas exchanges heat with cold TES
    for j = 1:n_slices
        % update gas temp
        T_gas_equ = (1000*dm*cp_gas*T_in_slice + mg_slice*cp_gas*1000*Cold_TES(j))/((mg_slice+dm)*cp_gas*1000);
        % update solid temp
        T_solid_equ = (ms_slice*HeatCap_TES*1000*Cold_TES(j) - h*A_slice*(Cold_TES(j)-T_gas_equ)*tau)/(ms_slice*HeatCap_TES*1000);
        % update gas temp
        T_gas_equ = ((mg_slice+dm)*cp_gas*1000*T_gas_equ + h*A_slice*(Cold_TES(j)-T_gas_equ)*tau)/((mg_slice+dm)*cp_gas*1000);      
        % update temp in next slice
        Cold_TES(j) = T_solid_equ;
        T_in_slice = T_gas_equ;
    end
    
    % Store the net work for each mass incrmement 
    W_net_charge = W_net_charge + Work_hot + Work_cold;
    % keep account of the total mass of gas circulated
    mass_circulated = mass_circulated+dm;
    
end

kWh = 3600; % convert the kJ into kWh

fprintf('System has stored %.2f kWh of work \n', (W_net_charge/kWh))
fprintf('%.2f kg of gas circulated \n', mass_circulated)
figure;
x = linspace(0,height_hot_tank,n_slices+1);
x = x(2:end);
h1 = plot(x, Hot_TES);
hold on
h2 = plot(x, Cold_TES);
xlim([0 height_hot_tank])
ylim([200 550])
xlabel('length (m)')
ylabel('Temperature (K)')
legend([h1,h2], 'Hot store after charging', 'cold store after charging')

%%% DISCHARGING

% Since the fluid flow is reversed during teh discharging, then we flip the
% temperature profiles in the TES's for convenience
Hot_TES_dis = flip(Hot_TES);
Cold_TES_dis = flip(Cold_TES);

T_gas_dis = zeros(6,1);
p_gas_dis = zeros(6,1);

W_net_dis = 0;
mass_circulated_dis = 0;

%%%%%%%%%%%%%%%%%%%%%%
% IMPORTANT: How do we stop the discharge process? To stop the discharge,
% we can either use a temperature threshold (this would probably be the
% case in reality), or we can recirculate the same amount of mass. Actually
% we opt for just less.
%%%%%%%%%%%%%%%%%%%%%%

while mass_circulated_dis<mass_circulated - 1

    mass_cold_tank = mass_cold_tank - dm;

    T_gas_dis(1) = Cold_TES_dis(n_slices);
    p_gas_dis(1) = p_cold_tank;
    % new cold tank pressure and temp calculated at the end

    %%%%%%%%%%%%%% gas enters cold HEX %%%%%%%%%%%%%%
    T_inlet = T_gas_dis(1);
    p_inlet = p_gas_dis(1);
    [T_outlet, p_outlet] = Hex2(p_inlet, T_inlet);
    % set gas to the HEX outlet
    T_gas_dis(2) = T_outlet;
    p_gas_dis(2) = p_outlet;

    %%%%%%%%%%%%%% gas enters the cold machine %%%%%%%%%%%%%%
    T_inlet = T_gas_dis(2);
    p_inlet = p_gas_dis(2);
    [T_actual, T_outlet, p_outlet] = ColdMachine(p_inlet, p_hot_tank, T_inlet, gamma);
    Work_cold = dm*cp_gas*(T_actual-T_inlet);
    % set gas to the compressor outlet
    T_gas_dis(3) = T_actual;
    p_gas_dis(3) = p_outlet;

    %%%%%%%%%%%%%% gas enters the hot store %%%%%%%%%%%%%%
    mass_hot_tank = mass_hot_tank + dm;
    
    % check vol of dm is less than vol slice
    vol_dm = dm*R_gas*T_gas_dis(3)/p_gas(3);
    if vol_dm > gas_vol_slice
        % Hopefully we shouldn't  trigger this warning during discharge
        fprintf('Warning - gas volume greater than slice volume \n')
    end
    
    T_in_slice = T_gas_dis(3);
    
    % gas exchanges heat with hot TES
    for j = 1:n_slices
        % update gas temp
        T_gas_equ = (1000*dm*cp_gas*T_in_slice + mg_slice*cp_gas*1000*Hot_TES_dis(j))/((mg_slice+dm)*cp_gas*1000);
        % update solid temp
        T_solid_equ = (ms_slice*HeatCap_TES*1000*Hot_TES_dis(j) - h*A_slice*(Hot_TES_dis(j)-T_gas_equ)*tau)/(ms_slice*HeatCap_TES*1000);
        % update gas temp
        T_gas_equ = ((mg_slice+dm)*cp_gas*1000*T_gas_equ + h*A_slice*(Hot_TES_dis(j)-T_gas_equ)*tau)/((mg_slice+dm)*cp_gas*1000);      
        % update temp in next slice
        Hot_TES_dis(j) = T_solid_equ;
        T_in_slice = T_gas_equ;
    end
    
    % dm removed from hot tank
    mass_hot_tank = mass_hot_tank - dm;
    T_gas_dis(4) = T_in_slice;
    p_gas_dis(4) = p_gas_dis(3);
    
    %%%%%%%%%%%%%% gas enters the hot HEX %%%%%%%%%%%%%%
    T_inlet = T_gas_dis(4);
    p_inlet = p_gas_dis(4);
    [T_outlet, p_outlet] = Hex1(p_inlet, T_inlet);
    % set gas to the HEX outlet
    T_gas_dis(5) = T_outlet;
    p_gas_dis(5) = p_outlet;

    %%%%%%%%%%%%%% gas enters hot machine %%%%%%%%%%%%%%
    T_inlet = T_gas_dis(5);
    p_inlet = p_hot_tank;
    [T_actual, T_outlet, p_outlet] = HotMachine(p_inlet, p_cold_tank, T_inlet, gamma);
    Work_hot = dm*cp_gas*(T_actual-T_inlet);
    % set gas to the compressor outlet
    T_gas_dis(6) = T_actual;
    p_gas_dis(6) = p_outlet; 
    
    %%%%%%%%%%%%%% gas enters the cold store %%%%%%%%%%%%%%
    mass_cold_tank = mass_cold_tank + dm;
    
    T_in_slice = T_gas_dis(6);
    
    % gas exchanges heat with cold TES
    for j = 1:n_slices
        % update gas temp
        T_gas_equ = (1000*dm*cp_gas*T_in_slice + mg_slice*cp_gas*1000*Cold_TES_dis(j))/((mg_slice+dm)*cp_gas*1000);
        % update solid temp
        T_solid_equ = (ms_slice*HeatCap_TES*1000*Cold_TES_dis(j) - h*A_slice*(Cold_TES_dis(j)-T_gas_equ)*tau)/(ms_slice*HeatCap_TES*1000);
        % update gas temp
        T_gas_equ = ((mg_slice+dm)*cp_gas*1000*T_gas_equ + h*A_slice*(Cold_TES_dis(j)-T_gas_equ)*tau)/((mg_slice+dm)*cp_gas*1000);      
        % update temp in next slice
        Cold_TES_dis(j) = T_solid_equ;
        T_in_slice = T_gas_equ;
    end
    
W_net_dis = W_net_dis + Work_hot + Work_cold;
mass_circulated_dis = mass_circulated_dis+dm;
    
end

fprintf('System has released %.2f kWh of work \n', (-W_net_dis/kWh))
fprintf('System efficiency: %.2f \n', (-100*W_net_dis/W_net_charge))
fprintf('%.2f kg of gas circulated \n', mass_circulated_dis)
figure;
h1 = plot(x, Hot_TES_dis);
hold on
h2 = plot(x, Cold_TES_dis);
xlim([0 height_hot_tank])
ylim([200 550])
xlabel('length (m)')
ylabel('Temperature (K)')
legend([h1,h2], 'Hot store after discharging', 'cold store after discharging')



