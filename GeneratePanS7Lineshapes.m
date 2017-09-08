%% Written by David Feng, PhD candidate at Princeton University c. July 2017.

%%

clear; clc;
close all; 

k_B = 1.3806503e-23; % Boltzmann constant

lambda_0 = 532e-9; % m
theta = 90*(pi/180); % Angle between scatterer and observer / radians

% Comment in only one set of parameters (either for temperature or pressure)!
% For range of temperatures
T = [300 400]; % K
P_atm = ones(1,length(T))*1;
P_pa = 101325*P_atm;
n_den = P_pa.*(T*k_B).^(-1); % m^(-3)
T_norm = T(1); % K
n_norm = P_pa(1).*(T_norm*k_B).^(-1); % m^(-3)

% % For range of pressures
% P_atm = [1 2 3 4 5];
% T = ones(1,length(P_atm))*293; % K
% P_pa = 101325*P_atm; % Pa
% n_den = P_pa.*(T*k_B).^(-1); % m^(-3)
% P_norm = P_pa(end); % K
% n_norm = P_norm.*(T(1)*k_B).^(-1); % m^(-3)

N = 1000; % No. of points for half of the lineshape

gastypeinput = 1; % 1: Air, 2: CO2, 3: Ar, 4: C3H8

% Plot lineshapes normalized to n_norm
figure('DefaultAxesFontSize',18,'defaulttextfontsize',18,'defaultLineLineWidth',2)
for i = 1:length(T)
    
    P_in = P_atm(i);
    T_in = T(i);
    
    [x_Hz(:,i), lineshp(:,i), y] = Mikhail_S7_ver3_plot_S7_lineshapes...
        (gastypeinput, P_in, T_in, lambda_0, theta, N);
    
    x_GHz(:,i) = x_Hz(:,i)*1e-9;
    lineshp_cal(:,i) = (n_den(i) / n_norm) * lineshp(:,i); % Lineshape calibrated with respect to n_norm    
    lineshp_norm(:,i) = msnorm(x_Hz(:,i),lineshp_cal(:,i),'Max',1); % Lineshape normalized to unity
    
    plot(x_GHz(:,i), lineshp_cal(:,i), 'DisplayName', ['y = ' num2str(y, '%.2f'), ', T = ' num2str(T(i)), ' K, ' ...
        'P = ' num2str(P_atm(i), '%.1f'), ' atm'])
    
    hold on
    
end
xlabel('Detuning / GHz'); ylabel('Amplitude / A.U.');
legend('-DynamicLegend');

% Plot normalized to unity lineshapes
figure('DefaultAxesFontSize',18,'defaulttextfontsize',18,'defaultLineLineWidth',2)
for i = 1:length(T)
    
    plot(x_GHz(:,i), lineshp_norm(:,i), 'DisplayName', ['y = ' num2str(y, '%.2f'), ', T = ' num2str(T(i)), ' K, ' ...
        'P = ' num2str(P_atm(i), '%.1f'), ' atm'])
    
    hold on
    
end

xlabel('Detuning / GHz'); ylabel('Amplitude / A.U.');
legend('-DynamicLegend');
ylim([0 1])