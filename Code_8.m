%% Valve Lift and Flow CharacteristicsAnalysis

clc; clear; close all;

 
% Constants
gamma = 1.4;                %Specific heat ratio
R = 287;           % J/kg·K (gas constant)
Tt = 300;                  % K (throat temperature)
Pt = 105e3;                 %Pa (throatpressure)

Pb = 101e3;                 %Pa (backpressure)
A = 0.01;                  % m² (throat area - assumed)
P_critical = Pb / 0.528;    % Critical pressure (~191.3 kPa)
 

% Valve parameters (example values)
theta_vo = 10;              % Valve opening start angle (degrees)
theta_vc = 240;             % Valve closing angle (degrees)
theta_dur = theta_vc - theta_vo; % Valve duration (degrees)
l_max = 0.012;              % Maximum valve lift (m)
Cd = 0.7;                
  % Discharge coefficient

 
% Crankshaft angle range (0 to 720 degrees for 2 revolutions)
theta = linspace(0, 720, 1000);
 
% Valve lift calculation
l_theta = zeros(size(theta));
for i = 1:length(theta)
    if theta(i) >= theta_vo&& theta(i) <= theta_vc
        l_theta(i) = l_max * sin(pi * (theta(i) - theta_vo) / theta_dur);
    end
end
 
% Upstream pressure range (50 kPa to 300 kPa)
Pu = linspace(50e3, 300e3, 500);  % Pa
 
% Initialize arrays
Vt = zeros(size(Pu));
mdot = zeros(size(Pu));
 
% Calculate throat velocity and mass flow rate
for i = 1:length(Pu)
    if Pu(i) >= P_critical
        % Sonic flow conditions
        Vt(i) = sqrt(gamma * R * Tt);
        term =gamma * (2/(gamma+1))^((gamma+1)/(gamma-1));
        mdot(i) = (A * Pu(i) / sqrt(R*Tt)) * sqrt(term);
    else
        % Subsonic flow conditions
        Vt(i) = sqrt((2*gamma*R*Tt)/(gamma-1) * (1 - (Pt/Pu(i))^((gamma-1)/gamma)));
        term = (2*gamma/(gamma-1)) * ((Pt/Pu(i))^(2/gamma) - (Pt/Pu(i))^((gamma+1)/gamma));
        mdot(i) = (A * Pu(i) / sqrt(R*Tt)) * sqrt(term);
    end

end
 
% Convert pressure to kPa
Pu_kPa = Pu / 1e3;
P_critical_kPa = P_critical / 1e3;
 
%% Plot Valve Lift vs Crankshaft Angle
figure;
plot(theta, l_theta*1000, 'b-', 'LineWidth', 2);
xlabel('Crankshaft Angle(degrees)', 'FontSize', 12);
ylabel('Valve Lift (mm)', 'FontSize', 12);
title('Intake/Exhaust ValveLift vs Crankshaft Angle', 'FontSize', 14);
grid on;
xlim([0 720]);
set(gca,'XTick',0:90:720); 
legend('Valve Lift', 'Location', 'northeast');
 
%% Plot Throat Velocity vs UpstreamPressure
figure;
plot(Pu_kPa, Vt, 'r-', 'LineWidth', 2);
hold on;
line([P_critical_kPa P_critical_kPa], [min(Vt) max(Vt)], 'Color', 'k', 'LineStyle', '--');
xlabel('Upstream Pressure(kPa)', 'FontSize', 12);
ylabel('Throat Velocity (m/s)', 'FontSize', 12);
title('Throat Velocity vs.Upstream Pressure', 'FontSize', 14);
legend('Throat Velocity', 'Sonic/Subsonic Threshold', 'Location', 'southeast');
grid on;
xlim([min(Pu_kPa) max(Pu_kPa)]);
 
%% Plot Mass Flow Rate vs Upstream
Pressure
figure;
plot(Pu_kPa, mdot, 'g-', 'LineWidth', 2);
hold on;
line([P_critical_kPa P_critical_kPa], [min(mdot) max(mdot)], 'Color', 'k', 'LineStyle', '--');
xlabel('Upstream Pressure(kPa)', 'FontSize', 12);
ylabel('Mass Flow Rate (kg/s)', 'FontSize', 12);
title('Mass Flow Rate vs.Upstream Pressure', 'FontSize', 14);
legend('Mass Flow Rate', 'Sonic/Subsonic Threshold', 'Location', 'northwest');
grid on;
xlim([min(Pu_kPa) max(Pu_kPa)]);
 
%% Display critical pressure
fprintf('Critical Pressure(Sonic/Subsonic Threshold): %.2f kPa\n', P_critical_kPa);
%% Plot Throat Velocity vs Upstream
Pressure
figure;
plot(Pu_kPa, Vt, 'r-', 'LineWidth', 2);
hold on;
line([P_critical_kPaP_critical_kPa], [min(Vt) max(Vt)], 'Color', 'k', 'LineStyle', '--');
xlabel('Upstream Pressure(kPa)', 'FontSize', 12);
ylabel('Throat Velocity (m/s)', 'FontSize', 12);
title('Throat Velocity vs.Upstream Pressure', 'FontSize', 14);
legend('Throat Velocity', 'Sonic/Subsonic Threshold', 'Location', 'southeast');
grid on;
xlim([min(Pu_kPa) max(Pu_kPa)]);
 
disp('(Vt):');
disp(Vt(1:10));