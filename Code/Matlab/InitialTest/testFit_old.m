clearvars
clc
close all
clear all;


dataPosition = '../../../Data/IV-T_dependence_20250404_180103/';

temps = 10
%temps = 10:1:70;

offsets = [-1, 5]; % V
n_pulls = 2;
pulls = 0:n_pulls

flag_return_T = true;
flag_return_V = true;

voltages = [];
temperatures = [];
currents = [];
%voltages = zeroes( 50 * 2 * 3 * 61 * 2, 1); % 50 voltages, 2 directions (V_go and V_return), 3 pulls, 61 temperatures, 2 directions (T_go and T_return)
%temperatures = zeroes( 50 * 2 * 3 * 61 * 2, 1);
%currents = zeroes( 50 * 2 * 3 * 61 * 2, 1);


if flag_return_T
    ar_T_max = 2;
else
    ar_T_max = 1;
end

counter = 0;
first = true;
for T_direction = 1:ar_T_max
    if T_direction == 2
        temps = flip(temps);
    end
    for T = temps
        for pull = pulls
            filename = getFileName(T, offsets, pull, T_direction);

            raw_data = readmatrix(strcat(dataPosition, filename, '.txt'));

            Td_go = raw_data(:, 1);
            Td_return = raw_data(:, 2);
            Vcc = raw_data(:, 3);
            Vd_go = raw_data(:, 4);
            Vd_return = raw_data(:, 5);
            Err_Vd_go = raw_data(:, 6);
            Err_Vd_return = raw_data(:, 7);
            Id_go = raw_data(:, 12);
            Id_return = raw_data(:, 13);
            Err_Id_go = raw_data(:, 14);
            Err_Id_return = raw_data(:, 15);


            % add the data to the accumulation arrays
            

        

            voltages = [voltages ; Vd_go ; Vd_return];
            temperatures = [temperatures ; Td_go ; Td_return];
            currents = [currents ; Id_go ; Id_return];


%            plot3(Td_go, Vd_go, Id_go, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'black')
            if first
                hold on
                first = false;
            end
%            plot3(Td_return, Vd_return, Id_return, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'black')



            counter = counter +1;
        end


    end
end
counter * 50 * 2 % 50 voltages, 2 directions (V_go and V_return)


errorbar(Vd_go, Id_go, -Err_Id_go/2, Err_Id_go/2, -Err_Vd_go/2, Err_Vd_go/2, 'o', 'MarkerSize', 5, 'MarkerEdgeColor', 'blue')









p0 = [33, 6e-9];



I_exp = shockley(p0, Vd_go); % A
plot(Vd_go, I_exp, '--r', 'LineWidth', 2)


beta = nlinfit(Vd_go, Id_go, @shockley, p0)

I_exp = shockley(beta, Vd_go); % A
plot(Vd_go, I_exp, '--g', 'LineWidth', 2)


beta













function filename = getFileName(Temp, Vcc, pull, T_direction)
    if T_direction == 1
        T_direction = 'Go';
    elseif T_direction == 2
        T_direction = 'Return';
    end
    filename = sprintf('IV_T%.2f_V%.2f_%.2f_%dPull_T%s', Temp, Vcc(1), Vcc(2), pull, T_direction);
end


%{
function [I] = shockley(V, Vt, Is)
    % Shockley diode equation
    % Vt = k*T/q
    % Is = saturation current
    % V = voltage across the diode
    % T = temperature in Kelvin
    k = 1.380649e-23; % J/K
    q = 1.602176634e-19; % C
%    Vt = k*T/q; % V
    n = 1; % ideality factor (assumed to be 1 for simplicity)
    I = Is .* (exp(V./n.*Vt) - 1); % A
end
%}

function [I] = shockley(params, V)
    n = 2 % ideality factor
    I = params(2) .* ( exp(V./n*params(1)) - 1); % A
end




function [I] = tempDependance(Td, Eg, alpha)
    % Temperature dependence of saturation current
    k = 1.380649e-23; % J/K
%    Eg = 1.117; % eV
%    alpha = 2.5e-3; % A/K
    Is = alpha * Td.^2 .* exp(-Eg/(k.*Td)); % A
    I = Is; % return saturation current
end


function [I] = shockley3D(params, X)
    vv = X(1);
    tt = X(2);
    Is = tempDependance(tt, params(1), params(2));

    I = shockley([Is, n], vv);

end


% calculate the saturation current for each temperature


%Is = tempDependance(temperatures + 273.15, 2); % A

%I_expected = shockley(voltages, 0.7, temperatures + 273.15, Is(1)); % A




%for i = 1:length(temperatures)
%    Is(i) = tempDependance(temperatures(i) + 273.15, 2); % A
%end

% calculate the expected current for each voltage and temperature
%I_expected = zeros(length(voltages), length(temperatures));
%for i = 1:length(voltages)
%    for j = 1:length(temperatures)
%        I_expected(i, j) = shockley(voltages(i), 0.7, temperatures(j) + 273.15, Is(j)); % A
%    end
%end

%{
% plot the expected current for each temperature
figure
hold on
for i = 1:length(temperatures)
    plot3(temperatures, voltages, I_expected(:, i), 'DisplayName', sprintf('T = %.2f °C', temperatures(i)), 'LineWidth', 1.5)
end
xlabel('Temperature (C)')
ylabel('Voltage (V)')
zlabel('Current (A)')
title('Epected Current vs Voltage for Different Temperatures')
%legend('show')
grid on
hold off
%}



