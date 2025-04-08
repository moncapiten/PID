clearvars
clc
close all
clear all;


dataPosition = '../../../Data/IV-T_dependence_20250408_110423/';

temps = 16:1:70;

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
%            if first
%                hold on
%                first = false;
%            end
%            plot3(Td_return, Vd_return, Id_return, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'black')



            counter = counter +1;
        end


    end
end
counter * 50 * 2 % 50 voltages, 2 directions (V_go and V_return)


%errorbar(Vd_go, Id_go, -Err_Id_go/2, Err_Id_go/2, -Err_Vd_go/2, Err_Vd_go/2, 'o', 'MarkerSize', 5, 'MarkerEdgeColor', 'blue')


%plot3()

plot3(temperatures, voltages, abs(currents), 'ob', 'DisplayName', 'Raw Data')
grid on
grid minor
hold on

% set z log scale
set(gca, 'ZScale', 'log')

%{
X = [voltages, temperatures];
size(X)
p0 = [2.5e-3, 1.117, 0.7]; % [A, Eg, Vth]
%p0 = [33, 6e-9]; % [Vt, Is]

I_sym = complete_Shockley(p0, X);
%I_sym = shockley3D(p0, X);
%plot3(X(:, 2), X(:, 1), I_sym, 'vg', 'DisplayName', 'p0')

beta = nlinfit(X, currents, @complete_Shockley, p0);
%beta
I_fit = complete_Shockley(beta, X);

plot3(X(:, 2), X(:, 1), I_fit, 'or', 'DisplayName', 'Fitted Model')

legend()

xlabel('Temperature (C)')
ylabel('Voltage (V)')
zlabel('Current (A)')
title('Expected Current vs Voltage for Different Temperatures')




%}










function filename = getFileName(Temp, Vcc, pull, T_direction)
    if T_direction == 1
        T_direction = 'Go';
    elseif T_direction == 2
        T_direction = 'Return';
    end
    filename = sprintf('IV_T%.2f_V%.2f_%.2f_%dPull_T%s', Temp, Vcc(1), Vcc(2), pull, T_direction);
end




function [I] = shockley(params, V)
    n = 2; % ideality factor
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
    vv = X(:, 1);
    tt = X(:, 2);
    Is = tempDependance(tt, params(1), params(2));

    k = 1.380649e-23; % J/K
    q = 1.602176634e-19; % C
    Vt = k/q .* tt; % V

    I = shockley([Vt, Is], vv);

end



function [I] = complete_Shockley(params, X) % params = [A, Eg, Vth]
    vv = X(:, 1);
    TT = X(:, 2);
    n = 2;
    k = 1.380649e-23; % J/K
    q = 1.602176634e-19; % C
    Is = params(1) .* TT.^2 .* exp(-params(2)/(n*k.*TT)); % A
    Vt = k/q .* TT; % V
    I = Is .* (exp(vv./(n.*params(3))) - 1); % A
end

%beta(2) = Eg viene esatto 1.117? WTF???????????????????????


