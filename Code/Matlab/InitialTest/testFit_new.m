clearvars
clc
close all
clear all;

chooseDiode = 5; % 1: PN diode, 2: Schottky diode, 3: Zener diode 1, 4: Zener diode 2, 5: PN diode long take


names = {'PN diode', 'Schottky diode', 'Zener diode 1', 'Zener diode 2', 'PN diode long take'};




switch chooseDiode
    case 1
        offsets = [-1, 5]; % PN diode
        temps = 16:1:70;
        dataPosition = '../../../Data/IV-T_dependence_20250404_100546/';  % PN diode
    case 2
        offsets = [-1, 5]; % Schottky diode
        temps = 10:1:70;
        dataPosition = '../../../Data/IV-T_dependence_20250408_110423/'; % Schottky diode
    case 3
        offsets = [-3, 5]; % Zener diode 1
        temps = 10:1:70;
        dataPosition = '../../../Data/IV-T_dependence_20250408_124236/'; % Zener diode 1
    case 4
        offsets = [-5, 5]; % Zener diode 2
        temps = 10:1:70;
        dataPosition = '../../../Data/IV-T_dependence_20250411_094932/'; % Zener diode 2
    case 5
        offsets = [0, 5]; % PN diode long take
        temps = 10:0.5:70;
        dataPosition = '../../../Data/IV-T_dependence_20250414_092701/';  % PN diode long take
end
%offsets = [-1, 5]; % V
n_pulls = 2;
pulls = 0:n_pulls

flag_return_T = true;
flag_return_V = true;

voltages = [];
temperatures = [];
currents = [];


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






            counter = counter +1;
        end


    end
end
counter * 50 * 2 % 50 voltages, 2 directions (V_go and V_return)


set(gca, 'YScale', 'log')
legend('show')



temperatures = temperatures + 273.15; % C
currents = currents * 1e6; % A


X = [voltages'; temperatures'];

%p0 = { { 2e-8, 3/2, 1.8, 2 }, { 2e-8, 2, 0.3, 2 }, {}, {} }; % [I0, gam, Eg, n] PN diode

%p0 = [2e-8, 3/2, 1.1, 2 ]; % [I0, gam, Eg, n] PN diode
%p0 = [2e-8, 2, 0.3, 2 ]; % [I0, gam, Esh, n] Schottky diode
%p0 = [2e-8, 1.8, 1.8, 4 ]; % [I0, gam, Eg, n] Zener diode 1



%beta = nlinfit(X, currents', @PN_model, p0, 'Options', statset('MaxIter', 1000));
%beta = nlinfit(X, currents', @Schottky_model, p0, 'Options', statset('MaxIter', 1000));
%beta
switch chooseDiode
    case 1
        global p0
        p0 = [2e-8, 3/2, 1.1, 2 ]; % [I0, gam, Eg, n] PN diode
        beta = nlinfit(X, currents', @PN_model, p0, 'Options', statset('MaxIter', 1000));
        I_fit = PN_model(beta, X);
    case 2
        global p0
        p0 = [2e-8, 2, 0.3, 2 ]; % [I0, gam, Esh, n] Schottky diode
        beta = nlinfit(X, currents', @Schottky_model, p0, 'Options', statset('MaxIter', 1000));
        I_fit = Schottky_model(beta, X);
    case 3
        global p0
        p0 = [2e-8, 1.8, 1.8, 4 ]; % [I0, gam, Eg, n] Zener diode 1
        beta = nlinfit(X, currents', @Zener_model, p0, 'Options', statset('MaxIter', 1000));
        I_fit = Zener_model(beta, X);
    case 4
        global p0
        p0 = [2e-8, 1.8, 1.8, 4 ]; % [I0, gam, Eg, n] Zener diode 2
        beta = nlinfit(X, currents', @Zener_model, p0, 'Options', statset('MaxIter', 1000));
        Zener_model
    case 5
        global p0
        p0 = [2e-8, 3/2, 1.1, 2 ]; % [I0, gam, Eg, n] PN diode long take
        beta = nlinfit(X, currents', @PN_model, p0, 'Options', statset('MaxIter', 1000));
        I_fit = PN_model(beta, X);
end


%I_fit = PN_model(beta, X);
%I_fit = Schottky_model(beta, X);






t = tiledlayout(1, 2, "TileSpacing", "Tight", "Padding", "Compact");
t1 = nexttile(t);

%errorbar(Vd_go, Id_go, -Err_Id_go/2, Err_Id_go/2, -Err_Vd_go/2, Err_Vd_go/2, 'o', 'MarkerSize', 5, 'MarkerEdgeColor', 'blue')
plot3(temperatures, voltages, (currents'), 'ob', 'DisplayName', 'Raw Data')             % Raw Data
hold on

%plot3(X(2, :), X(1, :), PN_model(p0, X), 'og', 'DisplayName', 'PN p0')                 % PN p0
%plot3(X(2, :), X(1, :), Schottky_model(p0, X), 'vg', 'DisplayName', 'Schottky p0')     % Schottky p0
%plot3(X(2, :), X(1, :), Zener_model(p0, X), 'k', 'DisplayName', 'Zener p0')            % Zener p0

plot3(X(2, :), X(1, :), (I_fit'), 'vr', 'DisplayName', 'Fitted Model')                  % Fitted Model


hold off
grid on
grid minor

xlabel('Temperature [K]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('Voltage [V]', 'Interpreter', 'latex', 'FontSize', 14)
zlabel('Current ($ \mathrm{ \mu A } $) ', 'Interpreter', 'latex', 'FontSize', 14)
title('Linear Plot', 'Interpreter', 'latex', 'FontSize', 16)
legend()



t2 = nexttile(t);
plot3(temperatures, voltages, abs(currents'), 'ob', 'DisplayName', 'Raw Data')
hold on
plot3(X(2, :), X(1, :), abs(I_fit'), 'vr', 'DisplayName', 'Fitted Model')


hold off
grid on
grid minor

xlabel('Temperature [K]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('Voltage [V]', 'Interpreter', 'latex', 'FontSize', 14)
zlabel('Current ($ \mathrm{ \mu A } $) ', 'Interpreter', 'latex', 'FontSize', 14)
title('Semilog Plot', 'Interpreter', 'latex', 'FontSize', 16)
legend()
set(gca, 'ZScale', 'log')


title(t, strcat("IV-T curve for ", names{chooseDiode}), 'Interpreter', 'latex', 'FontSize', 18)









function filename = getFileName(Temp, Vcc, pull, T_direction)
    if T_direction == 1
        T_direction = 'Go';
    elseif T_direction == 2
        T_direction = 'Return';
    end
    filename = sprintf('IV_T%.2f_V%.2f_%.2f_%dPull_T%s', Temp, Vcc(1), Vcc(2), pull, T_direction);
end





function [I] = PN_model(params, X) % params = [I0, gam, Eg, n]
    vv = X(1, :);
    tt = X(2, :);

    k = 8.617333e-5; % eV/K
    Vt = 13e-3; % V
%    Vt = k * tt / 1.602176634e-19 % V
    Is = params(1) .* tt.^params(2) .* exp(-params(3)./(params(4).*k.*tt)); % A
    I = Is .* (exp(vv./(params(4).*Vt)) - 1) * 1e6; % A
%    I = swap_I * 1e9; % A
end

function [I] = Schottky_model(params, X) % params = [I0, gam, Esh, n]
    vv = X(1, :);
    tt = X(2, :);

    k = 8.617333e-5; % eV/K
    Vt = 13e-3; % V
    Is = params(1) .* tt.^params(2) .* exp(-params(3)./(k.*tt)); % A
    I = Is .* (exp(vv./(params(4).*Vt)) - 1) * 1e6; % A
%    I = swap_I * 1e9; % A
end

function [I] = Zener_model(params, X) % params = [I0, gam, Eg, n]
    vv = X(1, :);
    tt = X(2, :);

    k = 8.617333e-5; % eV/K
    Vt = 13e-3; % V
    Is = params(1) .* tt.^params(2) .* exp(-params(3)./(params(4).*k.*tt)); % A
    I = Is .* (exp(vv./(params(4).*Vt)) - 1) * 1e6; % A
end


