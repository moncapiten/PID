clearvars
clc
close all
clear all;

chooseDiode = 6; % 1: PN diode, 2: Schottky diode, 3: Zener diode 1, 4: Zener diode 2, 5: PN diode long take, 6: Schottky diode long take


names = {'PN diode', 'Schottky diode', 'Zener diode 1', 'Zener diode 2', 'PN diode long take', 'Schottky diode long take'};
    


threshold = 0.35; % V
switch chooseDiode
    case 1
        offsets = [-1, 5]; % PN diode
        temps = 16:1:70;
        dataPosition = '../../Data/IV-T_dependence_20250404_100546/';  % PN diode
    case 2
        offsets = [-1, 5]; % Schottky diode
        temps = 10:1:70;
        dataPosition = '../../Data/IV-T_dependence_20250408_110423/'; % Schottky diode
    case 3
        offsets = [-3, 5]; % Zener diode 1
        temps = 10:1:70;
        dataPosition = '../../Data/IV-T_dependence_20250408_124236/'; % Zener diode 1
    case 4
        offsets = [-5, 5]; % Zener diode 2
        temps = 10:1:70;
        dataPosition = '../../Data/IV-T_dependence_20250411_094932/'; % Zener diode 2
    case 5
        offsets = [0, 5]; % PN diode long take
        temps = 10:0.5:70;
        dataPosition = '../../Data/IV-T_dependence_20250414_092701/';  % PN diode long take
    case 6
        offsets = [-1, 5]; % Schottky diode long take
        temps = 10:0.5:70;
        dataPosition = '../../Data/IV-T_dependence_20250415_113544/';  % Schottky diode long take
        threshold = 0.15; % V
end

n_pulls = 2;
pulls = 0:n_pulls

flag_return_T = true;
flag_return_V = true;

voltages = [];
temperatures = [];
currents = [];

voltages2fit = [];
temperatures2fit = [];
currents2fit = [];

if chooseDiode == 3 || chooseDiode == 4
    negthreshold = -1; % V
    negvoltages2fit = [];
    negtemperatures2fit = [];
    negcurrents2fit = [];
end


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





            mask = voltages > threshold;
            voltages2fit = [voltages2fit ; voltages(mask)];
            temperatures2fit = [temperatures2fit ; temperatures(mask)];
            currents2fit = [currents2fit ; currents(mask)];

            if chooseDiode == 3 || chooseDiode == 4
                mask = voltages < negthreshold;
                negvoltages2fit = [negvoltages2fit ; voltages(mask)];
                negtemperatures2fit = [negtemperatures2fit ; temperatures(mask)];
                negcurrents2fit = [negcurrents2fit ; currents(mask)];
            end



            counter = counter +1;
        end


    end
end
counter * 50 * 2 % 50 voltages, 2 directions (V_go and V_return)


%set(gca, 'YScale', 'log')
%legend('show')










temperatures = temperatures + 273.15; % C
currents = currents * 1e6; % A

temperatures2fit = temperatures2fit + 273.15; % C
currents2fit = currents2fit * 1e6; % A

if chooseDiode == 3 || chooseDiode == 4
    negtemperatures2fit = negtemperatures2fit + 273.15; % C
    negcurrents2fit = negcurrents2fit * 1e6; % A
end


%X = [voltages'; temperatures'];
X = [voltages2fit'; temperatures2fit'];



switch chooseDiode
    case 1
        global p0
        p0 = [2e-8, 3/2, 1.1, 2 ]; % [I0, gam, Eg, n] PN diode
        lb = [0, 0, 0, 1];
        ub = [1e-6, 3, 2, 2];
        beta = lsqcurvefit(@PN_model, p0, X, currents2fit', lb, ub, optimset('MaxIter', 10000, 'MaxFunEvals', 10000));
        %[beta, R, ~, CovB] = nlinfit(X, currents2fit', @PN_model, p0, 'Options', statset('MaxIter', 10000));
        I_fit = PN_model(beta, X);
    case 2
        global p0
        p0 = [2e-8, 2, 0.3, 2 ]; % [I0, gam, Esh, n] Schottky diode
        [beta, R, ~, CovB] = nlinfit(X, currents2fit', @Schottky_model, p0, 'Options', statset('MaxIter', 1000));
        I_fit = Schottky_model(beta, X);
    case 3
        global p0
        p0 = [2e-8, 1.8, 1.8, 4 ]; % [I0, gam, Eg, n] Zener diode 1
%        [beta, R, ~, CovB] = nlinfit(X, currents', @Zener_model, p0, 'Options', statset('MaxIter', 1000));
        [beta, R, ~, CovB] = nlinfit(X, currents2fit', @Zener_model, p0, 'Options', statset('MaxIter', 1000));
        I_fit = Zener_model(beta, X);
    case 4
        global p0
        p0 = [2e-8, 1.8, 1.8, 4 ]; % [I0, gam, Eg, n] Zener diode 2
        [beta, R, ~, CovB] = nlinfit(X, currents2fit', @Zener_model, p0, 'Options', statset('MaxIter', 1000));
        I_fit = Zener_model(beta, X);
    case 5
        global p0
        p0 = [2e-8, 3/2, 1.1, 2 ]; % [I0, gam, Eg, n] PN diode long take
        [beta, R, ~, CovB] = nlinfit(X, currents2fit', @PN_model, p0, 'Options', statset('MaxIter', 1000));
        I_fit = PN_model(beta, X);
    case 6
        global p0
        p0 = [0.0154, 2, 0.2181, 3.2533 ];
        % 0.0002    2.6306    0.2013    3.2521  
        %        p0 = [5.0998, 0.7876, 1 ]; % [I0, Esh, n] Schottky diode long take
%        p0 = [2e-8, 2, 0.3, 2 ]; % [I0, gam, Esh, n] Schottky diode long take
%        [beta, R, ~, CovB] = nlinfit(X, currents', @Schottky_model, p0, 'Options', statset('MaxIter', 1000));
        [beta, R, ~, CovB] = nlinfit(X, currents2fit', @Schottky_model, p0, 'Options', statset('MaxIter', 1000));
        I_fit = Schottky_model(beta, X);
end

beta
%CovB





%{
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
%}
plot3(temperatures, voltages, abs(currents'), 'ob', 'DisplayName', 'Raw Data')
hold on

plot3(temperatures2fit, voltages2fit, abs(currents2fit'), 'xg', 'DisplayName', 'Fitted Data') % Fitted Data
if chooseDiode == 3 || chooseDiode == 4
    plot3(negtemperatures2fit, negvoltages2fit, abs(negcurrents2fit'), 'xg', 'DisplayName', 'Fitted Data') % Fitted Data
end

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


%title(t, strcat("IV-T curve for ", names{chooseDiode}), 'Interpreter', 'latex', 'FontSize', 18)









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
%    q = 1.602176634e-19; % C
    Vt = 13e-3; % V
%    Vt = k * tt / 1.602176634e-19 % V
    Is = params(1) .* tt.^params(2) .* exp(-params(3)./(params(4).*k.*tt)); % A
%    Is = params(1) .* tt.^params(2) .* exp(-params(3)./(k.*tt)); % A
    I = Is .* (exp(vv./(params(4).*Vt)) - 1) * 1e6; % A
%    I = Is .* ( exp( (vv*q)./(params(4).*k.*tt) ) - 1) * 1e6; % A
%    I = swap_I * 1e9; % A
end

function [I] = Schottky_model(params, X) % params = [I0, gam, Esh, n]
    vv = X(1, :);
    tt = X(2, :);

    k = 8.617333e-5; % eV/K
    Vt = 13e-3; % V
%    Is = params(1) .* tt.^2 .* exp(-params(2)./(k.*tt)); % A
%    I = Is .* (exp(vv./( params(3) .*Vt)) - 1) ; % A
    Is = params(1) .* tt.^params(2) .* exp(-params(3)./(k.*tt)); % A
    I = Is .* (exp(vv./(params(4).*Vt)) - 1); % A
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


