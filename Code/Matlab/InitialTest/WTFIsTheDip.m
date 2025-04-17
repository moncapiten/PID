clearvars
clc
close all
clear all;

chooseDiode = 1; % 1: PN diode, 2: Schottky diode, 3: Zener diode 1, 4: Zener diode 2


names = {'PN diode Short Take', 'Schottky diode', 'Zener diode 1', 'Zener diode 2', 'PN diode Long Take'};
filenames = { 'IV-T_dependence_20250404_100546/', 'IV-T_dependence_20250408_110423/', 'IV-T_dependence_20250408_124236/', 'IV-T_dependence_20250411_094932/', 'IV-T_dependence_20250414_092701/' };
dataPosition = strcat('../../../Data/', filenames{chooseDiode});

ranges = [ [-1, 5], [-1, 5], [-3, 5], [-5, 5], [0, 5] ]; % offsets for each diode
offsets = ranges(chooseDiode, :); % offsets for the chosen diode

tempRanges = [ [16,1,70], [10,1,70], [10,1,70], [10,1,70], [10,0.5,70] ]; % temperature ranges for each diode
temps = tempRanges(chooseDiode, 1):tempRanges(chooseDiode, 2):tempRanges(chooseDiode, 3); % temperatures for the chosen diode


log = true;

n_pulls = 2;
pulls = 0:n_pulls

flag_return_T = true;
flag_return_V = true;

voltages = [];
temperatures = [];
currents = [];


if flag_return_T ar_T_max = 2; else ar_T_max = 1; end

i0 = 0;
f = 0;
red_colors = ["#FFFF00", "#FFCC00", "#FF6600", "#CC0000", "#660000"];
blue_colors = ["#00FFFF", "#00CCFF", "#0099FF", "#0033CC", "#000080"];
legends = [];
legendsNames = {};




t = tiledlayout(2, 2, "TileSpacing", "Tight", "Padding", "Compact");
t1 = nexttile(t);
hold on

minI = 1e5;
minT = -1;
minIsGo = [];
minIsRe = [];
Ts = [];

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


            if(f == 0)
                f = 1;
                i0 = T;
            end

            if(pull == 0 )
                if(T_direction == 1)
                    if log
                        h = plot(Vd_go, abs(Id_go)*1e6, 'o-', 'Color', red_colors( 1+mod(T-i0+1, 5) ), 'Parent', t1);
                    else
                        h = plot(Vd_go, (Id_go)*1e6, 'o-', 'Color', red_colors( 1+mod(T-i0+1, 5) ), 'Parent', t1);
                    end
                    legendsNames{end+1} = strcat( 'T: ', num2str(T), ' C, pull: ', num2str(pull), ' (return) - T Go');
                else
                    if log
                        h = plot(Vd_go, abs(Id_go)*1e6, 'o-', 'Color', red_colors( 1+mod(T-i0+1, 5) ), 'Parent', t1);
                    else
                        h = plot(Vd_go, (Id_go)*1e6, 'o-', 'Color', red_colors( 1+mod(T-i0+1, 5) ), 'Parent', t1);
                    end
                    legendsNames{end+1} = strcat( 'T: ', num2str(T), ' C, pull: ', num2str(pull), ' (return) - T Return');
                end
                legends(end+1) = h;
            else
                if(T_direction == 1)
                    if log
                        plot(Vd_go, abs(Id_go)*1e6, 'o-', 'Color', red_colors( 1+mod(T-i0+1, 5) ), 'Parent', t1);
                        plot(Vd_return, abs(Id_return)*1e6, 'o-', 'Color', red_colors( 1+mod(T-i0+1, 5) ), 'Parent', t1);
                    else
                        plot(Vd_go, (Id_go)*1e6, 'o-', 'Color', red_colors( 1+mod(T-i0+1, 5) ), 'Parent', t1);
                        plot(Vd_return, (Id_return)*1e6, 'o-', 'Color', red_colors( 1+mod(T-i0+1, 5) ), 'Parent', t1);
                    end
                else
                    if log
                        plot(Vd_go, abs(Id_go)*1e6, 'o-', 'Color', red_colors( 1+mod(T-i0+1, 5) ), 'Parent', t1);
                        plot(Vd_return, abs(Id_return)*1e6, 'o-', 'Color', red_colors( 1+mod(T-i0+1, 5) ), 'Parent', t1);
                    else
                        plot(Vd_go, (Id_go)*1e6, 'o-', 'Color', red_colors( 1+mod(T-i0+1, 5) ), 'Parent', t1);
                        plot(Vd_return, (Id_return)*1e6, 'o-', 'Color', red_colors( 1+mod(T-i0+1, 5) ), 'Parent', t1);
                    end
                end
            end

                
            
            if min(abs(Id_go)) < minI || min(abs(Id_return)) < minI
                minI = min([min(abs(Id_go)), min(abs(Id_return))]);
                minT = T;
            end

            minIsGo = [minIsGo, min(abs(Id_go))];
            l = find(abs(Id_go) == min(abs(Id_go)), 1, 'first');
            minIsRe = [minIsRe, min(abs(Id_return))];
            m = find(abs(Id_return) == min(abs(Id_return)), 1, 'first');
            Ts = [Ts; Td_go(l), Td_return(m)];

%            minIsGo = [minIsGo, min(abs(Id_go))];
%            minIsRe = [minIsRe, min(abs(Id_return))];
%            Ts = [Ts; T, T];

%            end
            


            % add the data to the accumulation arrays
            

        

%            voltages = [voltages ; Vd_go ; Vd_return];
%            temperatures = [temperatures ; Td_go ; Td_return];
%            currents = [currents ; Id_go ; Id_return];






            counter = counter +1;
        end


    end
end
counter * 50 * 2 % 50 voltages, 2 directions (V_go and V_return)


%minI
%minT
grid on;
grid minor;
xlabel(t1, 'Voltage [V]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel(t1, 'Current [$ \mathrm{ \mu A} $]', 'Interpreter', 'latex', 'FontSize', 14)
%title('IV characteristics', 'Interpreter', 'latex', 'FontSize', 18)
if log set(t1, 'YScale', 'log'); end
xlim(t1, [0, 0.3])
%legend( legends, legendsNames, 'Interpreter', 'latex', 'FontSize', 10, 'Location', 'nw')
title(t1, strcat("IV families in T - ", names{chooseDiode}), 'Interpreter', 'latex', 'FontSize', 16)

t2 = nexttile(t, [2, 1]);
t2Legends = [];
t2LegendsNames = {};

hold on
h = plot(Ts(:, 1), minIsGo*1e6, 'o-', 'Color', red_colors(5));
plot(Ts(:, 2), minIsRe*1e6, 'o-', 'Color', red_colors(5));
t2Legends(end+1) = h;
t2LegendsNames{end+1} = names{chooseDiode};
grid on;
grid minor;
xlabel(t2, 'Temperature [K]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel(t2, 'Minimum Current [$ \mathrm{ \mu A} $]', 'Interpreter', 'latex', 'FontSize', 14)
set(t2, 'YScale', 'log')
xlim(t2, [30, 70])











































chooseDiode = 5; % 1: PN diode, 2: Schottky diode, 3: Zener diode 1, 4: Zener diode 2, 5: PN diode long take
dataPosition = strcat('../../../Data/', filenames{chooseDiode});
f = 0;

switch chooseDiode
    case 1
        offsets = [-1, 5]; % PN diode
        temps = 16:1:70;
    case 2
        offsets = [-1, 5]; % Schottky diode
        temps = 10:1:70;
    case 3
        offsets = [-3, 5]; % Zener diode 1
        temps = 10:1:70;
    case 4
        offsets = [-5, 5]; % Zener diode 2
        temps = 10:1:70;
    case 5
        offsets = [0, 5]; % PN diode long take
        temps = 10:0.5:70;
end


t3 = nexttile(t);
hold on

minI = 1e5;
minT = -1;
minIsGo = [];
minIsRe = [];
Ts = [];

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


            if(f == 0)
                f = 1;
                i0 = T;
            end

            if(pull == 0 )
                if(T_direction == 1)
                    if log
                        h = plot(Vd_go, abs(Id_go)*1e6, 'o-', 'Color', blue_colors( 1+floor(mod(T-i0+1,5)) ), 'Parent', t3);
                    else
                        h = plot(Vd_go, (Id_go)*1e6, 'o-', 'Color', blue_colors( 1+floor(mod(T-i0+1,5)) ), 'Parent', t3);
                    end
                    legendsNames{end+1} = strcat( 'T: ', num2str(T), ' C, pull: ', num2str(pull), ' (return) - T Go');
                else
                    if log
                            h = plot(Vd_go, abs(Id_go)*1e6, 'o-', 'Color', blue_colors( 1+floor(mod(T-i0+1,5)) ), 'Parent', t3);
                    else
                            h = plot(Vd_go, (Id_go)*1e6, 'o-', 'Color', blue_colors( 1+floor(mod(T-i0+1,5)) ), 'Parent', t3);
                    end
                    legendsNames{end+1} = strcat( 'T: ', num2str(T), ' C, pull: ', num2str(pull), ' (return) - T Return');
                end
                legends(end+1) = h;
            else
                if(T_direction == 1)
                    if log
                        plot(Vd_go, abs(Id_go)*1e6, 'o-', 'Color', blue_colors( 1+floor(mod(T-i0+1,5)) ), 'Parent', t3);
                        plot(Vd_return, abs(Id_return)*1e6, 'o-', 'Color', blue_colors( 1+floor(mod(T-i0+1,5)) ), 'Parent', t3);
                    else
                        plot(Vd_go, (Id_go)*1e6, 'o-', 'Color', blue_colors( 1+floor(mod(T-i0+1,5)) ), 'Parent', t3);
                        plot(Vd_return, (Id_return)*1e6, 'o-', 'Color', blue_colors( 1+floor(mod(T-i0+1,5)) ), 'Parent', t3);
                    end
                else
                    if log
                            plot(Vd_go, abs(Id_go)*1e6, 'o-', 'Color', blue_colors( 1+floor(mod(T-i0+1,5)) ), 'Parent', t3);
                            plot(Vd_return, abs(Id_return)*1e6, 'o-', 'Color', blue_colors( 1+floor(mod(T-i0+1,5)) ), 'Parent', t3);
                    else
                            plot(Vd_go, (Id_go)*1e6, 'o-', 'Color', blue_colors( 1+floor(mod(T-i0+1,5)) ), 'Parent', t3);
                            plot(Vd_return, (Id_return)*1e6, 'o-', 'Color', blue_colors( 1+floor(mod(T-i0+1,5)) ), 'Parent', t3);
                    end
                end
            end

                
            
                if min(abs(Id_go)) < minI || min(abs(Id_return)) < minI
                    minI = min([min(abs(Id_go)), min(abs(Id_return))]);
                    minT = T;
                end

                minIsGo = [minIsGo, min(abs(Id_go))];
                l = find(abs(Id_go) == min(abs(Id_go)), 1, 'first');
                minIsRe = [minIsRe, min(abs(Id_return))];
                m = find(abs(Id_return) == min(abs(Id_return)), 1, 'first');
                Ts = [Ts; Td_go(l), Td_return(m)];



            counter = counter +1;
        end


    end
end
counter * 100 * 2 % 50 voltages, 2 directions (V_go and V_return)



grid on;
grid minor;
xlabel(t3, 'Voltage [V]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel(t3, 'Current [$ \mathrm{ \mu A} $]', 'Interpreter', 'latex', 'FontSize', 14)
%title('IV characteristics', 'Interpreter', 'latex', 'FontSize', 18)
if log set(t3, 'YScale', 'log'); end
xlim(t3, [0, 0.3])


h = plot(t2, Ts(:, 1), minIsGo*1e6, 'o-', 'Color', blue_colors(5));
plot(Ts(:, 2), minIsRe*1e6, 'o-', 'Color', blue_colors(5), 'Parent', t2);
t2Legends(end+1) = h;
t2LegendsNames{end+1} = names{chooseDiode};

%legend(t2, 'PN diode', 's', 'b', 'PN diode long take', 'Interpreter', 'latex', 'FontSize', 10, 'Location', 'nw')
legend(t2, t2Legends, t2LegendsNames, 'Interpreter', 'latex', 'FontSize', 10, 'Location', 'nw')









title(t, strcat("Current Dip Analysis in PN Diode"), 'Interpreter', 'latex', 'FontSize', 18)

title(t3, strcat("IV families in T - ", names{chooseDiode}), 'Interpreter', 'latex', 'FontSize', 16)
title(t2, strcat("Current minimums on Temperature - "), 'Interpreter', 'latex', 'FontSize', 16)










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


