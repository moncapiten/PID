clearvars
clc
close all
clear all;

chooseDiode = 1; % 1: PN diode, 2: Schottky diode, 3: Zener diode 1, 4: Zener diode 2


names = {'PN diode', 'Schottky diode', 'Zener diode 1', 'Zener diode 2'};
filenames = { 'IV-T_dependence_20250404_100546/', 'IV-T_dependence_20250408_110423/', 'IV-T_dependence_20250408_124236/', 'IV-T_dependence_20250411_094932/' };
dataPosition = strcat('../../../Data/', filenames{chooseDiode});
%dataPosition = '../../../Data/IV-T_dependence_20250404_100546/';  % PN diode
%dataPosition = '../../../Data/IV-T_dependence_20250408_110423/'; %Schottky diode
%dataPosition = '../../../Data/IV-T_dependence_20250408_124236/'; % Zener diode 1
%dataPosition = '../../../Data/IV-T_dependence_20250411_094932/'; % Zener diode 2

temps = 16:1:70;

switch chooseDiode
    case 1
        offsets = [-1, 5]; % PN diode
    case 2
        offsets = [-1, 5]; % Schottky diode
    case 3
        offsets = [-3, 5]; % Zener diode 1
    case 4
        offsets = [-5, 5]; % Zener diode 2
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

i0 = 0;
f = 0;
red_colors = ["#FFFF00", "#FFCC00", "#FF6600", "#CC0000", "#660000"];
blue_colors = ["#00FFFF", "#00CCFF", "#0099FF", "#0033CC", "#000080"];
legends = [];
names = {};


t = tiledlayout(1, 2, "TileSpacing", "Tight", "Padding", "Compact");
t1 = nexttile(t);
%t2 = nexttile(t);
hold on

minI = 1e5;
minT = -1;


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



            Vd_go = Vd_go(10:end);
            Id_go = Id_go(10:end);
            Vd_return = Vd_return(10:end);
            Id_return = Id_return(10:end);



            if( T >= 50 && T <= 54 )
%            if( T >= 51 && T <= 53 )
%            if( T >= 52 && T <= 52 )
                if(f == 0)
                    f = 1;
                    i0 = T;
                end

                if(pull == 0 )
                    if(T_direction == 1)
                        h = plot(Vd_go, abs(Id_go)*1e6, '-', 'Color', red_colors(T-i0+1), 'Parent', t1);
                        names{end+1} = strcat( 'T: ', num2str(T), ' C, pull: ', num2str(pull), ' (return) - T Go');
                    else
                        h = plot(Vd_go, abs(Id_go)*1e6, '-', 'Color', blue_colors(T-i0+1), 'Parent', t1);
                        names{end+1} = strcat( 'T: ', num2str(T), ' C, pull: ', num2str(pull), ' (return) - T Return');
                    end
                    legends(end+1) = h;
                else
                    if(T_direction == 1)
                        plot(Vd_go, (Id_go)*1e6, '-', 'Color', red_colors(T-i0+1), 'Parent', t1);
                        plot(Vd_return, (Id_return)*1e6, '-', 'Color', red_colors(T-i0+1), 'Parent', t1);
                    else
                        plot(Vd_go, (Id_go)*1e6, '-', 'Color', blue_colors(T-i0+1), 'Parent', t1);
                        plot(Vd_return, (Id_return)*1e6, '-', 'Color', blue_colors(T-i0+1), 'Parent', t1);
                    end
                end

                
            
                if min(abs(Id_go)) < minI || min(abs(Id_return)) < minI
                    minI = min([min(abs(Id_go)), min(abs(Id_return))]);
                    minT = T;
                end


            end

%            min(abs(Id_go)) < minI || min(abs(Id_return)) < minI
            


            % add the data to the accumulation arrays
            

        

            voltages = [voltages ; Vd_go ; Vd_return];
            temperatures = [temperatures ; Td_go ; Td_return];
            currents = [currents ; Id_go ; Id_return];






            counter = counter +1;
        end


    end
end
counter * 50 * 2 % 50 voltages, 2 directions (V_go and V_return)


minI
minT
grid on;
grid minor;
xlabel(t1, 'Voltage [V]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel(t1, 'Current [$ \mathrm{ \mu A} $]', 'Interpreter', 'latex', 'FontSize', 14)
%title('IV characteristics', 'Interpreter', 'latex', 'FontSize', 18)
%set(t1, 'YScale', 'log')
%legend( legends, names, 'Interpreter', 'latex', 'FontSize', 10, 'Location', 'nw')





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


