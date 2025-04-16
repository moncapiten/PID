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


%temps = 16:1:70;
%{
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
%offsets = [-1, 5]; % V
%}
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




            counter = counter +1;
        end


    end
end
counter * 50 * 2 % 50 voltages, 2 directions (V_go and V_return)






















































































function filename = getFileName(Temp, Vcc, pull, T_direction)
    if T_direction == 1
        T_direction = 'Go';
    elseif T_direction == 2
        T_direction = 'Return';
    end
    filename = sprintf('IV_T%.2f_V%.2f_%.2f_%dPull_T%s', Temp, Vcc(1), Vcc(2), pull, T_direction);
end



