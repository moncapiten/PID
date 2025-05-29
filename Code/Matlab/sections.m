clearvars
clc
close all
clear all;

chooseDiode = 1; % 1: PN diode, 2: Schottky diode, 3: Zener diode 1, 4: Zener diode 2, 5: PN diode long take


names = {'PN diode short take', 'Schottky diode', 'Zener diode 1', 'Zener diode 2', 'PN diode long take', 'Schottky diode long take'};
filenames = { 'IV-T_dependence_20250404_100546/', 'IV-T_dependence_20250408_110423/', 'IV-T_dependence_20250408_124236/', 'IV-T_dependence_20250411_094932/', 'IV-T_dependence_20250414_092701/', 'IV-T_dependence_20250415_113544/' };
dataPosition = strcat('../../Data/', filenames{chooseDiode});

ranges = [ [-1, 5]; [-1, 5]; [-3, 5]; [-5, 5]; [0, 5]; [-1, 5] ]; % offsets for each diode
offsets = ranges(chooseDiode, :); % offsets for the chosen diode

tempRanges = [ [16,1,70]; [10,1,70]; [10,1,70]; [10,1,70]; [10,0.5,70]; [10, 0.5, 70]]; % temperature ranges for each diode
temps = tempRanges(chooseDiode, 1):tempRanges(chooseDiode, 2):tempRanges(chooseDiode, 3); % temperatures for the chosen diode

n_pulls = 2;
pulls = 0:n_pulls

flag_return_T = true;
flag_return_V = true;
log = true;


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

Vrange = [-1.2, 0.7];
Vstep = 0.2;

Trange = [0, 80];
Tstep = 10;

Irange = [-25, 500];
Istep = 50;

Ilogrange = [1e-4, 1e3];



t = tiledlayout(2, 2, "TileSpacing", "Tight", "Padding", "Compact");
title(t, sprintf('Projections of IV-T dependence of %s', names{chooseDiode}), 'Interpreter', 'latex', 'FontSize', 18)

t1 = nexttile(t);
hold on
grid on
grid minor
ylabel('$ I_d [ \mathrm{ \mu A }] $', 'Interpreter', 'latex', 'FontSize', 14)
title('Projections of IV-T as curves - constant V', 'Interpreter', 'latex', 'FontSize', 16)

t2 = nexttile(t);
hold on
grid on
grid minor
title('Projections of IV-T as curves - constant T', 'Interpreter', 'latex', 'FontSize', 16)

t3 = nexttile(t);
hold on
grid on
grid minor
set(t3, 'YScale', 'log')
xlabel('$ T_d [ \mathrm{C} ] $ ', 'Interpreter', 'latex', 'FontSize', 14)
ylabel(' $ \mid I_d \mid  [ \mathrm{ \mu A} ] $ ', 'Interpreter', 'latex', 'FontSize', 14)

t4 = nexttile(t);
hold on
grid on
grid minor
set(t4, 'YScale', 'log')
xlabel(' $ V_d [V] $ ', 'Interpreter', 'latex', 'FontSize', 14)


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

%            plot(t1, Td_go, Id_go*1e6, 'o-', 'Color', '#ff0000');
%            plot(t1, Td_return, Id_return*1e6, 'o-', 'Color', '#0027bd');

%            plot(t2, Vd_go, Id_go*1e6, 'o-', 'Color', '#ff0000');
%            plot(t2, Vd_return, Id_return*1e6, 'o-', 'Color', '#0027bd');

            plot(t3, Td_go, abs(Id_go)*1e6, 'o-', 'Color', '#ff0000');
            plot(t3, Td_return, abs(Id_return)*1e6, 'o-', 'Color', '#0027bd');

            plot(t4, Vd_go, abs(Id_go)*1e6, 'o-', 'Color', '#ff0000');
            plot(t4, Vd_return, abs(Id_return)*1e6, 'o-', 'Color', '#0027bd');

            errorbar(t1, Td_go, Id_go*1e6, Err_Id_go*1e6, 'o-', 'Color', '#ff0000');
            errorbar(t1, Td_return, Id_return*1e6, Err_Id_return*1e6, 'o-', 'Color', '#0027bd');

            errorbar(t2, Vd_go, Id_go*1e6, Err_Id_go*1e6, 'o-', 'Color', '#ff0000');
            errorbar(t2, Vd_return, Id_return*1e6, Err_Id_return*1e6, 'o-', 'Color', '#0027bd');

%            errorbar(t3, Td_go, abs(Id_go)*1e6, 'o-', 'Color', '#ff0000');
%            errorbar(t3, Td_return, abs(Id_return)*1e6, 'o-', 'Color', '#0027bd');

%            errorbar(t4, Vd_go, abs(Id_go)*1e6, 'o-', 'Color', '#ff0000');
%            errorbar(t4, Vd_return, abs(Id_return)*1e6, 'o-', 'Color', '#0027bd');



            counter = counter +1;
        end


    end
end
counter * 50 * 2 % 50 voltages, 2 directions (V_go and V_return)









%ylim(t1, [-25, 500])
%xlim(t1, [10, 80])

%ylim(t2, [-25, 500])
%xlim(t2, [-1.2, 0.6])

%ylim(t3, [1e-4, 1e3])
%xlim(t3, [10, 80])

%ylim(t4, [1e-4, 1e3])
%xlim(t4, [-1.2, 0.6])

xlim(t1, Trange)
ylim(t1, Irange)

xlim(t2, Vrange)
ylim(t2, Irange)

xlim(t3, Trange)
ylim(t3, Ilogrange)

xlim(t4, Vrange)
ylim(t4, Ilogrange)




%pause(20)

xticklabels(t1, '')
xticks(t1, Trange(1):Tstep:Trange(2))
yticks(t2, Irange(1):Istep:Irange(2))
%xticks(t1, 10:10:80)
%yticks(t1, 0:50:500)

xticklabels(t2, '')
yticklabels(t2, '')
xticks(t2, Vrange(1):Vstep:Vrange(2))
yticks(t2, 0:50:500)
%xticks(t2, -1.2:0.2:0.6)
%yticks(t2, 0:50:500)

xticks(t3, Trange(1):Tstep:Trange(2))
yticks(t3, logspace(log10(Ilogrange(1)), log10(Ilogrange(2)), 8))
%xticks(t3, 10:10:80)
%yticks(t3, logspace(-4, 3, 8))

yticklabels(t4, '')
xticks(t4, Vrange(1):Vstep:Vrange(2))
yticks(t4, logspace(log10(Ilogrange(1)), log10(Ilogrange(2)), 8))
%xticks(t4, -1.2:0.2:0.6)
%yticks(t4, logspace(-4, 3, 8))





a = yticks(t1);
yticks(t2, a)
yticklabels(t2, "")






















































function filename = getFileName(Temp, Vcc, pull, T_direction)
    if T_direction == 1
        T_direction = 'Go';
    elseif T_direction == 2
        T_direction = 'Return';
    end
    filename = sprintf('IV_T%.2f_V%.2f_%.2f_%dPull_T%s', Temp, Vcc(1), Vcc(2), pull, T_direction);
end



