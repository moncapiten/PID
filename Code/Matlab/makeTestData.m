dataPosition = '../../../Data/FakeData/';
filename = 'IV_T10.00_V-5.00_5.00_0Pull';

function filename = getFilename(Temp, Vcc, pull)
    filename = sprintf('IV_T%.2f_V%s_%s_%sPull', Temp, Vcc(1), Vcc(2), pull);
end
%raw_data = readmatrix(strcat(dataPosition, filename, '.txt'));


function Id = shockley(Vd, Vth, T, n)
    Is = tempDependance(T, n); % A/K
    Id = Is * (exp(Vd / (n * Vth)) - 1);
end

function Is = tempDependance(Td, n)
    k = 1.380649e-23; % J/K
    Eg = 1.117; % eV
    alpha = 2.5e-3; % A/K
    Is = alpha * Td * Td * exp(- Eg / n * k * Td);
end

tempRange = [10, 70]; % °C
ntemps = 70;
temperatures = linspace(tempRange(1), tempRange(2), ntemps); % °C


Vcc = [-5, 5]; % V
nvolts = 50;
voltages = linspace(Vcc(1), Vcc(2), nvolts); % V





%R = 1e4 % Ohm

if flag_return_T
    ar_T_max = 2
else
    ar_T_max = 1
end

if flag_return_V
    ar_V_max = 2
else
    ar_V_max = 1
end


for ar_T = 1:ar_T_max
    if ar_T == 2
        temperatures.reverse()
    end
    for t = temperatures

        for ar_V = 1:ar_V_max
            data_T = [];
            data_V = [];
            data_I = [];
            if ar_V == 2
                voltages.reverse()
            end
            for v = voltages
                I = shockley(v, 0.7, t + 273.15, 2); % A

                if ar_V ==1
                    data_T_go = [data_T, t+random('Normal', 0, 0.1)];
                    data_V_go = [data_V, v+random('Normal', 0, 0.1)];
                    data_I_go = shockley(data_V, 0.7, data_T + 273.15, 2) + random('Normal', 0, 0.1); % A

                    plot3(data_T_go, data_V_go, data_I_go, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'black')
                else
                    data_T_return = [data_T, t+random('Normal', 0, 0.1)];
                    data_V_return = [data_V, v+random('Normal', 0, 0.1)];
                    data_I_return = shockley(data_V, 0.7, data_T + 273.15, 2) + random('Normal', 0, 0.1); % A

                    plot3(data_T_return, data_V_return, data_I_return, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'black')
                end


                

                if t == temperatures(1) && v == voltages(1)
                    hold on
                end


            end

            matrix_to_save = [ data_T', data_V', data_I'];
        end

    end
end




