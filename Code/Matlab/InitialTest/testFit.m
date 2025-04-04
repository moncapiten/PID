dataPosition = '../../../Data/IV-T_dependence_20250403_101055/';
filename = 'IV_T10.00_V-5.00_5.00_0Pull';
raw_data = readmatrix(strcat(dataPosition, filename, '.txt'));

% Td_go[�C]	Td_return[�C]	Vcc	Vd_go[V]	Vd_return[V]	ErrVd_go[V]	ErrVd_return[V]	Vr_go[V]	Vr_return[V]	ErrVr_go[V]	ErrVr_return[V]	Id_go[A]	Id_return[A]	ErrId_go[A]	ErrId_return[A]

Td_go = raw_data(:, 1);
Td_return = raw_data(:, 2);
Vcc = raw_data(:, 3);
Vd_go = raw_data(:, 4);
Vd_return = raw_data(:, 5);
ErrVd_go = raw_data(:, 6);
ErrVd_return = raw_data(:, 7);
Id_go = raw_data(:, 12);
Id_return = raw_data(:, 13);


plot3(Td_go, Vd_go, Id_go, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'black')
hold on
plot3(Td_return, Vd_return, Id_return, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'black')
hold off
xlabel('Td [°C]')
ylabel('Vd [V]')
zlabel('Id [A]')
title('3D Scatter Plot of Vd and Id vs Td')
grid on
legend('Td go', 'Td return')
view(3)
axis equal
set(gca, 'XScale', 'linear', 'YScale', 'linear', 'ZScale', 'linear')





