%%
clc; 
% clear; 
close all;

format long

set(groot,'defaultLineLineWidth',2)  %sets graph line width as 2
set(groot,'defaultAxesFontSize',24)  %sets graph axes font size as 18
set(groot,'defaulttextfontsize',24)  %sets graph text font size as 18
set(groot,'defaultLineMarkerSize',8) %sets line marker size as 8
set(groot,'defaultAxesXGrid','on')   %sets X axis grid on 
set(groot,'defaultAxesYGrid','on')   %sets Y axis grid on
set(groot,'DefaultAxesBox', 'on')   %sets Axes boxes on

picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.65; % feel free to play with this ratio
%%

data = load("E:\Y2\computing CW\Final CPP CW\Output_4.txt");

%%
% Extract time and data for vessels
% Assuming each line has: t y1 v1 h1 y2 v2 h2 ...
time = data(:, 1); % Time column
num_vessels = (size(data, 2) - 1) / 3; % Determine the number of vessels

% Initialize storage for vessel data
positions = zeros(length(time), num_vessels);
velocities = zeros(length(time), num_vessels);
water_levels = zeros(length(time), num_vessels);

% Extract positions, velocities, and water levels for each vessel
for i = 1:num_vessels
    positions(:, i) = data(:, 3 * i - 1);
    velocities(:, i) = data(:, 3 * i);
    water_levels(:, i) = data(:, 3 * i + 1);
end

% Part (a): Plot displacement of vessel base (y) over time for each vessel
% plot Analytical solution
Fig_name_1 = "Displacement_4_FE";

Line_Style = [":","-.","-"];
Fig_name = figure;
%set(Analytical_solution,"WindowState","maximized");
set(findall(Fig_name,'-property','FontSize'),'FontSize',24);
set(findall(Fig_name,'-property','Interpreter'),'Interpreter','latex') 
set(findall(Fig_name,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(Fig_name,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(Fig_name,'Position');
set(Fig_name,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

hold on;
for i = 1:num_vessels
    plot(time, positions(:, i),'DisplayName', ['Vessel ' num2str(i)],LineStyle = Line_Style(i));
end
xlabel('Time (s)');
ylabel('Displacement (m)');
% title('Displacement of Vessel Base Over Time');
grid on;
legend(Location= "southeast");
hold off;

saveas(Fig_name,'E:\Y2\computing CW\Final CPP CW\Figures_Part_a'+Fig_name_1+'.svg');


%%

% %% Plot Velocities
% 
% figure;
% hold on;
% for i = 1:num_vessels
%     plot(time, velocities(:, i), 'DisplayName', ['Vessel ' num2str(i)]);
% end
% xlabel('Time (s)');
% ylabel('Velocity (m/s)');
% title('Velocity of Vessel Base Over Time');
% grid on;
% legend;
% hold off;
% 
% %% Plot Water Level
% 
% figure;
% 
% 
% hold on;
% for i = 1:num_vessels
%     plot(time, water_levels(:, i), 'DisplayName', ['Vessel ' num2str(i)]);
% end
% xlabel('Time (s)');
% ylabel('Water Level (m)');
% title('Water Level of Vessel Base Over Time');
% grid on;
% legend;
% % hold off;


%% Part (b): Maximum timestep for convergence
% For simplicity, this part would involve analyzing output for different dt
% and determining convergence. Here, it is assumed convergence has been
% ensured via proper selection of dt in C++.

% Table for final values
final_displacements = positions(end, :);
final_velocities = velocities(end, :);
final_water_levels = water_levels(end, :);

%Display the results in MATLAB's command window
disp('Final values for each vessel:');
disp('Vessel | Final Displacement (m) | Final Velocity (m/s) | Final Water Level (m)');
for i = 1:num_vessels
    fprintf('%6d | %22.6f | %20.6f | %20.6f\n', i, final_displacements(i), final_velocities(i), final_water_levels(i));
end


%% Part B (ii)

%mass graph

m1 = [0.0005, 0.001, 0.002, 0.003, 0.004, 0.005, 0.05:0.05:1, 1.25, 1.50, 1.75 ,2.00 ];

data_m = load("E:\Y2\computing CW\Final CPP CW\Figures_Part_b\m\m1.txt");

m1_freq = calculateAverageFrequency(data_m);

mvsfreq = figure;

set(findall(mvsfreq,'-property','FontSize'),'FontSize',24);
set(findall(mvsfreq,'-property','Interpreter'),'Interpreter','latex') 
set(findall(mvsfreq,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(mvsfreq,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(mvsfreq,'Position');
set(mvsfreq,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

hold on;
plot(m1,m1_freq,"k-x")


xlabel('Mass (kg)');
ylabel('Frequency (Hz)');

grid on;
% legend('Location','best');
hold off;

saveas(mvsfreq,'E:\Y2\computing CW\Final CPP CW\Figures_Part_b\m\mvsfreq.svg');


%% A vs  

A1 = [0.1, 0.5, 1:1:10, 15, 20:10:70];

data_A = load("E:\Y2\computing CW\Final CPP CW\Figures_Part_b\A\A1.txt");

A1_freq = calculateAverageFrequency(data_A);

Avsfreq = figure;

set(findall(Avsfreq,'-property','FontSize'),'FontSize',24);
set(findall(Avsfreq,'-property','Interpreter'),'Interpreter','latex') 
set(findall(Avsfreq,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(Avsfreq,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(Avsfreq,'Position');
set(Avsfreq,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

hold on;

plot(A1,A1_freq,"r-x")


xlabel('Area (m^2)');
ylabel('Frequency (Hz)');

grid on;
% legend('Location','best');
hold off;

saveas(Avsfreq,'E:\Y2\computing CW\Final CPP CW\Figures_Part_b\A\Avsfreq.svg');


%% K vs  

K1 = [0:5:35,40:20:400];
data_K = load("E:\Y2\computing CW\Final CPP CW\Figures_Part_b\K\K1.txt");

K1_freq = calculateAverageFrequency(data_K);

Kvsfreq = figure;

set(findall(Kvsfreq,'-property','FontSize'),'FontSize',24);
set(findall(Kvsfreq,'-property','Interpreter'),'Interpreter','latex') 
set(findall(Kvsfreq,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(Kvsfreq,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(Kvsfreq,'Position');
set(Kvsfreq,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

hold on;

plot(K1,K1_freq,"b-x")


xlabel('K');
ylabel('Frequency (Hz)');

grid on;
% legend('Location','best');
hold off;

saveas(Kvsfreq,'E:\Y2\computing CW\Final CPP CW\Figures_Part_b\K\Kvsfreq.svg');


%% H_0 vs  


h1 = [0:0.01:0.04, 0.05:0.02:0.21, 0.22:0.01:0.3 ];
data_h = load("E:\Y2\computing CW\Final CPP CW\Figures_Part_b\h_0\h_0.txt");

h1_freq = calculateAverageFrequency(data_h);

hvsfreq = figure;

set(findall(hvsfreq,'-property','FontSize'),'FontSize',24);
set(findall(hvsfreq,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hvsfreq,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hvsfreq,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hvsfreq,'Position');
set(hvsfreq,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

hold on;
plot(h1,h1_freq,'g-x')


xlabel('h_0 (m)');
ylabel('Frequency (Hz)');

grid on;
% legend('Location','best');
hold off;

saveas(hvsfreq,'E:\Y2\computing CW\Final CPP CW\Figures_Part_b\h_0\h_0vsfreq.svg');



%% y_0 vs  


y_01 = -1:0.2:1;
data_y_0 = load("E:\Y2\computing CW\Final CPP CW\Figures_Part_b\y_0\y_0.txt");

y_01_freq = calculateAverageFrequency(data_y_0);

y_0vsfreq = figure;

set(findall(y_0vsfreq,'-property','FontSize'),'FontSize',24);
set(findall(y_0vsfreq,'-property','Interpreter'),'Interpreter','latex') 
set(findall(y_0vsfreq,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(y_0vsfreq,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(y_0vsfreq,'Position');
set(y_0vsfreq,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

hold on;

plot(y_01,y_01_freq,Marker= "x", Color="#D95319")


xlabel('y_0 (m)');
ylabel('Frequency (Hz)');

grid on;
% legend('Location','best');
hold off;

saveas(y_0vsfreq,'E:\Y2\computing CW\Final CPP CW\Figures_Part_b\y_0\y_0vsfreq.svg');


%% y_dot_0 vs  


y_dot_01 = -1:0.2:1;
data_y_dot_0 = load("E:\Y2\computing CW\Final CPP CW\Figures_Part_b\y_dot_0\y_dot_0.txt");

y_01_freq = calculateAverageFrequency(data_y_dot_0);

y_dot_0vsfreq = figure;

set(findall(y_dot_0vsfreq,'-property','FontSize'),'FontSize',24);
set(findall(y_dot_0vsfreq,'-property','Interpreter'),'Interpreter','latex') 
set(findall(y_dot_0vsfreq,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(y_dot_0vsfreq,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(y_dot_0vsfreq,'Position');
set(y_dot_0vsfreq,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

hold on;

plot(y_dot_01,y_01_freq,marker = "x", Color= "#EDB120")


xlabel('Initial Velocity, v_0 (m)');
ylabel('Frequency (Hz)');

grid on;
% legend('Location','best');
hold off;

saveas(y_dot_0vsfreq,'E:\Y2\computing CW\Final CPP CW\Figures_Part_b\y_dot_0\y_dot_0vsfreq.svg');








%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [dominant_frequency] = calculateAverageFrequency(data)
    % Load the data from the file


    % Extract the time column (first column)
    time = data(:, 1);

    % Determine the number of vessels from the data
    num_vessels = (size(data, 2) - 1) / 3; % Assuming y, v, h per vessel

    % Initialize a variable to store the sum of dominant frequencies
    dominant_frequency = zeros(num_vessels,1);
    % Loop through each vessel to calculate its dominant frequency
    for i = 1:num_vessels
        % Extract displacement data for the vessel
        displacement = data(:, 3 * i - 1);

        % Perform Fourier Transform on the displacement data
        Y = fft(displacement);
        L = length(time);

        % Compute the two-sided spectrum and then the single-sided spectrum
        P2 = abs(Y / L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2 * P1(2:end-1);

        % Define the frequency domain
        f = (0:(L/2)) / max(time);

        % Find the frequency corresponding to the maximum amplitude
        [~, max_idx] = max(P1);
        dominant_frequency(i) = f(max_idx);

        % fprintf('The frequency of vessel %.1f is %.4f Hz\n',i ,dominant_frequency);
    end

    % Calculate the average frequency


    % Display the result
end