
clear all; close all; clc;


% Electrical Parameters
V_line = 400e3;                      % Line voltage (V) - 400 kV
V_phase = V_line / sqrt(3);          % Phase voltage (V)
I_phase = 1000;                      % Phase current (A)
frequency = 50;                      % Frequency (Hz)
omega = 2 * pi * frequency;          % Angular frequency (rad/s)\

% Physical Constants
eps_0 = 8.854e-12;                   % Permittivity of free space (F/m)
mu_0 = 4 * pi * 1e-7;                % Permeability of free space (H/m)

% Tower and Conductor Geometry (400 kV Double Circuit Tower)
tower_height = 50;                   % Tower height (m)
tower_base_width = 12;               % Tower base width (m)
tower_top_width = 3;                 % Tower top width (m)

% Conductor positions (x, y coordinates in meters)
h_conductor = 45;                    % Conductor height above ground (m)
h_lower = 40;                        % Lower cross-arm height (m)
h_upper = 45;                        % Upper cross-arm height (m)

% Phase conductor positions (Double circuit - vertical configuration)
spacing_horizontal = 8;              % Horizontal spacing (m)

% Upper cross-arm conductors (Circuit 1)
x_a1 = -spacing_horizontal; y_a1 = h_upper;      % Phase A - Circuit 1
x_b1 = 0;                   y_b1 = h_upper;      % Phase B - Circuit 1
x_c1 = spacing_horizontal;  y_c1 = h_upper;      % Phase C - Circuit 1

% Lower cross-arm conductors (Circuit 2)
x_a2 = -spacing_horizontal; y_a2 = h_lower;      % Phase A - Circuit 2
x_b2 = 0;                   y_b2 = h_lower;      % Phase B - Circuit 2
x_c2 = spacing_horizontal;  y_c2 = h_lower;      % Phase C - Circuit 2

r_conductor = 0.015;                 % Conductor radius (m)

% Ground wire position
h_ground_wire = tower_height;        % Ground wire at tower top (m)
x_ground_wire = 0;                   % Center position (m)



% Define calculation grid
x_min = -60; x_max = 60;             % Horizontal range (m)
y_min = 0; y_max = 60;               % Vertical range (m)
grid_spacing = 1;                    % Grid spacing (m)

x_range = x_min:grid_spacing:x_max;
y_range = y_min:grid_spacing:y_max;
[X, Y] = meshgrid(x_range, y_range);

% Initialize field arrays
E_x = zeros(size(X));
E_y = zeros(size(Y));
B_x = zeros(size(X));
B_y = zeros(size(Y));


fprintf('Calculating Electric Field...\n');

% Phase voltages at time t = 0 (instantaneous values)
V_a = V_phase * cos(0);              % Phase A at 0 degrees
V_b = V_phase * cos(-2*pi/3);        % Phase B at -120 degrees
V_c = V_phase * cos(2*pi/3);         % Phase C at 120 degrees

% Simplified charge calculation (line charge density)
% For infinite line charge: C ≈ 2*pi*eps_0 / ln(2h/r)
q_a1 = V_a * 2*pi*eps_0 / log(2*y_a1/r_conductor);
q_b1 = V_b * 2*pi*eps_0 / log(2*y_b1/r_conductor);
q_c1 = V_c * 2*pi*eps_0 / log(2*y_c1/r_conductor);

q_a2 = V_a * 2*pi*eps_0 / log(2*y_a2/r_conductor);
q_b2 = V_b * 2*pi*eps_0 / log(2*y_b2/r_conductor);
q_c2 = V_c * 2*pi*eps_0 / log(2*y_c2/r_conductor);

% Calculate electric field from each conductor and its image
for i = 1:numel(X)
    x_point = X(i);
    y_point = Y(i);
    
    min_dist = 0.5;     % Minimum distance threshold (m)
    
    % Circuit 1 - Phase A
    r_a1 = sqrt((x_point - x_a1)^2 + (y_point - y_a1)^2);
    r_a1_img = sqrt((x_point - x_a1)^2 + (y_point + y_a1)^2);
    if r_a1 > min_dist
        E_x(i) = E_x(i) + q_a1/(2*pi*eps_0) * ((x_point - x_a1)/r_a1^2 - (x_point - x_a1)/r_a1_img^2);
        E_y(i) = E_y(i) + q_a1/(2*pi*eps_0) * ((y_point - y_a1)/r_a1^2 - (y_point + y_a1)/r_a1_img^2);
    end
    
    % Circuit 1 - Phase B
    r_b1 = sqrt((x_point - x_b1)^2 + (y_point - y_b1)^2);
    r_b1_img = sqrt((x_point - x_b1)^2 + (y_point + y_b1)^2);
    if r_b1 > min_dist
        E_x(i) = E_x(i) + q_b1/(2*pi*eps_0) * ((x_point - x_b1)/r_b1^2 - (x_point - x_b1)/r_b1_img^2);
        E_y(i) = E_y(i) + q_b1/(2*pi*eps_0) * ((y_point - y_b1)/r_b1^2 - (y_point + y_b1)/r_b1_img^2);
    end
    
    % Circuit 1 - Phase C
    r_c1 = sqrt((x_point - x_c1)^2 + (y_point - y_c1)^2);
    r_c1_img = sqrt((x_point - x_c1)^2 + (y_point + y_c1)^2);
    if r_c1 > min_dist
        E_x(i) = E_x(i) + q_c1/(2*pi*eps_0) * ((x_point - x_c1)/r_c1^2 - (x_point - x_c1)/r_c1_img^2);
        E_y(i) = E_y(i) + q_c1/(2*pi*eps_0) * ((y_point - y_c1)/r_c1^2 - (y_point + y_c1)/r_c1_img^2);
    end
    
    % Circuit 2 - Phase A
    r_a2 = sqrt((x_point - x_a2)^2 + (y_point - y_a2)^2);
    r_a2_img = sqrt((x_point - x_a2)^2 + (y_point + y_a2)^2);
    if r_a2 > min_dist
        E_x(i) = E_x(i) + q_a2/(2*pi*eps_0) * ((x_point - x_a2)/r_a2^2 - (x_point - x_a2)/r_a2_img^2);
        E_y(i) = E_y(i) + q_a2/(2*pi*eps_0) * ((y_point - y_a2)/r_a2^2 - (y_point + y_a2)/r_a2_img^2);
    end
    
    % Circuit 2 - Phase B
    r_b2 = sqrt((x_point - x_b2)^2 + (y_point - y_b2)^2);
    r_b2_img = sqrt((x_point - x_b2)^2 + (y_point + y_b2)^2);
    if r_b2 > min_dist
        E_x(i) = E_x(i) + q_b2/(2*pi*eps_0) * ((x_point - x_b2)/r_b2^2 - (x_point - x_b2)/r_b2_img^2);
        E_y(i) = E_y(i) + q_b2/(2*pi*eps_0) * ((y_point - y_b2)/r_b2^2 - (y_point + y_b2)/r_b2_img^2);
    end
    
    % Circuit 2 - Phase C
    r_c2 = sqrt((x_point - x_c2)^2 + (y_point - y_c2)^2);
    r_c2_img = sqrt((x_point - x_c2)^2 + (y_point + y_c2)^2);
    if r_c2 > min_dist
        E_x(i) = E_x(i) + q_c2/(2*pi*eps_0) * ((x_point - x_c2)/r_c2^2 - (x_point - x_c2)/r_c2_img^2);
        E_y(i) = E_y(i) + q_c2/(2*pi*eps_0) * ((y_point - y_c2)/r_c2^2 - (y_point + y_c2)/r_c2_img^2);
    end
end

% Calculate electric field magnitude
E_magnitude = sqrt(E_x.^2 + E_y.^2);
E_kV_per_m = E_magnitude / 1000;     % Convert to kV/m

fprintf('Electric Field Calculation Complete.\n');
fprintf('Maximum Electric Field: %.2f kV/m\n', max(E_kV_per_m(:)));


fprintf('Calculating Magnetic Field...\n');

% Phase currents at time t = 0 (instantaneous values)
I_a = I_phase * cos(0);              % Phase A at 0 degrees
I_b = I_phase * cos(-2*pi/3);        % Phase B at -120 degrees
I_c = I_phase * cos(2*pi/3);         % Phase C at 120 degrees

% Calculate magnetic field from each conductor
for i = 1:numel(X)
    x_point = X(i);
    y_point = Y(i);
    
    min_dist = 0.5; % Minimum distance threshold (m)
    
    % Circuit 1 - Phase A
    r_a1 = sqrt((x_point - x_a1)^2 + (y_point - y_a1)^2);
    if r_a1 > min_dist
        B_mag_a1 = mu_0 * abs(I_a) / (2*pi*r_a1);
        B_x(i) = B_x(i) + B_mag_a1 * (-(y_point - y_a1)/r_a1) * sign(I_a);
        B_y(i) = B_y(i) + B_mag_a1 * ((x_point - x_a1)/r_a1) * sign(I_a);
    end
    
    % Circuit 1 - Phase B
    r_b1 = sqrt((x_point - x_b1)^2 + (y_point - y_b1)^2);
    if r_b1 > min_dist
        B_mag_b1 = mu_0 * abs(I_b) / (2*pi*r_b1);
        B_x(i) = B_x(i) + B_mag_b1 * (-(y_point - y_b1)/r_b1) * sign(I_b);
        B_y(i) = B_y(i) + B_mag_b1 * ((x_point - x_b1)/r_b1) * sign(I_b);
    end
    
    % Circuit 1 - Phase C
    r_c1 = sqrt((x_point - x_c1)^2 + (y_point - y_c1)^2);
    if r_c1 > min_dist
        B_mag_c1 = mu_0 * abs(I_c) / (2*pi*r_c1);
        B_x(i) = B_x(i) + B_mag_c1 * (-(y_point - y_c1)/r_c1) * sign(I_c);
        B_y(i) = B_y(i) + B_mag_c1 * ((x_point - x_c1)/r_c1) * sign(I_c);
    end
    
    % Circuit 2 - Phase A
    r_a2 = sqrt((x_point - x_a2)^2 + (y_point - y_a2)^2);
    if r_a2 > min_dist
        B_mag_a2 = mu_0 * abs(I_a) / (2*pi*r_a2);
        B_x(i) = B_x(i) + B_mag_a2 * (-(y_point - y_a2)/r_a2) * sign(I_a);
        B_y(i) = B_y(i) + B_mag_a2 * ((x_point - x_a2)/r_a2) * sign(I_a);
    end
    
    % Circuit 2 - Phase B
    r_b2 = sqrt((x_point - x_b2)^2 + (y_point - y_b2)^2);
    if r_b2 > min_dist
        B_mag_b2 = mu_0 * abs(I_b) / (2*pi*r_b2);
        B_x(i) = B_x(i) + B_mag_b2 * (-(y_point - y_b2)/r_b2) * sign(I_b);
        B_y(i) = B_y(i) + B_mag_b2 * ((x_point - x_b2)/r_b2) * sign(I_b);
    end
    
    % Circuit 2 - Phase C
    r_c2 = sqrt((x_point - x_c2)^2 + (y_point - y_c2)^2);
    if r_c2 > min_dist
        B_mag_c2 = mu_0 * abs(I_c) / (2*pi*r_c2);
        B_x(i) = B_x(i) + B_mag_c2 * (-(y_point - y_c2)/r_c2) * sign(I_c);
        B_y(i) = B_y(i) + B_mag_c2 * ((x_point - x_c2)/r_c2) * sign(I_c);
    end
end

% Calculate magnetic flux density magnitude
B_magnitude = sqrt(B_x.^2 + B_y.^2);
B_microTesla = B_magnitude * 1e6;    % Convert to µT

fprintf('Magnetic Field Calculation Complete.\n');
fprintf('Maximum Magnetic Field: %.2f µT\n', max(B_microTesla(:)));


fprintf('Generating visualizations...\n');

% Create figure with subplots
figure('Position', [100, 100, 1400, 900]);


subplot(2,2,1);
hold on; grid on; box on;

% Plot electric field contour
contourf(X, Y, E_kV_per_m, 20, 'LineStyle', 'none');
colormap(gca, 'jet');
c1 = colorbar;
ylabel(c1, 'Electric Field (kV/m)', 'FontSize', 11);
caxis([0 min(6, max(E_kV_per_m(:)))]);

% Draw tower structure
draw_tower(tower_height, tower_base_width, tower_top_width, 'k', 2);

% Draw conductors
plot([x_a1, x_b1, x_c1], [y_a1, y_b1, y_c1], 'ro', 'MarkerSize', 8, ...
     'MarkerFaceColor', 'g', 'LineWidth', 2);
plot([x_a2, x_b2, x_c2], [y_a2, y_b2, y_c2], 'bo', 'MarkerSize', 8, ...
     'MarkerFaceColor', 'r', 'LineWidth', 2);
plot(x_ground_wire, h_ground_wire, 'ko', 'MarkerSize', 6, ...
     'MarkerFaceColor', 'k', 'LineWidth', 2);

% Add safety limit line
plot([x_min, x_max], [1, 1], 'w--', 'LineWidth', 2);
text(-55, 1.5, 'Human Height (1m)', 'Color', 'w', 'FontSize', 9, 'FontWeight', 'bold');

% Labels and formatting
xlabel('Horizontal Distance (m)', 'FontSize', 12);
ylabel('Height Above Ground (m)', 'FontSize', 12);
title('Electric Field Distribution (400kV Double Circuit)', 'FontSize', 13, 'FontWeight', 'bold');
xlim([x_min, x_max]);
ylim([y_min, y_max]);
legend('', 'Tower', 'Upper Circuit', 'Lower Circuit', 'Ground Wire', 'Location', 'northeast');

subplot(2,2,2);
hold on; grid on; box on;

% Plot magnetic field contour
contourf(X, Y, B_microTesla, 20, 'LineStyle', 'none');
colormap(gca, 'jet');
c2 = colorbar;
ylabel(c2, 'Magnetic Flux Density (µT)', 'FontSize', 11);
caxis([0 min(15, max(B_microTesla(:)))]);

% Draw tower structure
draw_tower(tower_height, tower_base_width, tower_top_width, 'k', 2);

% Draw conductors
plot([x_a1, x_b1, x_c1], [y_a1, y_b1, y_c1], 'ro', 'MarkerSize', 8, ...
     'MarkerFaceColor', 'g', 'LineWidth', 2);
plot([x_a2, x_b2, x_c2], [y_a2, y_b2, y_c2], 'bo', 'MarkerSize', 8, ...
     'MarkerFaceColor', 'r', 'LineWidth', 2);
plot(x_ground_wire, h_ground_wire, 'ko', 'MarkerSize', 6, ...
     'MarkerFaceColor', 'k', 'LineWidth', 2);

% Add safety limit line
plot([x_min, x_max], [1, 1], 'w--', 'LineWidth', 2);
text(-55, 1.5, 'Human Height (1m)', 'Color', 'w', 'FontSize', 9, 'FontWeight', 'bold');

% Labels and formatting
xlabel('Horizontal Distance (m)', 'FontSize', 12);
ylabel('Height Above Ground (m)', 'FontSize', 12);
title('Magnetic Field Distribution (400kV Double Circuit)', 'FontSize', 13, 'FontWeight', 'bold');
xlim([x_min, x_max]);
ylim([y_min, y_max]);
legend('', 'Tower', 'Upper Circuit', 'Lower Circuit', 'Ground Wire', 'Location', 'northeast');


subplot(2,2,3);
hold on; grid on; box on;

% Extract horizontal profile at 15 m height
y_profile = 15; % meters (human exposure height)
[~, y_idx] = min(abs(y_range - y_profile));
E_horizontal = E_kV_per_m(y_idx, :);
B_horizontal = B_microTesla(y_idx, :);

% Plot electric field
yyaxis left
plot(x_range, E_horizontal, 'b-', 'LineWidth', 2);
ylabel('Electric Field (kV/m)', 'FontSize', 11);
ylim([0 max(E_horizontal)*1.2]);
% ICNIRP safety limit
yline(5, 'b--', 'LineWidth', 1.5);
text(30, 5.3, 'ICNIRP Limit (5 kV/m)', 'Color', 'b', 'FontSize', 9);

% Plot magnetic field
yyaxis right
plot(x_range, B_horizontal, 'r-', 'LineWidth', 2);
ylabel('Magnetic Field (µT)', 'FontSize', 11);
ylim([0 max(B_horizontal)*1.2]);
% ICNIRP safety limit
yline(100, 'r--', 'LineWidth', 1.5);
text(30, 105, 'ICNIRP Limit (100 µT)', 'Color', 'g', 'FontSize', 9);

% Tower position indicator
xline(0, 'k:', 'LineWidth', 1.5);
text(0, max(B_horizontal)*0.5, 'Tower', 'FontSize', 9, 'Rotation', 90);

xlabel('Horizontal Distance from Tower Center (m)', 'FontSize', 12);
title('Field Strength at 15 m Above Ground (Horizontal Profile)', 'FontSize', 13, 'FontWeight', 'bold');
legend('Electric Field', 'Magnetic Field', 'Location', 'northeast');


subplot(2,2,4);
hold on; grid on; box on;

% Extract vertical profile at x = 0
x_profile = 0; % meters (tower center)
[~, x_idx] = min(abs(x_range - x_profile));
E_vertical = E_kV_per_m(:, x_idx);
B_vertical = B_microTesla(:, x_idx);

% Plot electric field
yyaxis left
plot(E_vertical, y_range, 'b-', 'LineWidth', 2);
xlabel('Electric Field (kV/m)', 'FontSize', 11);
xlim([0 max(E_vertical)*1.2]);
% ICNIRP safety limit
xline(5, 'b--', 'LineWidth', 1.5);

% Plot magnetic field
yyaxis right
plot(B_vertical, y_range, 'r-', 'LineWidth', 2);
xlabel('Magnetic Field (µT)', 'FontSize', 11);
xlim([0 max(B_vertical)*1.2]);
% ICNIRP safety limit
xline(100, 'r--', 'LineWidth', 1.5);

% Conductor height indicators
yline(h_lower, 'b:', 'LineWidth', 1);
text(max(B_vertical)*0.3, h_lower+1, 'Lower Conductors', 'FontSize', 9);
yline(h_upper, 'r:', 'LineWidth', 1);
text(max(B_vertical)*0.3, h_upper+1, 'Upper Conductors', 'FontSize', 9);

ylabel('Height Above Ground (m)', 'FontSize', 12);
title('Field Strength at Tower Center (Vertical Profile)', 'FontSize', 13, 'FontWeight', 'bold');
legend('Electric Field', 'Magnetic Field', 'Location', 'northeast');
grid on;


fprintf('\n========================================\n');
fprintf('ELECTROMAGNETIC FIELD ANALYSIS RESULTS\n');
fprintf('========================================\n');
fprintf('Transmission Line: 400 kV Double Circuit\n');
fprintf('Tower Height: %.1f m\n', tower_height);
fprintf('Conductor Height (Upper): %.1f m\n', h_upper);
fprintf('Conductor Height (Lower): %.1f m\n', h_lower);
fprintf('Phase Current: %.0f A\n', I_phase);
fprintf('----------------------------------------\n');
fprintf('ELECTRIC FIELD:\n');
fprintf('  Maximum: %.2f kV/m\n', max(E_kV_per_m(:)));
fprintf('  At ground (center): %.2f kV/m\n', E_kV_per_m(1, round(length(x_range)/2)));
fprintf('  At 15m height (center): %.2f kV/m\n', E_horizontal(round(length(x_range)/2)));
fprintf('  ICNIRP Limit: 5.0 kV/m\n');
fprintf('----------------------------------------\n');
fprintf('MAGNETIC FIELD:\n');
fprintf('  Maximum: %.2f µT\n', max(B_microTesla(:)));
fprintf('  At ground (center): %.2f µT\n', B_microTesla(1, round(length(x_range)/2)));
fprintf('  At 15m height (center): %.2f µT\n', B_horizontal(round(length(x_range)/2)));
fprintf('  ICNIRP Limit: 100 µT\n');
fprintf('========================================\n\n');

fprintf('Program completed successfully.\n');
fprintf('All visualizations have been generated.\n');


function draw_tower(height, base_width, top_width, color, linewidth)
    % Draws a simplified 2D representation of a lattice transmission tower
    
    % Tower legs (tapered from base to top)
    base_half = base_width / 2;
    top_half = top_width / 2;
    
    % Left leg
    plot([-base_half, -top_half], [0, height], color, 'LineWidth', linewidth);
    % Right leg
    plot([base_half, top_half], [0, height], color, 'LineWidth', linewidth);
    
    % Cross bracing (X-pattern at multiple levels)
    num_sections = 5;
    for i = 0:(num_sections-1)
        h1 = i * height / num_sections;
        h2 = (i + 1) * height / num_sections;
        
        w1_left = -base_half + (base_half - top_half) * i / num_sections;
        w1_right = base_half - (base_half - top_half) * i / num_sections;
        w2_left = -base_half + (base_half - top_half) * (i+1) / num_sections;
        w2_right = base_half - (base_half - top_half) * (i+1) / num_sections;
        
        % Diagonal 1
        plot([w1_left, w2_right], [h1, h2], color, 'LineWidth', linewidth*0.6);
        % Diagonal 2
        plot([w1_right, w2_left], [h1, h2], color, 'LineWidth', linewidth*0.6);
        % Horizontal
        plot([w2_left, w2_right], [h2, h2], color, 'LineWidth', linewidth*0.6);
    end
    
    % Cross-arms (horizontal beams for conductors)
    arm_length = base_width * 1.2;
    
    % Upper cross-arm (around 90% height)
    h_upper_arm = 0.90 * height;
    plot([-arm_length/2, arm_length/2], [h_upper_arm, h_upper_arm], color, 'LineWidth', linewidth);
    
    % Lower cross-arm (around 80% height)
    h_lower_arm = 0.80 * height;
    plot([-arm_length/2, arm_length/2], [h_lower_arm, h_lower_arm], color, 'LineWidth', linewidth);
    
    % Top peak
    plot([0, 0], [height-2, height], color, 'LineWidth', linewidth*1.2);
end
