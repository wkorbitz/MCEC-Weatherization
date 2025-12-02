% =========================================================================
% SOFC ANODE 1D TRANSIENT SIMULATION
% DISCLAIMER: On line 522, the "length of the control volume" params.cv_length is wonky and
% "reasonable" results depend on it
% - Solves the Dusty-Gas Model with SMR and WGS chemistry.
% - Uses method of lines with ode15s routine.
% BUG: CH4 values should go down!
% BUG: 1 - psi_eq < 0, reaction may be moving in the opposite direction
% =========================================================================
close all; clc; clearvars;
tic

%% 1. DEFINE PARAMETERS
% --- Species Information
% ORDER: 1=H2, 2=H2O, 3=CO2, 4=CO, 5=CH4
% Fuel cell vs electrolyzer mode depends on the reactant and product
% concentrations

%% DIFFERENT INLET CONFIGURATIONS
% Mastropasqua et al. configuration
x_inlet = [0.5932, 0.2226, 0.1095, 0.0662, 0.00846]';
%x_inlet = [0.2, 0.2226, 0.1095, 0.0662, 0.4]';
% Electrolyzer
%x_inlet = [0.41, 0.04, 0.32, 0.22, 0.01]'; % Electrolyzer init

%% Other parameters
current_ramp = 2900; % say run at 3000 A/m^2

params.species = {'H_2', 'H_2O', 'CO_2', 'CO', 'CH_4'};
params.MM = [2.016, 18.015, 44.009, 28.009, 16.043]';

% --- Operating Conditions
params.T = 800 + 273.15; % [K]
params.p_inlet = 101325; % [Pa]

% --- Physical Constants & Properties
params.A = 80e-4; % [m^2] Electrode active area
params.w = 15 * sqrt(params.A); % [m] Width of fuel channel or entire electrode
params.cv_length = sqrt(params.A) / 15; % length of control volume is wonky

params.R = 8.314; % [J/mol-K]
params.F = 96485.33; % [C/mol]
params.epsilon = 0.45;
params.tau = 3.50;
params.r_p = 1.25e-6; % [m]
params.L = 240e-6; % [m] Thickness of the anode
params.d_p = 2 * params.r_p; % [m]
D_binary_raw = [0,1.58e-3,1.12e-3,1.34e-3,1.19e-3;
                  1.58e-3,0,1.40e-4,1.79e-4,2.26e-4;
                  1.12e-3,1.40e-4,0,8.29e-5,7.34e-5;
                  1.34e-3,1.79e-4,8.29e-5,0,1.09e-4;
                  1.19e-3,2.26e-4,7.34e-5,1.09e-4,0];
params.D_binary = D_binary_raw; % [% m^2/s]
params.mu = 3.567e-05; % [Pa-s]

% --- Kinetic Parameters for Case 1
params.A_WGS = 3.11e-7; % mol/m^2-s
params.E_act_WGS = 103.1e3; % J/mol
params.A_SMR = 4.48e-6; % mol/m^2-s
params.E_act_SMR = 102e3; % J/mol

% --- Kinetic Parameters for Case 5
% params.A_WGS = 1.80e-7; % mol/m^2-s
% params.E_act_WGS = 103.1e3; % J/mol
% params.A_SMR = 1.51e-6; % mol/m^2-s
% params.E_act_SMR = 102e3; % J/mol

% --- Stoichiometry Vectors
% ORDER: 1=H2, 2=H2O, 3=CO2, 4=CO, 5=CH4
% Check units
% Methane should drop
params.nu_WGS = [1, -1, 1, -1, 0]';      % H2O + CO <=> H2 + CO2
params.nu_SMR = [3, -1, 0, 1, -1]';     % CH4 + H2O -> 3H2 + CO
params.nu_el  = [-1, 1, 0, 0, 0]';      % H2 + O^2- -> H2O


%% --- Discretization ---

% Fixed grid
% params.Nz = 50; % Total number of points
% params.dz = params.L / (params.Nz - 1);
% params.z = 0:params.dz:params.L;

%% Chebyshev
params.Nz = 100; % Total number of points
%1. Define the polynomial order (N+1 points, where N = Nz - 1)
N = params.Nz - 1;

% 2. Create the indices from N down to 0
% This generates the nodes in ascending order from -1 to 1.
k_indices = N:-1:0;

% 3. Calculate the standard Chebyshev nodes on the interval [-1, 1]
% x_k = cos(k*pi / N)
cheb_nodes_11 = cos(k_indices * pi / N);
ghost = 2*cheb_nodes_11(N+1) - cheb_nodes_11(N);

% 4. Map the nodes from [-1, 1] to the physical domain [0, L]
% The mapping is z = L * (x + 1) / 2
params.z = (params.L / 2) * (cheb_nodes_11 + 1);
params.zg = (params.L / 2) * (ghost + 1);

% 4. Ensure the start and end points are exact
params.z(1) = 0;
params.z(end) = params.L;

params.zplot = params.z;
params.z = params.z(2:end);
params.Nz=params.Nz-1;
% --- Calculate the non-uniform step sizes (for reference) ---
params.dz_vec = diff(params.z);

params.Ns = length(params.MM);
% --- Butler-Volmer Electrochemical Parameters (from Zhu et al.) ---
params.i0_ref = 8.5; % [A/cm^2] Reference Exchange Current Density [cite: 213]
params.beta = 1.5;   % Apparent transfer coefficient (beta_a + 1) [cite: 215]
params.a_V = 6 * (1-params.epsilon) / params.d_p * 0.15; % Rahmanipour et al.

params.D_eff_binary = params.D_binary * params.epsilon / params.tau;
params.D_Kn = (4/3*params.r_p*sqrt(8*params.R*params.T./(pi()*params.MM))) * params.epsilon / params.tau;
params.B_g = params.epsilon^3*(params.d_p)^2/(72*params.tau*(1-params.epsilon)^2);

%% --- Initial Conditions ---
% ORDER: 1=H2, 2=H2O, 3=CO2, 4=CO, 5=CH4
%x_inlet = [0.04, 0.66, 0.24, 0.04, 0.01]';

x_inlet = x_inlet / sum(x_inlet);
params.x_inlet = x_inlet;
C_tot_inlet = params.p_inlet / (params.R * params.T);
params.C_inlet = x_inlet * C_tot_inlet;
C_init = repmat(params.C_inlet, 1, params.Nz);

%% --- Start of Simulation ---
fprintf('==================================================\n');
fprintf('   MAIN SCRIPT: Simulation started.\n');
fprintf('==================================================\n');

% 1. Set initial conditions
min_diffusivity = min(params.D_eff_binary(params.D_eff_binary~=0));

t_start = 0;
t_char = params.L^2/min_diffusivity;
t_end = 5 * t_char;
tspan = [t_start, t_end];

y0 = C_init(:); % Must be a vector

% --- 2. Define the Mass Matrix ---
% PDE is M*dC/dt = F(C), M is related to porosity.

M_val = 1 / params.epsilon;
M = speye(params.Ns * params.Nz) * M_val;

% --- 3. Set the Solver Options ---
% This is where we give the solver our Jacobian function and Mass matrix.
options = odeset( ...
    'Mass', M, ...
    'RelTol', 1e-5, ... % Set solver tolerances
    'AbsTol', 1e-7, ...
    'NonNegative', 1:params.Nz*params.Ns, ...
    'Stats', 'on');   % Display solver statistics

% --- 4. Call the solver ---
[t_sol, C_sol_vectors] = ode15s(@(t,y) rhs_sofc(t, y, params), tspan, y0, options);
toc
params.Nt = numel(t_sol);
C_sol_inlet = repmat(params.C_inlet', params.Nt, 1);
%C_sol_vectors = [C_sol_inlet, C_sol_vectors];
%params.z = params.zplot;
%params.Nz=params.Nz+1;

% --- 5. Plot results ---

% ==================================================================
% FIGURE 1: STEADY-STATE PROFILES (Initial vs Final)
% Matches Zhu (2005) Fig. 8
% ==================================================================
if exist('C_sol_vectors','var') && size(C_sol_vectors,1) > 1
    
    Ns = params.Ns;
    Nz = params.Nz;

    % --- Initial & Final from solver ---
    C_init  = reshape(C_sol_vectors(1, :),  [Ns, Nz]);
    C_final = reshape(C_sol_vectors(end, :), [Ns, Nz]);

    % --- Mole fractions ---
    X_init  = C_init  ./ sum(C_init,  1);
    X_final = C_final ./ sum(C_final, 1);

    z_microns = params.z * 1e6;
    colors = {'b','r','g','c','k'};  % H2, H2O, CO2, CO, CH4

    % --- Figure 1 ---
    figure('Position', [100, 100, 800, 350], 'Color', 'white');

    % Initial
    subplot(1,2,1);
    for i = 1:Ns
        plot(z_microns, X_init(i,:), '--', 'LineWidth', 1.5, 'Color', colors{i});
        hold on;
    end
    xlabel('Anode Thickness (µm)'); ylabel('Mole Fraction');
    title('Initial Profile (t = 0)');
    xlim([params.z(1)*1e6, params.z(end)*1e6]); grid on; hold off;

    % Final
    subplot(1,2,2);
    for i = 1:Ns
        plot(z_microns, X_final(i,:), '-', 'LineWidth', 2, 'Color', colors{i}, ...
             'DisplayName', params.species{i});
        hold on;
    end
    xlabel('Anode Thickness (µm)'); ylabel('Mole Fraction');
    title('Final Profile (Steady State)');
    legend('Location', 'eastoutside');
    xlim([params.z(1)*1e6, params.z(end)*1e6]); grid on; hold off;

    sgtitle('SOFC Anode: Initial and Final Mole Fraction Profiles', ...
            'FontWeight', 'bold', 'FontSize', 12);
end

% ==================================================================
% FIGURE 2: TRANSIENT EVOLUTION
% Matches Mastropasqua (2019) Fig. 7
% ==================================================================
if exist('t_sol','var') && exist('C_sol_vectors','var')
    
    Ns = params.Ns;
    Nz = params.Nz;
    C_3D = reshape(C_sol_vectors', [Ns, Nz, length(t_sol)]);

    z_fuel = 1;

    z_target = 0.061;  % µm
    [~, idx_near] = min(abs(z_microns - z_target));
    z_near = z_microns(idx_near);

    z_mid  = max(2, min(Nz-1, round(Nz/2)));
    z_tpb  = Nz;

    X_fuel = squeeze(C_3D(:, z_fuel, :)) ./ sum(squeeze(C_3D(:, z_fuel, :)), 1);
    X_near = squeeze(C_3D(:, idx_near, :)) ./ sum(squeeze(C_3D(:, idx_near, :)), 1);
    X_mid  = squeeze(C_3D(:, z_mid , :)) ./ sum(squeeze(C_3D(:, z_mid , :)), 1);
    X_tpb  = squeeze(C_3D(:, z_tpb , :)) ./ sum(squeeze(C_3D(:, z_tpb , :)), 1);

    colors = {'b','r','g','c','k'};

    % --- Figure 2 ---
    figure('Position', [100, 100, 1000, 350], 'Color', 'white');

    % Fuel Channel
    subplot(1,4,1);
    for i = 1:Ns
        plot(t_sol, X_fuel(i,:), '-', 'LineWidth', 2, 'Color', colors{i});
        hold on;
    end
    xlabel('Time (s)'); ylabel('Mole Fraction');
    title(sprintf('Fuel Channel (z = 0 µm)'));
    xlim([0, t_sol(end)]); grid on; hold off;

    % Right after the inlet
    subplot(1,4,2)
    for i = 1:Ns
        plot(t_sol, X_near(i,:), '-', 'LineWidth', 2.5, 'Color', colors{i}, ...
             'DisplayName', params.species{i});
        hold on;
    end
    xlabel('Time (s)'); ylabel('Mole Fraction');
    title(sprintf('Transient at z = %.1f µm (Near Inlet)', z_near));
    xlim([0, t_sol(end)]);
    legend('Location', 'eastoutside'); grid on; hold off;

   
    % Middle
    subplot(1,4,3);
    for i = 1:Ns
        plot(t_sol, X_mid(i,:), '-', 'LineWidth', 2, 'Color', colors{i});
        hold on;
    end
    xlabel('Time (s)');
    title(sprintf('Middle (z = %.1f µm)', params.z(z_mid)*1e6));
    xlim([0, t_sol(end)]); grid on; hold off;

    % TPB
    subplot(1,4,4);
    for i = 1:Ns
        plot(t_sol, X_tpb(i,:), '-', 'LineWidth', 2, 'Color', colors{i});
        hold on;
    end
    xlabel('Time (s)');
    title(sprintf('TPB (z = %.1f µm)', params.z(end)*1e6));
    xlim([0, t_sol(end)]); grid on; hold off;

    sgtitle('SOFC Anode: Transient Mole Fraction Evolution', ...
            'FontWeight', 'bold', 'FontSize', 12);
end

% ==================================================================
% 2x3 SUBPLOT: Initial/Final Profiles + Transients
% CORRECT: C_init from C_sol_vectors(1,:), C_final from end
% Matches Zhu (2005) Fig. 8 & Mastropasqua (2019) Fig. 7
% ==================================================================
% if exist('t_sol','var') && exist('C_sol_vectors','var') && size(C_sol_vectors,1) > 1
% 
%     Ns = params.Ns;
%     Nz = params.Nz;
% 
%     % --- TRUE INITIAL & FINAL from solver output ---
%     C_init_vec  = C_sol_vectors(1, :);     % First time step
%     C_final_vec = C_sol_vectors(end, :);   % Last time step
% 
%     C_init  = reshape(C_init_vec,  [Ns, Nz]);
%     C_final = reshape(C_final_vec, [Ns, Nz]);
% 
%     % --- z in microns ---
%     z_microns = params.z * 1e6;
% 
%     % --- Mole fractions ---
%     X_init  = C_init  ./ sum(C_init,  1);
%     X_final = C_final ./ sum(C_final, 1);
% 
%     % --- Transients ---
%     C_3D = reshape(C_sol_vectors', [Ns, Nz, length(t_sol)]);
% 
%     z_fuel = 1;
%     z_mid  = max(2, min(Nz-1, round(Nz/2)));
%     z_tpb  = Nz;
% 
%     X_fuel = squeeze(C_3D(:, z_fuel, :)) ./ sum(squeeze(C_3D(:, z_fuel, :)), 1);
%     X_mid  = squeeze(C_3D(:, z_mid , :)) ./ sum(squeeze(C_3D(:, z_mid , :)), 1);
%     X_tpb  = squeeze(C_3D(:, z_tpb , :)) ./ sum(squeeze(C_3D(:, z_tpb , :)), 1);
% 
%     % --- Colors ---
%     colors = {'b','r','g','c','k'};  % H2, H2O, CO2, CO, CH4
% 
%     % === FIGURE ===
%     figure('Position', [100, 100, 1000, 650]);
% 
%     % --- 1. Initial Profile (from t=0) ---
%     subplot(2,3,1);
%     for i = 1:Ns
%         plot(z_microns, X_init(i,:), '--', 'LineWidth', 1.5, 'Color', colors{i});
%         hold on;
%     end
%     xlabel('Anode Thickness (µm)'); ylabel('Mole Fraction');
%     title('Initial Profile (t = 0)');
%     xlim([0, params.L*1e6]); grid on; hold off;
% 
%     % --- 2. Final Profile (steady state) ---
%     subplot(2,3,2);
%     for i = 1:Ns
%         plot(z_microns, X_final(i,:), '-', 'LineWidth', 2, 'Color', colors{i}, ...
%              'DisplayName', params.species{i});
%         hold on;
%     end
%     xlabel('Anode Thickness (µm)'); ylabel('Mole Fraction');
%     title('Final Profile');
%     legend('Location', 'eastoutside'); xlim([0, params.L*1e6]); grid on; hold off;
% 
%     % --- 3. Transient: Fuel Channel ---
%     subplot(2,3,4);
%     for i = 1:Ns
%         plot(t_sol, X_fuel(i,:), '-', 'LineWidth', 2, 'Color', colors{i});
%         hold on;
%     end
%     xlabel('Time (s)'); ylabel('Mole Fraction');
%     title('Fuel Channel (z = 0 µm)');
%     grid on; hold off;
% 
%     % --- 4. Transient: Middle ---
%     subplot(2,3,5);
%     for i = 1:Ns
%         plot(t_sol, X_mid(i,:), '-', 'LineWidth', 2, 'Color', colors{i});
%         hold on;
%     end
%     xlabel('Time (s)'); ylabel('Mole Fraction');
%     title(sprintf('Middle (z = %.1f µm)', params.z(z_mid)*1e6));
%     grid on; hold off;
% 
%     % --- 5. Transient: TPB ---
%     subplot(2,3,6);
%     for i = 1:Ns
%         plot(t_sol, X_tpb(i,:), '-', 'LineWidth', 2, 'Color', colors{i});
%         hold on;
%     end
%     xlabel('Time (s)'); ylabel('Mole Fraction');
%     title(sprintf('TPB (z = %.1f µm)', params.z(end)*1e6));
%     grid on; hold off;
% end

% ------------------------------------------------------------------
% OLD TRANSIENT SUBPLOT – ZHU & MASTROPASQUA STYLE
% ------------------------------------------------------------------
% if exist('t_sol','var') && exist('C_sol_vectors','var')
%     Nt = length(t_sol);
%     Ns = params.Ns;
%     Nz = params.Nz;
%     z_microns = params.z * 1e6;
% 
%     C_init  = reshape(C_init,  [Ns, Nz]);   % [5×50]
%     C_final = reshape(C_final, [Ns, Nz]);   % [5×50]
%     X_init  = C_init  ./ sum(C_init, 1);
%     X_final = C_final ./ sum(C_final, 1);
% 
%     % Reshape: [Nt x (Ns*Nz)] → [Ns x Nz x Nt]
%     C_3D = reshape(C_sol_vectors', [Ns, Nz, Nt]);
% 
%     % --- Safe node indices ---
%     z_fuel = 1;
%     z_mid  = max(2, min(Nz-1, round(Nz/2)));
%     z_tpb  = Nz;
% 
%     % --- Extract mole fractions [Ns x Nt] ---
%     X_fuel = squeeze(C_3D(:, z_fuel, :)) ./ sum(squeeze(C_3D(:, z_fuel, :)), 1);
%     X_mid  = squeeze(C_3D(:, z_mid , :)) ./ sum(squeeze(C_3D(:, z_mid , :)), 1);
%     X_tpb  = squeeze(C_3D(:, z_tpb , :)) ./ sum(squeeze(C_3D(:, z_tpb , :)), 1);
% 
%     % --- Create subplot ---
%     figure('Name', 'Transient Mole Fraction Profiles', 'Position', [100, 100, 900, 600]);
%     colors = {'b','r','g','c','k'};  % H2, H2O, CO2, CO, CH4
% 
%     % Subplot 1: Fuel Channel
%     subplot(3,1,1);
%     for i = 1:Ns
%         plot(t_sol, X_fuel(i,:), '-', 'LineWidth', 2, 'Color', colors{i});
%         hold on;
%     end
%     ylabel('Mole Fraction');
%     title(sprintf('Fuel Channel (z = 0 µm)'));
%     grid on; legend(params.species, 'Location', 'eastoutside'); hold off;
% 
%     % Subplot 2: Middle
%     subplot(3,1,2);
%     for i = 1:Ns
%         plot(t_sol, X_mid(i,:), '-', 'LineWidth', 2, 'Color', colors{i});
%         hold on;
%     end
%     ylabel('Mole Fraction');
%     title(sprintf('Middle (z = %.1f µm)', params.z(z_mid)*1e6));
%     grid on; hold off;
% 
%     % Subplot 3: TPB
%     subplot(3,1,3);
%     for i = 1:Ns
%         plot(t_sol, X_tpb(i,:), '-', 'LineWidth', 2, 'Color', colors{i});
%         hold on;
%     end
%     xlabel('Time (s)'); ylabel('Mole Fraction');
%     title(sprintf('TPB (z = %.1f µm)', params.z(end)*1e6));
%     grid on; hold off;
% 
%     % Link x-axes for zooming
%     linkaxes([subplot(3,1,1), subplot(3,1,2), subplot(3,1,3)], 'x');
% end

% =========================================================================
%                        HELPER FUNCTIONS
% =========================================================================

function J_f = compute_J_f(C_guess, pert, params)
    N_vars = numel(C_guess); J_f = zeros(N_vars, N_vars);
    f_base_vec = compute_F(C_guess, params); f_base_vec = f_base_vec(:);
    for k = 1:N_vars
        C_pert = C_guess; C_pert(k) = C_pert(k) + pert;
        f_pert_vec = compute_F(C_pert, params); f_pert_vec = f_pert_vec(:);
        J_f(:, k) = (f_pert_vec - f_base_vec) / pert;
    end
end
function R = compute_residual(C_guess, C_n, dt, params)
    F_guess = compute_F(C_guess, params);
    R_matrix = C_guess - C_n - dt * F_guess;
    R = R_matrix(:);
end

function [F_of_C, r_WGS, r_SMR] = compute_F(C, params)
    [Ns, Nz] = size(C);
    z = params.z;
    %J = zeros(Ns, Nz);

    % --- Transport Term (Flux Calculation) ---
    % This part remains the same as it correctly models the Dusty-Gas Model
    C_inlet = params.x_inlet * (params.p_inlet / (params.R * params.T));
    %% - INLET BOUNDARY CONDITIONS HAVE TO BE ENFORCED AT EVERY TIME STEP

    %C(:,1) = C_inlet;  % Force concentration at fuel channel
    %h = [z(end) - z(end-1), z(end) - z(end-2)];
    %b = 1 / (h(1)^2/h(2) - h(1));
    %c = -b*h(1)^2/h(2)^2;
    %a = -b-c;
    %C(:, end) = C(:,end-1); % 3-point 1st derivative

    %% - Continue Calculations
    C_tot = sum(C, 1);
    P = C_tot * params.R * params.T;
    P_inlet = sum(params.x_inlet,1) * params.R * params.T;
    X = C ./ C_tot;
    dJdz = zeros(Ns, Nz);
    
    % build one-sided stencil coefficients
    h1 = z(1);
    h2 = z(2)-z(1);
    hN1 = z(Nz) - z(Nz-1);
    %hN2 = hN1;
    %gRb = 1/(-hN1 + hN1^2/hN2); gRc = (-hN1^2/hN2^2)*gRb; gRa=-gRb-gRc;
    %gradR = [gRa,gRb,gRc];
    gLa = 1/(-h1 - h1^2/h2); gLc = (-h1^2/h2^2)*gLa; gLb=-gLa-gLc;
    gradL = [gLa,gLb,gLc];
    lLb = 1/(h1^2 - h1*h2); lLc = (-h1/h2)*lLb; lLa=-lLb-lLc;
    lapL = [lLa,lLb,lLc];
    %lRb = 1/(hN1^2 - hN1*hN2); lRc = (-hN1/hN2)*lRb; lRa=-lRb-lRc;
    %lapR = [lRa,lRb,lRc];
    
    % build left-hand (z_1) terms
    X_local = C(1) / C_tot(1);
    H = build_H_matrix(X_local, params.D_eff_binary, params.D_Kn);
    D_DGM = H \ eye(Ns);
    grad_X = gradL(1)*params.x_inlet + gradL(2)*X(:,1) + gradL(3)*X(:,2);
    grad_P = gradL(1)*P_inlet + gradL(2)*P(1) + gradL(3)*P(2);
    lap_X = lapL(1)*params.x_inlet + lapL(2)*X(:,1) + lapL(3)*X(:,2);
    lap_P = lapL(1)*P_inlet + lapL(2)*P(1) + lapL(3)*P(2);
    dJdz(:,1) = -D_DGM*lap_X - (params.B_g / params.mu) * ...
        (grad_X*grad_P + params.x_inlet*lap_P);

    % build right-hand (z_N) terms
    X_local = X(:, Nz);
    H = build_H_matrix(X_local, params.D_eff_binary, params.D_Kn);
    D_DGM = H \ eye(Ns);
    %grad_X = 0;
    %grad_P = gradR(1)*P(Nz) + gradR(2)*P(Nz-1) + gradR(3)*P(Nz-2);
    %P inherits neumann BC from X
    lap_X = 2*(X(:,Nz-1)-X(:,Nz))/hN1^2;
    lap_P = 2*(P(:,Nz-1)-P(:,Nz))/hN1^2;
    dJdz(:,Nz) = -D_DGM * lap_X - (params.B_g / params.mu) * X(:,Nz)*lap_P;

    for i = 2:Nz-1
        X_local = X(:, i);
        H = build_H_matrix(X_local, params.D_eff_binary, params.D_Kn);
        D_DGM = H \ eye(Ns);
        grad_X = (X(:, i+1) - X(:, i-1)) / ((z(i+1) - z(i-1)));
        hp = z(i+1) - z(i);
        hn = z(i) - z(i-1);
        kp = 2 / (hp*hn + hp^2);
        kn = hp/hn * kp;
        k = -kn -kp;
        lap_X = kp*X(:,i+1) + kn*X(:,i-1) + k*X(:,i);
        % grad_P = (P(:, i+1) - P(:, i-1)) / ((z(i+1) - z(i-1))); % scalar
        grad_P = (P(i+1) - P(i-1)) / (z(i+1) - z(i-1));
        lap_P = kp*P(i+1) + kn*P(i-1) + k*P(i);
        %kn_term = D_DGM * (X_local ./ params.D_Kn);
        dJdz(:,i) = -D_DGM * lap_X - (params.B_g / params.mu) * ...
            (grad_X*grad_P + X(:,i)*lap_P);
        %J(:,i) = -D_DGM * grad_X - (params.B_g / params.mu) * kn_term * grad_P;
    end
    
    

    %% - Continue Calculations
    % --- Transport Term (Flux Divergence) ---
    % This part also remains the same

    %for i = 2:Nz-1
    %    dJdz(:,i) = (J(:, i+1) - J(:, i-1)) / ((z(i+1) - z(i-1)));
    %end

    %% - OUTLET BOUNDARY CONDITIONS HAVE TO BE ENFORCED AT EVERY TIME STEP
    %dJdz(:, Nz) = -J(:, Nz-1)/(z(i)-z(i-1));

    %dJdz(:, end) = 0;

    % --- Source Term (Reaction Kinetics) ---

    % 1. Calculate Partial Pressures
    %C_tot = sum(C, 1);
    P_gas = (C * params.R * params.T);
    p_H2  = P_gas(1,:);
    p_H2O = P_gas(2,:);
    p_CO2 = P_gas(3,:);
    p_CO  = P_gas(4,:);
    p_CH4 = P_gas(5,:);

    % 2. Calculate Forward Rate Constants (k_f) mol/m^2-s, Arrhenius
    % formulation
    k_f_WGS = params.A_WGS * exp(-params.E_act_WGS / (params.R * params.T)); % mol/m^3-s Pa^???
    k_f_SMR = params.A_SMR * exp(-params.E_act_SMR / (params.R * params.T)); % mol/m^2-s Pa^???

    % JANAF Table Polynomial Fits
    deltaG_steam = -2e-10 * params.T^3 + 2e-6 * params.T^2 + 0.0512 * params.T - 244.82;
    deltaG_CO_2 = -3e-18 * params.T^5 + 6e-14 * params.T^4 -4e-10 * params.T^3 + 2e-6 * params.T^2 - 0.004 * params.T - 393.36;
    deltaG_CO = 2e-6 * params.T^2 - 0.0905 * params.T - 111.08;
    deltaG_CH_4 = -2e-10 * params.T^3 + 3e-6 * params.T^2 + 0.1018 * params.T - 82.935;

    deltaG_r_SMR = (deltaG_CO + 3*0) - (deltaG_CH_4 + deltaG_steam);
    deltaG_r_WGS = deltaG_CO_2 - (deltaG_CO + deltaG_steam);

    deltaG_r_SMR = deltaG_r_SMR * 1000; % convert kJ/mol to J/mol
    deltaG_r_WGS = deltaG_r_WGS * 1000; % convert kJ/mol to J/mol

    K_eq_SMR = exp(-deltaG_r_SMR / (params.R * params.T));
    K_eq_WGS = exp(-deltaG_r_WGS / (params.R * params.T));

    % 4. Calculate the Equilibrium Term (psi_eq)
    Kp_WGS = (p_H2 .* p_CO2) ./ (p_H2O .* p_CO);
    psi_eq_WGS = Kp_WGS ./ K_eq_WGS;
    
    Kp_SMR = ((p_H2.^3) .* p_CO) ./ (p_H2O .* p_CH4);
    psi_eq_SMR = Kp_SMR ./ K_eq_SMR;
    
    % 5. Calculate Final Net Reaction Rates (r) [mol/m2-s]
    % Using the paper's power-law form (Eq. 23 and 24), assuming reaction orders are 1

    % k_a_CO = 1;

    %% Kinetic Parameters for Case 1
    delta_H2O = 0.5;
    chi_CO2 = 1;
    alpha_CH4 = 1;
    beta_H2O = 0;
    gamma_CO = 1;

    %% Kinetic Parameters for Case 5
    % gamma_CO = 0.813;
    % delta_H2O = 0.492;
    % chi_CO2 = 0.810;
    % alpha_CH4 = 0.967;
    % beta_H2O = 0.005;

    driving_WGS = 1 - psi_eq_WGS;

    params.length_char = params.L * params.w * params.cv_length ./ params.A; % length of control volume is wonky
    r_WGS = k_f_WGS .* p_CO .^ gamma_CO .* p_H2O .^ delta_H2O .* p_CO2 .^ chi_CO2 .* driving_WGS .* params.length_char;

    driving_SMR = 1 - psi_eq_SMR;
    
    r_SMR = k_f_SMR .* p_CH4 .^ alpha_CH4 .* p_H2O .^ beta_H2O .* driving_SMR .* params.length_char;
    
    %sprintf("ksWGS = %d, kfSMR=%d",k_f_WGS,k_f_SMR)
    %sprintf("driving WGS = %d, SMR = %d",driving_WGS, driving_SMR)
    S_het = params.a_V * (params.nu_WGS * r_WGS + params.nu_SMR * r_SMR); % mol/m^3-s

    % % 2. Electrochemical Source Term (S_el)
    
    % --- COMBINE TERMS to get dC/dt ---
    F_of_C = -dJdz + S_het ;%+ S_el;
end

function H = build_H_matrix(x, D_binary, D_Kn)
    Ns = length(x); H = zeros(Ns, Ns);
    for k = 1:Ns
        sum_term = 0;
        for j = 1:Ns, if j ~= k, sum_term = sum_term + x(j) / D_binary(k, j); end, end
        H(k, k) = (1 / D_Kn(k)) + sum_term;
        for l = 1:Ns, if k ~= l, H(k, l) = -x(k) / D_binary(k, l); end, end
    end
end

% ODE15S Implementation

function F = rhs_sofc(t, y, params)
    % This function computes the right-hand side (RHS) of the ODE system dy/dt = F(t,y).
    % It's a wrapper for your existing compute_F function.
    
    % 1. Reshape the state vector 'y' into a matrix 'C'
    % ode15s uses a column vector, but our physics uses a [Ns x Nz] matrix.
    C = reshape(y, [params.Ns, params.Nz]);

    % 2. Call your physics function to get the time derivative
    % We only need the first output (F_of_C) from your compute_F function.
    [F_matrix, r_WGS, r_SMR] = compute_F(C, params);
    
    % 3. Reshape the resulting matrix back into a column vector for ode15s
    F = F_matrix(:);

    % Store for plotting
    persistent r_WGS_history r_SMR_history t_history
    if isempty(r_WGS_history)
        r_WGS_history = []; r_SMR_history = []; t_history = [];
    end
    r_WGS_history = [r_WGS_history; r_WGS];
    r_SMR_history = [r_SMR_history; r_SMR];
    t_history = [t_history; t];
    assignin('base', 'r_WGS_transient', r_WGS_history);
    assignin('base', 'r_SMR_transient', r_SMR_history);
    assignin('base', 't_transient', t_history);
end

function J = jacobian_sofc(t, y, params)
    % This function computes the Jacobian of the ODE system.
    % It's a wrapper for your existing compute_J_f function.
    
    % 1. Reshape the state vector 'y' into a matrix 'C'
    C = reshape(y, [params.Ns, params.Nz]);
    
    % 2. Call your function to compute the Jacobian matrix
    % For stiff problems, the perturbation 'pert' should be small.
    pert = 1e-8; 
    J = compute_J_f(C, pert, params);
end

function plot_transient(X_mat, t_vec, params, title_str)
    figure('Name',title_str);
    colors = {'b','r','g','c','k'};
    for i = 1:params.Ns
        plot(t_vec, X_mat(i,:), '-', 'LineWidth',2, ...
             'Color',colors{i}, 'DisplayName',params.species{i});
        hold on;
    end
    xlabel('Time (s)'); ylabel('Mole Fraction');
    title(title_str);
    legend('Location','best'); grid on;
    hold off;
end

C_final = reshape(C_sol_vectors(end,:), [params.Ns, params.Nz]);
[~, r_WGS_final, r_SMR_final] = compute_F(C_final, params);

z_microns = params.z * 1e6;

figure('Position', [100, 100, 800, 350]);
subplot(1,2,1);
plot(z_microns, r_WGS_final, 'r-', 'LineWidth', 2);
xlabel('Anode Thickness (µm)'); ylabel('r_{WGS} (mol/m²-s)');
title('WGS Reaction Rate'); grid on;

subplot(1,2,2);
plot(z_microns, r_SMR_final, 'b-', 'LineWidth', 2);
xlabel('Anode Thickness (µm)'); ylabel('r_{SMR} (mol/m²-s)');
title('SMR Reaction Rate'); grid on;

sgtitle('Reaction Rates at Steady State');
