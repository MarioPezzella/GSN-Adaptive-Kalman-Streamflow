% =========================================================================
% ADAPTIVE KALMAN FILTERING OF DAILY DISCHARGE - CAMASTRA RESERVOIR
%
% Software and code developed by: MARIO PEZZELLA
%
% This software is associated with the manuscript
%       "Adaptive Kalman Filtering for Streamflow Models with Application
%        to the Camastra Basin"
%
% Authors:
%       Mario Pezzella†, Luigi Cimorelli, Luisa D'Amore, Salvatore Cuomo
%                         † Corresponding author.
%
% Distributed under the terms of the GNU GPL v. 3 license (see LICENSE.md).
%
% =========================================================================
% This script applies the adaptive Kalman filter implemented in
% discharge_kalman.m to the daily inflow record of the Camastra reservoir
% (Basilicata, Italy), 1984-2020, described by the Generalized Shot Noise
% (GSN) model calibrated by Cimorelli et al. (2021).
%
% The filter is run over the complete record, with no training/test split
% and no rolling forecast, under two rainfall configurations:
%
%   SCENARIO A - a synthetic effective-rainfall sequence z_hat_d(n),
%                generated from the monthly PWNE model, is supplied as a
%                deterministic input. The instantaneous fraction alpha_0
%                of it reaches the output directly, outside the
%                measurement update.
%   SCENARIO B - no rainfall information is supplied. The unresolved
%                forcing is represented through the process-noise
%                covariance Qd^(B), initialized from the annual mean PWNE
%                spectral constant.
%
% Two filter modes are available:
%   USE_ADAPTIVE = true   Sage-Husa adaptive filter, Q and R estimated
%                         online (default, used for the manuscript)
%   USE_ADAPTIVE = false  Standard Kalman filter, Q and R calibrated
%                         offline by MultiStart so that mean NIS = 1
%                         (requires the Optimization and Global
%                         Optimization Toolboxes)
%
% The script then evaluates accuracy, consistency and predictive-band
% coverage of both scenarios, over the full record and after a burn-in
% period, and conditionally on the hydrological regime.
%
% -------------------------------------------------------------------------
% OUTPUTS
% -------------------------------------------------------------------------
%   Command Window:
%     full metric report, post-burn-in metrics, regime-conditional
%     metrics for Scenario B, and the ten largest NIS values with dates.
%
%   Figures (exported to FIG_OUTPUT_DIR = ./figures_out, vector PDF):
%     fig01_timeseries_data_scenA_scenB   full record, 3 stacked panels
%     fig02_hydrograph_zoom               representative years, daily
%     fig03_relerror_scenA_scenB          relative error, 2x2
%                                         (scenario x pre/post update)
%     fig04_monthly_RMSE_bias             monthly RMSE and bias, 2x2
%     fig05_adaptive_R_Qd                 online R_n and Qd_n (adaptive)
%     figS1_calibration_diagnostics       PIT | reliability | NIS QQ-plot
%                                         (only if PLOT_CALIB_DIAG = true)
%
%   LaTeX tables, printed to the Command Window (PRINT_TABLES = true)
%   and optionally written to FIG_OUTPUT_DIR (TABLES_TO_FILE = true):
%     tab_results   accuracy and consistency, both scenarios
%     tab_regime    regime-conditional metrics, Scenario B
%     tab_events    ten largest NIS values, Scenario B
%
% -------------------------------------------------------------------------
% DEPENDENCIES
% -------------------------------------------------------------------------
%   Functions, expected in the same folder as this script:
%     discharge_kalman.m              the Kalman filter
%     simulate_pwne_year_inverse.m    PWNE rainfall generator (Scenario A)
%   Data:
%     ./Dati/q_YYYY.txt   one file per year, one daily discharge value
%                         [m^3/s] per line, NaN for missing observations
%   Toolboxes:
%     Base MATLAB only in adaptive mode. The offline calibration branch
%     (USE_ADAPTIVE = false) additionally requires the Optimization and
%     Global Optimization Toolboxes.
%
% -------------------------------------------------------------------------
% REPRODUCIBILITY
% -------------------------------------------------------------------------
%   The Scenario A input is a Monte Carlo realization of the PWNE process:
%   it reproduces the calibrated monthly statistics of the effective
%   rainfall, but its individual events are not synchronized with those
%   that generated the observed discharge. RNG_SEED fixes the realization
%   used in the manuscript; change it to explore the sampling variability.
%   Scenario B does not depend on the seed.
% =========================================================================

clear; clc; close all;

RNG_SEED = 30;
rng(RNG_SEED);

%% 0. Configuration
% =========================================================================
% Run options, reporting choices and figure styling
% =========================================================================
USE_ADAPTIVE    = true;      % true = Sage-Husa, false = offline-calibrated
EXPORT_FIGURES  = true;      % export figures to FIG_OUTPUT_DIR
PRINT_TABLES    = true;      % LaTeX tables to the Command Window
% Distributional diagnostics (PIT histogram, reliability diagram, NIS
% QQ-plot against chi^2_1), exported as a supplementary figure. They go
% beyond the first-moment consistency E[NIS] = 1 targeted by the
% adaptation, and show how the full NIS distribution departs from chi^2_1.
PLOT_CALIB_DIAG = false;
TABLES_TO_FILE  = false;     % also write the tables as \input-able files
FIG_FORMATS     = {'pdf'};   % e.g. {'pdf','png'}

yr_start = 1984;   yr_end = 2020;

% Post-burn-in reporting. Qd^(B) is initialized from the climatological
% annual average of G_m and then adapted (fig05), so the first years
% include an initialization transient. Full-record and post-burn-in
% metrics are reported side by side.
BURN_IN_YEARS   = 5;

% Years shown at daily resolution in fig02: a wet year with a major event
% and a contrasting one.
ZOOM_YEARS      = [2013 2017];

% Seasonal split for the regime-conditional table. The dry months reflect
% lambda_m = 0 from July to September in the PWNE calibration.
WET_MONTHS      = [11 12 1 2 3 4];
DRY_MONTHS      = [6 7 8 9];

BAND_Z          = 1.96;      % uncertainty band multiplier used in figures
REL_ERR_FLOOR_Q = 1;       % [m^3/s] denominator floor for relative error

% -------------------------------------------------------------------------
% Global figure styling (applies to every figure created below)
% -------------------------------------------------------------------------
set(groot,'defaultAxesFontSize',11);
set(groot,'defaultAxesLineWidth',0.9);
set(groot,'defaultAxesBox','on');
set(groot,'defaultLegendFontSize',9);
set(groot,'defaultLineLineWidth',1.2);
set(groot,'defaultFigureColor','w');

% Colour scheme, defined once and reused throughout
clr_obs    = [0.12 0.12 0.12];
clr_A      = [0.13 0.36 0.71];   % blue  - Scenario A
clr_B      = [0.09 0.56 0.30];   % green - Scenario B
clr_A_filt = [0.00 0.08 0.45];   clr_A_pred = [0.25 0.55 0.95];
clr_B_filt = [0.00 0.30 0.05];   clr_B_pred = [0.40 0.80 0.25];
clr_mean   = [0.80 0.10 0.10];
clr_p90    = [0.45 0.45 0.45];

% Figure handles are collected as they are created, each paired with its
% export filename, so that the export does not depend on creation order.
FIGS = cell(0,2);

%% 1. GSN model parameters
% =========================================================================
% PCSTC calibration of the Camastra basin (Cimorelli et al. 2021):
% single-reservoir unit response and monthly PWNE parameters
% =========================================================================
alpha = [0.201, 0.799];    % [alpha_0, alpha_1]
k_vec = [38.431];          % reservoir characteristic time [days]
dt    = 1.0;               % time step [days]
ns    = numel(k_vec);      % number of reservoirs
c     = (alpha(2:end) ./ k_vec)';   % output vector, c_j = alpha_j/k_j

% Monthly Poisson rates lambda_m [1/day] and mean magnitudes eta_c,m
% (Jan..Dec). lambda_m = 0 from July to September.
lambda_month = [0.0707, 0.1195, 0.1625, 0.0949, 0.0223, 0.0013, ...
                0.0000, 0.0000, 0.0000, 0.0037, 0.0321, 0.0658];
eta_c_month  = [93.8420, 83.5171, 61.0562, 62.9453, 78.7795, 40.1870, ...
                 0.0000,  0.0000,  0.0000, 89.9786, 121.4535, 101.5331];

% Scenario B represents the seasonally varying spectral constant
% G_m = 2*eta_c,m^2*lambda_m by its annual mean, re-expressed in the
% unit-rate PWNE parametrization: lambda_B = 1, eta_c,B = sqrt(G/2).
G_annual = mean(2 .* eta_c_month.^2 .* lambda_month);
lambda_B = 1;
eta_c_B  = sqrt(G_annual / 2);

month_names = {'Jan','Feb','Mar','Apr','May','Jun', ...
               'Jul','Aug','Sep','Oct','Nov','Dec'};

%% 2. Data loading
% =========================================================================
% Reads the yearly discharge files, checks the number of days of each
% year against the calendar, and builds the daily time axis
% =========================================================================
base_dir       = fileparts(mfilename('fullpath'));
if isempty(base_dir), base_dir = pwd; end   % allows running by copy-paste
FIG_OUTPUT_DIR = fullfile(base_dir, 'figures_out');
data_dir       = fullfile(base_dir, 'Dati');

assert(exist('discharge_kalman','file') == 2, ...
       'discharge_kalman.m not found on the MATLAB path.');
assert(exist('simulate_pwne_year_inverse','file') == 2, ...
       'simulate_pwne_year_inverse.m not found on the MATLAB path.');
assert(isfolder(data_dir), 'Data folder not found: %s', data_dir);

years   = yr_start:yr_end;
n_years = numel(years);

y_cells = cell(n_years,1);
yr_days = zeros(n_years,1);
for i = 1:n_years
    fname = fullfile(data_dir, sprintf('q_%d.txt', years(i)));
    assert(isfile(fname), 'Missing data file: %s', fname);
    y_cells{i} = readmatrix(fname);
    yr_days(i) = numel(y_cells{i});
end

% -------------------------------------------------------------------------
% Leap-year check: an off-by-one in a yearly file would shift the whole
% record against the calendar, so the length of each file is verified.
% -------------------------------------------------------------------------
expected = zeros(n_years,1);
for i = 1:n_years
    expected(i) = days(datetime(years(i)+1,1,1) - datetime(years(i),1,1));
end
bad = find(yr_days ~= expected);
if isempty(bad)
    fprintf('Leap-year check: all %d year files correct.\n', n_years);
else
    for kb = bad(:)'
        warning('Year %d: expected %d days, found %d.', ...
                years(kb), expected(kb), yr_days(kb));
    end
end

% -------------------------------------------------------------------------
% Concatenated record and calendar
% -------------------------------------------------------------------------
y_meas   = vertcat(y_cells{:})';            % 1 x N row vector
cum_days = [0; cumsum(yr_days)];
N        = cum_days(end);
t        = (1:N)';                          % day index [1...N]
dates    = datetime(yr_start,1,1) + caldays(0:N-1);
month_d  = month(dates);
year_d   = year(dates);

fprintf('Dataset : %d years, %d days (%d-%d)\n', n_years, N, yr_start, yr_end);
fprintf('Missing : %d  (%.2f %%)\n\n', sum(isnan(y_meas)), 100*mean(isnan(y_meas)));

%% 3. PWNE effective-rainfall realization (Scenario A only)
% =========================================================================
% One Monte Carlo realization of the monthly PWNE process, year by year.
% Statistically representative of the effective rainfall, but not
% synchronized with the historical events (see the header).
% =========================================================================
Z = zeros(N,1);
for i = 1:n_years
    idx    = cum_days(i)+1 : cum_days(i+1);
    Z(idx) = simulate_pwne_year_inverse(lambda_month, eta_c_month, yr_days(i));
end
z_hat_d = Z';

%% 4. Filter options
% =========================================================================
% Scenario-specific options, Sage-Husa settings and, in standard mode,
% offline calibration of Q and R
% =========================================================================
% Preliminary measurement variance: 5% of the variance of the record
R0_heur = max(var(y_meas,'omitnan') * 0.05, 1e-8);

opts_A = struct('scenario','A', 'z_hat_d',z_hat_d, 'Q_diag',1, ...
                'R',R0_heur, 'gain_tol',1e-13, 'gain_patience',5);

opts_B = struct('scenario','B', 'lambda',lambda_B, 'eta_c',eta_c_B, ...
                'R',R0_heur, 'gain_tol',1e-13, 'gain_patience',5);

% Sage-Husa settings shared by both scenarios (NIS buffers of 180 days)
sage_p = struct('adaptive',true, 'sage_bootstrap',100, ...
                'sage_lambda_min',0.95, 'sage_window',30*6, ...
                'sage_R_min',1e-8, 'sage_R_max',Inf, ...
                'sage_Q_min',1e-10, 'sage_Q_max',Inf);

if USE_ADAPTIVE
    mode_str = 'Sage-Husa';
    fprintf('Mode: Sage-Husa adaptive\n\n');
    opts_A = merge_s(opts_A, sage_p);
    opts_B = merge_s(opts_B, sage_p);
    % The bootstrap estimator of R relies on a small prior variance over
    % the bootstrap window. Under Scenario B, Qd^(B) encodes the full
    % rainfall uncertainty, the open-loop prior variance is large and the
    % subtraction in the estimator becomes negative. The bootstrap is
    % therefore disabled for Scenario B, and R is transferred from
    % Scenario A in Section 5. R is a property of the gauge, independent of
    % the process-noise model, which is what makes the transfer consistent.
    opts_B.sage_bootstrap = 0;
else
    mode_str = 'Standard Kalman';
    fprintf('Mode: standard Kalman (MultiStart calibration)\n\n');
    fopt = optimoptions('fmincon','Algorithm','interior-point','Display','off', ...
                        'MaxIterations',500,'OptimalityTolerance',1e-9);

    % Scenario A: joint calibration of log10(q) and log10(R)
    fprintf('  Calibrating Scenario A ...\n');
    objA = @(x) nis_obj2D(x, t, y_meas, alpha, k_vec, opts_A, 1e-13);
    xA   = multistart_2D(objA, [0,-6], [6,1], fopt);
    opts_A.Q_diag = 10^xA(1);  opts_A.R = 10^xA(2);
    fprintf('  A: q=%.3e  R=%.3e\n', 10^xA(1), 10^xA(2));

    % Scenario B: Qd^(B) fixed by the PWNE parameters, calibration of R
    fprintf('  Calibrating Scenario B (R only) ...\n');
    objB = @(r) nis_obj1D(r, t, y_meas, alpha, k_vec, opts_B);
    rB   = multistart_1D(objB, -6, 1, fopt);
    opts_B.R = 10^rB;
    fprintf('  B: R=%.3e\n\n', 10^rB);
end

%% 5. Filtering of the full record
% =========================================================================
% Scenario A is run first; in adaptive mode its bootstrap estimate R_boot
% initializes the measurement variance of Scenario B.
% =========================================================================
fprintf('Running Scenario A ...\n');
out_A = discharge_kalman(t', y_meas, alpha, k_vec, opts_A);

R_boot_A = NaN;
if USE_ADAPTIVE && isfield(out_A,'R_boot') && ~isnan(out_A.R_boot) ...
                && out_A.R_boot > sage_p.sage_R_min + eps
    R_boot_A = out_A.R_boot;
    opts_B.R = R_boot_A;
    fprintf('Scenario B initialised with R_boot_A = %.4e\n', R_boot_A);
end

fprintf('Running Scenario B ...\n');
out_B = discharge_kalman(t', y_meas, alpha, k_vec, opts_B);
fprintf('\n');

%% 6. Metrics
% =========================================================================
% Two pointwise errors are distinguished throughout:
%   filtered  (post-update) : q_filt(n) - y(n)       uses y(n)
%   predicted (pre-update)  : q_pred(n) - y(n) = -innov(n)
% Since out.innov(n) = y(n) - c'*x_{n|n-1} is the one-step-ahead error,
% q_pred = y - innov, and no further filter output is needed. The gap
% between the two quantifies the contribution of the measurement update.
% =========================================================================

% Validity masks: full record, and record after the burn-in years
v_all  = ~isnan(y_meas) & ~isnan(out_A.NIS) & ~isnan(out_B.NIS);
n_burn = sum(yr_days(1:min(BURN_IN_YEARS, n_years)));
v_conv = v_all;  v_conv(1:n_burn) = false;

M_A = full_metrics(out_A, y_meas, v_all,  REL_ERR_FLOOR_Q);
M_B = full_metrics(out_B, y_meas, v_all,  REL_ERR_FLOOR_Q);
C_A = full_metrics(out_A, y_meas, v_conv, REL_ERR_FLOOR_Q);
C_B = full_metrics(out_B, y_meas, v_conv, REL_ERR_FLOOR_Q);

clim_std = std(y_meas(v_all));

% Pointwise error series, used by the figures below
err_A         = out_A.q_filt - y_meas;
err_B         = out_B.q_filt - y_meas;
relerr_A      = abs(err_A) ./ max(y_meas, REL_ERR_FLOOR_Q);
relerr_B      = abs(err_B) ./ max(y_meas, REL_ERR_FLOOR_Q);
pred_err_A    = -out_A.innov;
pred_err_B    = -out_B.innov;
pred_relerr_A = abs(pred_err_A) ./ max(y_meas, REL_ERR_FLOOR_Q);
pred_relerr_B = abs(pred_err_B) ./ max(y_meas, REL_ERR_FLOOR_Q);

% -------------------------------------------------------------------------
% Console report
% -------------------------------------------------------------------------
hr = repmat('=',1,64);
fprintf('%s\n', hr);
fprintf('  %s -- Camastra %d-%d   (N_valid = %d of %d)\n', ...
        mode_str, yr_start, yr_end, sum(v_all), N);
fprintf('%s\n', repmat('-',1,64));
fprintf('  %-36s | %10s | %10s\n', 'Metric', 'Scen. A', 'Scen. B');
fprintf('%s\n', repmat('-',1,64));
fprintf('  ACCURACY\n');
fprintf('  %-36s | %10.4f | %10.4f\n', 'RMSE, filtered  [m^3/s]',    M_A.rmse_f,     M_B.rmse_f);
fprintf('  %-36s | %10.4f | %10.4f\n', 'RMSE, predicted [m^3/s]',    M_A.rmse_p,     M_B.rmse_p);
fprintf('  %-36s | %10.4f | %10.4f\n', 'MAE,  filtered  [m^3/s]',    M_A.mae_f,      M_B.mae_f);
fprintf('  %-36s | %10.3f | %10.3f\n', 'NSE,  filtered  (target 1)', M_A.nse_f,      M_B.nse_f);
fprintf('  %-36s | %10.3f | %10.3f\n', 'NSE,  predicted (target 1)', M_A.nse_p,      M_B.nse_p);
fprintf('  %-36s | %10.3f | %10.3f\n', 'Rel. err, median (filtered)',M_A.rel_med_f,  M_B.rel_med_f);
fprintf('  %-36s | %10.3f | %10.3f\n', 'Rel. err, mean   (filtered)',M_A.rel_mean_f, M_B.rel_mean_f);
fprintf('  %-36s | %10.3f | %10.3f\n', 'Rel. err, P90    (filtered)',M_A.rel_p90_f,  M_B.rel_p90_f);
fprintf('  %-36s | %10.3f | %10.3f\n', 'Rel. err, median (predicted)',M_A.rel_med_p, M_B.rel_med_p);
fprintf('  %-36s | %10.3f | %10.3f\n', 'Rel. err, P90    (predicted)',M_A.rel_p90_p, M_B.rel_p90_p);
fprintf('  CALIBRATION\n');
fprintf('  %-36s | %10.3f | %10.3f\n', 'NIS, mean   (target 1)',     M_A.nis_mean,   M_B.nis_mean);
fprintf('  %-36s | %10.3f | %10.3f\n', 'NIS, median (target 0.455)', M_A.nis_med,    M_B.nis_med);
fprintf('  %-36s | %10.1f | %10.1f\n', 'NIS, max',                   M_A.nis_max,    M_B.nis_max);
fprintf('  %-36s | %9.1f%% | %9.1f%%\n','Coverage, nominal 50%',        100*M_A.cov50, 100*M_B.cov50);
fprintf('  %-36s | %9.1f%% | %9.1f%%\n','Coverage, nominal 90%',        100*M_A.cov90, 100*M_B.cov90);
fprintf('  %-36s | %9.1f%% | %9.1f%%\n','Coverage, nominal 95%',        100*M_A.cov95, 100*M_B.cov95);
fprintf('  %-36s | %9.1f%% | %9.1f%%\n','Coverage, nominal 99%',        100*M_A.cov99, 100*M_B.cov99);
fprintf('  %-36s | %9.1f%% | %9.1f%%\n','  [identity check, 95%]',      100*M_A.cov95_alt, 100*M_B.cov95_alt);
fprintf('  %-36s | %9.1f%% | %9.1f%%\n','  [state-sigma 95%, mis-scaled]',100*M_A.cov95_state,100*M_B.cov95_state);
fprintf('  %-36s | %10.2f | %10.2f\n', '  [state-sigma inflation factor]',M_A.state_inflation,M_B.state_inflation);
fprintf('  %-36s | %9.1f%% | %9.1f%%\n','Pr{NIS > 3.84}    (target 5%)', M_A.f95,     M_B.f95);
fprintf('  %-36s | %9.1f%% | %9.1f%%\n','Pr{NIS > 6.63}    (target 1%)', M_A.f99,     M_B.f99);
fprintf('  %-36s | %10.3f | %10.3f\n', 'PIT spread ratio (target 1)',M_A.pit_spread,M_B.pit_spread);
fprintf('%s\n', repmat('-',1,64));
fprintf('  Climatological std of the record : %.4f m^3/s\n', clim_std);
if ~isnan(R_boot_A)
    fprintf('  R_boot transferred A -> B       : %.4e (m^3/s)^2\n', R_boot_A);
end
fprintf('  Update gain, Scenario B  : RMSE %.4f -> %.4f  (%.1f%% reduction)\n', ...
        M_B.rmse_p, M_B.rmse_f, 100*(1 - M_B.rmse_f/M_B.rmse_p));
fprintf('  Update gain, Scenario B  : NSE  %.3f -> %.3f\n', M_B.nse_p, M_B.nse_f);
fprintf('  Update gain, Scenario A  : RMSE %.4f -> %.4f  (%.1f%% reduction)\n', ...
        M_A.rmse_p, M_A.rmse_f, 100*(1 - M_A.rmse_f/M_A.rmse_p));
fprintf('%s\n\n', hr);

fprintf('Post-burn-in (first %d years excluded, N_valid = %d):\n', ...
        BURN_IN_YEARS, sum(v_conv));
fprintf('  RMSE, filtered   A: %8.4f   B: %8.4f  [m^3/s]\n', C_A.rmse_f, C_B.rmse_f);
fprintf('  NSE,  filtered   A: %8.3f   B: %8.3f\n',          C_A.nse_f,  C_B.nse_f);
fprintf('  NIS,  mean       A: %8.3f   B: %8.3f\n',          C_A.nis_mean, C_B.nis_mean);
fprintf('  Coverage 95%%      A: %7.1f%%   B: %7.1f%%\n\n',   100*C_A.cov95, 100*C_B.cov95);

% -------------------------------------------------------------------------
% Regime-conditional metrics, Scenario B
% -------------------------------------------------------------------------
% The record is not a homogeneous population. With lambda_m = 0 from July
% to September, the dry-season GSN response reduces to a deterministic
% exponential recession: innovations are very small and S_n is bounded
% below by R_n, so NIS -> 0 for structural reasons. Pooling these days
% with the storm-driven wet season lowers the NIS median well below the
% chi^2_1 value 0.455. Conditioning on season and on flow level separates
% this mixing effect from the behaviour during high-flow events.
% -------------------------------------------------------------------------
y_p50 = pctl_(y_meas(v_all), 50);
y_p90 = pctl_(y_meas(v_all), 90);

reg_names = {sprintf('Wet (%s)', month_list_str(WET_MONTHS, month_names)), ...
             sprintf('Dry (%s)', month_list_str(DRY_MONTHS, month_names)), ...
             'Low flow (< P50)', 'High flow (> P90)'};
reg_masks = {v_all &  ismember(month_d, WET_MONTHS), ...
             v_all &  ismember(month_d, DRY_MONTHS), ...
             v_all & (y_meas <= y_p50), ...
             v_all & (y_meas >  y_p90)};

reg = struct('name', reg_names, 'M', cell(1,numel(reg_names)));
fprintf('Regime-conditional accuracy and calibration, Scenario B:\n');
fprintf('  %-22s | %6s | %9s | %9s | %9s | %9s\n', ...
        'Regime','N','NIS mean','NIS med','Cov. 95%','RMSE');
fprintf('  %s\n', repmat('-',1,78));
for r = 1:numel(reg_names)
    Mr = full_metrics(out_B, y_meas, reg_masks{r}, REL_ERR_FLOOR_Q);
    reg(r).M = Mr;
    fprintf('  %-22s | %6d | %9.3f | %9.3f | %8.1f%% | %9.3f\n', ...
            reg_names{r}, Mr.n, Mr.nis_mean, Mr.nis_med, ...
            100*Mr.cov95, Mr.rmse_f);
end
fprintf('\n');

% -------------------------------------------------------------------------
% Ten largest NIS values, Scenario B
% -------------------------------------------------------------------------
% Dating the largest NIS values shows whether they are scattered over the
% record or concentrated on storm-response days.
% -------------------------------------------------------------------------
NIS_B_masked = out_B.NIS;  NIS_B_masked(~v_all) = NaN;
[~, iSort]   = sort(NIS_B_masked, 'descend', 'MissingPlacement','last');
n_top        = min(10, sum(v_all));
top_idx      = iSort(1:n_top);

top_tbl = table(dates(top_idx)', out_B.NIS(top_idx)', y_meas(top_idx)', ...
                out_B.q_filt(top_idx)', ...
                out_B.q_filt(top_idx)' - y_meas(top_idx)', ...
                'VariableNames', {'Date','NIS','q_obs','q_filt','error'});
fprintf('Ten largest NIS values, Scenario B:\n');
disp(top_tbl);

%% 7. Figures
% =========================================================================
% Year tick labels, computed once and reused by every figure. yr_ticks
% and yr_labels are already paired (position, label); they must be sliced
% together, otherwise ticks and years fall out of step.
% =========================================================================
tick_stride = 2;
tick_idx    = 1:tick_stride:n_years;
yr_ticks    = cum_days(tick_idx + 1)';
yr_labels   = arrayfun(@(y) sprintf('%d',y), years(tick_idx), 'UniformOutput',false);

% -------------------------------------------------------------------------
% Fig 1: full record, three stacked panels sharing the x-axis
% -------------------------------------------------------------------------
% Observations and the two filtered series are kept in separate panels:
% at this point density an overlay is unreadable, while the shared x-axis
% still allows a comparison at the same time index. Over 37 years this
% figure shows the envelope; individual events are shown in fig02.
% -------------------------------------------------------------------------
fig1 = figure('Name','Fig1 - Full time series (data, Scen. A, Scen. B)');
tl1  = tiledlayout(fig1, 3, 1, 'TileSpacing','compact', 'Padding','compact');

% Top panel: observed discharge
nexttile(tl1); hold on; box on;
plot(t, y_meas, 'Color', clr_obs, 'LineWidth', 0.7, 'DisplayName', 'Observed');
set(gca,'XTick',yr_ticks,'XTickLabel',yr_labels,'FontSize',9,'XTickLabelRotation',45);
legend('Location','northeast','FontSize',8);
ylabel('q  [m^3/s]');
title('Observed discharge','FontSize',11,'FontWeight','bold');
ylim([-1, max(y_meas)*1.05]);

% Middle and bottom panels: filtered discharge with +-BAND_Z*sigma band
ts_panels = {out_A, clr_A, M_A.rmse_f, 'A'; out_B, clr_B, M_B.rmse_f, 'B'};
for sc = 1:2
    [out, clr, rmse_sc, sname] = ts_panels{sc,:};
    nexttile(tl1); hold on; box on;
    fill_band(t, out.q_filt, out.q_sigma, clr, BAND_Z);
    plot(t, out.q_filt, 'Color', clr, 'LineWidth', 1.1, ...
         'DisplayName', sprintf('Filtered, RMSE = %.3f m^3/s', rmse_sc));
    set(gca,'XTick',yr_ticks,'XTickLabel',yr_labels,'FontSize',9,'XTickLabelRotation',45);
    legend('Location','northeast','FontSize',8);
    ylabel('q  [m^3/s]');
    title(sprintf('Scenario %s - filtered', sname),'FontSize',11,'FontWeight','bold');
    ylim([-1, max(y_meas)*1.05]);
end
xlabel(tl1, 'Year', 'FontSize', 11);
title(tl1, sprintf('Camastra %d-%d  [%s]', yr_start, yr_end, mode_str), ...
      'FontSize', 13, 'FontWeight','bold');
FIGS(end+1,:) = {fig1, 'fig01_timeseries_data_scenA_scenB'};

% -------------------------------------------------------------------------
% Fig 2: daily hydrograph for representative years
% -------------------------------------------------------------------------
% Shows the filter tracking the hydrograph at daily resolution, and the
% uncertainty band narrowing on recession limbs and widening at peaks.
% -------------------------------------------------------------------------
fig2 = figure('Name','Fig2 - Hydrograph zoom, representative years');
tl2  = tiledlayout(fig2, numel(ZOOM_YEARS), 1, 'TileSpacing','compact', ...
                   'Padding','compact');
for iz = 1:numel(ZOOM_YEARS)
    yz  = ZOOM_YEARS(iz);
    idx = find(year_d == yz);
    if isempty(idx)
        warning('Zoom year %d is not in the record; skipped.', yz);
        continue;
    end
    nexttile(tl2); hold on; box on;
    % Plotted against a numeric day-of-year index rather than a datetime
    % axis, since fill with datetime x-data is not supported uniformly
    % across releases. Month boundaries are set as explicit ticks below.
    dz = 1:numel(idx);
    fill_band(dz, out_B.q_filt(idx), out_B.q_sigma(idx), clr_B, BAND_Z);
    plot(dz, y_meas(idx),       '-',  'Color',clr_obs, 'LineWidth',1.3, ...
         'DisplayName','Observed');
    plot(dz, out_B.q_filt(idx), '-',  'Color',clr_B,   'LineWidth',1.5, ...
         'DisplayName',sprintf('Scen. B - filtered (\\pm%.2f\\sigma)', BAND_Z));
    plot(dz, out_A.q_filt(idx), '--', 'Color',clr_A,   'LineWidth',1.2, ...
         'DisplayName','Scen. A - filtered');
    % First day of each month within this year
    m_of_day = month_d(idx);
    m_ticks  = find([true, diff(m_of_day) ~= 0]);
    set(gca,'XTick',m_ticks,'XTickLabel',month_names(m_of_day(m_ticks)),'FontSize',9);
    ylabel('q  [m^3/s]','FontSize',10);
    title(sprintf('%d', yz),'FontSize',11,'FontWeight','bold');
    if iz == 1, legend('Location','northwest','FontSize',8,'Box','off'); end
    xlim([dz(1) dz(end)]);
end
title(tl2, sprintf('Daily hydrograph, representative years  [%s]', mode_str), ...
      'FontSize',13,'FontWeight','bold');
FIGS(end+1,:) = {fig2, 'fig02_hydrograph_zoom'};

% -------------------------------------------------------------------------
% Fig 3: pointwise relative error, both scenarios, pre and post update
% -------------------------------------------------------------------------
% Rows are scenarios, columns are the estimate after (filtered) and before
% (predicted) the measurement update.
%
% Rows have independent y-ranges, since Scenario A errors are about an
% order of magnitude larger; the two columns of a row share the same
% range, which is where the pre/post comparison is made.
%
% The error |q - y|/max(y, REL_ERR_FLOOR_Q) is a dimensionless ratio on a
% logarithmic axis; the floor prevents near-zero summer flows from
% producing disproportionately large values.
%
% Each panel reports the median, the mean and the 90th percentile. Median
% and mean differ by about an order of magnitude: the error distribution
% is sharply peaked with a heavy tail, so the median describes the
% typical day while the mean reflects the large errors that drive the
% RMSE.
% -------------------------------------------------------------------------
relerr_floor_plot = min(min([relerr_A;  pred_relerr_A; ...
                    relerr_B; pred_relerr_B]));   

rel_rows = {relerr_A, pred_relerr_A, clr_A_filt, clr_A_pred, 'A'; ...
            relerr_B, pred_relerr_B, clr_B_filt, clr_B_pred, 'B'};

fig3 = figure('Name','Fig3 - Pointwise relative error, both scenarios');
tl3  = tiledlayout(fig3, 2, 2, 'TileSpacing','compact', 'Padding','compact');

for sc = 1:2
    [re_f, re_p, c_f, c_p, sname] = rel_rows{sc,:};

    % Summary statistics over the valid days
    stat_f = struct('med', median(re_f(v_all)), 'avg', mean(re_f(v_all)), ...
                    'p90', pctl_(re_f(v_all), 90));
    stat_p = struct('med', median(re_p(v_all)), 'avg', mean(re_p(v_all)), ...
                    'p90', pctl_(re_p(v_all), 90));

    % Common y-range for the two columns of this row
    all_re = [re_f(v_all), re_p(v_all)];
    y_lim  = [relerr_floor_plot, max(pctl_(all_re, 99.5)*2, 0.1)];

    row_panels = {max(re_f, relerr_floor_plot), c_f, '.',  8, stat_f, 'Filtered (post-update)'; ...
                  max(re_p, relerr_floor_plot), c_p, 'x', 10, stat_p, 'Predicted (pre-update)'};

    for pp = 1:2
        [dat, clr, mk, msz, st, ttl] = row_panels{pp,:};
        nexttile(tl3); hold on; box on;
        if strcmp(mk,'.')
            scatter(t(v_all), dat(v_all), msz, clr, mk, ...
                    'MarkerEdgeAlpha',0.55, 'HandleVisibility','off');
        else
            scatter(t(v_all), dat(v_all), msz, clr, mk, 'LineWidth',0.8, ...
                    'MarkerEdgeAlpha',0.55, 'HandleVisibility','off');
        end
        % Reference lines: median, mean, 90th percentile
        yline(st.med, '--', 'Color','k',      'LineWidth',1.5, ...
              'DisplayName', sprintf('Median = %.3f', st.med));
        yline(st.avg, '-.', 'Color',clr_mean, 'LineWidth',1.6, ...
              'DisplayName', sprintf('Mean = %.3f',   st.avg));
        yline(st.p90, ':',  'Color',clr_p90,  'LineWidth',1.5, ...
              'DisplayName', sprintf('P90 = %.3f',    st.p90));
        set(gca,'YScale','log','XTick',yr_ticks,'XTickLabel',yr_labels, ...
                'FontSize',9,'XTickLabelRotation',45,'LineWidth',0.9);
        legend('Location','southwest','FontSize',8,'Box','off');
        if pp == 1
            ylabel(sprintf('Scenario %s\nrel. error  (ratio)', sname), ...
                   'FontSize',10,'FontWeight','bold');
        end
        if sc == 1, title(ttl,'FontSize',12,'FontWeight','bold'); end
        ylim(y_lim);
    end
end

xlabel(tl3, 'Year', 'FontSize', 11);
title(tl3, sprintf('Pointwise relative error |q-y|/max(y,1)  [%s]', mode_str), ...
      'FontSize',13,'FontWeight','bold');
FIGS(end+1,:) = {fig3, 'fig03_relerror_scenA_scenB'};

% -------------------------------------------------------------------------
% Fig S1: distributional diagnostics  [PLOT_CALIB_DIAG only]
% -------------------------------------------------------------------------
% Three standard tools of probabilistic forecast verification, each
% answering a distinct question:
%
%   (a) PIT histogram. p_n = Phi(eps_n/sqrt(S_n)) is Uniform(0,1) for a
%       consistent filter, so the reference is a flat histogram. Mass
%       concentrated near 0.5 indicates an over-dispersed filter (band
%       wider than the realized error); mass pushed towards 0 and 1
%       indicates an over-confident one.
%
%   (b) Reliability diagram. Empirical against nominal coverage of the
%       predictive band q_pred +- z_p*sqrt(S_n). Points above the
%       diagonal are conservative, points below are over-confident; the
%       vertical distance measures the miscalibration at each level.
%
%   (c) NIS QQ-plot against chi^2_1 on log-log axes. It locates the
%       departure from chi^2_1 along the quantile range, and shows the
%       bulk and the tail of the distribution departing in opposite
%       directions.
% -------------------------------------------------------------------------
if PLOT_CALIB_DIAG
    pit_A = pit_signed(out_A.innov(v_all), out_A.S(v_all));
    pit_B = pit_signed(out_B.innov(v_all), out_B.S(v_all));

    figS1 = figure('Name','Fig S1 - Calibration diagnostics');
    tlS1 = tiledlayout(figS1, 1, 3, 'TileSpacing','compact', 'Padding','compact');

    % (a) PIT histogram
    nexttile(tlS1); hold on; box on;
    edges_pit = linspace(0,1,21);
    histogram(pit_A, edges_pit, 'Normalization','pdf', 'FaceColor',clr_A, ...
              'FaceAlpha',0.45, 'EdgeColor','none', ...
              'DisplayName', sprintf('Scenario A (spread %.2f)', M_A.pit_spread));
    histogram(pit_B, edges_pit, 'Normalization','pdf', 'FaceColor',clr_B, ...
              'FaceAlpha',0.55, 'EdgeColor','none', ...
              'DisplayName', sprintf('Scenario B (spread %.2f)', M_B.pit_spread));
    yline(1, 'k-', 'LineWidth',2.0, 'DisplayName','U(0,1)  (calibrated)');
    xlabel('PIT   p_n = \Phi(\epsilon_n / \surdS_n)'); ylabel('Density');
    title('(a) PIT histogram','FontSize',12,'FontWeight','bold');
    legend('Location','north','FontSize',8,'Box','off');
    set(gca,'FontSize',9); xlim([0 1]);

    % (b) Reliability diagram
    nexttile(tlS1); hold on; box on;
    p_nom = [0.50 0.60 0.70 0.80 0.90 0.95 0.975 0.99];
    plot([0.4 1],[0.4 1],'k-','LineWidth',1.6,'DisplayName','Perfect calibration');
    plot(p_nom, coverage_curve(out_A, y_meas, v_all, p_nom), '-o', 'Color',clr_A, ...
         'LineWidth',1.8,'MarkerFaceColor',clr_A,'MarkerSize',6,'DisplayName','Scenario A');
    plot(p_nom, coverage_curve(out_B, y_meas, v_all, p_nom), '-s', 'Color',clr_B, ...
         'LineWidth',1.8,'MarkerFaceColor',clr_B,'MarkerSize',6,'DisplayName','Scenario B');
    xlabel('Nominal coverage'); ylabel('Empirical predictive coverage');
    title('(b) Reliability diagram (predictive band)', ...
          'FontSize',12,'FontWeight','bold');
    legend('Location','southeast','FontSize',8,'Box','off');
    set(gca,'FontSize',9); axis([0.45 1 0.3 1.02]); grid on;

    % (c) NIS QQ-plot against chi^2_1
    nexttile(tlS1); hold on; box on;
    qq_pairs = {out_A.NIS(v_all), clr_A, 'Scenario A'; ...
                out_B.NIS(v_all), clr_B, 'Scenario B'};
    for sc = 1:2
        [Nv, clr, nm] = qq_pairs{sc,:};
        Nv = sort(Nv(~isnan(Nv)));
        nn = numel(Nv);
        q_theo = chi2inv_1dof(((1:nn) - 0.5) / nn);
        ok = q_theo > 0 & Nv > 0;
        plot(q_theo(ok), Nv(ok), '.', 'Color',clr, 'MarkerSize',6, 'DisplayName',nm);
    end
    lim_qq = [1e-5 1e4];
    plot(lim_qq, lim_qq, 'k-', 'LineWidth',1.6, 'DisplayName','1:1  (\chi^2_1)');
    set(gca,'XScale','log','YScale','log','FontSize',9);
    xlabel('\chi^2_1 theoretical quantile'); ylabel('Empirical NIS quantile');
    title('(c) NIS QQ-plot','FontSize',12,'FontWeight','bold');
    legend('Location','southeast','FontSize',8,'Box','off');
    axis([lim_qq lim_qq]);

    title(tlS1, sprintf('Filter calibration diagnostics - %s', mode_str), ...
          'FontSize',13,'FontWeight','bold');
    FIGS(end+1,:) = {figS1, 'figS1_calibration_diagnostics'};
end


% -------------------------------------------------------------------------
% Fig 4: monthly RMSE and bias, filtered and predicted
% -------------------------------------------------------------------------
% All years are pooled by calendar month. The filtered bias is a weak
% check of the seasonal adequacy of the model, because the update pulls
% q_filt towards y at every step and damps any systematic mismatch. The
% predicted (pre-update) bias, not yet informed by the current
% observation, is the more informative check; it is the model-side
% analogue of the monthly first-moment comparison of Morlando et al.
% (2016). Both are shown side by side.
% -------------------------------------------------------------------------
[rmse_mA_f, bias_mA_f] = monthly_err_stats(err_A,      month_d);
[rmse_mB_f, bias_mB_f] = monthly_err_stats(err_B,      month_d);
[rmse_mA_p, bias_mA_p] = monthly_err_stats(pred_err_A, month_d);
[rmse_mB_p, bias_mB_p] = monthly_err_stats(pred_err_B, month_d);

fig4 = figure('Name','Fig4 - Monthly RMSE and bias (filtered vs. predicted)');
tl4  = tiledlayout(fig4, 2, 2, 'TileSpacing','compact', 'Padding','compact');

% Top row: monthly RMSE, filtered (left) and predicted (right)
rmse_panels = {rmse_mA_f, rmse_mB_f, 'Filtered (post-update)', true; ...
               rmse_mA_p, rmse_mB_p, 'Predicted (pre-update)', false};
for pp = 1:2
    [rA, rB, ttl, do_ylab] = rmse_panels{pp,:};
    nexttile(tl4); hold on; box on;
    bh = bar(1:12, [rA; rB]', 'grouped');
    bh(1).FaceColor = clr_A;  bh(1).FaceAlpha = 0.85;
    bh(2).FaceColor = clr_B;  bh(2).FaceAlpha = 0.85;
    set(gca,'XTick',1:12,'XTickLabel',month_names,'FontSize',9);
    legend({'Scen. A','Scen. B'},'FontSize',8,'Location','northwest');
    if do_ylab, ylabel('RMSE  [m^3/s]'); end
    title(ttl,'FontSize',11,'FontWeight','bold');
end

% Bottom row: monthly bias, filtered (left) and predicted (right)
x_ = 1:12;  w_ = 0.35;
bias_panels = {bias_mA_f, bias_mB_f, true; bias_mA_p, bias_mB_p, false};
for pp = 1:2
    [bA, bB, do_ylab] = bias_panels{pp,:};
    nexttile(tl4); hold on; box on;
    bar(x_-w_/2, bA, w_, 'FaceColor',clr_A,'FaceAlpha',0.85,'DisplayName','Scen. A');
    bar(x_+w_/2, bB, w_, 'FaceColor',clr_B,'FaceAlpha',0.85,'DisplayName','Scen. B');
    yline(0,'k-','LineWidth',1.1,'HandleVisibility','off');
    set(gca,'XTick',1:12,'XTickLabel',month_names,'FontSize',9);
    legend('FontSize',8,'Location','southeast');
    xlabel('Month');
    if do_ylab, ylabel('Bias [m^3/s]  (+: overestimate)'); end
end

title(tl4, sprintf('Monthly RMSE and bias - %s', mode_str), ...
      'FontSize',13,'FontWeight','bold');
FIGS(end+1,:) = {fig4, 'fig04_monthly_RMSE_bias'};

% -------------------------------------------------------------------------
% Fig 5: online adaptation of R and Qd  [adaptive mode only]
% -------------------------------------------------------------------------
% Stacked panels sharing the year axis, logarithmic ordinates since both
% quantities span several decades. Values at days without an observation
% are carried forward from the last adaptive step. The dashed line marks
% the bootstrap estimate R_boot,A transferred to Scenario B.
% -------------------------------------------------------------------------
if USE_ADAPTIVE && isfield(out_A,'R_hist') && isfield(out_B,'R_hist') ...
                && isfield(out_A,'Qd_hist') && isfield(out_B,'Qd_hist')
    fig5 = figure('Name','Fig5 - Online adaptation of R and Q_d');
    tl5  = tiledlayout(fig5, 2, 1, 'TileSpacing','compact', 'Padding','compact');

    % (a) Measurement-noise variance R_n
    nexttile(tl5); hold on; box on;
    plot(t, fillmissing(out_A.R_hist,'previous'), 'Color',clr_A, ...
         'LineWidth',1.2, 'DisplayName','Scenario A');
    plot(t, fillmissing(out_B.R_hist,'previous'), 'Color',clr_B, ...
         'LineWidth',1.2, 'DisplayName','Scenario B');
    if ~isnan(R_boot_A)
        yline(R_boot_A,'k--','LineWidth',1.2, ...
              'DisplayName',sprintf('R_{boot,A} = %.3e', R_boot_A));
    end
    set(gca,'YScale','log','XTick',yr_ticks,'XTickLabel',yr_labels, ...
            'FontSize',9,'XTickLabelRotation',45);
    legend('FontSize',8,'Location','best');
    ylabel('R  [(m^3/s)^2]');
    title('(a) Measurement-noise variance R_n','FontSize',11,'FontWeight','bold');

    % (b) Process-noise diagonal Qd_n
    nexttile(tl5); hold on; box on;
    plot(t, fillmissing(out_A.Qd_hist,'previous'), 'Color',clr_A, ...
         'LineWidth',1.2, 'DisplayName','Scenario A');
    plot(t, fillmissing(out_B.Qd_hist,'previous'), 'Color',clr_B, ...
         'LineWidth',1.2, 'DisplayName','Scenario B');
    set(gca,'YScale','log','XTick',yr_ticks,'XTickLabel',yr_labels, ...
            'FontSize',9,'XTickLabelRotation',45);
    legend('FontSize',8,'Location','best');
    ylabel('Q_d  [(m^3/s)^2]');
    title('(b) Process-noise diagonal Q_{d,n}','FontSize',11,'FontWeight','bold');

    xlabel(tl5, 'Year', 'FontSize', 11);
    title(tl5, sprintf('Online noise-covariance adaptation - Scenario A vs B  [%s]', ...
          mode_str), 'FontSize',13,'FontWeight','bold');
    FIGS(end+1,:) = {fig5, 'fig05_adaptive_R_Qd'};
end

%% 8. Export
% =========================================================================
% LaTeX tables to the Command Window and/or to file, and figure export
% =========================================================================
if PRINT_TABLES
    % Printed as ready-to-paste LaTeX source. Set TABLES_TO_FILE = true
    % to also write them as \input-able fragments.
    fprintf('\n');
    fprintf('%% ===== COPY FROM HERE INTO THE MANUSCRIPT =====================\n\n');
    write_results_table(1, M_A, M_B, C_A, C_B, clim_std, BURN_IN_YEARS, ...
                        mode_str, yr_start, yr_end, sum(v_all));
    fprintf('\n');
    write_regime_table(1, reg);
    fprintf('\n');
    write_events_table(1, top_tbl);
    fprintf('\n%% ===== COPY TO HERE ===========================================\n\n');
end

if TABLES_TO_FILE
    if ~isfolder(FIG_OUTPUT_DIR), mkdir(FIG_OUTPUT_DIR); end
    tbl_specs = {'tab_results.tex', @(f) write_results_table(f, M_A, M_B, C_A, C_B, ...
                     clim_std, BURN_IN_YEARS, mode_str, yr_start, yr_end, sum(v_all)); ...
                 'tab_regime.tex',  @(f) write_regime_table(f, reg); ...
                 'tab_events.tex',  @(f) write_events_table(f, top_tbl)};
    for i = 1:size(tbl_specs,1)
        fn  = fullfile(FIG_OUTPUT_DIR, tbl_specs{i,1});
        fid = fopen(fn, 'w');
        assert(fid > 0, 'Cannot open %s for writing.', fn);
        tbl_specs{i,2}(fid);
        fclose(fid);
        fprintf('  %s\n', fn);
    end
end

if EXPORT_FIGURES
    fprintf('Exporting figures to %s ...\n', FIG_OUTPUT_DIR);
    export_named_figures(FIG_OUTPUT_DIR, FIG_FORMATS, FIGS);
end
fprintf('Done.\n');


% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================
% MATLAB requires the local functions of a script to be defined at the
% end of the file.
% =========================================================================

% =========================================================================
% Accuracy and consistency metrics
% =========================================================================
% Metrics of one scenario over one validity mask: RMSE, MAE, bias and
% Nash-Sutcliffe efficiency of the filtered and predicted estimates,
% relative-error statistics, NIS statistics, predictive-band coverage at
% four nominal levels, and PIT spread ratio.
%
% NSE = 1 - MSE/var(y) over the evaluation subset, so NSE > 0 means that
% the estimator improves on the climatological mean. Computing it for both
% the filtered and the predicted estimate expresses the contribution of
% the measurement update as a gain in skill.
% =========================================================================
function M = full_metrics(out, y_ref, mask, q_floor)
    m  = mask & ~isnan(y_ref) & ~isnan(out.NIS);
    yv = y_ref(m);
    qf = out.q_filt(m);
    qp = y_ref(m) - out.innov(m);      % q_pred(n) = y(n) - innov(n)
    sg = out.q_sigma(m);
    Nv = out.NIS(m);
    Sn = out.S(m);

    % Per-step measurement variance, needed for the posterior band.
    % Adaptive mode stores it in R_hist (NaN where no observation was
    % assimilated, filled here); standard mode uses the fixed value.
    if isfield(out,'R_hist')
        Rn_all = fillmissing(fillmissing(out.R_hist,'previous'),'next');
    else
        Rn_all = out.filter_params.R * ones(size(out.S));
    end
    Rn = Rn_all(m);

    var_clim = var(yv);

    % ---------------------------------------------------------------------
    % Accuracy
    % ---------------------------------------------------------------------
    M.n      = sum(m);
    M.rmse_f = sqrt(mean((qf - yv).^2));
    M.rmse_p = sqrt(mean((qp - yv).^2));
    M.mae_f  = mean(abs(qf - yv));
    M.bias_f = mean(qf - yv);
    M.nse_f  = 1 - mean((qf - yv).^2) / var_clim;
    M.nse_p  = 1 - mean((qp - yv).^2) / var_clim;

    % Relative error with floored denominator
    den = max(yv, q_floor);
    rf  = abs(qf - yv) ./ den;
    rp  = abs(qp - yv) ./ den;
    M.rel_med_f  = median(rf);  M.rel_mean_f = mean(rf);  M.rel_p90_f = pctl_(rf,90);
    M.rel_med_p  = median(rp);  M.rel_mean_p = mean(rp);  M.rel_p90_p = pctl_(rp,90);

    % ---------------------------------------------------------------------
    % NIS statistics
    % ---------------------------------------------------------------------
    M.nis_mean = mean(Nv);
    M.nis_med  = median(Nv);
    M.nis_max  = max(Nv);
    M.f95      = 100*mean(Nv > 3.8415);    % chi^2_1 95% critical value
    M.f99      = 100*mean(Nv > 6.6349);    % chi^2_1 99% critical value

    % ---------------------------------------------------------------------
    % Uncertainty-band coverage
    % ---------------------------------------------------------------------
    % For a scalar observation there is a single coverage. The predictive
    % band y_n in q_pred +- z*sqrt(S_n) is equivalent to NIS_n <= z^2,
    % since S_n is the predictive variance of the observation. The
    % posterior band uses the scale R_n/sqrt(S_n), and defines the same
    % event: delta_n = (R_n/S_n)*eps_n exactly, so
    %     |delta_n| <= z*R_n/sqrt(S_n)  <=>  eps_n^2/S_n <= z^2.
    % Coverage is therefore reported at four nominal levels.
    M.cov50 = mean(Nv <= 0.6745^2);
    M.cov90 = mean(Nv <= 1.6449^2);
    M.cov95 = mean(Nv <= 1.9600^2);
    M.cov99 = mean(Nv <= 2.5758^2);

    % Numerical check of the identity above; agrees to rounding
    sd_post    = Rn ./ sqrt(Sn);
    M.cov95_alt = mean(abs(yv - qf) <= 1.9600*sd_post);

    % Coverage built on the state standard deviation q_sigma =
    % sqrt(c'*P_{n|n}*c). This is the uncertainty of the state estimate,
    % not the predictive spread of an observation: it exceeds the
    % posterior scale by the factor sqrt(S_n*(S_n-R_n))/R_n, so the
    % resulting coverage is conservative by construction. Reported for
    % reference only.
    M.cov95_state     = mean(abs(yv - qf) <= 1.9600*sg);
    M.state_inflation = median(sqrt(Sn .* max(Sn - Rn, 0)) ./ Rn);
    M.R_med  = median(Rn);
    M.S_med  = median(Sn);

    % ---------------------------------------------------------------------
    % PIT spread ratio
    % ---------------------------------------------------------------------
    % Standard deviation of the signed PIT divided by its Uniform(0,1)
    % value 1/sqrt(12). Equal to 1 for a calibrated filter, below 1 for an
    % over-dispersed filter (PIT concentrated near 0.5), above 1 for an
    % over-confident one (PIT pushed towards both ends).
    pit          = pit_signed(out.innov(m), out.S(m));
    M.pit_spread = std(pit) / (1/sqrt(12));
    M.pit_mean   = mean(pit);
end


% =========================================================================
% Signed probability integral transform
% =========================================================================
% p_n = Phi(eps_n/sqrt(S_n)), Uniform(0,1) for a consistent filter,
% computed through erf (no Statistics Toolbox needed).
%
% The signed innovation is used instead of NIS: folding maps both an
% over-dispersed and an over-confident filter onto mass at one end of
% [0,1], so that the two failures could not be distinguished. With the
% signed transform, mass near 0.5 means over-dispersed and mass at both
% ends means over-confident.
% =========================================================================
function p = pit_signed(innov, S)
    ok = ~isnan(innov) & S > 0;
    p  = 0.5 * (1 + erf(innov(ok) ./ sqrt(2*S(ok))));
end


% =========================================================================
% Inverse chi-square CDF with one degree of freedom
% =========================================================================
% Inverting erf(sqrt(x/2)) = p gives x = 2*erfinv(p)^2. Toolbox-free
% equivalent of chi2inv(p,1).
% =========================================================================
function q = chi2inv_1dof(p)
    q = 2 * erfinv(min(max(p, 0), 1 - eps)).^2;
end


% =========================================================================
% Empirical predictive coverage curve
% =========================================================================
% Fraction of observations inside q_pred +- z_p*sqrt(S_n) at each nominal
% two-sided level p_nom, computed as P(NIS_n <= z_p^2) with
% z_p = sqrt(2)*erfinv(p). Toolbox-free.
% =========================================================================
function cv = coverage_curve(out, y_ref, mask, p_nom)
    m  = mask & ~isnan(y_ref) & ~isnan(out.NIS);
    Nv = out.NIS(m);
    cv = arrayfun(@(p) mean(Nv <= 2*erfinv(p).^2), p_nom);
end


% =========================================================================
% Percentile
% =========================================================================
% Linear-interpolation percentile following the prctile convention
% (sample points at (i-0.5)/n, clamped outside that range). Implemented
% locally so that no Statistics Toolbox is needed. NaN values are ignored.
% =========================================================================
function q = pctl_(x, p)
    x = sort(x(~isnan(x)));
    n = numel(x);
    if n == 0, q = NaN; return; end
    if n == 1, q = x*ones(size(p)); return; end
    pos = 100 * ((1:n) - 0.5) / n;
    q   = interp1(pos, x, p, 'linear');
    q(p <  pos(1))   = x(1);
    q(p >  pos(end)) = x(end);
end


% =========================================================================
% Monthly error statistics
% =========================================================================
% Per-calendar-month RMSE and bias, pooling all years. Used for both
% filtered errors (q_filt - y) and predicted errors (q_pred - y = -innov).
% Months with fewer than five valid days are left at zero.
% =========================================================================
function [rmse_m, bias_m] = monthly_err_stats(err, month_vec)
    rmse_m = zeros(1,12); bias_m = zeros(1,12);
    for m = 1:12
        idx = (month_vec == m) & ~isnan(err);
        if sum(idx) < 5, continue; end
        e = err(idx);
        rmse_m(m) = sqrt(mean(e.^2));
        bias_m(m) = mean(e);
    end
end


% =========================================================================
% Shaded uncertainty band
% =========================================================================
% Draws the band qf +- z*qs as a transparent patch, excluded from the
% legend. z is passed in so that figures and captions stay consistent.
% =========================================================================
function fill_band(tx, qf, qs, clr, z)
    fill([tx(:); flipud(tx(:))], ...
         [qf(:) + z*qs(:); flipud(qf(:) - z*qs(:))], ...
         clr, 'FaceAlpha',0.18, 'EdgeColor','none', 'HandleVisibility','off');
end


% =========================================================================
% Month range label
% =========================================================================
% Compact label such as "Nov-Apr", built from the first and last month of
% the set.
% =========================================================================
function s = month_list_str(months, names)
    if isempty(months), s = ''; return; end
    if isscalar(months), s = names{months}; return; end
    s = sprintf('%s-%s', names{months(1)}, names{months(end)});
end


% =========================================================================
% LaTeX table writers
% =========================================================================
% The tables are generated from the same variables used by the console
% report. Each writer prints to the file identifier fid; fid = 1 is the
% Command Window.
% =========================================================================

% =========================================================================
% Table of results: accuracy, consistency and coverage, both scenarios
% =========================================================================
% Reports the main accuracy metrics, the mean NIS and the predictive-band
% coverage at the 50% and 95% nominal levels, over the full record and
% after the burn-in years. Further diagnostics (90% and 99% coverage,
% chi^2_1 exceedance fractions, PIT spread ratio) are available in the
% console report and in figS1.
% =========================================================================
function write_results_table(fid, M_A, M_B, C_A, C_B, clim_std, nburn, ...
                             mode_str, y0, y1, nval)
    fprintf(fid, '%% Auto-generated by Main_GSN_ADPT_Kalman_Camastra.m\n');
    fprintf(fid, '\\begin{table}[h!]\n\\centering\n');
    fprintf(fid, '\\renewcommand{\\arraystretch}{1.15}\n');
    fprintf(fid, '\\begin{tabular}{lrrc}\n\\toprule\n');
    fprintf(fid, 'Metric & Scenario~A & Scenario~B & Target \\\\\n\\midrule\n');
    fprintf(fid, '\\multicolumn{4}{l}{\\emph{Accuracy}} \\\\\n');
    fprintf(fid, '\\quad RMSE, filtered [m$^3$/s]    & %.3f & %.3f & --- \\\\\n', M_A.rmse_f, M_B.rmse_f);
    fprintf(fid, '\\quad RMSE, predicted [m$^3$/s]   & %.3f & %.3f & --- \\\\\n', M_A.rmse_p, M_B.rmse_p);
    fprintf(fid, '\\quad MAE, filtered [m$^3$/s]     & %.3f & %.3f & --- \\\\\n', M_A.mae_f, M_B.mae_f);
    fprintf(fid, '\\quad NSE, filtered               & $%.3f$ & $%.3f$ & $1$ \\\\\n', M_A.nse_f, M_B.nse_f);
    fprintf(fid, '\\quad NSE, predicted              & $%.3f$ & $%.3f$ & $1$ \\\\\n', M_A.nse_p, M_B.nse_p);
    fprintf(fid, '\\quad Rel.\\ error, median (filt.) & %.3f & %.3f & $0$ \\\\\n', M_A.rel_med_f, M_B.rel_med_f);
    fprintf(fid, '\\quad Rel.\\ error, mean (filt.)   & %.3f & %.3f & $0$ \\\\\n', M_A.rel_mean_f, M_B.rel_mean_f);
    fprintf(fid, '\\quad Rel.\\ error, P90 (filt.)    & %.3f & %.3f & $0$ \\\\\n', M_A.rel_p90_f, M_B.rel_p90_f);
    fprintf(fid, '\\quad Rel.\\ error, median (pred.) & %.3f & %.3f & $0$ \\\\\n', M_A.rel_med_p, M_B.rel_med_p);
    fprintf(fid, '\\addlinespace\n\\multicolumn{4}{l}{\\emph{Consistency and band coverage}} \\\\\n');
    fprintf(fid, '\\quad $\\NIS$, mean                & %.3f & %.3f & $1$ \\\\\n', M_A.nis_mean, M_B.nis_mean);
    % Predictive-band coverage at the 50% and 95% nominal levels
    fprintf(fid, '\\quad Coverage, nominal $50\\%%$     & %.1f\\%% & %.1f\\%% & $50\\%%$ \\\\\n', 100*M_A.cov50, 100*M_B.cov50);
    fprintf(fid, '\\quad Coverage, nominal $95\\%%$     & %.1f\\%% & %.1f\\%% & $95\\%%$ \\\\\n', 100*M_A.cov95, 100*M_B.cov95);
    % Post-burn-in block
    fprintf(fid, '\\addlinespace\n\\multicolumn{4}{l}{\\emph{Post-burn-in (first %d years excluded)}} \\\\\n', nburn);
    fprintf(fid, '\\quad RMSE, filtered [m$^3$/s]    & %.3f & %.3f & --- \\\\\n', C_A.rmse_f, C_B.rmse_f);
    fprintf(fid, '\\quad NSE, filtered               & $%.3f$ & $%.3f$ & $1$ \\\\\n', C_A.nse_f, C_B.nse_f);
    fprintf(fid, '\\quad $\\NIS$, mean                & %.3f & %.3f & $1$ \\\\\n', C_A.nis_mean, C_B.nis_mean);
    fprintf(fid, '\\quad Coverage, nominal $95\\%%$     & %.1f\\%% & %.1f\\%% & $95\\%%$ \\\\\n', 100*C_A.cov95, 100*C_B.cov95);
    fprintf(fid, '\\midrule\n');
    fprintf(fid, '\\multicolumn{4}{l}{Climatological standard deviation of the record: $%.3f$ m$^3$/s} \\\\\n', clim_std);
    fprintf(fid, '\\bottomrule\n\\end{tabular}\n');
    fprintf(fid, '\\caption{Accuracy and calibration of the %s adaptive filter on the Camastra daily discharge record (%d--%d, $N_{\\mathrm{valid}}=%d$). ', mode_str, y0, y1, nval);
    fprintf(fid, 'Filtered quantities are post-update ($\\hat q_{n|n}$); predicted quantities are the one-step-ahead forecast before the update ($\\hat q_{n|n-1}$), so the difference between the two isolates the contribution of the measurement update. ');
    fprintf(fid, 'NSE denotes the Nash--Sutcliffe efficiency, defined in~\\eqref{eq:NSE}. The $\\NIS$ entry is the empirical mean of~\\eqref{eq:NIS}, the first-moment consistency diagnostic established by Proposition~\\ref{prop:Q_consistent}. Coverage is the fraction of observations inside $\\hat q_{n|n-1}\\pm z_p\\sqrt{S_n}$, equivalently $\\mathbb{P}(\\NIS_n\\le z_p^2)$ by~\\eqref{eq:NIS}. By the shrinkage identity~\\eqref{eq:q_shrink} the posterior band $\\hat q_{n|n}\\pm z_p R_n/\\sqrt{S_n}$, whose scale follows from Theorem~\\ref{prop:R_unbiased}, defines exactly the same event, so a single coverage is reported. A single coverage is therefore reported at each of the two nominal levels. ');
    fprintf(fid, 'Scenario~A is run with a statistically calibrated but event-desynchronised rainfall input and therefore functions as a stress test rather than as an assessment of the deterministic-input formulation; see Section~\\ref{subsec:design}.}\n');
    fprintf(fid, '\\label{tab:results}\n\\end{table}\n');
end


% =========================================================================
% Table of regime-conditional metrics, Scenario B
% =========================================================================
function write_regime_table(fid, reg)
    fprintf(fid, '%% Auto-generated by Main_GSN_ADPT_Kalman_Camastra.m\n');
    fprintf(fid, '\\begin{table}[h!]\n\\centering\n\\renewcommand{\\arraystretch}{1.15}\n');
    fprintf(fid, '\\begin{tabular}{lrrrrr}\n\\toprule\n');
    fprintf(fid, 'Regime & $N$ & $\\NIS$ (mean) & $\\NIS$ (median) & Coverage $95\\%%$ & RMSE [m$^3$/s] \\\\\n\\midrule\n');
    for r = 1:numel(reg)
        % Escape < and > for LaTeX math mode
        lbl = strrep(strrep(reg(r).name, '<', '$<$'), '>', '$>$');
        fprintf(fid, '%s & %d & %.3f & %.3f & %.1f\\%% & %.3f \\\\\n', ...
                lbl, reg(r).M.n, reg(r).M.nis_mean, reg(r).M.nis_med, ...
                100*reg(r).M.cov95, reg(r).M.rmse_f);
    end
    fprintf(fid, '\\bottomrule\n\\end{tabular}\n');
    fprintf(fid, '\\caption{Regime-conditional accuracy and calibration under Scenario~\\ref{Scenario_B}. ');
    fprintf(fid, 'The seasonal split reflects the vanishing of $\\lambda_m$ from July to September in the PWNE calibration (Table~\\ref{tab:monthly_pwne}); the flow-based split uses the median and the ninetieth percentile of the observed record. ');
    fprintf(fid, 'The near-unit aggregate $\\NIS$ of Table~\\ref{tab:results} is shown here to be an average over regimes of very different innovation energy rather than evidence of uniform agreement with $\\chi^2_1$.}\n');
    fprintf(fid, '\\label{tab:regime}\n\\end{table}\n');
end


% =========================================================================
% Table of the ten largest NIS values, Scenario B
% =========================================================================
function write_events_table(fid, T)
    fprintf(fid, '%% Auto-generated by Main_GSN_ADPT_Kalman_Camastra.m\n');
    fprintf(fid, '\\begin{table}[h!]\n\\centering\n\\renewcommand{\\arraystretch}{1.15}\n');
    fprintf(fid, '\\begin{tabular}{lrrrr}\n\\toprule\n');
    fprintf(fid, 'Date & $\\NIS_n$ & $q_n^{\\mathrm{obs}}$ & $\\hat q_{n|n}$ & Error \\\\\n');
    fprintf(fid, ' & & [m$^3$/s] & [m$^3$/s] & [m$^3$/s] \\\\\n\\midrule\n');
    for i = 1:height(T)
        fprintf(fid, '%s & %.1f & %.2f & %.2f & $%+.2f$ \\\\\n', ...
                datestr(T.Date(i), 'dd mmm yyyy'), T.NIS(i), ...
                T.q_obs(i), T.q_filt(i), T.error(i));  %#ok<DATST>
    end
    fprintf(fid, '\\bottomrule\n\\end{tabular}\n');
    fprintf(fid, '\\caption{The ten largest normalised innovation squared values under Scenario~\\ref{Scenario_B}, with the corresponding observed and filtered discharge. Dating the outlying values, rather than reporting only their maximum, shows both that every one of them is an underestimate of a rising limb and where in the record they fall.}\n');
    fprintf(fid, '\\label{tab:events}\n\\end{table}\n');
end


% =========================================================================
% Figure export
% =========================================================================
% Writes every collected figure under its paired filename, in each of the
% requested formats. Uses exportgraphics (R2020a or later, vector output)
% and falls back to print/saveas on older releases.
% =========================================================================
function export_named_figures(outdir, formats, FIGS)
    if ~isfolder(outdir), mkdir(outdir); end
    for i = 1:size(FIGS,1)
        h = FIGS{i,1};  base = FIGS{i,2};
        if ~isgraphics(h, 'figure'), continue; end
        for f = 1:numel(formats)
            fn = fullfile(outdir, sprintf('%s.%s', base, formats{f}));
            try
                exportgraphics(h, fn, 'ContentType','vector', ...
                               'BackgroundColor','white');
            catch
                % Fallback for releases without exportgraphics
                switch lower(formats{f})
                    case 'pdf', print(h, fn, '-dpdf',  '-vector', '-bestfit');
                    case 'eps', print(h, fn, '-depsc', '-vector');
                    case 'png', print(h, fn, '-dpng',  '-r300');
                    otherwise,  saveas(h, fn);
                end
            end
            fprintf('  %s\n', fn);
        end
    end
end


% =========================================================================
% Offline calibration helpers  [USE_ADAPTIVE = false only]
% =========================================================================
% Objective: drive the mean NIS to one, (mean(NIS) - 1)^2, plus a very
% small penalty on mean(log S) that only breaks ties among parameter sets
% already satisfying the NIS condition. Parameters are searched in log10
% scale; runs with fewer than 10 valid innovations are rejected.
% =========================================================================

% Scenario A: x = [log10(q), log10(R)]
function J = nis_obj2D(x, t, y, alpha, k, opts, lam)
    opts.Q_diag = 10^x(1)*ones(1,numel(k));
    opts.R      = 10^x(2);
    opts.adaptive = false;
    ws = warning('off','all');
    out = discharge_kalman(t, y, alpha, k, opts);
    warning(ws);
    S = out.S(:); e = out.innov(:); v = ~isnan(e) & S > 0;
    if sum(v) < 10, J = 1e12; return; end
    J = lam*mean(log(S(v))) + (mean(e(v).^2./S(v)) - 1)^2;
end


% Scenario B: lr = log10(R), Qd^(B) fixed
function J = nis_obj1D(lr, t, y, alpha, k, opts)
    opts.R = 10^lr;
    opts.adaptive = false;
    ws = warning('off','all');
    out = discharge_kalman(t, y, alpha, k, opts);
    warning(ws);
    S = out.S(:); e = out.innov(:); v = ~isnan(e) & S > 0;
    if sum(v) < 10, J = 1e12; return; end
    J = 1e-13*mean(log(S(v))) + (mean(e(v).^2./S(v)) - 1)^2;
end


% MultiStart fmincon in 2D: 60 points along the diagonal of the box plus
% 80 uniformly random points
function [x_opt, fval] = multistart_2D(obj, lb, ub, fopt)
    n  = 60;
    rp = [linspace(lb(1),ub(1),n)', linspace(lb(2),ub(2),n)'];
    rp = [rp; lb(1)+(ub(1)-lb(1))*rand(80,1), lb(2)+(ub(2)-lb(2))*rand(80,1)];
    spts = CustomStartPointSet(unique(rp,'rows'));
    prob = createOptimProblem('fmincon','objective',obj,'x0',[0.5,-1], ...
                              'lb',lb,'ub',ub,'options',fopt);
    ms = MultiStart('Display','off','UseParallel',true,'StartPointsToRun','bounds');
    [x_opt, fval] = run(ms, prob, spts);
end


% MultiStart fmincon in 1D: 200 equally spaced start points
function [x_opt, fval] = multistart_1D(obj, lb, ub, fopt)
    spts = CustomStartPointSet(linspace(lb,ub,200)');
    prob = createOptimProblem('fmincon','objective',obj,'x0',0, ...
                              'lb',lb,'ub',ub,'options',fopt);
    ms = MultiStart('Display','off','UseParallel',true,'StartPointsToRun','bounds');
    [x_opt, fval] = run(ms, prob, spts);
end


% =========================================================================
% Struct merge
% =========================================================================
% Field-wise merge, with the fields of s2 overriding those of s1.
% =========================================================================
function s = merge_s(s1, s2)
    s = s1;
    fn = fieldnames(s2);
    for k = 1:numel(fn), s.(fn{k}) = s2.(fn{k}); end
end
