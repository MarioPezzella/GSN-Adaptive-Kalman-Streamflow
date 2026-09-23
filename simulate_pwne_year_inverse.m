function z = simulate_pwne_year_inverse(lambda, eta_c, year_days)
% =========================================================================
% SIMULATE_PWNE_YEAR_INVERSE
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
% =========================================================================
% This function generates one year of daily effective rainfall from the
% monthly Poisson White Noise Exponential (PWNE) model of the Generalized
% Shot Noise framework.
%
% Within each day of month m, effective rainfall events arrive as a
% Poisson process of rate lambda(m) and carry independent exponential
% magnitudes of mean eta_c(m). The daily total is therefore a compound
% Poisson-exponential random variable, whose CDF has an atom at zero,
%
%       F(0) = exp(-lambda*dt)                      (probability of a dry day)
%
% and, for z > 0,
%
%       F(z) = exp(-lambda*dt)
%              + sum_{nu>=1} Pois(nu; lambda*dt) * P(nu, z/eta_c),
%
% where P(nu,.) is the regularized lower incomplete gamma function, i.e.
% the CDF of the sum of nu exponential magnitudes.
%
% Each daily value is drawn by inverse-transform sampling: a uniform
% number r ~ U(0,1) is generated; if r <= F(0) the day is dry, otherwise
% the equation F(z) = r is solved numerically with fzero.
%
% The function is called once per year by Main_GSN_ADPT_Kalman_Camastra.m
% to build the synthetic effective-rainfall input of Scenario A. The
% resulting sequence reproduces the calibrated monthly statistics of the
% PWNE process, but its individual events are not synchronized with those
% that generated the observed discharge.
%
% INPUTS:
%   lambda    - [12x1] or [1x12] monthly Poisson rates [1/day]
%   eta_c     - [12x1] or [1x12] monthly mean pulse magnitudes. Values of
%               months with lambda = 0 are never used.
%   year_days - Number of days of the year, 365 (default) or 366. With
%               366 days February is given 29 days.
%
% OUTPUT:
%   z         - [year_days x 1] simulated daily effective rainfall
% =========================================================================

    if nargin < 3
        year_days = 365;
    end

    % =====================================================================
    % Input checks
    % =====================================================================
    lambda = lambda(:)';
    eta_c  = eta_c(:)';
    if length(lambda) ~= 12 || length(eta_c) ~= 12
        error('lambda and eta_c must contain 12 values.');
    end

    % =====================================================================
    % Days per month (leap years handled through year_days)
    % =====================================================================
    days_in_month = [31 28 31 30 31 30 31 31 30 31 30 31];
    if year_days == 366
        days_in_month(2) = 29;
    end
    if sum(days_in_month) ~= year_days
        error('Inconsistent number of days.');
    end

    % =====================================================================
    % Initialization
    % =====================================================================
    z = zeros(year_days,1);
    dt = 1;                         % daily time step [days]
    day_idx = 1;

    % =====================================================================
    % Month-by-month, day-by-day sampling
    % =====================================================================
    for m = 1:12
        lam = lambda(m);
        eta = eta_c(m);
        p0 = exp(-lam*dt);          % probability of no event in one day
        for d = 1:days_in_month(m)
            % -------------------------------------------------------------
            % Uniform random number
            % -------------------------------------------------------------
            r = rand;
            % -------------------------------------------------------------
            % Dry day: r falls inside the atom at zero. For lambda = 0,
            % p0 = 1 and every day of the month is dry.
            % -------------------------------------------------------------
            if r <= p0
                z(day_idx) = 0;
            else
                % ---------------------------------------------------------
                % Wet day: solve F(z) = r
                % ---------------------------------------------------------
                F = @(x) compound_poisson_cdf(x,lam,eta,dt) - r;
                % Initial bracketing interval
                zmax = 100*eta;
                % Double the upper end until the root is bracketed
                while F(zmax) < 0
                    zmax = 2*zmax;
                end
                % Numerical inversion of the CDF
                z(day_idx) = fzero(F,[0,zmax]);

            end
            day_idx = day_idx + 1;
        end
    end
end

% =========================================================================
% Compound Poisson-exponential CDF
% =========================================================================
% Evaluates F(z) for the daily total of a Poisson(lambda*dt) number of
% independent exponential magnitudes of mean eta. The series over the
% number of events nu is truncated when its terms fall below tol, with a
% hard stop at nu = 500.
% =========================================================================
function val = compound_poisson_cdf(z,lambda,eta,dt)
    if z <= 0
        val = exp(-lambda*dt);
        return;
    end
    % Mean number of events in one time step
    mu = lambda*dt;
    % Zero-event term (atom at z = 0)
    val = exp(-mu);
    % Series truncation
    tol = 1e-12;
    term = inf;
    nu = 1;
    while term > tol
        % Poisson probability of exactly nu events
        w = exp(-mu) * mu^nu / factorial(nu);
        % CDF of the sum of nu exponentials: Gamma(nu, eta), i.e. the
        % regularized lower incomplete gamma function
        P = gammainc(z/eta,nu,'lower');
        term = w * P;
        val = val + term;
        nu = nu + 1;
        % Safety stop
        if nu > 500
            break;
        end
    end
end
