function out = discharge_kalman(t, y_meas, alpha, k_vec, opts)
% =========================================================================
% DISCHARGE_KALMAN
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
% This function applies a Kalman filter, standard or Sage-Husa adaptive,
% to a daily discharge record described by the Generalized Shot Noise
% (GSN) model (Murrone et al. 1997; Morlando et al. 2016; Cimorelli et
% al. 2021), and optionally produces a pure-prediction forecast beyond
% the end of the record.
%
% The number of linear reservoirs ns = numel(k_vec) is inferred from the
% input. For the PCSTC calibration of the Camastra basin, ns = 1.
%
% -------------------------------------------------------------------------
% PHYSICAL MODEL
% -------------------------------------------------------------------------
%   Unit response:  h(t) = a0*delta(t) + sum_j (alpha_j/k_j)*exp(-t/k_j),
%                   with sum(alpha) = 1
%   State j:        x_j(t) = int_0^t exp(-(t-s)/k_j) z(s) ds
%   Stationarity:   E[x_j] = k_j*E[z] = k_j*E[q]     (used for x0)
%   Discharge:      q(t) = c'*x(t) + a0*z(t),  c_j = alpha(j+1)/k_vec(j)
%
% -------------------------------------------------------------------------
% DISCRETE STATE-SPACE FORM (exact, time step dt)
% -------------------------------------------------------------------------
%   x_{n+1} = Ad*x_n + u_n + w_n,   w_n ~ N(0,Qd)
%   y_n     = c'*x_n + v_n,         v_n ~ N(0,R)
%
%   Ad_jj   = exp(-dt/k_j)                              exact transition
%
%   SCENARIO A - effective rainfall z_hat_d(n) supplied as a known input
%     u_n(j)  = z_hat_d(n)*k_j*(1 - exp(-dt/k_j))/dt    rectangular pulse
%     Qd_jj   = q_j*k_j/2*(1 - exp(-2*dt/k_j))          diagonal Qd^(A)
%     The instantaneous fraction a0 bypasses the reservoirs, so it is
%     removed from the data before filtering and restored afterwards:
%       y~_n = y_n - a0*z_hat_d(n)
%       q^_n = c'*x_{n|n} + a0*z_hat_d(n)
%
%   SCENARIO B - effective rainfall unknown, absorbed as process noise
%     u_n     = 0
%     Qd_ij   = G*k_i*k_j/(k_i+k_j)*(1 - exp(-(1/k_i+1/k_j)*dt))
%     with the PWNE spectral constant G = 2*eta_c^2*lambda  (Qd^(B)).
%
% -------------------------------------------------------------------------
% FILTER LOOP: PREDICT -> UPDATE
% -------------------------------------------------------------------------
%   At the start of step n, x_cur holds x_{n-1|n-1}.
%     PREDICT:  x_{n|n-1} = Ad*x_{n-1|n-1} + u_n
%     UPDATE:   x_{n|n}   = x_{n|n-1} + K_n*eps_n
%   With this ordering a run can be warm-started exactly: passing
%   x0 = x_{N|N} makes the first prediction compute Ad*x_{N|N} + u_{N+1}.
%   Missing observations (NaN) trigger the prediction step only.
%
% -------------------------------------------------------------------------
% STANDARD MODE  (opts.adaptive = false, default)
% -------------------------------------------------------------------------
%   Q and R are fixed. The steady-state gain K_inf is computed from the
%   discrete algebraic Riccati equation (DARE), and the filter switches to
%   it permanently once ||K_n - K_inf||/||K_inf|| < gain_tol holds for
%   gain_patience consecutive steps.
%
% -------------------------------------------------------------------------
% SAGE-HUSA ADAPTIVE MODE  (opts.adaptive = true)
% -------------------------------------------------------------------------
%   Q and R are estimated online, following the structure of Algorithm 1
%   in Li et al. (Appl. Sci. 2025, 15, 1731), with the adaptations for
%   the GSN model listed further below. At each step with a valid
%   observation, in order:
%
%   PRE-PREDICT
%   (I)    Save the scalars of P_{n-1|n-1}, overwritten by the prediction.
%   (II)   Forgetting weight d1 from the prior-innovation buffer buf_e:
%            L1 = -mean(buf_e)
%            b1 = lambda_min + (1 - lambda_min)*2^L1
%            d1 = (1 - b1)/(1 - b1^n_adapt)
%   (III)  Update R with the POSTERIOR residual delta_{n-1}:
%            R_n = (1-d1)*R_{n-1} + d1*(delta_{n-1}^2 + c'*P_{n-1|n-1}*c)
%          Both terms are non-negative, so R stays positive.
%
%   PREDICT
%   (IV)   x_{n|n-1} = Ad*x_{n-1|n-1} + u_n
%          P_{n|n-1} = Ad*P_{n-1|n-1}*Ad' + Qd     (full Riccati step)
%
%   UPDATE
%   (V)    Prior innovation and its variance:
%            eps_n = y~_n - c'*x_{n|n-1},   S_n = c'*P_{n|n-1}*c + R_n
%   (VI)   Gain and posterior, Joseph stabilized form:
%            K_n     = P_{n|n-1}*c/S_n
%            P_{n|n} = (I-K_n*c')*P_{n|n-1}*(I-K_n*c')' + K_n*R_n*K_n'
%            x_{n|n} = x_{n|n-1} + K_n*eps_n
%   (VII)  Posterior residual:
%            delta_n = y~_n - c'*x_{n|n}   (used by the next R update)
%
%   POST-UPDATE
%   (VIII) Forgetting weight d2 from the posterior-residual buffer buf_d,
%          evaluated before delta_n is appended (causal):
%            L2 = -mean(buf_d),  b2 and d2 as in (II).
%   (IX)   Update of the Qd diagonal, with propagation correction:
%            Qd_jj = (1-d2)*Qd_jj
%                    + d2*(K_j^2*eps_n^2 + P_{n|n,jj} - Ad_jj^2*P_{n-1|n-1,jj})
%          The PRIOR innovation eps_n enters, as in the term K*e*e'*K' of
%          the Sage-Husa process-noise estimator.
%   (X)    Append eps_n^2/S_n to buf_e and delta_n^2/S_n to buf_d.
%
% -------------------------------------------------------------------------
% ADAPTATIONS FOR THE GSN MODEL
% -------------------------------------------------------------------------
%   1. NIS-normalized forgetting weights (scale invariance).
%      With L = -mean(e^2), as for residuals of order one, discharge
%      residuals in m^3/s give 2^L ~ 0, hence b ~ lambda_min at every step
%      and no adaptivity. Both buffers therefore store NIS-normalized
%      values, L = -mean(e^2/S), which is about -1 for a consistent filter.
%
%   2. Propagation correction in the Q update (slow dynamics).
%      Without the term -Ad_jj^2*P_{n-1|n-1,jj}, the covariance propagated
%      through the state dynamics is counted as new process noise. For
%      k = 38.4 d (Ad_jj = 0.974), P_ss,jj ~ Qd_jj/(1 - Ad_jj^2) ~ 20*Qd_jj
%      and the recursion has no finite stationary point. With it, using
%      the DARE identity P_ss,jj*(1 - Ad_jj^2) = Qd_jj - K_j^2*S, the
%      steady-state increment is K_j^2*S*(NIS - 1), which vanishes for a
%      self-consistent filter (NIS = 1).
%
%   3. No covariance matching.
%      Scaling R by min(NIS,alpha_max) whenever eps_n^2 > S_n inflates R
%      geometrically for a consistent filter, since P(NIS > 1) ~ 31.7%,
%      and drives K towards 0. This step is not used; the posterior-
%      residual R update (III) already corrects a systematic
%      underestimation of R.
%
%   4. Bootstrap R estimate and equilibrium initialization.
%      When R0 >> R_true, K ~ 0 and delta ~ eps, so R0 is a fixed point of
%      the R update. A bootstrap of sage_bootstrap standard Kalman steps
%      is run first, and
%        R_boot = max(mean(eps^2) - mean(c'*P_pred*c), R_min),
%      which is unbiased for any gain because E[eps^2] = c'*P_pred*c + R.
%      The filter then starts from x0 = mean(y)*k_vec and P0 = P_ss(R_boot),
%      so that K_1 = K_inf and no cold-start transient occurs. Both buffers
%      are pre-filled with ones (NIS = 1), which gives a small, finite
%      forgetting weight from the first step.
%
% -------------------------------------------------------------------------
% INPUTS:
%   t        - [1xN] observation times [days], uniform, dt = t(2) - t(1)
%   y_meas   - [1xN] observed discharge [m^3/s], NaN for missing days
%   alpha    - [1x(ns+1)] routing fractions [a0, a1, ..., a_ns], all >= 0,
%              summing to 1
%   k_vec    - [1xns] reservoir characteristic times [days], all > 0
%   opts     - (optional) struct of options, all fields optional:
%
%     General
%     .scenario          'A' (default) or 'B'
%     .z_hat_d   [1xN]   effective rainfall input (Scenario A)
%     .Q_diag    [1xns]  spectral densities q_j of Qd^(A) (Scenario A,
%                        default R*0.01*ones)
%     .lambda, .eta_c    PWNE parameters defining G (Scenario B, required)
%     .R         scalar  initial (adaptive) or fixed (standard) measurement
%                        noise variance (default 5% of var(y_meas))
%     .x0        [nsx1]  initial state (overrides the default)
%     .P0        [nsxns] initial covariance (overrides the default and
%                        disables the bootstrap)
%     .gain_tol          steady-state gain tolerance (default 1e-3)
%     .gain_patience     steps below tolerance before switching to K_inf
%                        (default 3)
%
%     Sage-Husa adaptation (.adaptive = true)
%     .adaptive          false (default) | true
%     .sage_bootstrap    bootstrap steps for the R estimate (default 100;
%                        0 disables it)
%     .sage_lambda_min   lower bound lambda_min of b (default 0.95)
%     .sage_window       width of both NIS buffers (default 30)
%     .sage_R_min/max    bounds on R (defaults 1e-8 / Inf)
%     .sage_Q_min/max    bounds on the Qd diagonal (defaults 1e-10 / Inf)
%
%     Forecast (active when .t_fore is non-empty)
%     .t_fore    [1xF]   forecast times
%     .fore_option       drift: 0 = free decay, 1 = z_hat_d_fore,
%                        2 = constant PWNE mean, 3 = seasonal PWNE mean
%     .z_hat_d_fore [1xF]                               (option 1)
%     .lambda_fore, .eta_c_fore                         (option 2)
%     .lambda_monthly [1x12], .eta_c_monthly [1x12],
%     .months_fore [1xF]                                (option 3)
%
% OUTPUT (struct):
%   .q_filt    [1xN]      filtered discharge q^_n = c'*x_{n|n} (+ a0*z_hat_d)
%   .q_sigma   [1xN]      state standard deviation sqrt(c'*P_{n|n}*c)
%   .x_filt    [nsxN]     filtered states x_{n|n}
%   .P_filt    [nsxnsxN]  posterior covariances P_{n|n}
%   .innov     [1xN]      prior innovations eps_n (NaN if missing)
%   .S         [1xN]      innovation variances S_n
%   .K_hist    [nsxN]     Kalman gains K_n
%   .NIS       [1xN]      normalized innovation squared eps_n^2/S_n
%   .n_ss      scalar     step at which K_inf was activated (NaN = never)
%   .K_ss      [nsx1]     steady-state gain K_inf from the DARE
%   .filter_params        filter settings at the end of the run, retained
%                         for warm-starting a subsequent batch
%   .t_fore, .q_fore, .q_fore_lo, .q_fore_hi, .x_fore, .P_fore
%                         forecast and its 95% band (empty if no forecast)
%   Adaptive mode only:
%   .R_boot    scalar     bootstrap R estimate (NaN if not run)
%   .R_hist    [1xN]      R at each adaptive step (NaN otherwise)
%   .Qd_hist   [nsxN]     Qd diagonal at each adaptive step (NaN otherwise)
% =========================================================================
if nargin < 5, opts = struct(); end

% =========================================================================
% 0. Options
% =========================================================================
scenario      = upper(getopt(opts,'scenario',       'A'));
R             = getopt(opts,'R', max(var(y_meas,'omitnan')*0.05, 1e-8));
x0_in         = getopt(opts,'x0',                   []);
P0_in         = getopt(opts,'P0',                   []);
gain_tol      = getopt(opts,'gain_tol',             1e-3);
gain_patience = getopt(opts,'gain_patience',        3);
t_fore        = getopt(opts,'t_fore',               []);
fore_option   = getopt(opts,'fore_option',          0);

do_adaptive  = getopt(opts,'adaptive',          false);
sh_bootstrap = getopt(opts,'sage_bootstrap',    100);
sh_lam_min   = getopt(opts,'sage_lambda_min',   0.95);
sh_S         = getopt(opts,'sage_window',       30);
sh_R_min     = getopt(opts,'sage_R_min',        1e-8);
sh_R_max     = getopt(opts,'sage_R_max',        Inf);
sh_Q_min     = getopt(opts,'sage_Q_min',        1e-10);
sh_Q_max     = getopt(opts,'sage_Q_max',        Inf);

% =========================================================================
% 1. System matrices
% =========================================================================
N  = numel(t);
dt = t(2) - t(1);
ns = numel(k_vec);   % number of reservoirs, inferred from the model
In = eye(ns);

Ad      = diag(exp(-dt ./ k_vec));   % [nsxns] exact state transition
c       = (alpha(2:end) ./ k_vec)';  % [nsx1]  c_j = alpha(j+1)/k_vec(j)
a0      = alpha(1);                  % instantaneous (feed-through) fraction
Ad_diag = diag(Ad);                  % [nsx1]  for the element-wise Q update

% =========================================================================
% 2. Initial process noise Qd
% =========================================================================
Q_diag_used=[]; lambda_used=[]; eta_c_used=[];
switch scenario
    case 'A'
        % Diagonal Qd^(A) from the spectral densities q_j
        Q_diag_used = getopt(opts,'Q_diag', R*0.01*ones(1,ns));
        Qd = diag(Q_diag_used .* k_vec/2 .* (1-exp(-2*dt./k_vec)));
    case 'B'
        % Full Qd^(B) from the PWNE spectral constant G = 2*eta_c^2*lambda
        lambda_used = getopt(opts,'lambda',[]);
        eta_c_used  = getopt(opts,'eta_c', []);
        if isempty(lambda_used)||isempty(eta_c_used)
            error('discharge_kalman: Scenario B requires opts.lambda and opts.eta_c');
        end
        Qd = build_Qd_B(2*eta_c_used^2*lambda_used, k_vec, dt, ns);
    otherwise
        error('discharge_kalman: opts.scenario must be ''A'' or ''B''');
end

% =========================================================================
% 3. Input sequence and corrected measurements
% =========================================================================
% Scenario A: deterministic input U and removal of the feed-through a0*z.
% Scenario B: no input, measurements used as they are.
if strcmp(scenario,'A')
    z_hat_d = getopt(opts,'z_hat_d', zeros(1,N));
    if numel(z_hat_d)~=N
        error('discharge_kalman: opts.z_hat_d must have length %d', N);
    end
    U      = build_U_rect(z_hat_d, k_vec, dt, ns, N);
    y_corr = y_meas(:)' - a0*z_hat_d(:)';
else
    z_hat_d=[]; U=zeros(ns,N); y_corr=y_meas(:)';
end

% =========================================================================
% 4. Bootstrap pass (adaptive mode only; skipped when P0 is provided)
% =========================================================================
% Adaptation 4: when R0 >> R_true the Sage-Husa R update has a stable
% fixed point at R0, and R never converges to R_true. A short standard
% Kalman pass with the preliminary R provides the unbiased estimate
%       R_boot = max( mean(eps_i^2) - mean(c'*P_pred,i*c), R_min ),
% since E[eps_i^2] = c'*P_pred,i*c + R_true for any gain K. The estimate
% is accepted only if at least 10 valid observations were processed.
% =========================================================================
R_boot_estimate = NaN;

if do_adaptive && isempty(P0_in)
    % Steady-state posterior covariance for the preliminary R
    P_rough = dare_sym(Ad, c, Qd, R);
    S_rough = c'*P_rough*c + R;
    K_rough = P_rough*c / S_rough;
    P_ss_r  = (In-K_rough*c')*P_rough*(In-K_rough*c')' + K_rough*R*K_rough';

    % Start from the stationary mean E[x_j] = k_j*E[q]
    x_boot = mean(y_meas,'omitnan') * k_vec(:);
    P_boot = P_ss_r;
    inn_sq=0; cPc_s=0; nb_valid=0;

    % Standard Kalman steps until sh_bootstrap valid observations are used
    for nb = 1:N
        if nb_valid >= sh_bootstrap, break; end
        xp=Ad*x_boot+U(:,nb); Pp=Ad*P_boot*Ad'+Qd;
        if ~isnan(y_corr(nb))
            inn_nb = y_corr(nb)-c'*xp; cPc_nb=c'*Pp*c;
            inn_sq = inn_sq + inn_nb^2;
            cPc_s  = cPc_s  + cPc_nb;
            S_nb=cPc_nb+R; K_nb=Pp*c/S_nb;
            x_boot=xp+K_nb*inn_nb;
            P_boot=(In-K_nb*c')*Pp*(In-K_nb*c')'+K_nb*R*K_nb';
            nb_valid=nb_valid+1;
        else
            x_boot=xp; P_boot=Pp;
        end
    end

    % Bootstrap estimate, projected onto [R_min, R_max]
    if nb_valid >= 10
        R_raw           = inn_sq/nb_valid - cPc_s/nb_valid;
        R_boot_estimate = min(max(R_raw, sh_R_min), sh_R_max);
        R               = R_boot_estimate;
    end
end

% =========================================================================
% 5. DARE - steady-state gain K_inf
% =========================================================================
% Solves  P_inf = Ad*(P_inf - P_inf*c*(c'*P_inf*c + R)^(-1)*c'*P_inf)*Ad' + Qd
% and sets
%       K_inf = P_inf*c/(c'*P_inf*c + R)
%       P_ss  = (I - K_inf*c')*P_inf*(I - K_inf*c')' + K_inf*R*K_inf'
%
% DARE identity used by Adaptation 2:
%       Ad*P_ss*Ad' + Qd = P_inf
%   =>  P_ss,jj*(1 - Ad_jj^2) = Qd_jj - K_inf,j^2*S_inf
% =========================================================================
P_inf = dare_sym(Ad, c, Qd, R);
S_inf = c'*P_inf*c + R;
K_inf = P_inf*c / S_inf;
P_ss  = (In-K_inf*c')*P_inf*(In-K_inf*c')' + K_inf*R*K_inf';

% =========================================================================
% 6. Initial state and covariance
% =========================================================================
% Adaptive default (Adaptation 4):
%   x0 = mean(y)*k_vec   stationary mean E[x_j] = k_j*E[q]
%   P0 = P_ss            so that K_1 = K_inf and no cold-start transient
% Standard default:
%   x0 = 0,  P0 = 1e4*(Qd + R*I)   diffuse prior
% Both are overridden when opts.x0 / opts.P0 are supplied (warm start).
% =========================================================================
if ~isempty(x0_in)
    if numel(x0_in)~=ns, error('opts.x0 must be [%dx1]',ns); end
    x_cur = x0_in(:);
elseif do_adaptive
    x_cur = mean(y_meas,'omitnan') * k_vec(:);
else
    x_cur = zeros(ns,1);
end

if ~isempty(P0_in)
    P_cur = P0_in;
elseif do_adaptive
    P_cur = P_ss;
else
    P_cur = 1e4*(Qd+R*In);
end

% =========================================================================
% 7. Filter loop (PREDICT -> UPDATE)
% =========================================================================
x_filt=zeros(ns,N); P_filt=zeros(ns,ns,N);
innov =zeros(1,N);  S_out =zeros(1,N);
K_hist=zeros(ns,N);
use_ss=false; n_ss=NaN; below_count=0;

if do_adaptive
    % Pre-fill both NIS buffers with ones (Adaptation 4): mean(buf) = 1
    % gives L = -1 and a small, finite forgetting weight from the first
    % step, instead of the complete replacement d = 1.
    n_adapt = sh_S;
    buf_e   = ones(1,sh_S);   % prior-innovation NIS     -> b1, d1 -> R
    buf_d   = ones(1,sh_S);   % posterior-residual NIS   -> b2, d2 -> Q

    sh_d_prev = NaN;           % delta_{n-1}: NaN at n = 1 (no previous step)
    sh_cPc    = c'*P_cur*c;
    sh_Pd     = diag(P_cur);

    R_hist  = NaN(1,N);
    Qd_hist = NaN(ns,N);
end

for n = 1:N

    valid_obs   = ~isnan(y_corr(n));
    in_adaptive = do_adaptive && valid_obs;

    % ---------------------------------------------------------------------
    % PRE-PREDICT, steps (I)-(III)
    % ---------------------------------------------------------------------
    if in_adaptive
        % (I) Save the P_{n-1|n-1} scalars before the prediction
        sh_cPc = c'*P_cur*c;
        sh_Pd  = diag(P_cur);

        % (II) b1, d1 from the prior-innovation buffer buf_e, which
        %      stores eps_i^2/S_i (Adaptation 1). b1 is clamped to
        %      [lambda_min, 1 - 1e-9] to keep d1 finite.
        L1 = -mean(buf_e);
        b1 = min(max(sh_lam_min + (1-sh_lam_min)*2^L1, sh_lam_min), 1-1e-9);
        d1 = (1-b1) / (1-b1^n_adapt);

        % (III) R update with the posterior residual delta_{n-1}.
        %       Consistent in expectation:
        %       E[delta^2 + c'*P*c] = (R - c'*P*c) + c'*P*c = R.
        %       Skipped when delta_{n-1} is unavailable (first step, or
        %       previous observation missing).
        if ~isnan(sh_d_prev)
            R_new = (1-d1)*R + d1*(sh_d_prev^2 + sh_cPc);
            R = min(max(R_new, sh_R_min), sh_R_max);
        end
    end

    % ---------------------------------------------------------------------
    % PREDICT, step (IV)
    % ---------------------------------------------------------------------
    x_pred = Ad*x_cur + U(:,n);
    if do_adaptive || ~valid_obs || ~use_ss
        P_pred = Ad*P_cur*Ad' + Qd;   % full Riccati step
    else
        P_pred = P_inf;               % steady-state shortcut (standard mode)
    end

    % ---------------------------------------------------------------------
    % UPDATE
    % ---------------------------------------------------------------------
    if ~valid_obs
        % Missing observation: prediction only, no gain
        inn=NaN; S_n=c'*P_pred*c+R; K_n=zeros(ns,1);
        x_cur=x_pred; P_cur=P_pred;
        if do_adaptive, sh_d_prev=NaN; end

    elseif in_adaptive

        % (V) Prior innovation and innovation variance
        inn = y_corr(n) - c'*x_pred;
        S_n = c'*P_pred*c + R;

        % (VI) Kalman gain and posterior covariance (Joseph form)
        K_n   = P_pred*c / S_n;
        IKC   = In - K_n*c';
        P_cur = IKC*P_pred*IKC' + K_n*R*K_n';
        x_cur = x_pred + K_n*inn;

        % (VII) Posterior residual delta_n
        sh_d_n = y_corr(n) - c'*x_cur;

        % (VIII) b2, d2 from the posterior-residual buffer buf_d,
        %        evaluated BEFORE appending delta_n (past values only)
        L2 = -mean(buf_d);
        b2 = min(max(sh_lam_min + (1-sh_lam_min)*2^L2, sh_lam_min), 1-1e-9);
        d2 = (1-b2) / (1-b2^n_adapt);

        % (IX) Q update with propagation correction (Adaptation 2):
        %
        %   Qd_jj = (1-d2)*Qd_jj
        %           + d2*(K_j^2*eps_n^2 + P_{n|n,jj} - Ad_jj^2*P_{n-1|n-1,jj})
        %
        %   At steady state P_{n|n} ~ P_{n-1|n-1} ~ P_ss, and by the DARE
        %   identity of Section 5 the increment reduces to
        %       K_j^2*eps_n^2 + P_ss,jj*(1 - Ad_jj^2) - Qd_jj
        %     = K_j^2*S*(NIS - 1),
        %   which vanishes for a self-consistent filter (NIS = 1).
        %   The result is projected onto [Q_min, Q_max].
        if d2 > 0
             q_new = (1-d2)*diag(Qd) + d2*(K_n.^2 * inn^2 ...
                     + diag(P_cur) - Ad_diag.^2 .* sh_Pd);
            % Uncorrected update, without the propagation term. It has no
            % finite stationary point for slow reservoirs and is kept only
            % for comparison:
            % q_new = (1-d2)*diag(Qd) + d2*(K_n.^2 * inn^2 ...
                    % + diag(P_cur) );
            Qd = diag(min(max(q_new, sh_Q_min), sh_Q_max));
        end

        % (X) Append to both circular buffers, after b1 and b2 have been
        %     computed (causality)
        n_adapt = n_adapt + 1;
        idx     = 1 + mod(n_adapt-1, sh_S);
        buf_e(idx) = inn^2    / S_n;   % prior NIS      -> drives R
        buf_d(idx) = sh_d_n^2 / S_n;   % posterior NIS  -> drives Q

        sh_d_prev    = sh_d_n;
        R_hist(n)    = R;
        Qd_hist(:,n) = diag(Qd);

    elseif use_ss
        % Standard mode, steady-state gain already activated
        inn=y_corr(n)-c'*x_pred; K_n=K_inf; S_n=S_inf;
        P_cur=P_ss; x_cur=x_pred+K_n*inn;

    else
        % Standard mode, time-varying gain; monitor convergence to K_inf
        inn=y_corr(n)-c'*x_pred; S_n=c'*P_pred*c+R;
        K_n=P_pred*c/S_n; IKC=In-K_n*c';
        P_cur=IKC*P_pred*IKC'+K_n*R*K_n';
        x_cur=x_pred+K_n*inn;
        rho=norm(K_n-K_inf)/(norm(K_inf)+eps);
        if rho<gain_tol, below_count=below_count+1; else, below_count=0; end
        if below_count>=gain_patience, use_ss=true; n_ss=n; end
    end

    % Store the step
    x_filt(:,n)=x_cur; P_filt(:,:,n)=P_cur;
    innov(n)=inn; S_out(n)=S_n; K_hist(:,n)=K_n;

end  % filter loop

% =========================================================================
% 8. Filtered discharge, uncertainty and NIS
% =========================================================================
% In Scenario A the feed-through a0*z_hat_d, removed in Section 3, is
% restored here. A warning is issued when the mean NIS departs from its
% self-consistent value 1 by more than 0.3.
% =========================================================================
if strcmp(scenario,'A') && a0~=0
    q_filt = c'*x_filt + a0*z_hat_d(:)';
else
    q_filt = c'*x_filt;
end
q_var=zeros(1,N);
for n=1:N, q_var(n)=c'*P_filt(:,:,n)*c; end
q_sigma=sqrt(max(q_var,0));

NIS=innov.^2./S_out;
if abs(mean(NIS,'omitnan')-1)>0.3
    warning('discharge_kalman: mean NIS=%.3f (expected ~1).', mean(NIS,'omitnan'));
end

% =========================================================================
% 9. Forecast (pure prediction, no update)
% =========================================================================
% Starting from x_{N|N}, P_{N|N}, the state is propagated with the drift
% selected by fore_option; the covariance grows with the last Qd. The band
% q_fore +- 1.96*sigma is returned in q_fore_lo / q_fore_hi.
% =========================================================================
if ~isempty(t_fore)
    Nf=numel(t_fore);
    lam_fore=getopt(opts,'lambda_fore',lambda_used);
    ec_fore =getopt(opts,'eta_c_fore', eta_c_used);
    switch fore_option
        case 0,  U_fore=zeros(ns,Nf);                      % free decay
        case 1                                             % known input
            zhf=getopt(opts,'z_hat_d_fore',zeros(1,Nf));
            if numel(zhf)~=Nf, error('z_hat_d_fore must have length %d',Nf); end
            U_fore=build_U_rect(zhf,k_vec,dt,ns,Nf);
        case 2                                             % constant mean
            if isempty(lam_fore)||isempty(ec_fore)
                error('fore_option 2: set lambda_fore and eta_c_fore');
            end
            U_fore=build_U_rect(lam_fore*ec_fore*dt*ones(1,Nf),k_vec,dt,ns,Nf);
        case 3                                             % seasonal mean
            lam_m=getopt(opts,'lambda_monthly',[]);
            ec_m =getopt(opts,'eta_c_monthly', []);
            mon  =getopt(opts,'months_fore',   []);
            if isempty(lam_m)||isempty(ec_m)||isempty(mon)||numel(mon)~=Nf
                error('fore_option 3: set lambda_monthly, eta_c_monthly, months_fore');
            end
            U_fore=build_U_rect(lam_m(mon(:)').*ec_m(mon(:)')*dt,k_vec,dt,ns,Nf);
        otherwise, error('fore_option must be 0..3');
    end
    xf=x_filt(:,N); Pf=P_filt(:,:,N);
    q_fore=zeros(1,Nf); q_fv=zeros(1,Nf);
    x_fore=zeros(ns,Nf); P_fore=zeros(ns,ns,Nf);
    for f=1:Nf
        xf=Ad*xf+U_fore(:,f); Pf=Ad*Pf*Ad'+Qd;
        x_fore(:,f)=xf; P_fore(:,:,f)=Pf; q_fore(f)=c'*xf; q_fv(f)=c'*Pf*c;
    end
    q_fs=sqrt(max(q_fv,0));
    out.t_fore=t_fore; out.fore_option=fore_option;
    out.q_fore=q_fore;
    out.q_fore_lo=q_fore-1.96*q_fs; out.q_fore_hi=q_fore+1.96*q_fs;
    out.x_fore=x_fore; out.P_fore=P_fore;
else
    out.t_fore=[]; out.fore_option=fore_option;
    out.q_fore=[]; out.q_fore_lo=[]; out.q_fore_hi=[];
    out.x_fore=[]; out.P_fore=[];
end

% =========================================================================
% 10. Filter parameters at the end of the run (for warm starts)
% =========================================================================
% In adaptive mode the final Qd diagonal is converted back to spectral
% densities q_j, and the bootstrap is disabled, since R is already
% adapted.
% =========================================================================
fp.scenario=scenario; fp.gain_tol=gain_tol; fp.gain_patience=gain_patience;
fp.adaptive=do_adaptive;
if do_adaptive
    scale_j  = k_vec/2.*(1-exp(-2*dt./k_vec));
    fp.Q_diag = max(diag(Qd)',1e-20)./scale_j;
    fp.R      = R;
    fp.sage_bootstrap  = 0;
    fp.sage_lambda_min = sh_lam_min;
    fp.sage_window     = sh_S;
    fp.sage_R_min=sh_R_min; fp.sage_R_max=sh_R_max;
    fp.sage_Q_min=sh_Q_min; fp.sage_Q_max=sh_Q_max;
else
    fp.R=R;
    switch scenario
        case 'A', fp.Q_diag=Q_diag_used;
        case 'B', fp.lambda=lambda_used; fp.eta_c=eta_c_used;
    end
end

% =========================================================================
% 11. Output
% =========================================================================
out.t=t; out.y_meas=y_meas; out.scenario=scenario;
out.q_filt=q_filt; out.q_sigma=q_sigma;
out.x_filt=x_filt; out.P_filt=P_filt;
out.innov=innov; out.S=S_out; out.K_hist=K_hist;
out.NIS=NIS; out.n_ss=n_ss; out.K_ss=K_inf;
out.filter_params=fp;
if do_adaptive
    out.R_hist=R_hist; out.Qd_hist=Qd_hist; out.R_boot=R_boot_estimate;
end

end  % discharge_kalman


% =========================================================================
% Local functions
% =========================================================================

% =========================================================================
% Process-noise covariance, Scenario B
% =========================================================================
% Exact discretization of the rank-one continuous covariance G*b*b' over
% one time step, for the PWNE spectral constant G = 2*eta_c^2*lambda:
%       Qd_ij = G*k_i*k_j/(k_i+k_j)*(1 - exp(-(1/k_i+1/k_j)*dt))
% =========================================================================
function Qd = build_Qd_B(G, k_vec, dt, ns)
    Qd=zeros(ns);
    for i=1:ns
        for j=1:ns
        ki=k_vec(i); kj=k_vec(j);
        Qd(i,j)=G*ki*kj/(ki+kj)*(1-exp(-(1/ki+1/kj)*dt));
        end
    end
end

% =========================================================================
% Deterministic input, Scenario A
% =========================================================================
% Exact integration of a rectangular rainfall pulse of rate z(n)/dt over
% one time step:  u_n(j) = z(n)*k_j*(1 - exp(-dt/k_j))/dt.
% Returns U of size [ns x N].
% =========================================================================
function U = build_U_rect(z, k_vec, dt, ns, N)
    scale=(k_vec.*(1-exp(-dt./k_vec))/dt)';
    U=scale*z(:)';
end

% =========================================================================
% Discrete algebraic Riccati equation
% =========================================================================
% Uses dare (Control System Toolbox) when available; otherwise falls back
% to a fixed-point iteration of the Riccati recursion, stopped at a
% relative change below 1e-12 or after 5000 iterations.
% =========================================================================
function P = dare_sym(Ad, c, Qd, R)
    try
        P=dare(Ad,c,Qd,R);
    catch
        P=Qd+eye(size(Qd,1));
        for it_=1:5000
            S_=c'*P*c+R; Pn=Ad*(P-P*c*(S_\(c'*P)))*Ad'+Qd;
            if norm(Pn-P,'fro')/(norm(P,'fro')+eps)<1e-12, break; end
            P=Pn;
        end
        P=Pn; %#ok<UNRCH>
    end
end

% =========================================================================
% Option reader
% =========================================================================
% Returns s.(field) if present and non-empty, the default otherwise.
% =========================================================================
function v = getopt(s, field, default)
    if isfield(s,field)&&~isempty(s.(field)), v=s.(field); else, v=default; end
end
