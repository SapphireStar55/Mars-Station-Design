%% SHIELD_DESIGN_FROM_FUSION.m
%
% Sizes a LiH + W shield around the D-3He fusion source defined by
% fusion_power_design.m, using the same 1-D slab Monte Carlo transport
% approach as the original, but now driven by actual neutron
% source rates (n/s) at the two physical D-3He neutron energies (2.45
% MeV from D-D, 14.1 MeV from secondary D-T burn) instead of a generic
% fission spectrum, and iterating shield thickness to meet a dose target
% at a specified crew standoff distance.
%
% GEOMETRY MODEL: thin-shell approximation. The 1-D slab MC gives the
% shield's neutron/gamma TRANSMISSION FRACTION (fraction of source
% particles that punch through the shield). We then apply standard
% inverse-square spreading from the reactor core to the crew standoff
% distance, i.e.:
%     flux_at_crew = S_total * T_shield / (4*pi*standoff_cm^2)
% This is accurate when shield thickness << standoff distance (true
% here: shield ~tens of cm, standoff ~10 m) -- it is NOT a full 3-D
% spherical-shell transport solution.
% DOSE TARGET: NASA's current astronaut CAREER limit (NASA-STD-3001) is
% 600 mSv total effective dose, with a separate nuclear-technology-
% specific sub-limit (requirement V1 4032) whose exact numeric value I
% do not have verified -- you should look that up directly before
% finalizing a design. In the absence of that number, this script uses
% a placeholder ANNUAL reactor dose budget (user input below) that is a
% small fraction of the career limit, consistent with the standard
% ALARA (As Low As Reasonably Achievable) approach of budgeting only a
% portion of total allowable dose to any one source (GCR, SPE, and
% reactor all draw from the same career budget).
%
% Style: user inputs at top, helper functions at bottom, self-contained.

clear; clc; close all;

%% =====================================================================
%  USER INPUTS
%  =====================================================================
load('fusion_source_terms.mat', 'fusion_source');

S_n_2p45 = fusion_source.S_n_2p45_MeV;   % n/s at 2.45 MeV
S_n_14p1 = fusion_source.S_n_14p1_MeV;   % n/s at 14.1 MeV
S_n_total = S_n_2p45 + S_n_14p1;
standoff_cm = fusion_source.standoff_m * 100;

fprintf('Loaded fusion source: %.3e n/s @ 2.45 MeV, %.3e n/s @ 14.1 MeV\n', S_n_2p45, S_n_14p1);
fprintf('Crew standoff distance: %.1f m\n', fusion_source.standoff_m);
fprintf('Operating point: %s\n\n', fusion_source.operating_point_note);

% --- Dose target (PLACEHOLDER -- verify against NASA-STD-3001 V1 4032) ---
annual_dose_budget_mSv = 5.0;      % mSv/year allocated to this reactor
                                     % (ALARA-style budget fraction of the
                                     % 600 mSv career limit; NOT a
                                     % verified regulatory number for
                                     % nuclear-source-specific limits --
                                     % update once you have the real V1
                                     % 4032 figure)
hours_per_year_exposed = 8760;      % assume continuous exposure (worst case)
target_dose_rate_mrem_hr = annual_dose_budget_mSv * 100 / hours_per_year_exposed; % mSv->mrem: x100

fprintf('Target dose rate at standoff: %.4f mrem/hr (%.2f mSv/yr budget, continuous exposure)\n\n', ...
    target_dose_rate_mrem_hr, annual_dose_budget_mSv);

% --- Material properties (same as mc_shield_transport.m) ---
rho_LiH = 0.78; M_LiH = 7.95; N_LiH = rho_LiH*6.02214e23/M_LiH;
mat_data.LiH = struct('N', N_LiH, 'frac_Li6', 0.075, 'frac_Li7', 0.925, 'rho', rho_LiH);
rho_W = 19.3; M_W = 183.84; N_W = rho_W*6.02214e23/M_W;
mat_data.W = struct('N', N_W, 'rho', rho_W, 'Z', 74);

% --- Search grid for shield thickness optimization ---
% Mass-minimization proxy: search combos of LiH/W thickness, keep the
% lightest (by areal density) combo that meets the dose target.
LiH_search_cm = 10:5:60;
W_search_cm   = 2:2:30;
N_search      = 15000;   % reduced particle count for the search grid
N_final       = 300000;  % high-fidelity run at the selected design point

% --- Dose conversion factors (ANSI/ANS-6.1.1-style, approximate) ---
neutron_dose_E   = [0.01 0.1 0.5 1 2 5 10 14 20];
neutron_dose_fac = [3.5e-7 3.2e-7 1.1e-6 4.5e-6 8.0e-6 1.0e-5 1.5e-5 2.0e-5 2.4e-5]; % rem/hr per n/cm2-s
gamma_dose_E     = [0.01 0.1 0.5 1 2 5 10];
gamma_dose_fac   = [1.0e-7 3.0e-7 6.0e-7 1.0e-6 1.7e-6 3.0e-6 4.5e-6]; % rem/hr per gamma/cm2-s

%% =====================================================================
%  SHIELD THICKNESS SEARCH -- ANALYTICAL (point-kernel), NOT MONTE CARLO
%  =====================================================================
% IMPORTANT LESSON FROM DEVELOPING THIS SCRIPT: an earlier version of
% this search ran the analog MC transport at reduced particle count
% (N~15000) for every grid point. That is WRONG and gives silently bad
% answers: at the attenuation depths this source requires (source is
% ~1e16 n/s; target dose needs transmission suppressed by roughly
% 8-9 orders of magnitude), the true transmission probability per
% history is far below 1/N, so most search-grid MC runs return exactly
% ZERO transmitted particles not because the shield is adequate but
% because the rare event was never sampled. That false zero looks like
% "dose met" and the search will happily pick a shield that is actually
% thousands of times over the dose target once you check it with a
% larger N (this is exactly what happened in testing: a 20 cm LiH + 2 cm
% W combination "passed" the N=3000 search, then measured ~35,000x over
% target at N=20000). Analog MC cannot resolve attenuation this deep
% without variance reduction (importance sampling / weight windows),
% which is beyond a hand-rolled script -- that's what real shielding
% codes (MCNP, Serpent, OpenMC with weight-window generators) are for.
%
% So: the SEARCH below uses a fast analytical point-kernel exponential-
% attenuation model with a flat buildup factor (standard first-pass
% shielding-engineering practice), which has no such statistical floor.
% The Monte Carlo is then run ONCE at the end at the selected design
% point purely as a diagnostic (spectrum shape, fate breakdown) -- NOT
% to re-derive the absolute dose number, which is outside what N_final
% histories can statistically resolve at this attenuation depth.

buildup_factor = 5; % flat conservative buildup factor (multiplies the
                     % pure-exponential dose to approximate scattered/
                     % secondary dose not captured by uncollided-flux
                     % attenuation alone). Real buildup factors are
                     % energy- and thickness-dependent (see ANSI/ANS-
                     % 6.4.3 tables) -- 5x is a commonly-used
                     % conservative placeholder for few-mean-free-path
                     % shields, not a substitute for the real tables.

fprintf('Searching shield thickness grid (%d x %d combinations, analytical attenuation)...\n', ...
    numel(LiH_search_cm), numel(W_search_cm));

best_areal_density = Inf;
best_LiH = NaN; best_W = NaN;
results_table = [];

for iL = 1:numel(LiH_search_cm)
    for iW = 1:numel(W_search_cm)
        tL = LiH_search_cm(iL);
        tW = W_search_cm(iW);
        areal_density = tL*rho_LiH + tW*rho_W; % g/cm^2, mass proxy

        if areal_density >= best_areal_density
            continue; % can't possibly beat current best, skip
        end

        dose_mrem_hr = analytical_dose_estimate(tL, tW, S_n_2p45, S_n_14p1, ...
            standoff_cm, mat_data, buildup_factor, neutron_dose_E, neutron_dose_fac);

        results_table = [results_table; tL, tW, areal_density, dose_mrem_hr]; %#ok<AGROW>

        if dose_mrem_hr <= target_dose_rate_mrem_hr && areal_density < best_areal_density
            best_areal_density = areal_density;
            best_LiH = tL; best_W = tW;
        end
    end
end

if isnan(best_LiH)
    warning(['No combination in the search grid met the dose target. ' ...
             'Widen LiH_search_cm / W_search_cm ranges and re-run -- or ' ...
             'increase standoff distance, which is usually cheaper than ' ...
             'more shield mass for a source this intense.']);
    [~, idx] = min(results_table(:,4));
    best_LiH = results_table(idx,1);
    best_W   = results_table(idx,2);
    fprintf('Falling back to lowest-dose combination found: LiH=%.0f cm, W=%.0f cm\n', best_LiH, best_W);
else
    fprintf('\nLightest shield meeting dose target (analytical estimate): LiH = %.0f cm, W = %.0f cm (areal density %.1f g/cm^2)\n', ...
        best_LiH, best_W, best_areal_density);
end

%% =====================================================================
%  FINAL MONTE CARLO RUN -- DIAGNOSTIC ONLY, NOT A DOSE VERIFICATION
%  =====================================================================
% This run reports the transmitted-neutron energy spectrum and the
% neutron fate breakdown (useful for sanity-checking the physics and for
% seeing which layer is doing the work). It does NOT re-verify the
% analytical dose number above -- see the note in the search section:
% at this attenuation depth, N_final histories will very likely show
% ZERO transmitted neutrons, which is the CORRECT and EXPECTED MC
% outcome (it means true transmission is below the ~1/N statistical
% floor), not evidence the shield is either adequate or inadequate.
fprintf('\nRunning diagnostic Monte Carlo (N=%d) at the selected design point...\n', N_final);
[~, diag, gdiag] = run_shield_transport(N_final, best_LiH, best_W, ...
    S_n_2p45, S_n_14p1, S_n_total, standoff_cm, mat_data, ...
    neutron_dose_E, neutron_dose_fac, gamma_dose_E, gamma_dose_fac);

analytical_dose_final = analytical_dose_estimate(best_LiH, best_W, S_n_2p45, S_n_14p1, ...
    standoff_cm, mat_data, buildup_factor, neutron_dose_E, neutron_dose_fac);

fprintf('\n=== FINAL SHIELD DESIGN ===\n');
fprintf('LiH thickness: %.0f cm\n', best_LiH);
fprintf('W thickness:   %.0f cm\n', best_W);
fprintf('Total shield mass (per unit area): %.1f g/cm^2\n', best_LiH*rho_LiH + best_W*rho_W);
fprintf('MC neutron transmission in N=%d histories: %d transmitted (%.4f%%)\n', ...
    N_final, numel(diag.transmitted_n_E), 100*diag.transmitted_frac);
fprintf('  -> statistical floor at this N is ~%.1e; a zero count here is\n', 1/N_final);
fprintf('     consistent with (not proof of) the analytical estimate below.\n');
fprintf('Analytical (point-kernel + buildup-factor-%.0fx) dose estimate at %.1f m: %.3e mrem/hr (target: %.3e mrem/hr)\n', ...
    buildup_factor, fusion_source.standoff_m, analytical_dose_final, target_dose_rate_mrem_hr);
fprintf('Implied annual dose at continuous exposure: %.3e mSv/yr\n', ...
    analytical_dose_final/100*hours_per_year_exposed);
fprintf(['\nCAVEAT: the analytical point-kernel model uses the SAME simplified, non-ENDF\n' ...
         'parametric cross sections as the MC (see get_neutron_xs/get_gamma_xs), plus a flat\n' ...
         'buildup factor. This is a first-pass sizing estimate, not a certifiable dose number.\n' ...
         'A real design needs ENDF/B cross sections and a production transport code (MCNP,\n' ...
         'Serpent, OpenMC) run with proper variance reduction at this attenuation depth.\n']);

%% =====================================================================
%  PLOTS
%  =====================================================================
figure('Name','Shield Design Results','Position',[100 100 1100 450]);

subplot(1,3,1);
if ~isempty(results_table)
    scatter(results_table(:,1), results_table(:,2), 60, log10(results_table(:,4)+1e-6), 'filled');
    colorbar; hold on;
    plot(best_LiH, best_W, 'rp', 'MarkerSize', 18, 'MarkerFaceColor','y');
    xlabel('LiH thickness (cm)'); ylabel('W thickness (cm)');
    title('log_{10}(dose rate) search grid');
end

subplot(1,3,2);
histogram(diag.transmitted_n_E, 30);
xlabel('Exit neutron energy (MeV)'); ylabel('Count');
title(sprintf('Transmitted neutrons (N=%d)', N_final));

subplot(1,3,3);
categories = categorical({'Transmitted','Absorbed LiH','Absorbed W'});
bar(categories, [numel(diag.transmitted_n_E), diag.absorbed_LiH, diag.absorbed_W]);
ylabel('Count'); title('Neutron fate (final design)');

sgtitle(sprintf('D-3He Fusion Shield: %.0f cm LiH + %.0f cm W, dose %.3f mrem/hr @ %.0f m', ...
    best_LiH, best_W, dose_mrem_hr_final, fusion_source.standoff_m));


%% =====================================================================
%  HELPER FUNCTIONS
%  =====================================================================

function [dose_mrem_hr, diag, gdiag] = run_shield_transport(N, t_LiH, t_W, ...
    S_n_2p45, S_n_14p1, S_n_total, standoff_cm, mat_data, ...
    neutron_dose_E, neutron_dose_fac, gamma_dose_E, gamma_dose_fac)
    % Runs the neutron + gamma MC transport through a 2-layer LiH/W slab
    % for N source histories, then scales the resulting transmission
    % fraction by the actual source rate and 1/(4*pi*r^2) spreading to
    % get a real dose rate at the crew standoff distance.

    layers = struct('name', {'LiH','W'}, 'thickness_cm', {t_LiH, t_W});
    n_layers = 2;
    bounds = [0, t_LiH, t_LiH+t_W];

    p_14p1 = S_n_14p1 / S_n_total; % probability a source neutron is the 14.1 MeV line

    transmitted_n_E = [];
    absorbed_LiH = 0; absorbed_W = 0;
    absorbed_positions = [];

    for ip = 1:N
        if rand() < p_14p1
            E = 14.1;
        else
            E = 2.45;
        end
        x = 0; mu = 1.0; alive = true;

        while alive
            layer_idx = get_layer_index(x, bounds);
            if layer_idx == 0
                alive = false; break; % reflected
            elseif layer_idx > n_layers
                transmitted_n_E(end+1) = E; %#ok<AGROW>
                alive = false; break;
            end
            mat = layers(layer_idx).name;
            [SigT, ~, SigA] = get_neutron_xs(mat, E, mat_data);
            s = -log(rand())/SigT;
            x_new = x + s*mu;
            [x_new, mu, layer_idx, crossed] = handle_boundary(x, x_new, mu, bounds, layer_idx);
            x = x_new;
            if crossed, continue; end

            if rand() < SigA/SigT
                absorbed_positions(end+1) = x; %#ok<AGROW>
                if layer_idx == 1, absorbed_LiH = absorbed_LiH+1; else, absorbed_W = absorbed_W+1; end
                alive = false;
            else
                A = get_effective_mass_number(mat);
                alpha = ((A-1)/(A+1))^2;
                E = E*(alpha + (1-alpha)*rand());
                mu = 2*rand()-1;
                if E < 1e-11
                    absorbed_positions(end+1) = x; %#ok<AGROW>
                    if layer_idx==1, absorbed_LiH=absorbed_LiH+1; else, absorbed_W=absorbed_W+1; end
                    alive = false;
                end
            end
        end
    end

    % --- capture gammas from W captures (2/capture, 0.5-3 MeV, as before) ---
    capture_gamma_E = []; capture_gamma_x = [];
    for k = 1:numel(absorbed_positions)
        li = get_layer_index(absorbed_positions(k), bounds);
        if li == 2
            for g = 1:2
                capture_gamma_E(end+1) = 0.5+2.5*rand(); %#ok<AGROW>
                capture_gamma_x(end+1) = absorbed_positions(k); %#ok<AGROW>
            end
        end
    end
    N_gamma_total = numel(capture_gamma_E);
    transmitted_g_E = [];
    for ig = 1:N_gamma_total
        E = capture_gamma_E(ig); x = capture_gamma_x(ig); mu = 1.0; alive = true;
        while alive
            layer_idx = get_layer_index(x, bounds);
            if layer_idx==0, alive=false; break;
            elseif layer_idx>n_layers, transmitted_g_E(end+1)=E; alive=false; break; %#ok<AGROW>
            end
            mat = layers(layer_idx).name;
            mu_att = get_gamma_xs(mat, E, mat_data);
            s = -log(rand())/mu_att;
            x_new = x + s*mu;
            [x_new, mu, layer_idx, crossed] = handle_boundary(x, x_new, mu, bounds, layer_idx);
            x = x_new;
            if crossed, continue; end
            r = rand();
            if E < 0.1, p_pe=0.7; p_cs=0.3; elseif E<1.022, p_pe=0.05; p_cs=0.95; else, p_pe=0.02; p_cs=0.68; end
            if r < p_pe
                alive = false;
            elseif r < p_pe+p_cs
                costheta = 1-2*rand();
                E = E/(1+(E/0.511)*(1-costheta));
                mu = costheta;
                if E < 0.01, alive = false; end
            else
                alive = false;
            end
        end
    end

    transmitted_frac = numel(transmitted_n_E)/N;

    % --- dose scaling: real source rate + inverse-square to standoff ---
    n_dose_factor = interp1(neutron_dose_E, neutron_dose_fac, transmitted_n_E, 'linear', 'extrap');
    n_dose_factor(transmitted_n_E > neutron_dose_E(end)) = neutron_dose_fac(end); % cap extrapolation
    n_flux_at_standoff_per_source = sum(n_dose_factor)/N * S_n_total / (4*pi*standoff_cm^2);

    if N_gamma_total > 0
        g_dose_factor = interp1(gamma_dose_E, gamma_dose_fac, transmitted_g_E, 'linear', 'extrap');
        % gamma "source rate" here is per-neutron-source-history capture-gamma yield;
        % scale by same S_n_total since capture gammas are generated per source neutron
        g_flux_at_standoff = sum(g_dose_factor)/N * S_n_total / (4*pi*standoff_cm^2);
    else
        g_flux_at_standoff = 0;
    end

    dose_mrem_hr = (n_flux_at_standoff_per_source + g_flux_at_standoff) * 1000; % rem/hr -> mrem/hr

    diag.transmitted_n_E = transmitted_n_E;
    diag.transmitted_frac = transmitted_frac;
    diag.absorbed_LiH = absorbed_LiH;
    diag.absorbed_W = absorbed_W;
    gdiag.transmitted_g_E = transmitted_g_E;
    gdiag.N_gamma_total = N_gamma_total;
end

function dose_mrem_hr = analytical_dose_estimate(t_LiH, t_W, S_n_2p45, S_n_14p1, ...
    standoff_cm, mat_data, buildup_factor, neutron_dose_E, neutron_dose_fac)
    % Point-kernel exponential attenuation estimate:
    %   flux_at_standoff(E) = S(E) * exp(-SigT_LiH(E)*t_LiH - SigT_W(E)*t_W)
    %                          * buildup_factor / (4*pi*standoff_cm^2)
    % applied separately to each source line (2.45 MeV, 14.1 MeV), then
    % converted to dose via the same flux-to-dose factors used elsewhere
    % and summed. This has no statistical noise floor (unlike MC), which
    % is exactly why it's used for the search -- see the note above.
    E_lines = [2.45, 14.1];
    S_lines = [S_n_2p45, S_n_14p1];
    dose_rem_hr = 0;
    for i = 1:2
        E = E_lines(i);
        [SigT_LiH, ~, ~] = get_neutron_xs('LiH', E, mat_data);
        [SigT_W, ~, ~]   = get_neutron_xs('W', E, mat_data);
        transmission = exp(-SigT_LiH*t_LiH - SigT_W*t_W);
        flux = S_lines(i) * transmission * buildup_factor / (4*pi*standoff_cm^2);
        dose_factor = interp1(neutron_dose_E, neutron_dose_fac, E, 'linear', 'extrap');
        dose_rem_hr = dose_rem_hr + flux*dose_factor;
    end
    dose_mrem_hr = dose_rem_hr * 1000;
end

function idx = get_layer_index(x, bounds)
    n_layers = numel(bounds)-1;
    if x < 0, idx = 0;
    elseif x >= bounds(end), idx = n_layers+1;
    else, idx = find(x >= bounds(1:end-1) & x < bounds(2:end), 1, 'first');
    end
end

function [x_out, mu_out, layer_idx_out, crossed] = handle_boundary(x_start, x_end, mu, bounds, layer_idx) %#ok<INUSL>
    crossed = false;
    n_layers = numel(bounds)-1;
    if x_end > bounds(min(layer_idx+1, n_layers+1)) && mu > 0 && layer_idx <= n_layers
        x_out = bounds(layer_idx+1) + 1e-9; mu_out = mu; layer_idx_out = layer_idx+1; crossed = true;
    elseif x_end < bounds(layer_idx) && mu < 0 && layer_idx >= 1
        x_out = bounds(layer_idx) - 1e-9; mu_out = mu; layer_idx_out = layer_idx-1; crossed = true;
    else
        x_out = x_end; mu_out = mu; layer_idx_out = layer_idx; crossed = false;
    end
end

function [SigT, SigS, SigA] = get_neutron_xs(mat, E, mat_data)
    % Same parametric cross sections as mc_shield_transport.m, extended
    % to the higher (14.1 MeV) energy range typical of D-T neutrons.
    % NOTE: at 14.1 MeV, (n,2n) and inelastic scattering in W become
    % significant in reality and are NOT modeled here (this model only
    % has elastic scatter + capture) -- real design needs ENDF/B data
    % and a code that handles (n,2n) multiplication.
    barn = 1e-24;
    if strcmp(mat, 'LiH')
        d = mat_data.LiH;
        sigma_s_H  = 20/sqrt(max(E,0.01))*barn;
        sigma_s_Li = 1.4*barn;
        sigma_a_Li6 = 940*sqrt(0.025e-6/max(E,1e-9))*barn;
        SigS = d.N*sigma_s_H + d.N*sigma_s_Li;
        SigA = d.N*d.frac_Li6*sigma_a_Li6;
        SigT = SigS+SigA;
    elseif strcmp(mat, 'W')
        d = mat_data.W;
        sigma_s = 5.0*barn;
        sigma_a = 2.0*sqrt(0.025e-6/max(E,1e-9))*barn + 0.5*barn;
        SigS = d.N*sigma_s; SigA = d.N*sigma_a; SigT = SigS+SigA;
    else
        error('Unknown material: %s', mat);
    end
    SigT = max(SigT, 1e-6);
end

function A = get_effective_mass_number(mat)
    if strcmp(mat,'LiH'), A = 1.0; else, A = 184.0; end
end

function mu_att = get_gamma_xs(mat, E, mat_data)
    if strcmp(mat,'LiH'), Zeff = 1.9; rho = mat_data.LiH.rho;
    elseif strcmp(mat,'W'), Zeff = 74; rho = mat_data.W.rho;
    else, error('Unknown material: %s', mat);
    end
    E = max(E, 0.01);
    mu_pe = 1e-3*(Zeff/20)^4.5/E^3;
    mu_cs = 0.15*(Zeff/20)/(1+E);
    if E > 1.022, mu_pp = 0.02*(Zeff/74)*(E-1.022); else, mu_pp = 0; end
    mu_att = rho*(mu_pe+mu_cs+mu_pp);
    mu_att = max(mu_att, 1e-5);
end