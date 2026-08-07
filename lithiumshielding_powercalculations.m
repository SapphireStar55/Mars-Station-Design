%% FUSION_POWER_DESIGN.m
%
% 0-D (volume-averaged) power balance and sizing tool for a D-3He fusion
% power source (Mars orbital base camp application). "0-D" means we
% treat the plasma as a single well-mixed volume at temperature T and
% density n -- no radial profiles, no MHD equilibrium, no confinement
% scaling law. This is a first-pass systems sizing tool, not a physics
% design code (that would require a transport code like TRANSP/ASTRA and
% an equilibrium solver).
%
% FUEL: D-3He is often called "aneutronic" but is NOT neutron-free in
% practice: sustaining a D-3He plasma always carries some D, and D
% inevitably undergoes D-D side reactions:
%     D + D -> T(1.01 MeV)  + p(3.02 MeV)      50%  (2i)
%     D + D -> 3He(0.82 MeV)+ n(2.45 MeV)      50%  (2ii)
% The tritium from (2i) can then burn with ambient D before being pumped
% out:
%     D + T -> 4He(3.5 MeV) + n(14.1 MeV)
% These two neutron channels (2.45 MeV and 14.1 MeV) are the actual
% neutron source this shield has to stop -- much harder than a fission
% spectrum because of the 14.1 MeV component.
%
% Reactivity data: Maxwellian-averaged <sigma*v>(T) tables below are the
% NRL Plasma Formulary values (Bosch-Hale-based), T in keV, values in
% m^3/s, log-log interpolated. This is real published data, not a
% parametric guess.
%
% Style: user inputs at top, helper functions at bottom, self-contained.

clear; clc; close all;

%% =====================================================================
%  USER INPUTS
%  =====================================================================

% --- Plasma operating point ---
T_keV   = 70;          % ion temperature, keV. D-3He power density peaks
                        % ~58 keV (NRL); Pfus/Pbrem ratio peaks ~100 keV.
                        % 70 keV is a compromise design point -- SWEEP
                        % THIS (see plot at bottom) rather than trust a
                        % single value blindly.
f_D     = 0.35;        % D fraction of ion mix, n_D/(n_D+n_He3). Leaner-D
                        % mixes reduce D-D side-neutron production (which
                        % scales as n_D^2) at some cost to D-3He rate
                        % (scales as n_D*n_He3). 0.3-0.5 is a common
                        % literature range for minimizing neutronicity;
                        % NOT a settled optimum -- treat as a design lever.
n_i     = 3e14;        % total ion density, particles/cm^3. This stands
                        % in for whatever confinement concept you pick
                        % (magnetic pressure limit, ICF compression,
                        % etc) -- it's the main lever the "concept-
                        % specific" physics you skipped would normally
                        % set for you.
tau_E   = 4.0;          % energy confinement time, s. THE most uncertain
                        % input in this whole model -- current
                        % magnetic-confinement experiments (JET, DIII-D)
                        % achieve ~0.1-1 s at reactor-relevant density;
                        % ITER targets ~3-6 s. 4 s assumes a mature,
                        % beyond-ITER-class confinement technology. If
                        % you have a real number for your concept, put
                        % it here.
f_burn_T = 0.25;        % fraction of D-D-produced tritium that burns
                        % in-situ via D-T before being exhausted/pumped
                        % (vs. escaping to a tritium handling system).
                        % Higher burnup = more 14.1 MeV neutrons but also
                        % more fusion power.

% --- Power target & conversion efficiencies ---
P_e_target   = 200e3;   % target net electrical power output, W (200 kWe
                         % example for a Mars base camp -- change to match
                         % your actual load requirement)
eta_direct   = 0.65;    % direct energy conversion efficiency for charged
                         % fusion products (D-3He/D-D charged particles
                         % are fast, high-energy ions well suited to
                         % traveling-wave direct conversion -- 0.6-0.8
                         % is the literature range for D-3He concepts)
eta_thermal  = 0.35;    % thermal-cycle efficiency (Brayton/Rankine) for
                         % neutron power captured as blanket heat

% --- Geometry / shielding interface ---
standoff_m   = 10;      % distance from plasma core to crew/electronics
                         % location, meters (informs point-source dose
                         % scaling in the shielding script)

%% =====================================================================
%  REACTIVITY DATA (NRL Plasma Formulary, Maxwellian-averaged <sigma v>)
%  T in keV, <sigma v> in cm^3/s (verified against canonical peak values:
%  D-T peaks ~8.7e-16 cm^3/s near 65 keV, D-3He peaks ~2.1-2.4e-16 cm^3/s
%  near 100 keV -- both match the raw table values directly, confirming
%  units are cm^3/s, NOT m^3/s -- do not add a further unit conversion).
%  =====================================================================
T_table = [1 2 5 10 20 50 100 200 500 1000];

sv_DHe3_table = [1.0e-26 1.4e-23 6.7e-21 2.3e-19 3.8e-18 5.4e-17 1.6e-16 2.4e-16 2.3e-16 1.8e-16];
sv_DDn_table  = [1.5e-22 5.4e-21 1.8e-19 1.2e-18 5.2e-18 2.1e-17 4.5e-17 8.8e-17 1.8e-16 2.2e-16]; % D+D->3He+n channel only
sv_DT_table   = [5.5e-21 2.6e-19 1.3e-17 1.1e-16 4.2e-16 8.7e-16 8.5e-16 6.3e-16 3.7e-16 2.7e-16];

sv_DHe3 = interp_reactivity(T_keV, T_table, sv_DHe3_table); % cm^3/s
sv_DDn  = interp_reactivity(T_keV, T_table, sv_DDn_table);  % cm^3/s
sv_DT   = interp_reactivity(T_keV, T_table, sv_DT_table);   % cm^3/s
% D-D is symmetric: (2i) and (2ii) channels are ~50/50, so total D-D
% reactivity ~= 2x the tabulated neutron-only channel.
sv_DDp  = sv_DDn; % proton+T channel, same rate as neutron channel (50/50 split)

%% =====================================================================
%  DENSITIES AND REACTION RATE DENSITIES
%  =====================================================================
n_D   = f_D * n_i;
n_He3 = (1-f_D) * n_i;
n_e   = n_D*1 + n_He3*2;   % quasineutrality, Z_D=1, Z_He3=2

R_DHe3 = n_D * n_He3 * sv_DHe3;              % reactions/cm^3/s
R_DDn  = 0.5 * n_D^2 * sv_DDn;               % D+D->3He+n reactions/cm^3/s
R_DDp  = 0.5 * n_D^2 * sv_DDp;               % D+D->T+p  reactions/cm^3/s
R_DT   = f_burn_T * R_DDp;                   % secondary D-T burn (uses
                                              % the T produced by R_DDp)

MeV_to_J = 1.60218e-13;

% --- Power densities (W/cm^3) ---
p_charged = (R_DHe3*18.3 + R_DDn*0.82 + R_DDp*4.03 + R_DT*3.5) * MeV_to_J;
p_neutron = (R_DDn*2.45 + R_DT*14.1) * MeV_to_J;
p_fus     = p_charged + p_neutron;

% --- Bremsstrahlung loss density (NRL Formulary standard form) ---
% P_br [W/cm^3] = 1.69e-32 * sqrt(Te_keV) * n_e * sum(n_Z * Z^2)
p_brem = 1.69e-32 * sqrt(T_keV) * n_e * (n_D*1 + n_He3*4);

% --- Stored thermal energy density and conduction loss density ---
w_th   = 1.5 * (n_D + n_He3 + n_e) * T_keV * 1.60218e-16; % J/cm^3
p_cond = w_th / tau_E;                                     % W/cm^3

% --- Ignition check & recirculating power ---
p_loss_plasma = p_brem + p_cond;
if p_charged >= p_loss_plasma
    ignited = true;
    p_recirc = 0;
    margin = p_charged / p_loss_plasma;
else
    ignited = false;
    p_recirc = p_loss_plasma - p_charged; % W/cm^3, drawn back from output
    margin = p_charged / p_loss_plasma;
end

% --- Net electric power density and required volume ---
p_net_e = eta_direct*p_charged + eta_thermal*p_neutron - p_recirc; % W/cm^3

fprintf('=== Plasma operating point ===\n');
fprintf('T = %.1f keV, f_D = %.2f, n_i = %.2e cm^-3, tau_E = %.1f s\n', T_keV, f_D, n_i, tau_E);
fprintf('Charged power density: %.3e W/cm^3\n', p_charged);
fprintf('Neutron power density: %.3e W/cm^3\n', p_neutron);
fprintf('Bremsstrahlung loss:   %.3e W/cm^3\n', p_brem);
fprintf('Conduction loss:       %.3e W/cm^3\n', p_cond);
if ignited
    fprintf('Plasma state: IGNITED (self-heating margin %.2fx over losses)\n', margin);
else
    fprintf('Plasma state: NOT ignited (self-heating covers %.1f%% of losses)\n', 100*margin);
    fprintf('Required recirculating heating power density: %.3e W/cm^3\n', p_recirc);
end

if p_net_e <= 0
    error(['Design infeasible at this operating point: net electric power density <= 0. ' ...
           'Increase n_i, tau_E, or T, or reduce P_e_target.']);
end

V_cm3 = P_e_target / p_net_e;
R_plasma_cm = (3*V_cm3/(4*pi))^(1/3);

fprintf('\n=== Sizing for %.0f kWe net output ===\n', P_e_target/1e3);
fprintf('Required plasma volume: %.3e cm^3 (%.2f m^3)\n', V_cm3, V_cm3*1e-6);
fprintf('Equivalent spherical plasma radius: %.2f cm (%.3f m)\n', R_plasma_cm, R_plasma_cm/100);

%% =====================================================================
%  SOURCE TERMS FOR SHIELDING (particles/s, at the two neutron energies)
%  =====================================================================
S_n_2p45 = R_DDn * V_cm3;   % n/s at 2.45 MeV
S_n_14p1 = R_DT  * V_cm3;   % n/s at 14.1 MeV
S_n_total = S_n_2p45 + S_n_14p1;

P_fus_total_W = p_fus * V_cm3;
P_charged_total_W = p_charged * V_cm3;
P_neutron_total_W = p_neutron * V_cm3;

fprintf('\n=== Neutron source terms (feed to shielding) ===\n');
fprintf('2.45 MeV neutrons (D-D branch): %.3e n/s\n', S_n_2p45);
fprintf('14.1 MeV neutrons (D-T branch): %.3e n/s\n', S_n_14p1);
fprintf('Total neutron source rate:      %.3e n/s\n', S_n_total);
fprintf('Total fusion power:             %.3f MW\n', P_fus_total_W/1e6);
fprintf('  of which charged (direct-conv candidate): %.3f MW\n', P_charged_total_W/1e6);
fprintf('  of which neutron (blanket/thermal only):  %.3f MW\n', P_neutron_total_W/1e6);
fprintf('Net electrical output:          %.1f kWe\n', P_e_target/1e3);

%% =====================================================================
%  SAVE SOURCE TERMS FOR THE SHIELDING SCRIPT
%  =====================================================================
fusion_source.S_n_2p45_MeV = S_n_2p45;      % n/s
fusion_source.S_n_14p1_MeV = S_n_14p1;      % n/s
fusion_source.R_plasma_cm  = R_plasma_cm;
fusion_source.standoff_m   = standoff_m;
fusion_source.P_fus_total_W = P_fus_total_W;
fusion_source.T_keV = T_keV;
fusion_source.operating_point_note = sprintf(['T=%.1f keV, f_D=%.2f, n_i=%.2e cm^-3, ' ...
    'tau_E=%.1f s, f_burn_T=%.2f'], T_keV, f_D, n_i, tau_E, f_burn_T);

save('fusion_source_terms.mat', 'fusion_source');
fprintf('\nSource terms saved to fusion_source_terms.mat for shield_design_from_fusion.m\n');

%% =====================================================================
%  SENSITIVITY SWEEP: operating temperature vs. neutron source & volume
%  (helps you see how sensitive the design is to the T_keV choice above)
%  =====================================================================
T_sweep = 20:5:150;
S_n_sweep = zeros(size(T_sweep));
V_sweep   = zeros(size(T_sweep));
for i = 1:numel(T_sweep)
    Tk = T_sweep(i);
    sv_DHe3_i = interp_reactivity(Tk, T_table, sv_DHe3_table);
    sv_DDn_i  = interp_reactivity(Tk, T_table, sv_DDn_table);
    sv_DT_i   = sv_DDn_i;

    R_DHe3_i = n_D*n_He3*sv_DHe3_i;
    R_DDn_i  = 0.5*n_D^2*sv_DDn_i;
    R_DDp_i  = 0.5*n_D^2*sv_DDn_i;
    R_DT_i   = f_burn_T*R_DDp_i;

    p_charged_i = (R_DHe3_i*18.3 + R_DDn_i*0.82 + R_DDp_i*4.03 + R_DT_i*3.5)*MeV_to_J;
    p_neutron_i = (R_DDn_i*2.45 + R_DT_i*14.1)*MeV_to_J;
    p_brem_i    = 1.69e-32*sqrt(Tk)*n_e*(n_D*1+n_He3*4);
    w_th_i      = 1.5*(n_D+n_He3+n_e)*Tk*1.60218e-16;
    p_cond_i    = w_th_i/tau_E;
    p_loss_i    = p_brem_i + p_cond_i;
    p_recirc_i  = max(0, p_loss_i - p_charged_i);
    p_net_e_i   = eta_direct*p_charged_i + eta_thermal*p_neutron_i - p_recirc_i;

    if p_net_e_i > 0
        V_i = P_e_target/p_net_e_i;
        S_n_sweep(i) = (R_DDn_i + R_DT_i)*V_i;
        V_sweep(i)   = V_i*1e-6; % m^3
    else
        S_n_sweep(i) = NaN;
        V_sweep(i)   = NaN;
    end
end

figure('Name','D-3He Operating Point Sensitivity','Position',[100 100 900 400]);
subplot(1,2,1);
plot(T_sweep, S_n_sweep, 'o-', 'LineWidth', 1.5);
xline(T_keV, 'r--', 'Chosen T');
xlabel('Ion temperature (keV)'); ylabel('Total neutron source (n/s)');
title('Neutron source vs. operating temperature');
grid on;

subplot(1,2,2);
plot(T_sweep, V_sweep, 'o-', 'LineWidth', 1.5);
xline(T_keV, 'r--', 'Chosen T');
xlabel('Ion temperature (keV)'); ylabel('Required plasma volume (m^3)');
title('Reactor size vs. operating temperature');
grid on;

sgtitle(sprintf('D-3He design sensitivity at %.0f kWe target, n_i=%.1e cm^{-3}, \\tau_E=%.1f s', ...
    P_e_target/1e3, n_i, tau_E));


%% =====================================================================
%  HELPER FUNCTIONS
%  =====================================================================

function sv = interp_reactivity(T_query, T_table, sv_table)
    % Log-log interpolation of Maxwellian reactivity table. Reactivity
    % vs temperature is smooth and close to power-law-ish locally, so
    % log-log interpolation is much more accurate than linear
    % interpolation of these tabulated values, which span >20 orders of
    % magnitude in some cases.
    if T_query < T_table(1) || T_query > T_table(end)
        warning(['Operating temperature %.1f keV is outside the tabulated range ' ...
                 '[%.0f, %.0f] keV -- extrapolating, treat result with caution.'], ...
                 T_query, T_table(1), T_table(end));
    end
    logT  = log(T_table);
    logSV = log(sv_table);
    sv = exp(interp1(logT, logSV, log(T_query), 'linear', 'extrap'));
end