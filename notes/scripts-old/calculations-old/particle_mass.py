import math

# Updated PDG 2025 Masses (MeV; from pdg.lbl.gov/2025 sum-quarks)
PDG_MASSES = {
    'electron': 0.510998950,
    'muon': 105.6583755,
    'tau': 1776.86,
    'u': 2.16,
    'd': 4.67,
    'c': 1270,
    's': 93,
    't': 172690,
    'b': 4180,
    'proton': 938.272,
    'lambda': 1115.683,
    'sigma': 1189.37,
    'xi': 1314.86,
    'omega': 1672.45
}

# Constants
PHI = (1 + math.sqrt(5)) / 2
BETA = 1 / (2 * math.pi)

# Anchors
EPSILON_LEPTON = 0.0603
P_AVG_QUARK = 1.43
EPSILON_QUARK = 0.55
A_L = 2.734
KAPPA = 15.299
ZETA_FACTOR = 19.6
ZETA = KAPPA / (PHI**2 * ZETA_FACTOR)

LAMBDA_QCD = 250
H_BAR_MEV_S = 6.582e-22

# Instability etas (derived from Gamma/m or Lambda/m; higher for unstables)
ETA_INSTABILITY = {'t': 0.35, 's': -0.15, 'default': 0.0}  # Negative for s bound boost

def calculate_radius(n, p, epsilon):
    return (2 * n + 1) ** p * (1 + epsilon * n * (n - 1))

def apply_instability_correction(mass, particle):
    eta = ETA_INSTABILITY.get(particle, ETA_INSTABILITY['default'])
    return mass * (1 - eta * LAMBDA_QCD / (mass + 1e-6))

def predict_lepton_masses():
    preds = {}
    for n, part in enumerate(['electron', 'muon', 'tau']):
        a_n = calculate_radius(n, PHI, EPSILON_LEPTON)
        mass = PDG_MASSES['electron'] * a_n**3 if n > 0 else PDG_MASSES['electron']
        preds[part] = mass
    a_3 = calculate_radius(3, PHI, EPSILON_LEPTON)
    preds['fourth_lepton'] = PDG_MASSES['electron'] * a_3**3 / 1000  # GeV
    return preds

def predict_quark_masses():
    preds = {}
    p_up = P_AVG_QUARK + 0.5
    a_u = calculate_radius(0, p_up, EPSILON_QUARK)
    preds['u'] = PDG_MASSES['u']
    a_c = calculate_radius(1, p_up, EPSILON_QUARK)
    preds['c'] = apply_instability_correction(PDG_MASSES['u'] * (a_c / a_u)**3, 'c')
    a_t = calculate_radius(2, p_up, EPSILON_QUARK)
    preds['t'] = apply_instability_correction(PDG_MASSES['u'] * (a_t / a_u)**3, 't')
    a_4u = calculate_radius(3, p_up, EPSILON_QUARK)
    preds['fourth_up'] = PDG_MASSES['u'] * (a_4u / a_u)**3 / 1e6  # TeV

    p_down = P_AVG_QUARK - 0.5
    a_d = calculate_radius(0, p_down, EPSILON_QUARK)
    preds['d'] = PDG_MASSES['d']
    a_s = calculate_radius(1, p_down, EPSILON_QUARK)
    preds['s'] = apply_instability_correction(PDG_MASSES['d'] * (a_s / a_d)**3, 's')
    a_b = calculate_radius(2, p_down, EPSILON_QUARK)
    preds['b'] = apply_instability_correction(PDG_MASSES['d'] * (a_b / a_d)**3, 'b')
    a_4d = calculate_radius(3, p_down, EPSILON_QUARK)
    preds['fourth_down'] = PDG_MASSES['d'] * (a_4d / a_d)**3 / 1000  # GeV
    return preds

def predict_baryon_masses():
    a_s = PHI * A_L
    kappa_s = KAPPA / PHI**2
    eta = ZETA * PHI
    zeta_l = ZETA / PHI
    delta_int = BETA * math.log(a_s / A_L) if A_L != 0 else 0

    preds = {}
    preds['proton'] = 3 * KAPPA * A_L**3
    preds['lambda'] = 2 * KAPPA * A_L**3 + kappa_s * a_s**3 + 2 * zeta_l * a_s**3 * (1 + delta_int)
    preds['sigma'] = 2 * KAPPA * A_L**3 + kappa_s * a_s**3 + 2 * ZETA * a_s**3 * (1 + delta_int)
    preds['xi'] = KAPPA * A_L**3 + 2 * kappa_s * a_s**3 + 2 * ZETA * a_s**3 * (1 + delta_int) + eta * a_s**3 * (1 + delta_int)
    preds['omega'] = 3 * kappa_s * a_s**3 + 3 * eta * a_s**3 * (1 + delta_int)
    return preds

def estimate_free_quark_lifetimes():
    tau_min = H_BAR_MEV_S / 300
    tau_max = H_BAR_MEV_S / 200
    return tau_min, tau_max

def compute_percent_error(pred, actual):
    return 100 * abs(pred - actual) / actual if actual != 0 else float('inf')

def get_unit(mass):
    if mass > 1e6:
        return 'TeV', mass / 1e6
    elif mass > 1e3:
        return 'GeV', mass / 1e3
    return 'MeV', mass

def print_results(predictions, category):
    print(f"\n{category} Masses:")
    print(f"{'Particle':15} {'Predicted':12} {'Unit':5} {'PDG Actual':12} {'Unit':5} {'% Error':8}")
    print("-" * 60)
    for particle, pred in predictions.items():
        actual = PDG_MASSES.get(particle, None)
        if actual is not None:
            error = compute_percent_error(pred, actual)
            p_unit, p_val = get_unit(pred)
            a_unit, a_val = get_unit(actual)
            error_str = f"{error:.2f}" if error != float('inf') else '--'
            print(f"{particle:15} {p_val:.3f} {p_unit:5} {a_val:.3f} {a_unit:5} {error_str:8}")
        else:
            p_unit, p_val = get_unit(pred)
            print(f"{particle:15} {p_val:.3f} {p_unit:5} {'--':12} {'--':5} {'--':8}")

# Run
print("Improved Vortex Model Predictions with Instability Correction")
lepton_preds = predict_lepton_masses()
print_results(lepton_preds, "Leptons")

quark_preds = predict_quark_masses()
print_results(quark_preds, "Quarks")

baryon_preds = predict_baryon_masses()
print_results(baryon_preds, "Baryons")

tau_min, tau_max = estimate_free_quark_lifetimes()
print("\nFree Quark Lifetimes: {0:.2e} to {1:.2e} s".format(tau_min, tau_max))
