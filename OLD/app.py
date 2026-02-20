"""
Hypersaline pH Custom Buffer Calculator v1.3
Every formula traced cell-by-cell from the Excel file.

Formula source map:
  data!D10:D21   → user ion inputs (mg/L) and physical properties
  data!J15:J20   → titration setup
  data!E19,D19   → TDS unit and density (feed hidden!M10)
  data!E26,E27,E29,E32 → hydration state selections
  hidden!C3:C8   → ion mmol/kgw  = data!Dx / hidden!M10 / hidden!Qx
  hidden!M10     → water_mass = IF(N5="g/L",(1000*M4-M5)/1000,(1000*M4-M5*M4)/1000)
  hidden!C12:C18 → salt mmol/kgw (recipe)
  hidden!D12:D18 → salt MW (hydration-dependent)
  hidden!E12:E18 → salt g/kgw = D*C/1000
  hidden!I17     → remaining H2O = 1000 - J18*I18 - J19*I19 - J20*I20 - J21*I21
  hidden!M6      → H3BO3 total mmol = data!J15 * data!J17
  hidden!M7      → H3BO3 mmol/kgw  = M6 / (data!J16 * M10 / 1000)   → input!B41
  hidden!M9      → NaOH  total mmol = data!J18 * data!J19
  hidden!M8      → NaOH  mmol/kgw  = M9 / (data!J16 * M10 / 1000)   → input!B50
  input!C25:C31  → REACTION 1 amounts = hidden!C12:C18 / hidden!C23 (=200)
  input!C40      → H3BO3 H2O ratio = 55.5556 / data!J15
  input!C49      → NaOH  H2O ratio = 55.5556 / data!J18
  input!E50      → NaOH steps = data!J19 / data!J20 = always 20
  input!E21      → CO2 log = LOG(0.000426)
"""

import math
import traceback
from flask import Flask, render_template, request, jsonify

app = Flask(__name__)

# ── Atomic weights  (hidden!Q3:Q13) ──────────────────────────────────────────
Q3  = 22.989769   # Na
Q4  = 39.0983     # K
Q5  = 6.941       # Li
Q6  = 24.305      # Mg
Q7  = 40.078      # Ca
Q8  = 32.02       # S  (SO4 reported as mg-S/L)
Q9  = 35.453      # Cl
Q10 = 18.01528    # H2O
Q11 = 1.0         # H

# ── Salt molar masses  (hidden!J3:J13) ───────────────────────────────────────
J3  = 95.211                  # MgCl2  anhydrous
J4  = J3  + 6  * Q10         # MgCl2·6H2O  = 203.302
J6  = 110.984                 # CaCl2  anhydrous
J7  = J6  + 2  * Q10         # CaCl2·2H2O  = 147.015
J9  = 120.366                 # MgSO4  anhydrous
J10 = J9  + 7  * Q10         # MgSO4·7H2O  = 246.472
J12 = 142.04                  # Na2SO4 anhydrous
J13 = 322.2                   # Na2SO4·10H2O
D14 = 42.394                  # LiCl   (always anhydrous, hidden!D14 is constant)
D16 = 58.44                   # NaCl   (always anhydrous, hidden!D16 is constant)
D17 = 74.55                   # KCl    (always anhydrous, hidden!D17 is constant)

CO2_LOG = math.log10(0.000426)   # input!E21 = LOG(0.000426) = -3.370590401
RF = 200                          # hidden!C23  relative factor for REACTION amounts


# ─────────────────────────────────────────────────────────────────────────────
# WATER MASS  (hidden!M10)
# =IF(N5="g/L", (1000*M4-M5)/1000, (1000*M4-M5*M4)/1000)
#   M4 = density (data!D19, unit kgs/L)
#   M5 = TDS     (data!D20)
#   N5 = TDS unit (data!E20)  → "g/L" or "g/kgs"
# ─────────────────────────────────────────────────────────────────────────────
def calc_water_mass(density, tds, tds_unit):
    if tds_unit == 'g/L':
        return (1000 * density - tds) / 1000
    else:   # g/kgs
        return (1000 * density - tds * density) / 1000


# ─────────────────────────────────────────────────────────────────────────────
# ION → mmol/kgw  (hidden!C3:C8)
# hidden!C3 = data!D10 / $M$10 / hidden!Q3
#           = Na_mg_L  / water_mass / MW_Na
# NOTE: formula divides by M10 first, then by MW (order matters for clarity)
# ─────────────────────────────────────────────────────────────────────────────
def to_mmol_kgw(mg_per_L, MW, water_mass):
    return mg_per_L / water_mass / MW


# ─────────────────────────────────────────────────────────────────────────────
# BUFFER RECIPE  (hidden!B12:E18, data!E26:E32, hidden!I17)
#
# Flags:
#   A18 = IF(C8>C6, 1, 0)        SO4_mmol > Mg_mmol → use Na2SO4 path
#   A15 = IF(A18=1, 0, 1)        inverse → use MgSO4 path
#
# Salt mmol/kgw (hidden!C12:C18):
#   C12 MgCl2  = C6 - C8*A15          Mg minus SO4 if MgSO4 path, else just Mg
#   C13 CaCl2  = C7                    Ca
#   C14 LiCl   = C5                    Li
#   C15 MgSO4  = C8  [* A15 implicit]  SO4 (zero when Na2SO4 path via A15)
#   C16 NaCl   = C3 - 2*C18*A18       Na minus 2×SO4 if Na2SO4 path
#   C17 KCl    = C4                    K
#   C18 Na2SO4 = C8  [* A18 implicit]  SO4 (zero when MgSO4 path via A18)
#
# Salt MW (hidden!D12:D18) depends on data!E26:E32 (hydration dropdowns):
#   D12 = IF(data!I26="MgCl₂",  J3, J4)    anhydrous vs ·6H2O
#   D13 = IF(data!I27="CaCl₂",  J6, J7)    anhydrous vs ·2H2O
#   D15 = IF(data!I29="MgSO₄",  J9, J10)   anhydrous vs ·7H2O
#   D18 = IF(data!I31="Na₂SO₄", J12,J13)   anhydrous vs ·10H2O
#   (data!I26 = IF(E26="Anhydrous","MgCl₂","MgCl₂·6H₂O") etc.)
#
# g/kgw = D*C/1000  (hidden!E12:E18)
#
# Hydration water removed (hidden!I17):
#   I18 = C12*6*Q10/1000    water g from MgCl2·6H2O
#   J18 = IF(data!C26="MgCl₂",0,1)  0 if anhydrous (data!C26 shows the salt name)
#   I17 = 1000 - J18*I18 - J19*I19 - J20*I20 - J21*I21
# ─────────────────────────────────────────────────────────────────────────────
def build_recipe(ion_mmol_kgw, hydration):
    """
    ion_mmol_kgw : {Na, K, Li, Mg, Ca, SO4}  in mmol/kgw
    hydration    : {MgCl2:'Hexahydrate'/'Anhydrous',
                    CaCl2:'Anhydrous'/'Dihydrate',
                    MgSO4:'Heptahydrate'/'Anhydrous',
                    Na2SO4:'Anhydrous'/'Decahydrate'}
    """
    C3 = ion_mmol_kgw['Na']   # hidden!C3
    C4 = ion_mmol_kgw['K']    # hidden!C4
    C5 = ion_mmol_kgw['Li']   # hidden!C5
    C6 = ion_mmol_kgw['Mg']   # hidden!C6
    C7 = ion_mmol_kgw['Ca']   # hidden!C7
    C8 = ion_mmol_kgw['SO4']  # hidden!C8

    # Flags  (hidden!A18, A15)
    A18 = 1 if C8 > C6 else 0
    A15 = 0 if A18 == 1 else 1

    # Salt mmol/kgw  (hidden!C12:C18)
    C12 = max(0.0, C6 - C8 * A15)        # MgCl2
    C13 = max(0.0, C7)                    # CaCl2
    C14 = max(0.0, C5)                    # LiCl
    C15 = max(0.0, C8) * A15             # MgSO4  (0 when Na2SO4 path)
    C16 = max(0.0, C3 - 2 * C8 * A18)   # NaCl
    C17 = max(0.0, C4)                    # KCl
    C18 = max(0.0, C8) * A18             # Na2SO4 (0 when MgSO4 path)

    # Salt MW based on hydration choice  (hidden!D12:D18)
    # Exact match against the dropdown option values in the UI:
    #   MgCl2:  'Anhydrous' → J3=95.211   | 'Hexahydrate'  → J4=203.302
    #   CaCl2:  'Anhydrous' → J6=110.984  | 'Dihydrate'    → J7=147.015
    #   MgSO4:  'Anhydrous' → J9=120.366  | 'Heptahydrate' → J10=246.472
    #   Na2SO4: 'Anhydrous' → J12=142.04  | 'Heptahydrate' → J13=322.2 (10H2O in Excel)
    def is_anhydrous(salt_key):
        return hydration.get(salt_key, 'Anhydrous') == 'Anhydrous'

    D12 = J3  if is_anhydrous('MgCl2')  else J4    # anhydrous → J3, hexahydrate  → J4
    D13 = J6  if is_anhydrous('CaCl2')  else J7    # anhydrous → J6, dihydrate    → J7
    D15 = J9  if is_anhydrous('MgSO4')  else J10   # anhydrous → J9, heptahydrate → J10
    D18 = J12 if is_anhydrous('Na2SO4') else J13   # anhydrous → J12, heptahydrate → J13

    # g/kgw  (hidden!E12:E18)
    E12 = D12 * C12 / 1000
    E13 = D13 * C13 / 1000
    E14 = D14 * C14 / 1000
    E15 = D15 * C15 / 1000
    E16 = D16 * C16 / 1000
    E17 = D17 * C17 / 1000
    E18 = D18 * C18 / 1000

    # Hydration water (hidden!I18:I21, J18:J21, I17)
    # J18 = IF(data!C26="MgCl₂", 0, 1)  → 0 when anhydrous (name is plain MgCl₂)
    J18 = 0 if is_anhydrous('MgCl2')  else 1
    J19 = 0 if is_anhydrous('CaCl2')  else 1
    J20 = 0 if is_anhydrous('MgSO4')  else 1
    J21 = 0 if is_anhydrous('Na2SO4') else 1

    I18 = C12 * 6  * Q10 / 1000   # hidden!I18 = C12*G18*Q10/1000
    I19 = C13 * 2  * Q10 / 1000   # hidden!I19
    I20 = C15 * 7  * Q10 / 1000   # hidden!I20
    I21 = C18 * 10 * Q10 / 1000   # hidden!I21 = G21*C18*Q10/1000

    I16 = 1000.0   # hidden!I16 = 1000 g water base
    I17 = I16 - J18*I18 - J19*I19 - J20*I20 - J21*I21   # remaining H2O g

    # data!J33 = hidden!I17/hidden!Q10   (mmol H2O)
    # data!K33 = hidden!I17              (g H2O)
    water_g   = max(0.0, I17)
    water_mmol = water_g / Q10

    return {
        'MgCl₂':  {'mmol': round(C12, 4), 'g': round(E12, 4), 'mw': round(D12, 3), 'form': hydration.get('MgCl2',  'Anhydrous')},
        'CaCl₂':  {'mmol': round(C13, 4), 'g': round(E13, 4), 'mw': round(D13, 3), 'form': hydration.get('CaCl2',  'Anhydrous')},
        'LiCl':   {'mmol': round(C14, 4), 'g': round(E14, 4), 'mw': round(D14, 3), 'form': 'Anhydrous'},
        'MgSO₄':  {'mmol': round(C15, 4), 'g': round(E15, 4), 'mw': round(D15, 3), 'form': hydration.get('MgSO4',  'Anhydrous')},
        'NaCl':   {'mmol': round(C16, 4), 'g': round(E16, 4), 'mw': round(D16, 3), 'form': 'Anhydrous'},
        'KCl':    {'mmol': round(C17, 4), 'g': round(E17, 4), 'mw': round(D17, 3), 'form': 'Anhydrous'},
        'Na₂SO₄': {'mmol': round(C18, 4), 'g': round(E18, 4), 'mw': round(D18, 3), 'form': hydration.get('Na2SO4', 'Anhydrous')},
        'H₂O':    {'mmol': round(water_mmol, 2), 'g': round(water_g, 2), 'mw': round(Q10, 3), 'form': 'liquid'},
    }


# ─────────────────────────────────────────────────────────────────────────────
# PHREEQC INPUT STRING  (mirrors input sheet cell-by-cell)
#
# REACTION 1 species names  (input!B25:B31 = hidden!S3:S9):
#   MgCl2, CaCl2, LiCl, MgSO4, NaCl, KCl, Na2SO4
#
# REACTION 1 amounts  (input!C25:C31 = hidden!C12:C18 / input!B33):
#   = salt_mmol_kgw / RF   (RF = hidden!C23 = 200)
#
# REACTION 2 H3BO3 amount  (input!B41 = hidden!M7):
#   hidden!M6 = data!J15 * data!J17        = H3BO3_conc * H3BO3_vol  (mmol total)
#   hidden!M7 = M6 / (data!J16 * M10/1000) = mmol total / kgw_in_sample
#
# REACTION 2 NaOH amount  (input!B50 = hidden!M8):
#   hidden!M9 = data!J18 * data!J19        = NaOH_conc * NaOH_vol  (mmol total)
#   hidden!M8 = M9 / (data!J16 * M10/1000) = mmol total / kgw_in_sample
#
# NaOH steps  (input!E50 = data!J19/data!J20 = NaOH_vol/(NaOH_vol/20) = 20)
# H3BO3 H2O  (input!C40 = 55.5556/data!J15)
# NaOH  H2O  (input!C49 = 55.5556/data!J18)
# ─────────────────────────────────────────────────────────────────────────────
def build_phreeqc_input(ion_mmol_kgw, params, water_mass):
    C3 = ion_mmol_kgw['Na']
    C4 = ion_mmol_kgw['K']
    C5 = ion_mmol_kgw['Li']
    C6 = ion_mmol_kgw['Mg']
    C7 = ion_mmol_kgw['Ca']
    C8 = ion_mmol_kgw['SO4']

    A18 = 1 if C8 > C6 else 0
    A15 = 0 if A18 == 1 else 1

    C12 = max(0.0, C6 - C8 * A15)
    C13 = max(0.0, C7)
    C14 = max(0.0, C5)
    C15 = max(0.0, C8) * A15
    C16 = max(0.0, C3 - 2 * C8 * A18)
    C17 = max(0.0, C4)
    C18 = max(0.0, C8) * A18

    # data!J15:J20
    H3BO3_conc = params['H3BO3_conc']   # J15
    sample_vol = params['sample_vol']   # J16
    H3BO3_vol  = params['H3BO3_vol']    # J17
    NaOH_conc  = params['NaOH_conc']    # J18
    NaOH_vol   = params['NaOH_vol']     # J19
    # J20 = J19/20  (step size, used only in E50)

    # kgw in the sample aliquot
    kgw_sample = sample_vol * water_mass / 1000   # data!J16 * M10 / 1000

    # hidden!M6, M7, M9, M8
    M6 = H3BO3_conc * H3BO3_vol          # mmol H3BO3 total
    M7 = M6 / kgw_sample                  # mmol/kgw  → input!B41
    M9 = NaOH_conc  * NaOH_vol           # mmol NaOH total
    M8 = M9 / kgw_sample                  # mmol/kgw  → input!B50

    # input!C40, C49
    H3BO3_H2O = 55.5556 / H3BO3_conc
    NaOH_H2O  = 55.5556 / NaOH_conc

    # input!E50 = data!J19 / data!J20 = NaOH_vol / (NaOH_vol/20) = 20
    n_steps = 20

    def r(val):
        return f"{val / RF:.9f}"

    pqi = (
        "SELECTED_OUTPUT 1\n"
        "-molalities\tB(OH)3  B(OH)4- \n"
        "\tB3O3(OH)4-  B4O5(OH)4-2  MgB(OH)4+  CaB(OH)4+\n"
        "-ionic_strength       \ttrue\n"
        "-pH\ttrue\n"
        "-user_punch\ttrue\n"
        "-water\ttrue\n"
        "-alkalinity\ttrue\n"
        "USER_PUNCH \t1\n"
        "-headings \tVolume Density\n"
        "-start\n"
        "10\tPUNCH SOLN_VOL\n"
        "20\tPUNCH RHO\n"
        "30\tEND\n"
        "SOLUTION\t1\n"
        "\ttemp\t20\n"
        "\tpH\t7\n"
        "\tpe\t4\n"
        "\tredox\tpe\n"
        "\tunits\tmol/l\n"
        f"\tC(4)\t1\tCO2(g)\t{CO2_LOG:.9f}\n"
        "\t-water\t1\t#\tkg\n"
        "\n"
        "REACTION\t1\n"
        f"\tMgCl2\t{r(C12)}\n"   # hidden!S3 / (C12/B33)
        f"\tCaCl2\t{r(C13)}\n"   # hidden!S4
        f"\tLiCl\t{r(C14)}\n"    # hidden!S5
        f"\tMgSO4\t{r(C15)}\n"   # hidden!S6
        f"\tNaCl\t{r(C16)}\n"    # hidden!S7
        f"\tKCl\t{r(C17)}\n"     # hidden!S8
        f"\tNa2SO4\t{r(C18)}\n"  # hidden!S9
        "\n"
        f"\t{RF}\tmillimoles\tin \t1\tsteps\n"
        "SAVE\tSolution\t1\n"
        "END\n"
        "USE\tSolution\t1\n"
        "\n"
        "REACTION\t2\n"
        "\tH3BO3\t1\n"
        f"\tH2O\t{H3BO3_H2O:.3f}\n"
        f"\t{M7:.8f}\tmillimoles\tin \t1\tsteps\n"
        "\n"
        "SAVE\tSolution\t1\n"
        "END\n"
        "USE\tSolution\t1\n"
        "\n"
        "REACTION\t2\n"
        "\tNaOH\t1\n"
        f"\tH2O\t{NaOH_H2O:.4f}\n"
        f"\t{M8:.8f}\tmillimoles\tin \t{n_steps}\tsteps\n"
        "\n"
    )

    return pqi, n_steps


# ─────────────────────────────────────────────────────────────────────────────
# DATABASE FINDER
# ─────────────────────────────────────────────────────────────────────────────
def find_database(db_filename):
    import os, glob

    app_dir = os.path.dirname(os.path.abspath(__file__))
    p = os.path.join(app_dir, db_filename)
    if os.path.isfile(p):
        return p

    for root in [r"C:\Program Files\USGS", r"C:\Program Files (x86)\USGS"]:
        for m in glob.glob(os.path.join(root, '**', db_filename), recursive=True):
            return m

    try:
        import phreeqpy
        pkg = os.path.dirname(phreeqpy.__file__)
        for m in glob.glob(os.path.join(pkg, '**', db_filename), recursive=True):
            return m
    except Exception:
        pass

    return os.path.join(app_dir, db_filename)


# ─────────────────────────────────────────────────────────────────────────────
# IPHREEQC RUNNER
# ─────────────────────────────────────────────────────────────────────────────
def run_phreeqc(pqi):
    try:
        import phreeqpy.iphreeqc.phreeqc_dll as mod
        # Uncomment if DLL not found:
        # mod.IPhreeqc.IPHREEQC_DLL_PATH = r"C:\Program Files\USGS\IPhreeqcCOM 3.x.x\bin\IPhreeqc.dll"
        iphreeqc = mod.IPhreeqc()
    except ImportError:
        raise RuntimeError("phreeqpy not installed. Run: pip install phreeqpy")
    except Exception as e:
        raise RuntimeError(f"IPhreeqc init failed: {e}")

    db = find_database('pitzer.dat')
    iphreeqc.load_database(db)
    if iphreeqc.phc_error_count > 0:
        raise RuntimeError(
            f"Could not load {db}\n{iphreeqc.get_error_string()}\n"
            "Copy pitzer.dat from your IPhreeqcCOM installation's database folder."
        )

    print("\n" + "="*60 + "\nPHREEQC INPUT:\n" + "="*60 + "\n" + pqi + "="*60)

    iphreeqc.run_string(pqi)
    if iphreeqc.phc_error_count > 0:
        raise RuntimeError(f"PHREEQC error:\n{iphreeqc.get_error_string()}")

    n_rows = iphreeqc.row_count
    n_cols = iphreeqc.column_count
    arr    = iphreeqc.get_selected_output_array()

    headers = [str(arr[0][c]) for c in range(n_cols)]
    rows = []
    for r in range(1, n_rows):
        row = {}
        for c in range(n_cols):
            v = arr[r][c]
            try:
                row[headers[c]] = float(v)
            except (TypeError, ValueError):
                row[headers[c]] = str(v)
        rows.append(row)
    return rows


# ─────────────────────────────────────────────────────────────────────────────
# B4/B3 RATIO  (data!Q10 formula)
# =SUM(output!M,P,Q,N,O,O) / SUM(output!L,N,N,O,O)
#
# Output column map (from output sheet row 1):
#   col12 L = m_B(OH)3(mol/kgw)
#   col13 M = m_B(OH)4-(mol/kgw)
#   col14 N = m_B3O3(OH)4-(mol/kgw)
#   col15 O = m_B4O5(OH)4-2(mol/kgw)
#   col16 P = m_MgB(OH)4+(mol/kgw)
#   col17 Q = m_CaB(OH)4+(mol/kgw)
# ─────────────────────────────────────────────────────────────────────────────
def b4b3_ratio(row):
    L = row.get('m_B(OH)3(mol/kgw)',     0) or 0
    M = row.get('m_B(OH)4-(mol/kgw)',    0) or 0
    N = row.get('m_B3O3(OH)4-(mol/kgw)', 0) or 0
    O = row.get('m_B4O5(OH)4-2(mol/kgw)',0) or 0
    P = row.get('m_MgB(OH)4+(mol/kgw)',  0) or 0
    Q = row.get('m_CaB(OH)4+(mol/kgw)',  0) or 0
    num = M + P + Q + N + O + O    # SUM(M,P,Q,N,O,O)
    den = L + N + N + O + O        # SUM(L,N,N,O,O)
    return num / den if den > 1e-30 else 0.0


# ─────────────────────────────────────────────────────────────────────────────
# SHARED PAYLOAD PARSER
# ─────────────────────────────────────────────────────────────────────────────
def parse_payload(d):
    density  = float(d['density'])
    tds      = float(d['tds'])
    tds_unit = d['tds_unit']

    # hidden!M10
    wm = calc_water_mass(density, tds, tds_unit)

    # hidden!C3:C8  ion mmol/kgw
    ion = {
        'Na':  to_mmol_kgw(float(d['Na']),  Q3, wm),
        'K':   to_mmol_kgw(float(d['K']),   Q4, wm),
        'Li':  to_mmol_kgw(float(d['Li']),  Q5, wm),
        'Mg':  to_mmol_kgw(float(d['Mg']),  Q6, wm),
        'Ca':  to_mmol_kgw(float(d['Ca']),  Q7, wm),
        'SO4': to_mmol_kgw(float(d['SO4']), Q8, wm),
    }

    hyd = {
        'MgCl2':  d.get('hyd_MgCl2',  'Hexahydrate'),
        'CaCl2':  d.get('hyd_CaCl2',  'Anhydrous'),
        'MgSO4':  d.get('hyd_MgSO4',  'Heptahydrate'),
        'Na2SO4': d.get('hyd_Na2SO4', 'Anhydrous'),
    }

    return wm, ion, hyd


# ─────────────────────────────────────────────────────────────────────────────
# FLASK ROUTES
# ─────────────────────────────────────────────────────────────────────────────
@app.route('/')
def index():
    return render_template('index.html')


@app.route('/recipe', methods=['POST'])
def recipe_route():
    """Instant recipe — no PHREEQC, just the hidden sheet math."""
    try:
        d = request.get_json()
        wm, ion, hyd = parse_payload(d)
        return jsonify({'recipe': build_recipe(ion, hyd), 'water_mass': round(wm, 6)})
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/calculate', methods=['POST'])
def calculate():
    try:
        d = request.get_json()
        wm, ion, hyd = parse_payload(d)

        # Add B and Br for display (not in recipe but shown in ion table)
        ion_full = dict(ion)
        ion_full['B']  = to_mmol_kgw(float(d['B']),  10.811, wm)
        ion_full['Br'] = to_mmol_kgw(float(d['Br']), 79.904, wm)

        params = {
            'H3BO3_conc': float(d['H3BO3_conc']),
            'H3BO3_vol':  float(d['H3BO3_vol']),
            'sample_vol': float(d['sample_vol']),
            'NaOH_conc':  float(d['NaOH_conc']),
            'NaOH_vol':   float(d['NaOH_vol']),
            'pH':         float(d.get('pH', 8.5)),
        }

        pqi, n_steps = build_phreeqc_input(ion, params, wm)
        rows   = run_phreeqc(pqi)
        recipe = build_recipe(ion, hyd)

        step_ml = float(d['NaOH_vol']) / 20
        titration = []
        naoh_vol  = 0.0
        for row in rows:
            state = row.get('state', '')
            pH    = row.get('pH')
            ratio = b4b3_ratio(row)
            if state == 'i_soln':
                v = 0.0
            elif state == 'react':
                naoh_vol = round(naoh_vol + step_ml, 8)
                v = naoh_vol
            else:
                continue
            titration.append({
                'state':  state,
                'V_NaOH': v,
                'pH':     round(pH, 5) if pH is not None else None,
                'B4B3':   round(ratio, 7),
            })

        return jsonify({
            'titration':    titration,
            'recipe':       recipe,
            'water_mass':   round(wm, 6),
            'ion_mmol_kgw': {k: round(v, 5) for k, v in ion_full.items()},
            'n_steps':      n_steps,
        })

    except RuntimeError as e:
        return jsonify({'error': str(e)}), 500
    except Exception as e:
        return jsonify({'error': str(e), 'trace': traceback.format_exc()}), 500


@app.route('/show_input')
def show_input():
    """Preview the PHREEQC input string at http://localhost:5000/show_input"""
    d = {k: request.args.get(k, v) for k, v in {
        'Na':'2000','K':'5000','Li':'500','Mg':'20','Ca':'250',
        'SO4':'50','density':'1.3','tds':'350','tds_unit':'g/L',
        'H3BO3_conc':'0.4','H3BO3_vol':'4','sample_vol':'50',
        'NaOH_conc':'1','NaOH_vol':'2',
        'hyd_MgCl2':'Hexahydrate','hyd_CaCl2':'Anhydrous',
        'hyd_MgSO4':'Heptahydrate','hyd_Na2SO4':'Anhydrous',
    }.items()}
    wm, ion, hyd = parse_payload(d)
    params = {k: float(d[k]) for k in
              ['H3BO3_conc','H3BO3_vol','sample_vol','NaOH_conc','NaOH_vol']}
    pqi, n = build_phreeqc_input(ion, params, wm)
    return (
        '<h2 style="font-family:monospace">PHREEQC Input Preview</h2>'
        f'<p style="font-family:monospace">water_mass = {wm:.6f} kgw/L &nbsp;|&nbsp; n_steps = {n}</p>'
        '<pre style="font-family:monospace;font-size:13px;background:#111;color:#0f0;padding:20px;border-radius:8px">'
        + pqi.replace('<','&lt;').replace('>','&gt;') + '</pre>'
    )


@app.route('/debug')
def debug():
    """IPhreeqc diagnostics at http://localhost:5000/debug"""
    import os
    info = {}
    try:
        import phreeqpy
        info['phreeqpy_version']  = getattr(phreeqpy, '__version__', 'unknown')
        info['phreeqpy_location'] = phreeqpy.__file__
    except ImportError as e:
        return f"<pre>phreeqpy not installed: {e}</pre>"
    try:
        import phreeqpy.iphreeqc.phreeqc_dll as mod
        info['dll_path'] = getattr(mod.IPhreeqc, 'IPHREEQC_DLL_PATH', 'not set')
        obj = mod.IPhreeqc()
        methods = sorted(m for m in dir(obj) if not m.startswith('_'))
        info['all_methods'] = methods
    except Exception as e:
        info['init_error'] = str(e)

    for fn in ['pitzer.dat', 'phreeqc.dat']:
        p = find_database(fn)
        info[fn] = p + (' ✓' if os.path.isfile(p) else ' ✗ NOT FOUND')

    html = ['<h2 style="font-family:monospace">IPhreeqc Debug</h2><pre style="font-family:monospace;font-size:13px;line-height:1.6">']
    for k, v in info.items():
        if isinstance(v, list):
            html.append(f"{k}:\n  " + "\n  ".join(str(x) for x in v) + "\n")
        else:
            html.append(f"{k}:  {v}\n")
    html.append('</pre>')
    return '\n'.join(html)


if __name__ == '__main__':
    print("=" * 62)
    print("  Hypersaline pH Custom Buffer Calculator v1.3")
    print("  App:         http://localhost:5000")
    print("  PHREEQC in:  http://localhost:5000/show_input")
    print("  Diagnostics: http://localhost:5000/debug")
    print("=" * 62)
    app.run(debug=True, port=5000)
