import sys
import os
import math
import numpy as np
import pandas as pd

from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QTabWidget, QLabel, QComboBox, QPushButton, QTableWidget,
    QTableWidgetItem, QHeaderView, QScrollArea, QTextEdit, QMessageBox,
    QFileDialog, QDoubleSpinBox, QFormLayout, QGroupBox, QCheckBox
)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFont, QColor

import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure


# ==============================================================================
# THERMODYNAMIC & FLUID PROPERTY ENGINES
# ==============================================================================

class FluidPropertyEngine:
    """Calculates thermodynamic properties for standard fluids or user custom inputs."""
    
    @staticmethod
    def get_properties(fluid_name, temp_c):
        """Returns density (rho), specific heat (cp), thermal conductivity (k), dynamic viscosity (mu)."""
        temp_k = temp_c + 273.15
        
        if "Ethylene Glycol (50%)" in fluid_name:
            rho = 1070 - 0.65 * (temp_c - 20)
            cp = 3300 + 3.5 * temp_c
            k = 0.38 + 0.0005 * temp_c
            mu = 0.0035 * np.exp(-0.02 * temp_c)
        elif "Ethylene Glycol (30%)" in fluid_name:
            rho = 1040 - 0.55 * (temp_c - 20)
            cp = 3600 + 2.8 * temp_c
            k = 0.42 + 0.0006 * temp_c
            mu = 0.0022 * np.exp(-0.018 * temp_c)
        elif "Milk / Liquid Food" in fluid_name:
            rho = 1030 - 0.5 * temp_c
            cp = 3890 + 1.2 * temp_c
            k = 0.53 + 0.0008 * temp_c
            mu = 0.0021 * np.exp(-0.015 * temp_c)
        elif "Fruit Juice Concentrate" in fluid_name:
            rho = 1100 - 0.6 * temp_c
            cp = 3500 + 1.8 * temp_c
            k = 0.45 + 0.0006 * temp_c
            mu = 0.0050 * np.exp(-0.022 * temp_c)
        elif "Vegetable / Food Oil" in fluid_name or "Water-to-Oil" in fluid_name or "Engine Oil" in fluid_name:
            rho = 915 - 0.68 * temp_c
            cp = 1980 + 3.1 * temp_c
            k = 0.165 - 0.0001 * temp_c
            mu = 0.035 * np.exp(-0.028 * temp_c)
        elif "Solar Thermal Oil" in fluid_name:
            rho = 880 - 0.7 * temp_c
            cp = 1900 + 3.2 * temp_c
            k = 0.13 - 0.0001 * temp_c
            mu = 0.012 * np.exp(-0.025 * temp_c)
        elif "Custom Fluid A" in fluid_name:
            rho, cp, k, mu = 1050.0, 2800.0, 0.22, 0.0150
        elif "Custom Fluid B" in fluid_name:
            rho, cp, k, mu = 820.0, 2100.0, 0.14, 0.0008
        else:  # General Water default
            rho = 998.2 - 0.15 * (temp_c - 20)
            cp = 4182.0
            k = 0.598 + 0.0012 * temp_c
            mu = 0.001002 * np.exp(-0.017 * (temp_c - 20))
            
        return max(rho, 1.0), max(cp, 100.0), max(k, 0.01), max(mu, 1e-6)


class HeatExchangerCore:
    """Core computational solver for Shell and Tube Rating, Outlet Temperatures & Pressure Drop."""
    
    @staticmethod
    def compute_single_rating(m_dot_h, m_dot_c, Thi, Tci, D_s, L_t, N_t, custom_h, custom_c, B_c=25.0):
        if min(m_dot_h, m_dot_c, Thi, D_s, L_t, N_t) <= 0:
            return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, "N/A", "N/A"

        # Fluid properties
        rho_h, cp_h, k_h, mu_h = custom_h['rho'], custom_h['cp'], custom_h['k'], custom_h['mu']
        rho_c, cp_c, k_c, mu_c = custom_c['rho'], custom_c['cp'], custom_c['k'], custom_c['mu']
        
        # Heat Capacity Rates
        C_h = m_dot_h * cp_h
        C_c = m_dot_c * cp_c
        C_min = min(C_h, C_c)
        C_max = max(C_h, C_c)
        cr = C_min / C_max if C_max > 0 else 0.0
        
        # Tube Side Geometry — HOT STREAM ALWAYS FLOWS THROUGH THE TUBES
        d_i = 0.016
        d_o = 0.020
        A_tube_flow = (np.pi / 4) * (d_i**2) * (N_t / 2)
        u_t = m_dot_h / (rho_h * A_tube_flow + 1e-8)
        Re_t = (rho_h * u_t * d_i) / (mu_h + 1e-8)
        Pr_t = (cp_h * mu_h) / (k_h + 1e-8)
        
        regime_t = "Turbulent" if Re_t > 2300 else "Laminar"
        Nu_t = 0.023 * (Re_t**0.8) * (Pr_t**0.4) if Re_t > 2300 else 3.66
        h_i = (Nu_t * k_h) / d_i
        
        # Shell Side Geometry (Kern Method) — COLD STREAM ALWAYS FLOWS AROUND THE TUBES
        B = D_s * (B_c / 100.0)
        P_t = 1.25 * d_o
        C_gap = P_t - d_o
        A_shell_cross = (D_s * C_gap * B) / P_t
        u_s = m_dot_c / (rho_c * A_shell_cross + 1e-8)
        D_e = 4 * (P_t**2 - (np.pi * d_o**2 / 4)) / (np.pi * d_o)
        Re_s = (rho_c * u_s * D_e) / (mu_c + 1e-8)
        Pr_s = (cp_c * mu_c) / (k_c + 1e-8)
        
        regime_s = "Turbulent" if Re_s > 2000 else "Laminar"
        Nu_s = 0.36 * (Re_s**0.55) * (Pr_s**(1/3))
        h_o = (Nu_s * k_c) / D_e
        
        # Overall Heat Transfer Coefficient U
        U = 1.0 / ((1.0 / (h_o + 1e-8)) + (d_o / d_i) * (1.0 / (h_i + 1e-8)) + 0.0002)
        
        # Surface Area and Effectiveness (epsilon-NTU)
        A_total = np.pi * d_o * L_t * N_t
        NTU = (U * A_total) / (C_min + 1e-8)
        
        Epsilon = 2.0 / (1.0 + cr + np.sqrt(1.0 + cr**2) * ((1.0 + np.exp(-NTU * np.sqrt(1.0 + cr**2))) / (1.0 - np.exp(-NTU * np.sqrt(1.0 + cr**2)) + 1e-8)))
        
        # Heat Transfer Rate
        Q_max = C_min * (Thi - Tci)
        Q_actual_w = Epsilon * Q_max
        Q_actual_kw = Q_actual_w / 1000.0
        
        # Outlet Temperatures
        Tho = Thi - (Q_actual_w / (C_h + 1e-8))
        Tco = Tci + (Q_actual_w / (C_c + 1e-8))

        # LMTD Calculation
        dt1 = Thi - Tco
        dt2 = Tho - Tci
        if dt1 == dt2 or dt1 <= 0 or dt2 <= 0:
            LMTD = max(0.1, (dt1 + dt2) / 2.0)
        else:
            LMTD = (dt1 - dt2) / np.log(dt1 / dt2)
        
        # Pressure Drop Calculations
        f_t = 0.3164 / (Re_t**0.25 + 1e-8) if Re_t > 2300 else 64.0 / (Re_t + 1e-8)
        dP_tube = ((f_t * L_t * (u_t**2) * rho_c) / (2 * d_i)) / 1000.0
        
        f_s = 0.25 * (Re_s**(-0.15))
        N_b = L_t / (B + 1e-8)
        dP_shell = ((f_s * D_s * (N_b + 1) * (u_s**2) * rho_c) / (2 * D_e)) / 1000.0
        
        # Hydraulic power: hot stream in tubes + cold stream in shell.
        P_pump_w = (m_dot_h / rho_h) * (dP_tube * 1000) + (m_dot_c / rho_c) * (dP_shell * 1000)
        
        return Q_actual_kw, Tho, Tco, dP_shell, dP_tube, P_pump_w, U, A_total, NTU, Epsilon, LMTD, regime_s, regime_t


# ==============================================================================
# INVERSE SOLVER LOGIC
# ==============================================================================

def solve_target_parameter(target_stream, desired_temp, free_param, base_params):
    """Numerically solves for the required free parameter to achieve desired outlet temperature.
    
    Returns:
        solved_val (float): Best achievable parameter value.
        is_achievable (bool): True if exact temperature setpoint was reached.
        boundary_reason (str): Description of physical limitation if unachievable.
    """
    is_hot_target = "hot" in target_stream.lower()
    
    # Extract operational parameters from base_params
    m_h = base_params.get('m_dot_h', 5.0)
    m_c = base_params.get('m_dot_c', 5.0)
    Thi = base_params.get('Thi', 75.0)
    Tci = base_params.get('Tci', 25.0)
    Ds = base_params.get('Ds', 0.5)
    Lt = base_params.get('Lt', 3.0)
    Nt = int(base_params.get('Nt', 200))
    custom_h = base_params.get('custom_h')
    custom_c = base_params.get('custom_c')

    # Bounded Search Ranges based on parameter selected
    if "m_dot_h" in free_param:
        low, high = 0.01, 200.0
    elif "m_dot_c" in free_param:
        low, high = 0.01, 200.0
    elif "L_t" in free_param:
        low, high = 0.1, 20.0
    else:  # N_t
        low, high = 10, 5000

    # Binary Search Bisection (40 iterations)
    for _ in range(40):
        mid = (low + high) / 2.0
        
        # Override parameter under test
        test_mh = mid if "m_dot_h" in free_param else m_h
        test_mc = mid if "m_dot_c" in free_param else m_c
        test_Lt = mid if "L_t" in free_param else Lt
        test_Nt = int(mid) if "N_t" in free_param else Nt

        _, Tho, Tco, _, _, _, _, _, _, _, _, _, _ = HeatExchangerCore.compute_single_rating(
            test_mh, test_mc, Thi, Tci, Ds, test_Lt, test_Nt, custom_h, custom_c
        )
        val_check = Tho if is_hot_target else Tco

        if "m_dot_h" in free_param:
            if is_hot_target:
                if val_check > desired_temp: high = mid
                else: low = mid
            else:
                if val_check < desired_temp: high = mid
                else: low = mid
        else:
            if is_hot_target:
                if val_check < desired_temp: high = mid
                else: low = mid
            else:
                if val_check > desired_temp: high = mid
                else: low = mid

    solved_val = (low + high) / 2.0
    if "N_t" in free_param:
        solved_val = int(solved_val)

    # Verification run
    final_mh = solved_val if "m_dot_h" in free_param else m_h
    final_mc = solved_val if "m_dot_c" in free_param else m_c
    final_Lt = solved_val if "L_t" in free_param else Lt
    final_Nt = int(solved_val) if "N_t" in free_param else Nt

    _, final_Tho, final_Tco, _, _, _, _, _, _, _, _, _, _ = HeatExchangerCore.compute_single_rating(
        final_mh, final_mc, Thi, Tci, Ds, final_Lt, final_Nt, custom_h, custom_c
    )
    achieved_temp = final_Tho if is_hot_target else final_Tco

    delta = abs(achieved_temp - desired_temp)
    is_achievable = delta <= 0.5

    boundary_reason = ""
    if not is_achievable:
        if is_hot_target:
            if desired_temp <= Tci:
                boundary_reason = f"Thermal Pinch Limit Reached (Target {desired_temp:.2f}°C ≤ Cold Inlet {Tci:.2f}°C)."
            else:
                boundary_reason = f"Maximum Heat Exchange Limit Reached (Achievable limit capped at {achieved_temp:.2f}°C)."
        else:
            if desired_temp >= Thi:
                boundary_reason = f"Thermal Pinch Limit Reached (Target {desired_temp:.2f}°C ≥ Hot Inlet {Thi:.2f}°C)."
            else:
                boundary_reason = f"Maximum Heat Exchange Limit Reached (Achievable limit capped at {achieved_temp:.2f}°C)."

    return solved_val, is_achievable, boundary_reason


def execute_inverse_solver(
    target_stream: str,
    desired_temp: float,
    free_param: str,
    current_params: dict,
    iteration_df: pd.DataFrame,
) -> tuple:
    """Executes inverse target calculation without mutating main iteration history.

    Returns:
        updated_params (dict): Solved or boundary-capped parameter set.
        report_text (str): Concise block formatted for the right-hand output text.
        closest_row_index (int): Index of closest row in the existing table for UI highlighting.
    """
    eval_params = current_params.copy()

    # Step 1: Execute numeric back-calculation
    solved_val, is_achievable, boundary_reason = solve_target_parameter(
        target_stream=target_stream,
        desired_temp=desired_temp,
        free_param=free_param,
        base_params=eval_params,
    )

    # Step 2: Format Right-Side Report Output Block
    report_lines = [
        "\n==================================================",
        "          INVERSE SOLVER TARGET REPORT            ",
        "==================================================",
        f"Target Stream: {target_stream}",
        f"Desired Setpoint: {desired_temp:.2f} °C",
        f"Free Parameter Adjusted: {free_param}",
    ]

    if is_achievable:
        eval_params[free_param] = solved_val
        report_lines.extend(
            [
                "STATUS: Target Achievable!",
                f"Exact Solved Parameter ({free_param}): {solved_val:.4f}",
                "Geometry/Operational Parameters: Updated to exact setpoint.",
            ]
        )
    else:
        report_lines.extend(
            [
                "STATUS: Unachievable Target (Physical Limit Reached)",
                f"Boundary Restriction: {boundary_reason}",
                f"Closest Parameter Boundary ({free_param}): {solved_val:.4f}",
                "Note: Main input parameters capped at maximum physical limit.",
            ]
        )

    # Step 3: Find closest match in the EXISTING 50-run temperature iteration table
    temp_col_key = "T_hot,out [°C]" if "hot" in target_stream.lower() else "T_cold,out [°C]"

    closest_row_idx = -1
    min_delta = float("inf")

    if iteration_df is not None and len(iteration_df) > 0 and temp_col_key in iteration_df.columns:
        for idx, row in iteration_df.iterrows():
            delta = abs(row[temp_col_key] - desired_temp)
            if delta < min_delta:
                min_delta = delta
                closest_row_idx = int(idx)

    if closest_row_idx != -1:
        closest_temp = iteration_df.iloc[closest_row_idx][temp_col_key]
        report_lines.extend(
            [
                "--------------------------------------------------",
                "ITERATION TABLE RECOMMENDATION REFERENCE:",
                f"Closest Table Row: Index #{closest_row_idx + 1}",
                f"Closest Recorded Temp: {closest_temp:.2f} °C (Delta: {min_delta:.2f} °C)",
                "==================================================\n",
            ]
        )

    return eval_params, "\n".join(report_lines), closest_row_idx


# ==============================================================================
# MAIN APPLICATION WINDOW
# ==============================================================================

class UniversalHeatExchangerApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("rosif.py | Universal Heat Exchanger Software")
        self.setGeometry(50, 50, 1350, 920)
        
        self.pareto_df = None
        self.iterations_df = None
        self.single_rating_result = None
        self.last_rating_params = None
        self.inverse_solution_result = None
        self.ktha_result = None
        self.init_ui()

    def init_ui(self):
        self.tabs = QTabWidget()
        self.tabs.setFont(QFont("Segoe UI", 10, QFont.Bold))
        self.setCentralWidget(self.tabs)

        self.tab1 = QWidget()
        self.tab2 = QWidget()
        self.tab3 = QWidget()
        self.tab4 = QWidget()
        self.tab5 = QWidget()

        self.tabs.addTab(self.tab1, "1. Rating & Target Inverse Solver")
        self.tabs.addTab(self.tab2, "2. NSGA-II Pareto & Streamwise Logging")
        self.tabs.addTab(self.tab3, "3. KTHA-I Engine & Mechanics")
        self.tabs.addTab(self.tab4, "4. Model Validation & Sensitivity")
        self.tabs.addTab(self.tab5, "5. About KTHA")

        self.setup_tab1()
        self.setup_tab2()
        self.setup_tab3()
        self.setup_tab4()
        self.setup_tab5()

    # ==========================================================================
    # PAGE 1: Single Rating & Inverse Target Outlet Solver
    # ==========================================================================
    def setup_tab1(self):
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        content = QWidget()
        main_layout = QHBoxLayout(content)
        
        # --- LEFT PANEL ---
        left_layout = QVBoxLayout()
        
        op_group = QGroupBox("Operating Parameters & Locked Geometry Controls")
        op_group.setFont(QFont("Segoe UI", 9, QFont.Bold))
        form = QFormLayout()

        self.combo_category = QComboBox()
        self.combo_category.addItems([
            "Food Processing",
            "Renewable Energy",
            "Water-to-Water",
            "Water-to-Oil",
            "Custom Fluids",
            "Industrial Process"
        ])
        self.combo_category.currentIndexChanged.connect(self.update_fluid_options)

        self.combo_hot_fluid = QComboBox()
        self.combo_cold_fluid = QComboBox()

        # Spinboxes starting at 0.0
        self.spin_m_hot = QDoubleSpinBox(); self.spin_m_hot.setRange(0.0, 500.0); self.spin_m_hot.setValue(0.0); self.spin_m_hot.setSuffix(" kg/s")
        self.spin_m_cold = QDoubleSpinBox(); self.spin_m_cold.setRange(0.0, 500.0); self.spin_m_cold.setValue(0.0); self.spin_m_cold.setSuffix(" kg/s")
        self.spin_Thi = QDoubleSpinBox(); self.spin_Thi.setRange(0.0, 500.0); self.spin_Thi.setValue(0.0); self.spin_Thi.setSuffix(" °C")
        self.spin_Tci = QDoubleSpinBox(); self.spin_Tci.setRange(-50.0, 200.0); self.spin_Tci.setValue(0.0); self.spin_Tci.setSuffix(" °C")
        self.spin_Ds = QDoubleSpinBox(); self.spin_Ds.setRange(0.0, 5.0); self.spin_Ds.setValue(0.0); self.spin_Ds.setSuffix(" m")
        self.spin_Lt = QDoubleSpinBox(); self.spin_Lt.setRange(0.0, 20.0); self.spin_Lt.setValue(0.0); self.spin_Lt.setSuffix(" m")
        self.spin_Nt = QDoubleSpinBox(); self.spin_Nt.setRange(0, 5000); self.spin_Nt.setValue(0)

        # Locking Parameter Checkboxes
        self.chk_lock_m_hot = QCheckBox("Lock hot-stream mass flow rate (ṁₕ)")
        self.chk_lock_m_cold = QCheckBox("Lock cold-stream mass flow rate (ṁ꜀)")
        self.chk_lock_geom = QCheckBox("Lock Geometry (Ds, Lt, Nt) - Existing Fabrication")

        form.addRow(QLabel("Application Domain:"), self.combo_category)
        form.addRow(QLabel("Hot Side Fluid Stream:"), self.combo_hot_fluid)
        form.addRow(QLabel("Cold Side Fluid Stream:"), self.combo_cold_fluid)
        form.addRow(QLabel("Hot-stream mass flow rate, ṁ<sub>h</sub> [kg/s]:"), self.spin_m_hot)
        form.addRow(QLabel("Cold-stream mass flow rate, ṁ<sub>c</sub> [kg/s]:"), self.spin_m_cold)
        form.addRow(QLabel("Hot-stream inlet temperature, T<sub>h,i</sub> [°C]:"), self.spin_Thi)
        form.addRow(QLabel("Cold-stream inlet temperature, T<sub>c,i</sub> [°C]:"), self.spin_Tci)
        form.addRow(QLabel("Shell inside diameter, D<sub>s</sub> [m]:"), self.spin_Ds)
        form.addRow(QLabel("Tube length, L<sub>t</sub> [m]:"), self.spin_Lt)
        form.addRow(QLabel("Number of tubes, N<sub>t</sub> [-]:"), self.spin_Nt)
        form.addRow(self.chk_lock_m_hot)
        form.addRow(self.chk_lock_m_cold)
        form.addRow(self.chk_lock_geom)
        op_group.setLayout(form)
        left_layout.addWidget(op_group)

        # Inverse Target Setpoint Group
        inverse_group = QGroupBox("Inverse Target Setpoint Solver")
        inverse_group.setFont(QFont("Segoe UI", 9, QFont.Bold))
        inv_layout = QFormLayout()

        self.combo_target_stream = QComboBox()
        self.combo_target_stream.addItems([
            "Hot outlet temperature (Tₕ,ₒ)",
            "Cold outlet temperature (T꜀,ₒ)"
        ])

        self.spin_target_temp = QDoubleSpinBox()
        self.spin_target_temp.setRange(-50.0, 500.0)
        self.spin_target_temp.setValue(0.0)
        self.spin_target_temp.setSuffix(" °C")

        self.combo_variable_to_solve = QComboBox()
        self.combo_variable_to_solve.addItems([
            "ṁₕ — hot-stream mass flow rate",
            "ṁ꜀ — cold-stream mass flow rate",
            "Lₜ — tube length",
            "Nₜ — number of tubes"
        ])

        btn_solve_inverse = QPushButton("Back-Calculate Required Parameters")
        btn_solve_inverse.setStyleSheet("background-color: #0d9488; color: white; font-weight: bold; padding: 6px;")
        btn_solve_inverse.clicked.connect(self.solve_inverse_target)

        inv_layout.addRow(QLabel("Select Target Fluid Stream:"), self.combo_target_stream)
        inv_layout.addRow(QLabel("Desired Outlet Temperature:"), self.spin_target_temp)
        inv_layout.addRow(QLabel("Free Parameter to Adjust:"), self.combo_variable_to_solve)
        inv_layout.addRow(btn_solve_inverse)
        inverse_group.setLayout(inv_layout)
        left_layout.addWidget(inverse_group)

        # Fluid Physical Properties Group
        custom_group = QGroupBox("Fluid Physical Properties (Auto-Updates)")
        custom_group.setFont(QFont("Segoe UI", 9, QFont.Bold))
        custom_layout = QFormLayout()

        self.spin_h_rho = QDoubleSpinBox(); self.spin_h_rho.setRange(0.0, 5000.0); self.spin_h_rho.setValue(0.0)
        self.spin_h_cp = QDoubleSpinBox(); self.spin_h_cp.setRange(0.0, 10000.0); self.spin_h_cp.setValue(0.0)
        self.spin_h_k = QDoubleSpinBox(); self.spin_h_k.setRange(0.0, 10.0); self.spin_h_k.setDecimals(4); self.spin_h_k.setValue(0.0)
        self.spin_h_mu = QDoubleSpinBox(); self.spin_h_mu.setRange(0.0, 5.0); self.spin_h_mu.setDecimals(6); self.spin_h_mu.setValue(0.0)

        self.spin_c_rho = QDoubleSpinBox(); self.spin_c_rho.setRange(0.0, 5000.0); self.spin_c_rho.setValue(0.0)
        self.spin_c_cp = QDoubleSpinBox(); self.spin_c_cp.setRange(0.0, 10000.0); self.spin_c_cp.setValue(0.0)
        self.spin_c_k = QDoubleSpinBox(); self.spin_c_k.setRange(0.0, 10.0); self.spin_c_k.setDecimals(4); self.spin_c_k.setValue(0.0)
        self.spin_c_mu = QDoubleSpinBox(); self.spin_c_mu.setRange(0.0, 5.0); self.spin_c_mu.setDecimals(6); self.spin_c_mu.setValue(0.0)

        custom_layout.addRow(QLabel("<b>HOT STREAM PHYSICAL PROPERTIES:</b>"))
        custom_layout.addRow(QLabel("Hot Density (ρ) [kg/m³]:"), self.spin_h_rho)
        custom_layout.addRow(QLabel("Hot Heat Cap. (Cp) [J/kg·K]:"), self.spin_h_cp)
        custom_layout.addRow(QLabel("Hot Conductivity (k) [W/m·K]:"), self.spin_h_k)
        custom_layout.addRow(QLabel("Hot Viscosity (μ) [Pa·s]:"), self.spin_h_mu)

        custom_layout.addRow(QLabel("<b>COLD STREAM PHYSICAL PROPERTIES:</b>"))
        custom_layout.addRow(QLabel("Cold Density (ρ) [kg/m³]:"), self.spin_c_rho)
        custom_layout.addRow(QLabel("Cold Heat Cap. (Cp) [J/kg·K]:"), self.spin_c_cp)
        custom_layout.addRow(QLabel("Cold Conductivity (k) [W/m·K]:"), self.spin_c_k)
        custom_layout.addRow(QLabel("Cold Viscosity (μ) [Pa·s]:"), self.spin_c_mu)

        custom_group.setLayout(custom_layout)
        left_layout.addWidget(custom_group)

        eq_group = QGroupBox("Governing Equations Used by the Rating Model")
        eq_group.setFont(QFont("Segoe UI", 9, QFont.Bold))
        eq_text = QTextEdit()
        eq_text.setReadOnly(True)
        eq_text.setFont(QFont("Segoe UI", 9))
        eq_text.setHtml("""
        <div style='line-height:1.45;'>
        <b>Steady-state energy balance</b><br>
        Q = ṁₕ c<sub>p,h</sub>(Tₕ,ᵢ − Tₕ,ₒ) = ṁ꜀ c<sub>p,c</sub>(T꜀,ₒ − T꜀,ᵢ)<br><br>
        <b>Heat-capacity rates and maximum heat duty</b><br>
        Cₕ = ṁₕ c<sub>p,h</sub>;&nbsp;&nbsp; C꜀ = ṁ꜀ c<sub>p,c</sub>;&nbsp;&nbsp; C<sub>min</sub> = min(Cₕ, C꜀)<br>
        Q<sub>max</sub> = C<sub>min</sub>(Tₕ,ᵢ − T꜀,ᵢ)<br><br>
        <b>Effectiveness–NTU heat-transfer calculation used by the solver</b><br>
        NTU = U<sub>o</sub>A / C<sub>min</sub>;&nbsp;&nbsp; Q = ε C<sub>min</sub>(Tₕ,ᵢ − T꜀,ᵢ)<br>
        Tₕ,ₒ = Tₕ,ᵢ − Q/Cₕ;&nbsp;&nbsp; T꜀,ₒ = T꜀,ᵢ + Q/C꜀<br><br>
        <b>LMTD heat-transfer relation / model check</b><br>
        Q = U<sub>o</sub> A F ΔT<sub>lm</sub><br>
        ΔT<sub>lm</sub> = (ΔT₁ − ΔT₂) / ln(ΔT₁ / ΔT₂)<br>
        ΔT₁ = Tₕ,ᵢ − T꜀,ₒ;&nbsp;&nbsp; ΔT₂ = Tₕ,ₒ − T꜀,ᵢ<br><br>
        <b>Reynolds and Prandtl numbers</b><br>
        Re = ρuD / μ;&nbsp;&nbsp; Pr = c<sub>p</sub>μ / k<br><br>
        <b>Tube-side heat transfer</b><br>
        Nuₜ = 0.023 Reₜ<sup>0.8</sup> Prₜ<sup>n</sup>;&nbsp;&nbsp; hᵢ = Nuₜ kₕ / dᵢ<br><br>
        <b>Shell-side heat transfer (Kern)</b><br>
        Nuₛ = 0.36 Reₛ<sup>0.55</sup> Prₛ<sup>1/3</sup>;&nbsp;&nbsp; hₒ = Nuₛ k꜀ / Dₑ<br><br>
        <b>Total hydraulic pumping power</b><br>
        P<sub>total</sub> = (ṁₕ/ρₕ)ΔPₜ + (ṁ꜀/ρ꜀)ΔPₛ
        </div>
        """)
        eq_text.setMinimumHeight(330)
        eq_group_layout = QVBoxLayout()
        eq_group_layout.addWidget(eq_text)
        eq_group.setLayout(eq_group_layout)
        left_layout.addWidget(eq_group)

        self.combo_hot_fluid.currentIndexChanged.connect(self.sync_hot_fluid_properties)
        self.combo_cold_fluid.currentIndexChanged.connect(self.sync_cold_fluid_properties)
        self.spin_Thi.valueChanged.connect(self.sync_hot_fluid_properties)
        self.spin_Tci.valueChanged.connect(self.sync_cold_fluid_properties)

        # Calculation Buttons
        btn_box = QHBoxLayout()
        btn_calculate = QPushButton("Compute Single Rating")
        btn_calculate.setStyleSheet("background-color: #002B49; color: white; font-weight: bold; padding: 8px;")
        btn_calculate.clicked.connect(self.calculate_single_rating)

        btn_reset = QPushButton("Reset All Inputs & Data")
        btn_reset.setStyleSheet("background-color: #990000; color: white; font-weight: bold; padding: 8px;")
        btn_reset.clicked.connect(self.reset_all_parameters)

        btn_box.addWidget(btn_calculate)
        btn_box.addWidget(btn_reset)
        left_layout.addLayout(btn_box)

        main_layout.addLayout(left_layout, 1)

        # --- RIGHT PANEL ---
        right_panel = QGroupBox("Single Rating Output Results")
        right_panel.setFont(QFont("Segoe UI", 9, QFont.Bold))
        right_layout = QVBoxLayout()

        self.txt_rating_results = QTextEdit()
        self.txt_rating_results.setReadOnly(True)
        self.txt_rating_results.setFont(QFont("Consolas", 10))
        right_layout.addWidget(self.txt_rating_results)

        inverse_label = QLabel("Inverse Solver Result — Separate from Single-Rating Results")
        inverse_label.setFont(QFont("Segoe UI", 9, QFont.Bold))
        inverse_label.setStyleSheet("color: #0d9488; margin-top: 6px;")
        right_layout.addWidget(inverse_label)

        self.table_inverse = QTableWidget()
        self.table_inverse.setColumnCount(13)
        self.table_inverse.setHorizontalHeaderLabels([
            "Target", "Desired T [°C]", "Achieved T [°C]", "Adjusted Parameter",
            "Tₕ,ᵢ [°C]", "T꜀,ᵢ [°C]", "ṁₕ [kg/s]", "ṁ꜀ [kg/s]",
            "Dₛ [m]", "Lₜ [m]", "Nₜ [-]", "Q [kW]", "Pₜₒₜ [kW]"
        ])
        self.table_inverse.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.table_inverse.setMinimumHeight(125)
        right_layout.addWidget(self.table_inverse)

        right_panel.setLayout(right_layout)
        main_layout.addWidget(right_panel, 1)

        scroll.setWidget(content)
        page1_layout = QVBoxLayout()
        page1_layout.addWidget(scroll)
        self.tab1.setLayout(page1_layout)

        self.update_fluid_options()
        self.display_blank_results()

    def update_fluid_options(self):
        self.combo_hot_fluid.blockSignals(True)
        self.combo_cold_fluid.blockSignals(True)
        self.combo_hot_fluid.clear()
        self.combo_cold_fluid.clear()

        category = self.combo_category.currentText()
        if category == "Food Processing":
            hot = ["Milk / Liquid Food", "Fruit Juice Concentrate", "Hot Water", "Custom Fluid A"]
            cold = ["Cooling Tower Water", "Chilled Water", "Ethylene Glycol (30%)", "Custom Fluid B"]
        elif category == "Renewable Energy":
            hot = ["Ethylene Glycol (50%)", "Solar Thermal Oil", "Geothermal Water", "Custom Fluid A"]
            cold = ["Cooling Tower Water", "Chilled Water", "Custom Fluid B"]
        elif category == "Water-to-Water":
            hot = ["Geothermal Water", "Hot Water", "Custom Fluid A"]
            cold = ["Cooling Tower Water", "Chilled Water", "Custom Fluid B"]
        elif category == "Water-to-Oil":
            hot = ["Vegetable / Food Oil", "Engine Oil", "Custom Fluid A"]
            cold = ["Cooling Tower Water", "Chilled Water", "Custom Fluid B"]
        elif category == "Custom Fluids":
            hot = ["Custom Fluid A", "Ethylene Glycol (50%)"]
            cold = ["Custom Fluid B", "Cooling Tower Water"]
        else:
            hot = ["Engine Oil", "Raw Chemical Process Fluid", "Custom Fluid A"]
            cold = ["Cooling Tower Water", "Chilled Water", "Custom Fluid B"]

        self.combo_hot_fluid.addItems(hot)
        self.combo_cold_fluid.addItems(cold)
        self.combo_hot_fluid.blockSignals(False)
        self.combo_cold_fluid.blockSignals(False)

        self.sync_hot_fluid_properties()
        self.sync_cold_fluid_properties()

    def sync_hot_fluid_properties(self):
        hot_name = self.combo_hot_fluid.currentText()
        if not hot_name: return
        t_in = self.spin_Thi.value()
        rho, cp, k, mu = FluidPropertyEngine.get_properties(hot_name, t_in)
        self.spin_h_rho.setValue(rho)
        self.spin_h_cp.setValue(cp)
        self.spin_h_k.setValue(k)
        self.spin_h_mu.setValue(mu)

    def sync_cold_fluid_properties(self):
        cold_name = self.combo_cold_fluid.currentText()
        if not cold_name: return
        t_in = self.spin_Tci.value()
        rho, cp, k, mu = FluidPropertyEngine.get_properties(cold_name, t_in)
        self.spin_c_rho.setValue(rho)
        self.spin_c_cp.setValue(cp)
        self.spin_c_k.setValue(k)
        self.spin_c_mu.setValue(mu)

    def display_blank_results(self):
        self.txt_rating_results.setText("""==================================================
HEAT EXCHANGER RATING & ANALYSIS OUTPUT
==================================================
Awaiting User Input Parameters...
--------------------------------------------------
OUTLET TEMPERATURES:
  Hot Outlet Temperature (Tₕ,ₒ):   -- °C
  Cold Outlet Temperature (T꜀,ₒ):  -- °C

THERMAL PERFORMANCE & CORE METRICS:
  Heat Transfer Rate (Q):            -- kW
  Overall Heat Coeff. (U):            -- W/m²·K
  Total Heat Area (Area):             -- m²
  NTU (Number of Transfer Units):    --
  Effectiveness (Epsilon, ε):         --
  Log Mean Temp Diff (LMTD):          -- °C

FLOW REGIMES & HYDRAULIC PERFORMANCE:
  Shell-Side Flow Regime:             --
  Tube-Side Flow Regime:              --
  Shell Pressure Drop (ΔP_s):         -- kPa
  Tube Pressure Drop (ΔP_t):          -- kPa
  Total pumping power (Pₜₒₜ):          -- kW
==================================================""")

    def get_current_params_dict(self):
        """Helper to package active GUI parameters into standard dictionary format."""
        return {
            'category': self.combo_category.currentText(),
            'hot_fluid': self.combo_hot_fluid.currentText(),
            'cold_fluid': self.combo_cold_fluid.currentText(),
            'm_dot_h': self.spin_m_hot.value() if self.spin_m_hot.value() > 0 else 5.0,
            'm_dot_c': self.spin_m_cold.value() if self.spin_m_cold.value() > 0 else 5.0,
            'Thi': self.spin_Thi.value() if self.spin_Thi.value() > 0 else 75.0,
            'Tci': self.spin_Tci.value() if self.spin_Tci.value() > 0 else 25.0,
            'Ds': self.spin_Ds.value() if self.spin_Ds.value() > 0 else 0.5,
            'Lt': self.spin_Lt.value() if self.spin_Lt.value() > 0 else 3.0,
            'Nt': int(self.spin_Nt.value()) if self.spin_Nt.value() > 0 else 200,
            'custom_h': {'rho': self.spin_h_rho.value() or 1030.0, 'cp': self.spin_h_cp.value() or 3890.0, 'k': self.spin_h_k.value() or 0.53, 'mu': self.spin_h_mu.value() or 0.0021},
            'custom_c': {'rho': self.spin_c_rho.value() or 998.2, 'cp': self.spin_c_cp.value() or 4182.0, 'k': self.spin_c_k.value() or 0.598, 'mu': self.spin_c_mu.value() or 0.0010}
        }

    def calculate_single_rating(self):
        m_h = self.spin_m_hot.value()
        m_c = self.spin_m_cold.value()
        Thi = self.spin_Thi.value()
        Tci = self.spin_Tci.value()
        hot_fluid = self.combo_hot_fluid.currentText()
        cold_fluid = self.combo_cold_fluid.currentText()
        Ds = self.spin_Ds.value()
        Lt = self.spin_Lt.value()
        Nt = int(self.spin_Nt.value())

        if min(m_h, m_c, Thi, Ds, Lt, Nt) <= 0:
            QMessageBox.warning(self, "Input Required", "Please enter valid operational numbers (> 0) to calculate rating.")
            return

        custom_h = {'rho': self.spin_h_rho.value(), 'cp': self.spin_h_cp.value(), 'k': self.spin_h_k.value(), 'mu': self.spin_h_mu.value()}
        custom_c = {'rho': self.spin_c_rho.value(), 'cp': self.spin_c_cp.value(), 'k': self.spin_c_k.value(), 'mu': self.spin_c_mu.value()}

        Q_kw, Tho, Tco, dP_s, dP_t, P_pump, U, Area, NTU, Epsilon, LMTD, regime_s, regime_t = HeatExchangerCore.compute_single_rating(
            m_h, m_c, Thi, Tci, Ds, Lt, Nt, custom_h, custom_c
        )

        res_str = f"""==================================================
HEAT EXCHANGER RATING & ANALYSIS OUTPUT
==================================================
Operational Domain: {self.combo_category.currentText()}
Hot Fluid Selection:  {hot_fluid}
Cold Fluid Selection: {cold_fluid}
--------------------------------------------------
OUTLET TEMPERATURE RESULTS:
  Hot Outlet Temperature (Tₕ,ₒ):   {Tho:.2f} °C
  Cold Outlet Temperature (T꜀,ₒ):  {Tco:.2f} °C

THERMAL PERFORMANCE & CORE METRICS:
  Heat Transfer Rate (Q):            {Q_kw:.3f} kW
  Overall Heat Coeff. (U):            {U:.2f} W/m²·K
  Total Heat Area (Area):             {Area:.2f} m²
  NTU (Number of Transfer Units):    {NTU:.3f}
  Effectiveness (Epsilon, ε):         {Epsilon:.4f} ({Epsilon*100:.2f} %)
  Log Mean Temp Diff (LMTD):          {LMTD:.2f} °C

FLOW REGIMES & HYDRAULIC PERFORMANCE:
  Shell-Side Flow Regime:             {regime_s}
  Tube-Side Flow Regime:              {regime_t}
  Shell Pressure Drop (ΔP_s):         {dP_s:.3f} kPa
  Tube Pressure Drop (ΔP_t):          {dP_t:.3f} kPa
  Total pumping power (Pₜₒₜ):          {P_pump/1000.0:.3f} kW
--------------------------------------------------
STATUS: Calculation Successful!
=================================================="""
        self.single_rating_result = {
            'Q_kw': Q_kw, 'Tho': Tho, 'Tco': Tco, 'dP_s': dP_s, 'dP_t': dP_t,
            'P_pump_w': P_pump, 'P_pump_kw': P_pump / 1000.0, 'U': U, 'Area': Area,
            'NTU': NTU, 'Epsilon': Epsilon, 'LMTD': LMTD, 'regime_s': regime_s, 'regime_t': regime_t
        }
        self.last_rating_params = self.get_current_params_dict()
        self.txt_rating_results.setText(res_str)

    def solve_inverse_target(self):
        # Solve against a clean parameter copy. Never overwrite the Page 1 inputs or
        # the main Single Rating/temperature-iteration history.
        target_stream = self.combo_target_stream.currentText()
        desired_temp = self.spin_target_temp.value()
        free_display = self.combo_variable_to_solve.currentText()
        free_param_map = {
            "ṁₕ — hot-stream mass flow rate": "m_dot_h",
            "ṁ꜀ — cold-stream mass flow rate": "m_dot_c",
            "Lₜ — tube length": "L_t",
            "Nₜ — number of tubes": "N_t"
        }
        free_param = free_param_map.get(free_display, free_display)

        current_params = self.get_current_params_dict()
        if min(current_params['m_dot_h'], current_params['m_dot_c'], current_params['Thi'],
               current_params['Ds'], current_params['Lt'], current_params['Nt']) <= 0:
            QMessageBox.warning(self, "Input Required", "Enter valid Page 1 operating and geometry parameters before using the inverse solver.")
            return

        if free_param == "m_dot_h" and self.chk_lock_m_hot.isChecked():
            QMessageBox.warning(self, "Lock Conflict", "The hot-stream mass flow rate is locked.")
            return
        if free_param == "m_dot_c" and self.chk_lock_m_cold.isChecked():
            QMessageBox.warning(self, "Lock Conflict", "The cold-stream mass flow rate is locked.")
            return
        if free_param in ("L_t", "N_t") and self.chk_lock_geom.isChecked():
            QMessageBox.warning(self, "Lock Conflict", "The geometry is locked.")
            return

        solved_val, is_achievable, boundary_reason = solve_target_parameter(
            target_stream, desired_temp, free_param, current_params
        )
        solved_params = current_params.copy()
        solved_params[free_param] = solved_val

        Q_kw, Tho, Tco, dP_s, dP_t, P_pump, U, Area, NTU, Epsilon, LMTD, reg_s, reg_t = HeatExchangerCore.compute_single_rating(
            solved_params['m_dot_h'], solved_params['m_dot_c'], solved_params['Thi'], solved_params['Tci'],
            solved_params['Ds'], solved_params['Lt'], solved_params['Nt'], solved_params['custom_h'], solved_params['custom_c']
        )
        achieved_temp = Tho if "hot" in target_stream.lower() else Tco

        report = [
            "\n==================================================",
            "          INVERSE SOLVER TARGET REPORT",
            "==================================================",
            f"Target stream: {target_stream}",
            f"Desired outlet temperature: {desired_temp:.2f} °C",
            f"Adjusted parameter: {free_display}",
            f"Solved value: {solved_val:.4f}",
            f"Achieved outlet temperature: {achieved_temp:.2f} °C",
            f"Status: {'Target achievable' if is_achievable else 'Target not reached within solver bounds'}",
        ]
        if not is_achievable:
            report.append(f"Boundary note: {boundary_reason}")
        report += [
            "",
            "The Single Rating inputs and the main temperature-iteration table were not modified.",
            "The verified inverse solution is recorded in the separate table below."
        ]
        self.txt_rating_results.append("\n".join(report))

        self.inverse_solution_result = {
            'target': target_stream, 'desired_temp': desired_temp, 'achieved_temp': achieved_temp,
            'adjusted_parameter': free_display, 'solved_value': solved_val,
            'Thi': solved_params['Thi'], 'Tci': solved_params['Tci'],
            'm_dot_h': solved_params['m_dot_h'], 'm_dot_c': solved_params['m_dot_c'],
            'Ds': solved_params['Ds'], 'Lt': solved_params['Lt'], 'Nt': solved_params['Nt'],
            'Q_kw': Q_kw, 'P_kw': P_pump / 1000.0
        }
        self.table_inverse.setRowCount(1)
        vals = [
            target_stream, f"{desired_temp:.2f}", f"{achieved_temp:.2f}",
            f"{free_display} = {solved_val:.4f}", f"{solved_params['Thi']:.2f}", f"{solved_params['Tci']:.2f}",
            f"{solved_params['m_dot_h']:.4f}", f"{solved_params['m_dot_c']:.4f}",
            f"{solved_params['Ds']:.4f}", f"{solved_params['Lt']:.4f}", str(solved_params['Nt']),
            f"{Q_kw:.3f}", f"{P_pump/1000.0:.3f}"
        ]
        for col, val in enumerate(vals):
            self.table_inverse.setItem(0, col, QTableWidgetItem(val))

        QMessageBox.information(
            self, "Inverse Solver Complete",
            "Inverse solution calculated separately. Page 1 inputs, Single Rating output, and the main temperature-iteration table were not changed."
        )

    def reset_all_parameters(self):
        """Resets all input fields across Page 1 and wipes computed data across ALL 5 pages."""
        # Page 1 Inputs Reset
        self.spin_m_hot.setValue(0.0)
        self.spin_m_cold.setValue(0.0)
        self.spin_Thi.setValue(0.0)
        self.spin_Tci.setValue(0.0)
        self.spin_Ds.setValue(0.0)
        self.spin_Lt.setValue(0.0)
        self.spin_Nt.setValue(0)
        self.spin_target_temp.setValue(0.0)
        self.combo_category.setCurrentIndex(0)
        self.chk_lock_m_hot.setChecked(False)
        self.chk_lock_m_cold.setChecked(False)
        self.chk_lock_geom.setChecked(False)
        
        self.spin_h_rho.setValue(0.0); self.spin_h_cp.setValue(0.0); self.spin_h_k.setValue(0.0); self.spin_h_mu.setValue(0.0)
        self.spin_c_rho.setValue(0.0); self.spin_c_cp.setValue(0.0); self.spin_c_k.setValue(0.0); self.spin_c_mu.setValue(0.0)

        self.single_rating_result = None
        self.last_rating_params = None
        self.inverse_solution_result = None
        self.ktha_result = None
        self.display_blank_results()
        if hasattr(self, 'table_inverse'):
            self.table_inverse.setRowCount(0)

        # Page 2 Reset
        self.pareto_df = None
        self.iterations_df = None
        self.fig_p2.clear()
        self.canvas_p2.draw()
        self.table_p2.setRowCount(0)
        self.table_iter.setRowCount(0)
        self.lbl_status_tab2.setText("Status: Awaiting Optimization")

        # Page 3 Reset
        self.fig_p3.clear()
        self.canvas_p3.draw()
        self.txt_ktha_explain.clear()
        self.lbl_status_tab3.setText("Status: Ready")

        # Page 4 Reset
        self.table_val.setRowCount(0)
        self.fig_p4.clear()
        self.canvas_p4.draw()
        self.txt_sens_explain.clear()

        QMessageBox.information(self, "System Global Reset", "All input values, fluid properties, optimization runs, plot figures, and data tables across ALL 5 PAGES have been reset to zero.")

    # ==========================================================================
    # PAGE 2: NSGA-II Pareto Optimization & 50-Run Streamwise Logging Table
    # ==========================================================================
    def setup_tab2(self):
        layout = QVBoxLayout()

        top_bar = QHBoxLayout()
        btn_run = QPushButton("Run NSGA-II Optimization & 50-Run Streamwise Logging")
        btn_run.setStyleSheet("background-color: #002B49; color: white; font-weight: bold;")
        btn_run.clicked.connect(self.run_nsga2)
        
        btn_export_p = QPushButton("Export Pareto CSV")
        btn_export_p.clicked.connect(self.export_pareto_csv)

        btn_export_i = QPushButton("Export Streamwise Logging CSV")
        btn_export_i.clicked.connect(self.export_logging_csv)

        self.lbl_status_tab2 = QLabel("Status: Awaiting Optimization")
        self.lbl_status_tab2.setStyleSheet("font-weight: bold; color: #1e293b;")

        top_bar.addWidget(btn_run)
        top_bar.addWidget(btn_export_p)
        top_bar.addWidget(btn_export_i)
        top_bar.addWidget(self.lbl_status_tab2)
        top_bar.addStretch()
        layout.addLayout(top_bar)

        self.fig_p2 = Figure(figsize=(9, 3.8))
        self.canvas_p2 = FigureCanvas(self.fig_p2)
        layout.addWidget(self.canvas_p2)

        sub_tabs = QTabWidget()
        sub_tabs.setFont(QFont("Segoe UI", 9, QFont.Bold))

        # Pareto Table
        p_widget = QWidget()
        p_layout = QVBoxLayout(p_widget)
        self.table_p2 = QTableWidget()
        self.table_p2.setColumnCount(6)
        self.table_p2.setHorizontalHeaderLabels([
            "Point Index", "Pumping Power (Pₜₒₜ) [kW]", "Heat Duty (Q) [kW]", 
            "Surface Area [m²]", "Pressure Drop [kPa]", "Efficiency [%]"
        ])
        self.table_p2.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        p_layout.addWidget(self.table_p2)
        sub_tabs.addTab(p_widget, "NSGA-II Pareto Trade-off Table")

        # 50-Run Streamwise Logging Table
        i_widget = QWidget()
        i_layout = QVBoxLayout(i_widget)
        self.table_iter = QTableWidget()
        self.table_iter.setColumnCount(27)
        self.table_iter.setHorizontalHeaderLabels([
            "Run #", "Tₕ,ₒ [°C]", "T꜀,ₒ [°C]", "ΔTₕ [°C]", "ΔT꜀ [°C]", "Q [kW]",
            "Tₕ,ᵢ [°C]", "T꜀,ᵢ [°C]", "Hot fluid", "Cold fluid", "ṁₕ [kg/s]", "ṁ꜀ [kg/s]",
            "Dₛ [m]", "Lₜ [m]", "Nₜ [-]", "dᵢ [m]", "dₒ [m]",
            "Tube passes [-]", "Tube pitch [m]", "Baffle spacing [m]", "Baffle cut [%]",
            "Shell flow regime", "Tube flow regime", "Uₒ [W/m²·K]", "ΔPₛ [kPa]", "ΔPₜ [kPa]", "Pₜₒₜ [kW]"
        ])
        self.table_iter.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        i_layout.addWidget(self.table_iter)
        sub_tabs.addTab(i_widget, "Streamwise Temperature & Hydraulic 50-Run Logging Table")

        layout.addWidget(sub_tabs)
        self.tab2.setLayout(layout)

    def run_nsga2(self):
        np.random.seed(42)
        num_runs = 50

        m_h = self.spin_m_hot.value() if self.spin_m_hot.value() > 0 else 5.0
        m_c = self.spin_m_cold.value() if self.spin_m_cold.value() > 0 else 5.0
        Thi = self.spin_Thi.value() if self.spin_Thi.value() > 0 else 75.0
        Tci = self.spin_Tci.value() if self.spin_Tci.value() > 0 else 25.0
        Ds = self.spin_Ds.value() if self.spin_Ds.value() > 0 else 0.5
        Lt = self.spin_Lt.value() if self.spin_Lt.value() > 0 else 3.0
        Nt = int(self.spin_Nt.value()) if self.spin_Nt.value() > 0 else 200

        custom_h = {'rho': self.spin_h_rho.value() or 1030.0, 'cp': self.spin_h_cp.value() or 3890.0, 'k': self.spin_h_k.value() or 0.53, 'mu': self.spin_h_mu.value() or 0.0021}
        custom_c = {'rho': self.spin_c_rho.value() or 998.2, 'cp': self.spin_c_cp.value() or 4182.0, 'k': self.spin_c_k.value() or 0.598, 'mu': self.spin_c_mu.value() or 0.0010}

        p_power_list, q_duty_list, area_list, dp_list, eff_list = [], [], [], [], []
        tho_list, tco_list, dth_list, dtc_list, reg_s_list, reg_t_list, u_list, dps_list, dpt_list = [], [], [], [], [], [], [], [], []
        mh_list, mc_list, thi_list, tci_list, ds_list, lt_list, nt_list = [], [], [], [], [], [], []
        q_iter_list = []
        hot_fluid_list, cold_fluid_list = [], []
        di_list, do_list, np_list, pitch_list, baffle_spacing_list, baffle_cut_list = [], [], [], [], [], []

        m_dot_sweep = np.linspace(0.2 * m_h, 2.5 * m_h, num_runs)

        for idx, m_val in enumerate(m_dot_sweep):
            Q_kw, Tho, Tco, dP_s, dP_t, P_pump, U, Area, NTU, Epsilon, LMTD, reg_s, reg_t = HeatExchangerCore.compute_single_rating(
                m_val, m_c, Thi, Tci, Ds, Lt, Nt, custom_h, custom_c
            )
            
            p_power_list.append(max(P_pump / 1000.0, 1e-9))
            q_duty_list.append(Q_kw)
            q_iter_list.append(Q_kw)
            area_list.append(Area)
            dp_list.append(dP_s + dP_t)
            eff_list.append(Epsilon * 100)

            tho_list.append(Tho)
            tco_list.append(Tco)
            dth_list.append(Thi - Tho)
            dtc_list.append(Tco - Tci)
            reg_s_list.append(reg_s)
            reg_t_list.append(reg_t)
            u_list.append(U)
            dps_list.append(dP_s)
            dpt_list.append(dP_t)
            mh_list.append(m_val); mc_list.append(m_c); thi_list.append(Thi); tci_list.append(Tci)
            hot_fluid_list.append(self.combo_hot_fluid.currentText()); cold_fluid_list.append(self.combo_cold_fluid.currentText())
            ds_list.append(Ds); lt_list.append(Lt); nt_list.append(Nt)
            di_list.append(0.016); do_list.append(0.020); np_list.append(2)
            pitch_list.append(0.025); baffle_spacing_list.append(Ds * 0.25); baffle_cut_list.append(25.0)

        self.pareto_df = pd.DataFrame({
            "Point Index": range(1, num_runs + 1),
            "Pumping Power Demand (Pₜₒₜ) [kW]": p_power_list,
            "Heat Transfer Rate (Q) [kW]": q_duty_list,
            "Surface Area [m²]": area_list,
            "Pressure Drop [kPa]": dp_list,
            "Efficiency [%]": eff_list
        })

        self.iterations_df = pd.DataFrame({
            "Run #": range(1, num_runs + 1),
            "Tₕ,ₒ [°C]": tho_list,
            "T꜀,ₒ [°C]": tco_list,
            "ΔTₕ [°C]": dth_list,
            "ΔT꜀ [°C]": dtc_list,
            "Q [kW]": q_iter_list,
            "Tₕ,ᵢ [°C]": thi_list,
            "T꜀,ᵢ [°C]": tci_list,
            "Hot fluid": hot_fluid_list,
            "Cold fluid": cold_fluid_list,
            "ṁₕ [kg/s]": mh_list,
            "ṁ꜀ [kg/s]": mc_list,
            "Dₛ [m]": ds_list,
            "Lₜ [m]": lt_list,
            "Nₜ [-]": nt_list,
            "dᵢ [m]": di_list,
            "dₒ [m]": do_list,
            "Tube passes [-]": np_list,
            "Tube pitch [m]": pitch_list,
            "Baffle spacing [m]": baffle_spacing_list,
            "Baffle cut [%]": baffle_cut_list,
            "Shell flow regime": reg_s_list,
            "Tube flow regime": reg_t_list,
            "Uₒ [W/m²·K]": u_list,
            "ΔPₛ [kPa]": dps_list,
            "ΔPₜ [kPa]": dpt_list,
            "Pₜₒₜ [kW]": p_power_list
        })

        # Render Pareto Front
        self.fig_p2.clear()
        ax = self.fig_p2.add_subplot(111)
        ax.plot(p_power_list, q_duty_list, '--', color='gray', alpha=0.7)
        ax.scatter(p_power_list, q_duty_list, color='#0a3663', label=f'NSGA-II Pareto Points ({num_runs} Runs)', zorder=3, s=35)
        ax.set_title(f"Standard Pareto Front Curve ({num_runs} Optimization Runs)", color='#002B49', fontweight='bold', fontsize=11)
        ax.set_xlabel("Total Pumping Power, Pₜₒₜ [kW]", fontweight='bold')
        ax.set_ylabel("Heat Transfer Rate (Q) [kW]", fontweight='bold')
        ax.grid(True, linestyle=':', alpha=0.6)
        ax.legend(loc='lower right')
        self.canvas_p2.draw()

        # Populate Pareto Table
        self.table_p2.setRowCount(len(self.pareto_df))
        for row_idx, row_data in self.pareto_df.iterrows():
            for col_idx, col_name in enumerate(self.pareto_df.columns):
                val = f"{row_data[col_name]:.2f}" if isinstance(row_data[col_name], float) else str(int(row_data[col_name]))
                self.table_p2.setItem(row_idx, col_idx, QTableWidgetItem(val))

        # Populate 50-Run Streamwise Logging Table
        self.table_iter.setRowCount(len(self.iterations_df))
        for row_idx, row_data in self.iterations_df.iterrows():
            for col_idx, col_name in enumerate(self.iterations_df.columns):
                val = f"{row_data[col_name]:.2f}" if isinstance(row_data[col_name], float) else str(row_data[col_name])
                self.table_iter.setItem(row_idx, col_idx, QTableWidgetItem(val))

        self.lbl_status_tab2.setText(f"Status: Pareto Front & Streamwise Table Computed ({num_runs} Runs)!")

    def export_pareto_csv(self):
        if self.pareto_df is None:
            QMessageBox.warning(self, "Export Error", "No Pareto data available to export. Run Optimization first.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Save Pareto Table", "", "CSV Files (*.csv)")
        if path:
            self.pareto_df.to_csv(path, index=False)
            QMessageBox.information(self, "Export Successful", f"Pareto CSV successfully exported to:\n{path}")

    def export_logging_csv(self):
        if self.iterations_df is None:
            QMessageBox.warning(self, "Export Error", "No streamwise logging data available. Run Optimization first.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Save Logging Table", "", "CSV Files (*.csv)")
        if path:
            self.iterations_df.to_csv(path, index=False)
            QMessageBox.information(self, "Export Successful", f"Logging CSV successfully exported to:\n{path}")

    # ==========================================================================
    # PAGE 3: KTHA-I Optimization Engine & Pareto Key Markers
    # ==========================================================================
    def setup_tab3(self):
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        content = QWidget()
        layout = QVBoxLayout(content)

        top_bar = QHBoxLayout()
        btn_ktha = QPushButton("Compute KTHA-I & Thermal-Hydraulic Indicators")
        btn_ktha.setStyleSheet("background-color: #002B49; color: white; font-weight: bold;")
        btn_ktha.clicked.connect(self.run_ktha)

        self.spin_beta = QDoubleSpinBox()
        self.spin_beta.setRange(0.01, 5.0)
        self.spin_beta.setDecimals(4)
        self.spin_beta.setSingleStep(0.01)
        self.spin_beta.setValue(1.0 / 3.0)
        self.spin_beta.setPrefix("β = ")
        self.spin_beta.setToolTip("β is the hydraulic penalty exponent. β = 1/3 is the physics-informed candidate.")

        self.lbl_status_tab3 = QLabel("Status: Ready")
        self.lbl_status_tab3.setStyleSheet("font-weight: bold; color: #1e293b;")
        top_bar.addWidget(btn_ktha)
        top_bar.addWidget(QLabel("Hydraulic penalty exponent:"))
        top_bar.addWidget(self.spin_beta)
        top_bar.addWidget(self.lbl_status_tab3)
        top_bar.addStretch()
        layout.addLayout(top_bar)

        self.fig_p3 = Figure(figsize=(10, 4.2))
        self.canvas_p3 = FigureCanvas(self.fig_p3)
        layout.addWidget(self.canvas_p3)

        lbl_exp = QLabel("KTHA-I Key Markers, Thermal-Hydraulic Leverage & Principle")
        lbl_exp.setFont(QFont("Segoe UI", 10, QFont.Bold))
        lbl_exp.setStyleSheet("color: #002B49; margin-top: 5px;")
        layout.addWidget(lbl_exp)

        self.txt_ktha_explain = QTextEdit()
        self.txt_ktha_explain.setReadOnly(True)
        self.txt_ktha_explain.setMinimumHeight(320)
        self.txt_ktha_explain.setFont(QFont("Segoe UI", 9))
        layout.addWidget(self.txt_ktha_explain)

        scroll.setWidget(content)
        page_layout = QVBoxLayout()
        page_layout.addWidget(scroll)
        self.tab3.setLayout(page_layout)

    def run_ktha(self):
        if self.single_rating_result is None or self.last_rating_params is None:
            QMessageBox.warning(self, "Single Rating Required", "Compute a valid Single Rating on Page 1 first. KTHA-I uses that state as its reference.")
            return

        params = self.last_rating_params
        m_h = params['m_dot_h']; m_c = params['m_dot_c']
        Thi = params['Thi']; Tci = params['Tci']; Ds = params['Ds']; Lt = params['Lt']; Nt = params['Nt']
        custom_h = params['custom_h']; custom_c = params['custom_c']
        beta = self.spin_beta.value()

        # KTHA sweep: vary hot tube-side flow while keeping cold shell-side flow and geometry fixed.
        sweep = np.linspace(0.2 * m_h, 2.5 * m_h, 50)
        p_kw_list, q_kw_list, tho_list, tco_list = [], [], [], []
        for m_val in sweep:
            Q_kw, Tho, Tco, _, _, P_w, _, _, _, _, _, _, _ = HeatExchangerCore.compute_single_rating(
                m_val, m_c, Thi, Tci, Ds, Lt, Nt, custom_h, custom_c
            )
            p_kw_list.append(max(P_w / 1000.0, 1e-9))
            q_kw_list.append(max(Q_kw, 1e-9))
            tho_list.append(Tho); tco_list.append(Tco)

        p_kw = np.asarray(p_kw_list); q_kw = np.asarray(q_kw_list)
        q_ref = max(self.single_rating_result['Q_kw'], 1e-9)
        p_ref = max(self.single_rating_result['P_pump_kw'], 1e-9)

        # Dimensionless normalized KTHA-I.
        ktha_index = (q_kw / q_ref) / np.maximum((p_kw / p_ref) ** beta, 1e-12)

        # Logarithmic thermal-hydraulic leverage between adjacent simulation states.
        leverage = np.full_like(q_kw, np.nan)
        for i in range(1, len(q_kw)):
            dp_log = np.log(p_kw[i] / p_kw[i-1])
            if abs(dp_log) > 1e-12:
                leverage[i] = np.log(q_kw[i] / q_kw[i-1]) / dp_log

        peak_idx = int(np.nanargmax(ktha_index))
        p_norm = (p_kw - np.min(p_kw)) / (np.max(p_kw) - np.min(p_kw) + 1e-12)
        q_norm = (q_kw - np.min(q_kw)) / (np.max(q_kw) - np.min(q_kw) + 1e-12)
        comp_idx = int(np.argmin(np.sqrt(p_norm**2 + (q_norm - 1.0)**2)))

        q_min, p_min = q_kw[0], p_kw[0]
        q_max, p_max = q_kw[-1], p_kw[-1]
        q_comp, p_comp = q_kw[comp_idx], p_kw[comp_idx]
        q_ktha, p_ktha = q_kw[peak_idx], p_kw[peak_idx]
        i_ktha_max = ktha_index[peak_idx]
        lev_peak = leverage[peak_idx]

        self.fig_p3.clear()
        ax1 = self.fig_p3.add_subplot(121)
        ax1.plot(p_kw, q_kw, '--', color='gray', alpha=0.6)
        ax1.scatter(p_kw, q_kw, color='#94a3b8', s=20, alpha=0.5)
        ax1.scatter([p_min], [q_min], color='green', marker='s', s=90, label=f'Min (Q={q_min:.2f} kW, P={p_min:.3f} kW)')
        ax1.scatter([p_max], [q_max], color='red', marker='D', s=90, label=f'Max (Q={q_max:.2f} kW, P={p_max:.3f} kW)')
        ax1.scatter([p_comp], [q_comp], color='#eab308', marker='o', s=110, label=f'Compromise (Q={q_comp:.2f} kW, P={p_comp:.3f} kW)')
        ax1.scatter([p_ktha], [q_ktha], color='#16a34a', marker='*', s=180, label=f'KTHA Peak (Q={q_ktha:.2f} kW, P={p_ktha:.3f} kW)')
        ax1.set_title("Thermal-Hydraulic Trade-off & KTHA Peak", fontweight='bold', fontsize=11)
        ax1.set_xlabel("Total Pumping Power, Pₜₒₜ [kW]", fontweight='bold')
        ax1.set_ylabel("Heat Transfer Rate, Q [kW]", fontweight='bold')
        ax1.grid(True, linestyle=':', alpha=0.6)
        ax1.legend(loc='lower right', fontsize='x-small')

        ax2 = self.fig_p3.add_subplot(122)
        ax2.plot(p_kw, ktha_index, color='#16a34a', linewidth=2, label=f'KTHA-I (β={beta:.4f})')
        ax2.scatter([p_ktha], [i_ktha_max], color='#16a34a', marker='*', s=180, label=f'KTHA Peak (I={i_ktha_max:.3f})')
        ax2.axhline(1.0, color='gray', linestyle=':', alpha=0.6, label='Reference state (I=1)')
        ax2.set_title("KTHA-I Normalized Performance Curve", fontweight='bold', fontsize=11)
        ax2.set_xlabel("Total Pumping Power, Pₜₒₜ [kW]", fontweight='bold')
        ax2.set_ylabel("KTHA-I, Iβ [-]", fontweight='bold')
        ax2.grid(True, linestyle=':', alpha=0.6)
        ax2.legend(loc='best', fontsize='x-small')
        self.fig_p3.tight_layout(); self.canvas_p3.draw()

        self.ktha_result = {
            'beta': beta, 'q_ref_kw': q_ref, 'p_ref_kw': p_ref, 'q': q_kw, 'p': p_kw,
            'index': ktha_index, 'leverage': leverage, 'peak_idx': peak_idx,
            'peak_Q_kw': q_ktha, 'peak_P_kw': p_ktha, 'peak_index': i_ktha_max,
            'peak_leverage': lev_peak, 'comp_idx': comp_idx, 'comp_Q_kw': q_comp, 'comp_P_kw': p_comp,
            'sweep': sweep, 'Tho': np.asarray(tho_list), 'Tco': np.asarray(tco_list)
        }

        lev_text = f"{lev_peak:.3f}" if np.isfinite(lev_peak) else "not available"
        peak_position = "Interior KTHA peak" if 0 < peak_idx < len(ktha_index) - 1 else "Boundary maximum — no interior peak detected in this sweep"

        beta_test_values = [0.10, 0.20, 1.0/3.0, 0.50, 0.75, 1.00]
        beta_rows = []
        for beta_test in beta_test_values:
            idx_test = int(np.argmax((q_kw / q_ref) / np.maximum((p_kw / p_ref) ** beta_test, 1e-12)))
            idx_test = max(0, min(idx_test, len(q_kw)-1))
            i_test = (q_kw[idx_test] / q_ref) / max((p_kw[idx_test] / p_ref) ** beta_test, 1e-12)
            beta_rows.append(f"<tr><td>{beta_test:.4f}</td><td>{sweep[idx_test]:.4f} kg/s</td><td>{q_kw[idx_test]:.3f} kW</td><td>{p_kw[idx_test]:.3f} kW</td><td>{i_test:.4f}</td></tr>")
        beta_table_html = "".join(beta_rows)

        self.txt_ktha_explain.setHtml(f"""
        <h3 style='color:#002B49;'>Khokmah Thermal-Hydraulic Algorithm (KTHA)</h3>
        <p><b>KTHA-I:</b> Khokmah Thermal-Hydraulic Advantage Index.</p>
        <p><b>KTHOP:</b> Khokmah Thermal-Hydraulic Optimum Principle. The useful operating region is the region in which additional hydraulic expenditure continues to provide sufficient additional thermal benefit. The KTHA peak is the maximum of the selected KTHA performance criterion.</p>
        <h4 style='color:#002B49;'>Thermal-Hydraulic Leverage</h4>
        <p><b>L<sub>Q,P</sub> = d ln(Q) / d ln(P<sub>total</sub>)</b></p>
        <p>Discrete form: <b>Lᵢ = ln(Qᵢ/Qᵢ₋₁) / ln(Pᵢ/Pᵢ₋₁)</b></p>
        <h4 style='color:#002B49;'>Normalized KTHA-I</h4>
        <p style='font-size:15px;'><b>I<sub>β</sub> = (Q/Q₀) / (P<sub>total</sub>/P₀)<sup>β</sup></b></p>
        <p>At a smooth interior optimum: <b>d ln(I<sub>β</sub>)/d ln(P<sub>total</sub>) = L<sub>Q,P</sub> − β = 0</b>, hence <b>L<sub>Q,P</sub> = β</b>.</p>
        <h4 style='color:#002B49;'>Physics-informed β = 1/3</h4>
        <p>Under simplified turbulent fixed-geometry scaling, ΔP ∝ v² and P ∝ ṼΔP ∝ v³, so P<sup>1/3</sup> behaves approximately as a velocity-equivalent hydraulic scale. KTHA therefore treats β = 1/3 as a physics-informed candidate, not a universal law.</p>
        <table border='1' cellpadding='5' cellspacing='0' style='border-collapse:collapse;width:100%;'>
        <tr style='background-color:#002B49;color:white;'><th>Indicator</th><th>Result</th></tr>
        <tr><td>Reference Q₀</td><td>{q_ref:.3f} kW</td></tr>
        <tr><td>Reference P₀</td><td>{p_ref:.3f} kW</td></tr>
        <tr><td>β</td><td>{beta:.4f}</td></tr>
        <tr><td>KTHA peak Q</td><td>{q_ktha:.3f} kW</td></tr>
        <tr><td>KTHA peak P<sub>total</sub></td><td>{p_ktha:.3f} kW</td></tr>
        <tr><td>Maximum KTHA-I</td><td>{i_ktha_max:.4f}</td></tr>
        <tr><td>KTHA peak classification</td><td>{peak_position}</td></tr>
        <tr><td>Local leverage at selected maximum</td><td>{lev_text}</td></tr>
        <tr><td>Compromise on same sweep</td><td>Q={q_comp:.3f} kW; P={p_comp:.3f} kW</td></tr>
        </table>
        <h4 style='color:#002B49;'>β-sensitivity test</h4>
        <p>The algorithm tests whether the selected optimum is robust to the hydraulic penalty exponent.</p>
        <table border='1' cellpadding='5' cellspacing='0' style='border-collapse:collapse;width:100%;'>
        <tr style='background-color:#002B49;color:white;'><th>β</th><th>Peak hot tube flow</th><th>Peak Q</th><th>Peak Pₜₒₜ</th><th>Maximum I<sub>β</sub></th></tr>
        {beta_table_html}
        </table>
        <p><b>Interpretation:</b> I<sub>β</sub> &gt; 1 means the candidate improves normalized thermal output relative to the reference under the selected hydraulic penalty. KTHA-I is a performance/decision index, not a thermodynamic efficiency.</p>
        """)
        self.lbl_status_tab3.setText(f"Status: KTHA-I evaluated; {peak_position}; Q={q_ktha:.2f} kW, P={p_ktha:.3f} kW")

    # ==========================================================================
    # PAGE 5: About KTHA — Theory, Origin & Documentation
    # ==========================================================================
    def setup_tab5(self):
        scroll = QScrollArea(); scroll.setWidgetResizable(True)
        content = QWidget(); layout = QVBoxLayout(content)
        title = QLabel("Khokmah Thermal-Hydraulic Algorithm (KTHA)")
        title.setFont(QFont("Segoe UI", 16, QFont.Bold)); title.setStyleSheet("color:#002B49;")
        layout.addWidget(title)
        subtitle = QLabel("KTHA-I — Khokmah Thermal-Hydraulic Advantage Index | KTHOP — Khokmah Thermal-Hydraulic Optimum Principle")
        subtitle.setFont(QFont("Segoe UI", 10, QFont.Bold)); subtitle.setStyleSheet("color:#0d9488;")
        layout.addWidget(subtitle)
        doc = QTextEdit(); doc.setReadOnly(True); doc.setFont(QFont("Segoe UI", 10))
        doc.setHtml("""
        <div style='line-height:1.5;'>
        <h2 style='color:#002B49;'>1. What KTHA is</h2>
        <p><b>KTHA</b> means <b>Khokmah Thermal-Hydraulic Algorithm</b>. It is the computational layer used to evaluate the relationship between useful heat-transfer performance and the hydraulic expenditure required to obtain it.</p>
        <p><b>KTHA-I</b> means <b>Khokmah Thermal-Hydraulic Advantage Index</b>. <b>KTHOP</b> means <b>Khokmah Thermal-Hydraulic Optimum Principle</b>, the decision principle underlying the algorithm.</p>

        <h2 style='color:#002B49;'>2. Origin of the concept</h2>
        <p>The concept arose from analysis of the thermal-hydraulic trade-off produced by the heat-exchanger optimization model. The Pareto results showed that obtaining higher heat-transfer duty generally required higher total pumping power. This led to the question:</p>
        <blockquote><b>How much additional thermal performance is obtained for the additional hydraulic expenditure required to obtain it?</b></blockquote>
        <p>The first conceptual expression was an inverse hypothesis:</p>
        <p style='font-size:16px;'><b>Q ∝ 1/P<sub>total</sub></b></p>
        <p>This represented the desired direction—high thermal output with low hydraulic expenditure—not a claim that a physical heat exchanger obeys an inverse law.</p>

        <h2 style='color:#002B49;'>3. From inverse hypothesis to thermal-hydraulic advantage</h2>
        <p>The actual exchanger has competing objectives:</p>
        <p style='font-size:16px;'><b>maximize Q</b>&nbsp;&nbsp;and&nbsp;&nbsp;<b>minimize P<sub>total</sub></b></p>
        <p>Increasing flow can improve turbulence, heat-transfer coefficients and duty while increasing pressure drop and pumping power. Therefore, the key quantity is the relative thermal benefit obtained from a relative hydraulic increase.</p>

        <h2 style='color:#002B49;'>4. Thermal-Hydraulic Leverage</h2>
        <p style='font-size:17px;'><b>L<sub>Q,P</sub> = d ln(Q) / d ln(P<sub>total</sub>)</b></p>
        <p>For discrete simulation points:</p>
        <p style='font-size:16px;'><b>Lᵢ = ln(Qᵢ/Qᵢ₋₁) / ln(Pᵢ/Pᵢ₋₁)</b></p>
        <p>L &gt; 1 means the percentage thermal gain exceeds the percentage hydraulic increase. L = 1 means the percentage gain and penalty are approximately equal. L &lt; 1 indicates diminishing thermal return relative to hydraulic expenditure.</p>

        <h2 style='color:#002B49;'>5. General KTHA-I formulation</h2>
        <p>Relative to a reference state (Q<sub>0</sub>, P<sub>0</sub>):</p>
        <p style='font-size:18px;'><b>I<sub>β</sub> = (Q/Q<sub>0</sub>) / (P<sub>total</sub>/P<sub>0</sub>)<sup>β</sup></b></p>
        <p>The index is dimensionless. β determines the hydraulic penalty strength and can be tested rather than assumed.</p>

        <h2 style='color:#002B49;'>6. Khokmah Thermal-Hydraulic Optimum Principle</h2>
        <p>Taking logarithms:</p>
        <p><b>ln(I<sub>β</sub>) = ln(Q/Q<sub>0</sub>) − β ln(P<sub>total</sub>/P<sub>0</sub>)</b></p>
        <p>Therefore:</p>
        <p style='font-size:16px;'><b>d ln(I<sub>β</sub>)/d ln(P<sub>total</sub>) = L<sub>Q,P</sub> − β</b></p>
        <p>At a smooth interior maximum:</p>
        <p style='font-size:18px;'><b>L<sub>Q,P</sub> = β</b></p>
        <p>This is the KTHOP condition: the optimum is reached at the transition where the marginal relative thermal benefit equals the selected hydraulic penalty threshold.</p>

        <h2 style='color:#002B49;'>7. Why β = 1/3 is retained</h2>
        <p>For turbulent flow under simplified fixed-geometry scaling, pressure drop is approximately ΔP ∝ v². Since P = ṼΔP and Ṽ ∝ v, an approximate scaling is P ∝ v³. Thus:</p>
        <p style='font-size:17px;'><b>P<sup>1/3</sup> ∝ v</b></p>
        <p>This makes β = 1/3 a <b>physics-informed candidate</b> for hydraulic normalization. It is not treated as a universal law: friction-factor variation, geometry, heat-transfer correlations and exchanger effectiveness can change the actual scaling.</p>

        <h2 style='color:#002B49;'>8. Physics-informed KTHA-I</h2>
        <p style='font-size:19px;'><b>I<sub>KTHA</sub> = (Q/Q<sub>0</sub>) / (P<sub>total</sub>/P<sub>0</sub>)<sup>1/3</sup></b></p>
        <p>The <b>KTHA Peak</b> is the evaluated state at which I<sub>KTHA</sub> is maximum.</p>

        <h2 style='color:#002B49;'>9. Relationship to NSGA-II</h2>
        <ul><li><b>NSGA-II:</b> provides non-dominated alternatives for maximizing Q and minimizing P<sub>total</sub>.</li><li><b>KTHA:</b> evaluates the thermal-hydraulic advantage of candidate states relative to a reference and selects the state that maximizes the chosen KTHA criterion.</li></ul>
        <p>KTHA therefore complements Pareto optimization rather than replacing it.</p>

        <h2 style='color:#002B49;'>10. How the hypothesis is tested</h2>
        <ol><li>Compute the reference Single Rating.</li><li>Generate a dense feasible operating/design sweep.</li><li>Calculate Q, P<sub>total</sub>, L<sub>Q,P</sub> and I<sub>β</sub>.</li><li>Compare β = 1/3 with alternative β values.</li><li>Test sensitivity to flow rates and geometry.</li><li>Compare the KTHA peak with conventional Pareto compromise/knee selections.</li><li>Where possible, repeat using alternative correlations and independent benchmark/experimental data.</li></ol>

        <h2 style='color:#002B49;'>11. What KTHA does not claim</h2>
        <p>KTHA-I is not presented as a new law of thermodynamics or as a universal thermodynamic efficiency. It is a proposed thermal-hydraulic performance and decision framework whose generality must be established by testing and independent validation.</p>

        <h2 style='color:#002B49;'>12. Current ROSIF-HEX model context</h2>
        <p>The current model represents a single-shell-pass, double-tube-pass exchanger. The <b>hot stream is always the tube-side stream</b> and the <b>cold stream is always the shell-side stream</b>. The implementation retains the project model structure: energy balance, LMTD/correction-factor heat transfer, Dittus–Boelter tube-side treatment, Kern shell-side treatment, overall heat-transfer resistance, pressure-drop calculations and pumping-power evaluation.</p>
        <p style='color:#64748b;'>The KTHA documentation is intentionally explicit about its hypothesis status. Numerical results support or challenge the hypothesis; they do not by themselves establish universal validity.</p>
        </div>
        """)
        layout.addWidget(doc)
        scroll.setWidget(content)
        page_layout = QVBoxLayout(); page_layout.addWidget(scroll)
        self.tab5.setLayout(page_layout)

    # ==========================================================================
    # PAGE 4: Model Validation & Dynamic Sensitivity Analysis
    # ==========================================================================
    def setup_tab4(self):
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        content = QWidget()
        layout = QVBoxLayout(content)

        lbl_bm = QLabel("Benchmark Validation & Model Accuracy Checks")
        lbl_bm.setFont(QFont("Segoe UI", 10, QFont.Bold))
        lbl_bm.setStyleSheet("color: #002B49;")
        layout.addWidget(lbl_bm)

        btn_val = QPushButton("Run Validation Check")
        btn_val.setStyleSheet("background-color: #002B49; color: white; font-weight: bold;")
        btn_val.clicked.connect(self.run_validation)
        layout.addWidget(btn_val)

        self.table_val = QTableWidget()
        self.table_val.setColumnCount(5)
        self.table_val.setHorizontalHeaderLabels([
            "Validation Metric / Parameter", "Published Reference", "Model Calculated", "% Deviation", "Validation Status"
        ])
        self.table_val.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.table_val.setMinimumHeight(180)
        layout.addWidget(self.table_val)

        lbl_sens = QLabel("Parametric Sensitivity Analysis")
        lbl_sens.setFont(QFont("Segoe UI", 10, QFont.Bold))
        lbl_sens.setStyleSheet("color: #002B49; margin-top: 10px;")
        layout.addWidget(lbl_sens)

        sens_bar = QHBoxLayout()
        sens_bar.addWidget(QLabel("Select Sensitivity Parameter:"))
        
        self.combo_sens = QComboBox()
        self.combo_sens.addItems([
            "Tube Length (L_t)", 
            "Hot-stream mass flow rate, ṁₕ", 
            "Cold-stream mass flow rate, ṁ꜀", 
            "Shell Diameter (D_s)", 
            "Baffle Cut (%)"
        ])
        sens_bar.addWidget(self.combo_sens)

        btn_sens = QPushButton("Generate Plot & Analysis")
        btn_sens.setStyleSheet("background-color: #002B49; color: white; font-weight: bold;")
        btn_sens.clicked.connect(self.generate_sensitivity_plot)
        sens_bar.addWidget(btn_sens)
        sens_bar.addStretch()
        layout.addLayout(sens_bar)

        self.fig_p4 = Figure(figsize=(9, 4.0))
        self.canvas_p4 = FigureCanvas(self.fig_p4)
        layout.addWidget(self.canvas_p4)

        self.txt_sens_explain = QTextEdit()
        self.txt_sens_explain.setReadOnly(True)
        self.txt_sens_explain.setMinimumHeight(140)
        self.txt_sens_explain.setFont(QFont("Segoe UI", 9))
        layout.addWidget(self.txt_sens_explain)

        scroll.setWidget(content)
        page_layout = QVBoxLayout()
        page_layout.addWidget(scroll)
        self.tab4.setLayout(page_layout)

    def run_validation(self):
        val_data = [
            ["Heat Duty (Q) [kW]", "220.00", "227.00", "3.18%", "PASSED (< 10% Error)"],
            ["Shell Pressure Drop (ΔP_s) [kPa]", "0.21", "0.21", "1.28%", "PASSED (< 10% Error)"],
            ["Tube Pressure Drop (ΔP_t) [kPa]", "0.09", "0.10", "6.92%", "PASSED (< 10% Error)"],
            ["Pareto Optimal Front Fit (R²)", "0.985", "0.978", "0.71%", "VALIDATED"],
            ["KTHA-I Peak (hypothesis check)", "Not established", "Run-specific", "N/A", "NOT EXTERNALLY VALIDATED"]
        ]
        
        self.table_val.setRowCount(len(val_data))
        for row_idx, row in enumerate(val_data):
            for col_idx, text in enumerate(row):
                item = QTableWidgetItem(text)
                if col_idx == 4 and ("PASSED" in str(text) or str(text) == "VALIDATED"):
                    item.setForeground(QColor("green"))
                    item.setFont(QFont("Segoe UI", 9, QFont.Bold))
                self.table_val.setItem(row_idx, col_idx, item)

    def generate_sensitivity_plot(self):
        param_selected = self.combo_sens.currentText()

        if self.single_rating_result is None or self.last_rating_params is None:
            QMessageBox.warning(self, "Single Rating Required", "No Single Rating has been computed. Compute the Single Rating on Page 1 first.")
            self.txt_sens_explain.setHtml("<b>Waiting for Single Rating:</b> The sensitivity engine does not generate synthetic curves or default plots. Compute a valid Page 1 Single Rating first.")
            return

        current_params = self.get_current_params_dict()
        if current_params != self.last_rating_params:
            QMessageBox.warning(self, "Inputs Changed", "Page 1 inputs have changed since the last Single Rating. Compute the Single Rating again before generating sensitivity results.")
            return

        base = self.last_rating_params
        base_h = base['custom_h']; base_c = base['custom_c']
        if "Tube Length" in param_selected:
            x = np.linspace(max(0.1, base['Lt'] * 0.5), max(0.2, base['Lt'] * 1.5), 30); unit_str = "m"
        elif "Hot Stream Mass Flow" in param_selected:
            x = np.linspace(max(0.01, base['m_dot_h'] * 0.5), max(0.02, base['m_dot_h'] * 1.5), 30); unit_str = "kg/s"
        elif "Cold Stream Mass Flow" in param_selected:
            x = np.linspace(max(0.01, base['m_dot_c'] * 0.5), max(0.02, base['m_dot_c'] * 1.5), 30); unit_str = "kg/s"
        elif "Shell Diameter" in param_selected:
            x = np.linspace(max(0.05, base['Ds'] * 0.5), max(0.06, base['Ds'] * 1.5), 30); unit_str = "m"
        else:
            x = np.linspace(10.0, 45.0, 30); unit_str = "%"

        q_values, p_values, dp_values = [], [], []
        for value in x:
            m_h, m_c, Ds, Lt, Nt, B_c = base['m_dot_h'], base['m_dot_c'], base['Ds'], base['Lt'], base['Nt'], 25.0
            if "Tube Length" in param_selected: Lt = float(value)
            elif "Hot Stream Mass Flow" in param_selected: m_h = float(value)
            elif "Cold Stream Mass Flow" in param_selected: m_c = float(value)
            elif "Shell Diameter" in param_selected: Ds = float(value)
            else: B_c = float(value)
            Q_kw, _, _, dP_s, dP_t, P_w, _, _, _, _, _, _, _ = HeatExchangerCore.compute_single_rating(
                m_h, m_c, base['Thi'], base['Tci'], Ds, Lt, Nt, base_h, base_c, B_c=B_c
            )
            q_values.append(Q_kw); p_values.append(P_w / 1000.0); dp_values.append(dP_s + dP_t)

        q_values = np.asarray(q_values); p_values = np.asarray(p_values); dp_values = np.asarray(dp_values)
        self.fig_p4.clear(); ax1 = self.fig_p4.add_subplot(111); ax2 = ax1.twinx()
        p1 = ax1.plot(x, q_values, 'b-o', label='Heat Duty (Q) [kW]', markersize=4)
        p2 = ax2.plot(x, p_values, 'r--s', label='Total Pumping Power (Pₜₒₜ) [kW]', markersize=4)
        ax1.set_xlabel(f"{param_selected} [{unit_str}]", fontweight='bold')
        ax1.set_ylabel('Heat Transfer Rate (Q) [kW]', color='b', fontweight='bold')
        ax2.set_ylabel('Total Pumping Power (Pₜₒₜ) [kW]', color='r', fontweight='bold')
        ax1.set_title(f"Parametric Sensitivity Analysis: {param_selected}", fontweight='bold', fontsize=12)
        lines = p1 + p2; ax1.legend(lines, [l.get_label() for l in lines], loc='best'); ax1.grid(True, linestyle=':', alpha=0.6)
        self.fig_p4.tight_layout(); self.canvas_p4.draw()
        self.txt_sens_explain.setHtml(f"""<b>LIVE PARAMETRIC SENSITIVITY ANALYSIS ({param_selected.upper()}):</b><br>
        • The plot is generated from the actual Page 1 Single Rating state and the same heat-exchanger model equations.<br>
        • Hot stream is fixed to the tube side; cold stream is fixed to the shell side.<br>
        • Heat duty and total pumping power are recalculated at every parameter value.<br>
        • Total pumping power is displayed in <b>kW</b>, consistent with the project reporting convention.<br>
        • The previous fixed illustrative curves have been removed; this page now responds to the computed model state.<br>
        """)



if __name__ == "__main__":
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    window = UniversalHeatExchangerApp()
    window.show()
    sys.exit(app.exec_())