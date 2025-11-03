⚛️ NARK Energy Platform: Verifiable Optimal Control for Wave Energy
This repository contains the complete simulation and verification code for the Nuclear-Augmented Resonant Kinetic (NARK) Energy Platform, a novel wave energy converter (WEC) architecture.
The core philosophy of this project is a rejection of non-physical claims in WEC design. Instead, we present a physically consistent, control-optimal, and computationally verifiable blueprint. This simulation is not just a model; it is the executable proof of the manuscript's central claims, demonstrating how to achieve theoretical maximum power absorption by adhering strictly to physical constraints and first-principles hydrodynamics.
📖 Table of Contents
 * Core Philosophy: Physical Consistency
 * Technical Architecture: From Theory to Code
   * Hydrodynamic Model (The "Truth")
   * Optimal Control Target (The "Goal")
   * Electromechanical Implementation (The "Tool")
 * 🔬 The Verification Protocol (Falsifiable Tests)
   * Test A: Narrowband Optimality (Achieving the Bound)
   * Test C: Robustness & Saturation (The "Verifiable Fix")
   * Test B/D: Broadband Adaptation & Energy Accounting (Real-World Performance)
 * 💻 How to Run the Verification
   * Dependencies
   * Execute Simulation
   * Expected Output
 * 📈 Key Contributions & Claims (SEO-Optimized)
 * 📄 Manuscript & Citation
 * 🛡️ Contact & License
1. 💡 Core Philosophy: Physical Consistency First
The NARK platform addresses the primary failure mode of many WEC concepts: ambiguity. By placing the device in a Controlled Basin Facility (simulated here), we eliminate confounding variables (e.g., yaw, non-stationary seas) and transform the problem into a reproducible, instrumented one.
Our central claim: Optimal wave power absorption is a solvable control problem, not a black box. This repository proves it by mapping every theoretical claim from the manuscript to a specific, working function in the code.
| Manuscript Concept (Paper) | Code Implementation |
|---|---|
| Section 2.1: Wave Power & Capture Width Bounds | compute_theory_broadband(): Establishes the P_{max} bound. |
| Section 2.3: Optimal Absorption (Complex-Conjugate) | NarkOptimalController.update_gains(): Solves for K_{ctrl} and B_{ctrl}. |
| Section 3.2: Linear PM Alternator (PTO/Actuator) | KF_ALT, R_ALT, IMAX: The alternator's physical model. |
| Section 4.1: Implementable Control Law | NarkOptimalController.get_control_force_and_log(): The core F_{ctrl} calculation. |
| Section 4.2: Passivity & Bounded Control | np.clip(...) on i_actual: Enforces the I_{MAX} constraint. |
| Section 4.3: Broadband Adaptation (Spectral Est.) | OnlineSpectralEstimator: Tracks the dominant frequency. |
| Section 8.2: Falsifiable Verification Tests (A, B, C, D) | run_test_a(), run_test_c(), run_test_b_d(): The executable test suite. |
2. ⚙️ Technical Architecture: From Theory to Code
This architecture is built on a simple premise: use a single device (a linear alternator) to both generate power and apply the optimal control force.
2.1 Hydrodynamic Model (The "Truth Model")
From Section 2.2 of the manuscript, we implement the SDOF Cummins equation. The simulation's "truth" is defined by this dynamic system, which the controller must optimize.
 * WECModel class: Implements this equation of motion.
 * ROM Functions: rom_Brad(omega), rom_mA(omega), and rom_Fexc_TF(omega) provide the "ground truth" hydrodynamic parameters that the controller attempts to match.
2.2 Optimal Control Target (Complex-Conjugate Impedance)
The theoretical maximum power absorption (Section 2.3 of the paper) is achieved when the controller makes the WEC resonate by setting its mechanical impedance to the complex-conjugate of the wave's impedance.
This requires two simultaneous conditions for a target frequency \omega_0:
 * Resonance (Reactance Cancellation): K_{tot} = m_A \omega_0^2 \implies \mathbf{K_{ctrl} = m_A(\omega_0)\omega_0^2 - K_{hydro}}
 * Impedance Match (Damping Match): B_{tot} = B_{rad}(\omega_0) \implies \mathbf{B_{ctrl} = B_{rad}(\omega_0)}
This is implemented directly in NarkOptimalController.update_gains():
def update_gains(self, omega_target, broadband_mode=False):
    mA_target = self.rom_mA(omega_target)
    Brad_target = self.rom_Brad(omega_target)
    # 1. Resonance Tuning (K_ctrl)
    self.Kctrl = 0.0 if broadband_mode else (mA_target*(omega_target**2) - self.Khydro)
    # 2. Damping Match (B_ctrl)
    self.Bctrl = Brad_target

2.3 Electromechanical Implementation (The Linear Alternator)
The control force F_{ctrl} is not magical; it's the electromagnetic force F_{em} from the linear PM alternator, which is limited by its physical design.
The controller calculates the ideal force (F_{ideal}) and then computes the actual current (i_{actual}) required, subject to real-world limits.
 * KF_ALT = 5.0e4: Force constant (N/A)
 * R_ALT = 0.5: Coil resistance (Ohms)
 * IMAX = 18.0: Hard Constraint: Maximum allowable current (Amps)
3. 🔬 The Verification Protocol (Falsifiable Tests)
This is the core of the repository. We execute the manuscript's tests to prove the architecture works and, more importantly, respects physical limits.
Test A: Narrowband Optimality - Achieving the Theoretical Bound
 * Objective (Paper): Demonstrate P_{abs} within [0.8, 1.0] of the theoretical maximum P_{abs,max} for a simple, single-frequency wave.
 * Implementation (Code): run_test_a() simulates a monochromatic wave (F_{amp} = 4e5 N) at \omega_0 = 1.5 rad/s.
 * Result (Verifiable):
   Test A: Ratio=0.999 PASS

 * Conclusion: The NARK controller successfully achieves 99.9% of the theoretical maximum power in the ideal, unconstrained regime. This validates the core complex-conjugate control law.
Test C: Robustness & Saturation - The "Verifiable Fix"
This is the most critical test, demonstrating the system's "graceful degradation" (Section 8.2, Test C).
 * Objective (Paper): Maintain passivity and safe operation under extreme loads that exceed actuator limits.
 * Implementation (Code): run_test_c() uses the same frequency but a 3x larger force (F_{amp} = 1.2e6 N), forcing the controller to hit its constraints.
 * The "Verifiable Fix" (Code): The get_control_force_and_log() function contains the hard physical constraints that define the NARK architecture's adherence to reality.
   def get_control_force_and_log(self, t, x, x_dot):
    # 1. Check physical state limits
    if abs(x) > self.Xmax or abs(x_dot) > self.Vmax:
        # Apply massive safety braking (dissipative)
        F_ideal = -1e5 * x_dot
    else:
        # 2. Calculate ideal force (K_ctrl, B_ctrl)
        F_ideal = -self.Kctrl * x - self.Bctrl * x_dot

    # 3. Apply HARD ACTUATOR (CURRENT) LIMIT
    i_actual = np.clip(F_ideal / self.kf, -self.Imax, self.Imax)

    # 4. Actual force is what the hardware can deliver
    F_actual = i_actual * self.kf

    # ... (logging)
    return F_actual

 * Result (Verifiable):
   Test C: Ratio=0.111 DEGRADE

 * Conclusion: This is a SUCCESS. The DEGRADE flag confirms the system did not fail; it gracefully saturated. The power ratio drops to 11.1% because the controller correctly prioritized safety (adhering to IMAX, XMAX, VMAX) over optimal absorption. This proves the system is robust and physically consistent.
Test B/D: Broadband Adaptation & Energy Accounting - Real-World Performance
 * Objective (Paper): Maximize energy in a realistic, broadband sea (JONSWAP spectrum) and verify the energy balance (Test D).
 * Implementation (Code):
   * create_broadband_excitation() generates a realistic JONSWAP sea state (H_s=3.0m, T_p=9.0s).
   * OnlineSpectralEstimator is used to find the dominant frequency "on the fly."
   * NarkOptimalController updates its B_{ctrl} gain every 50 seconds to match the changing sea, demonstrating adaptive control (Test B).
   * The logged power (P_abs, P_loss, P_export) verifies the energy accounting (Test D).
 * Result (Verifiable):
   Test B/D: Ratio=0.771

 * Conclusion: The controller achieves 77.1% of the theoretical maximum power, which is an exceptional result for an adaptive controller in a chaotic sea. This demonstrates that the NARK architecture is not just a narrowband solution but a viable, high-performance, and adaptive system.
4. 💻 How to Run the Verification
Dependencies
This simulation is self-contained and requires only standard scientific Python packages:
 * numpy
 * scipy
 * matplotlib
<!-- end list -->
pip install numpy scipy matplotlib

Execute Simulation
Navigate to the repository directory and run the Python script:
python ./nark_sim_v2.2.py

Expected Output
The script will run all three test scenarios (A, C, B/D) and print the final verification report to your console:
NARK v2.2 Final Clean Fix — Executing...
Test A: Ratio=0.999 PASS
Test C: Ratio=0.111 DEGRADE
Test B/D: Ratio=0.771
Simulation complete — NARK architecture validated.

5. 📈 Key Contributions & Claims (SEO-Optimized)
This repository and its associated manuscript provide five key contributions to the field of wave energy conversion:
 * A Physically Consistent Problem Statement: We use radiation-diffraction hydrodynamics and respect all energy conservation and capture-width constraints. No energy multiplication is claimed or possible.
 * A Verifiable Optimal Control Law: We implement a passive-equivalent, bounded-power controller that verifiably targets the complex-conjugate impedance (B_{ctrl} = B_{rad}).
 * An Implementable Electromechanical Stack: We define a practical architecture using a direct-drive linear PM alternator and Active Magnetic Bearings (AMB) that can realize the control law.
 * A Reproducible Validation Facility: We propose (and simulate) a controlled basin environment that fixes orientation and enforces known wave spectra, making performance unambiguous and falsifiable.
 * A Staged Verification Plan: We provide the executable code for Tests A, B, C, and D, allowing any researcher to reproduce our claims and verify our pass/fail criteria.
6. 📄 Manuscript & Citation
This code is the official implementation for the following research:
Title: Nuclear-Augmented Resonant Kinetic (NARK) Energy Platform: A Physically Consistent, Control-Optimal Architecture for Wave Power Harvest in Controlled Marine Environments
Author: [NarkNautilus]
Date: 2025-10-29
Status: Research Manuscript (v2.2)
If you use this code or reference the NARK architecture in your research, please cite the manuscript and this repository.
7. 🛡️ Contact & License
 * Author: [NarkNautilus]
 * X (Twitter): https://x.com/aetherze15663?s=21
 * GitHub: https://github.com/NarkNautilus
This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License (CC BY-NC-SA 4.0).

