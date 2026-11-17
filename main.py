# `sys` reads command-line arguments; NumPy handles the complex matrices and
# vectors used by the power-flow equations; `time` measures elapsed runtime.
import sys
import numpy as np
import time

def parse_fivebus_pq(filepath):
    # Keep system settings, bus records, and line records in separate containers.
    system_data, buses, lines = {}, [], []
    # The input file is a whitespace-delimited text format; each record begins
    # with a keyword such as SYSTEM, BUS, or LINE.
    with open(filepath, 'r') as f:
        for line in f:
            # Lines beginning with '%' are documentation/comments; blank lines
            # contain no record to parse.
            if line.startswith('%') or not line.strip():
                continue
            # Split the record into fields; numeric fields are converted below.
            parts = line.split()
            if parts[0] == 'SYSTEM':
                # The base MVA converts specified MW/MVAr to per-unit power;
                # tolerance is the convergence threshold for power mismatch.
                system_data = {'name': parts[1], 'baseMVA': float(parts[2]), 'tolerance': float(parts[3])}
            elif parts[0] == 'BUS':
                # Bus type determines which state variables and power equations
                # are unknown: SL (slack), PV (specified P and voltage), or
                # PQ (specified P and Q).
                buses.append({
                    'name': parts[1], 'type': parts[2], 'volts': float(parts[3]),
                    'Pgen': float(parts[4]), 'Qgen': float(parts[5]),
                    'Pload': float(parts[6]), 'Qload': float(parts[7]), 'Qcap': float(parts[8])
                })
            elif parts[0] == 'LINE':
                # Series resistance/reactance define the line impedance;
                # shunt conductance/susceptance define its shunt admittance.
                # Rating is retained from the input but is not used by this
                # power-flow calculation.
                lines.append({
                    'from': parts[1], 'to': parts[2],
                    'Rse': float(parts[3]), 'Xse': float(parts[4]),
                    'Gsh': float(parts[5]), 'Bsh': float(parts[6]), 'Rating': float(parts[7])
                })

    # Return the parsed records for network construction and solving.
    return system_data, buses, lines

def build_ybus(buses, lines):
    # The bus-admittance matrix Ybus relates bus voltages to injected currents:
    # I = Ybus * V. There is one row and column per bus.
    n = len(buses)
    # Map bus names used in the file to zero-based matrix/vector indices.
    bus_idx = {bus['name']: i for i, bus in enumerate(buses)}
    # Complex entries represent conductance (real part) and susceptance
    # (imaginary part); start with no connections or shunts.
    Ybus = np.zeros((n, n), dtype=complex)
    
    for line in lines:
        # Find the matrix positions for this line's two endpoint buses.
        i, j = bus_idx[line['from']], bus_idx[line['to']]
        # Series impedance is Z = R + jX, so its admittance is the reciprocal
        # Y = 1/Z. The imaginary unit `1j` represents j = sqrt(-1).
        Zse = line['Rse'] + 1j * line['Xse']
        Yse = 1.0 / Zse
        # Combine shunt conductance and susceptance into complex admittance.
        Ysh = line['Gsh'] + 1j * line['Bsh']
        
        # Each line contributes its series and shunt admittance to both
        # endpoint self-admittances (diagonal terms).
        Ybus[i, i] += Yse + Ysh
        Ybus[j, j] += Yse + Ysh
        # The current from one endpoint to the other depends on the voltage
        # difference; this produces negative series-admittance off-diagonals.
        Ybus[i, j] -= Yse
        Ybus[j, i] -= Yse
        
    # Return the assembled nodal admittance matrix.
    return Ybus

def newton_raphson_power_flow(filepath):
    # Start timing before parsing so the reported runtime includes file reading,
    # network construction, and the iterative solution.
    start_time = time.time()
    # Read system/bus/line data, then form the electrical network model.
    system_data, buses, lines = parse_fivebus_pq(filepath)
    Ybus = build_ybus(buses, lines)

    print(f"Ybus Matrix:\n{Ybus} \n \n")

    # Newton-Raphson state variables are voltage magnitudes V and phase angles
    # delta. The input magnitudes provide the initial guess; all angles start
    # at zero radians, making the slack-bus angle the reference.
    V = np.array([bus['volts'] for bus in buses], dtype=float)
    delta = np.zeros(len(buses), dtype=float)

    print(f"Initial Voltage Magnitudes: {V}")
    print(f"Initial Angles: {delta}")
    
    # Power specifications are converted to per unit by dividing by base MVA.
    base_mva = system_data['baseMVA']
    # Stop when every included active/reactive power equation has mismatch
    # magnitude below this tolerance.
    tolerance = system_data['tolerance']
    # Bound the number of Newton steps to avoid an unending solve if it fails
    # to converge.
    max_iter = 500
    # Count only iterations where a Newton correction is applied.
    total_iter = 0

    print(f"Starting Newton-Raphson Power Flow Analysis with base MVA = {base_mva}, tolerance = {tolerance}, max iterations = {max_iter} \n \n ")
    
    for iteration in range(max_iter):
        # Re-evaluate the network equations at the current voltage estimate.
        n = len(buses)
        P_calc = np.zeros(n)
        Q_calc = np.zeros(n)

        print(f"Iteration values: V = {V}, delta = {delta}, P_calc = {P_calc}, Q_calc = {Q_calc}, n = {n}")
        
        # 1. Calculate each bus's injected complex power from S = V * conj(I),
        # where I = Ybus * V. Expanding by bus pair gives:
        # P_i = sum_j Vi*Vj*|Yij|*cos(theta_ij - delta_i + delta_j)
        # Q_i = -sum_j Vi*Vj*|Yij|*sin(theta_ij - delta_i + delta_j).
        # `np.angle` gives theta_ij and the trig functions use radians.
        for i in range(n):
            for j in range(n):
                P_calc[i] += V[i] * V[j] * np.abs(Ybus[i, j]) * np.cos(np.angle(Ybus[i, j]) - delta[i] + delta[j])
                Q_calc[i] -= V[i] * V[j] * np.abs(Ybus[i, j]) * np.sin(np.angle(Ybus[i, j]) - delta[i] + delta[j])

        print(f"Calculated Power: P_calc = {P_calc}, Q_calc = {Q_calc}")
        
        # Mismatch arrays hold specified minus calculated power. Entries for
        # buses whose corresponding equation is not solved remain zero.
        dP = np.zeros(n)
        dQ = np.zeros(n)
        # Every non-slack bus has an active-power equation and unknown angle.
        pq_pv_indices = [] # Indices for delta updates (PQ and PV buses)
        # Only PQ buses have a reactive-power equation and unknown magnitude.
        pq_indices = []    # Indices for voltage updates (PQ buses only)
        
        # 2. Form the specified net injections: generation minus demand.
        # Slack-bus power is not specified; it balances the network losses and
        # is therefore omitted from the mismatch equations.
        for i, bus in enumerate(buses):
            if bus['type'] == 'SL':
                continue
            
            # Net active injection in per unit: (Pgen - Pload) / base MVA.
            P_spec = (bus['Pgen'] - bus['Pload']) / base_mva
            # Newton's right-hand side is the required correction to calculated
            # power: specified injection minus the current calculated injection.
            dP[i] = P_spec - P_calc[i]
            pq_pv_indices.append(i)
            
            if bus['type'] == 'PQ':
                # Capacitor MVAr is modeled as positive reactive injection, so
                # add Qcap to generation before subtracting reactive load.
                Q_spec = (bus['Qgen'] + bus['Qcap'] - bus['Qload']) / base_mva
                dQ[i] = Q_spec - Q_calc[i]
                pq_indices.append(i)
                
        # The nonlinear equations are ordered as all non-slack active-power
        # mismatches followed by all PQ reactive-power mismatches. Slack P/Q
        # and PV reactive power are excluded because they are not specified
        # equations in this formulation.
        mismatch_vector = np.concatenate([dP[pq_pv_indices], dQ[pq_indices]])
        print(f"Mismatch vector: {mismatch_vector}")
        # Infinity norm (largest absolute component) checks the worst equation;
        # convergence requires every mismatch to be smaller than tolerance.
        mismatch_max = np.max(np.abs(mismatch_vector))
        print(f"\n \nIteration {iteration + 1}: Max mismatch = {mismatch_max:.6f}")
        if mismatch_max < tolerance:
            print(f"\n \nConverged in {iteration} iterations.")
            break
            
        # 3. Build the Jacobian J = [dP/d(delta)  dP/dV;
        #                            dQ/d(delta)  dQ/dV].
        # Its columns follow the unknown ordering [delta for PV/PQ buses,
        # V for PQ buses], matching the mismatch vector's equation ordering.
        num_angles = len(pq_pv_indices)
        num_volts = len(pq_indices)
        # There is one equation per unknown: one P equation per non-slack bus,
        # and one Q equation per PQ bus.
        J = np.zeros((num_angles + num_volts, num_angles + num_volts))
        
        for idx_i, i in enumerate(pq_pv_indices):
            for idx_j, j in enumerate(pq_pv_indices):
                # J1 = dP/d(delta). Off-diagonal derivative:
                # dP_i/d(delta_j) = -Vi*Vj*|Yij|*sin(theta_ij-delta_i+delta_j).
                if i != j:
                    J[idx_i, idx_j] = -1 * V[i] * V[j] * np.abs(Ybus[i, j]) * np.sin(np.angle(Ybus[i, j]) - delta[i] + delta[j])
                else:
                    # The diagonal derivative collects all terms involving
                    # delta_i: dP_i/d(delta_i) = -Q_i - Vi^2*imag(Yii).
                    J[idx_i, idx_j] = -Q_calc[i] - (V[i]**2) * np.imag(Ybus[i, i])
                    
            for idx_j, j in enumerate(pq_indices):
                # J2 = dP/dV. Columns exist only for PQ-bus magnitudes because
                # slack and PV magnitudes are held fixed in this formulation.
                col_j = num_angles + idx_j
                if i != j:
                    # Off-diagonal: dP_i/dV_j = Vi*|Yij|*cos(angle difference).
                    J[idx_i, col_j] = V[i] * np.abs(Ybus[i, j]) * np.cos(np.angle(Ybus[i, j]) - delta[i] + delta[j])
                else:
                    # Diagonal: dP_i/dV_i = P_i/Vi + Vi*real(Yii).
                    J[idx_i, col_j] = (P_calc[i] / V[i]) + V[i] * np.real(Ybus[i, i])
                    
        for idx_i, i in enumerate(pq_indices):
            # Reactive-power equations exist only for PQ buses; their Jacobian
            # rows follow the active-power rows.
            row_i = num_angles + idx_i
            for idx_j, j in enumerate(pq_pv_indices):
                # J3 = dQ/d(delta). Off-diagonal:
                # dQ_i/d(delta_j) = -Vi*Vj*|Yij|*cos(angle difference).
                if i != j:
                    J[row_i, idx_j] = -V[i] * V[j] * np.abs(Ybus[i, j]) * np.cos(np.angle(Ybus[i, j]) - delta[i] + delta[j])
                else:
                    # Diagonal: dQ_i/d(delta_i) = P_i - Vi^2*real(Yii).
                    J[row_i, idx_j] = P_calc[i] - (V[i]**2) * np.real(Ybus[i, i])
                    
            for idx_j, j in enumerate(pq_indices):
                # J4 = dQ/dV. These columns are also limited to PQ buses.
                col_j = num_angles + idx_j
                if i != j:
                    # Off-diagonal: dQ_i/dV_j = -Vi*|Yij|*sin(angle difference).
                    J[row_i, col_j] = -V[i] * np.abs(Ybus[i, j]) * np.sin(np.angle(Ybus[i, j]) - delta[i] + delta[j])
                else:
                    # Diagonal: dQ_i/dV_i = Q_i/Vi - Vi*imag(Yii).
                    J[row_i, col_j] = (Q_calc[i] / V[i]) - V[i] * np.imag(Ybus[i, i])

        print(f"Jacobian Matrix:\n{J}")

        # 4. Newton's linearization is J * correction = mismatch. Solve the
        # system directly (rather than explicitly forming inv(J)) for better
        # numerical stability and efficiency. The correction order is
        # [delta updates, then PQ voltage-magnitude updates].
        updates = np.linalg.solve(J, mismatch_vector)
        
        # Apply the Newton correction: update angles at PV/PQ buses, while the
        # slack angle remains the reference; update magnitudes only at PQ buses,
        # keeping slack and PV voltage magnitudes fixed.
        delta[pq_pv_indices] += updates[:num_angles]
        V[pq_indices] += updates[num_angles:]

        print(f"Iteration {iteration + 1}: Max mismatch = {mismatch_max:.6f}")
        # Count this iteration because the Newton state correction was applied.
        total_iter += 1
            
    # End the timer after convergence or after reaching the iteration limit.
    end_time = time.time()
    print(f"Processing time: {end_time - start_time:.4f} seconds.")
    print(f"Total iterations: {total_iter}")
    # Return solved voltage magnitudes and phase angles, in bus input order.
    return V, delta

if __name__ == "__main__":
    # When run as a script, use the first command-line argument as the input
    # path; otherwise use the provided five-bus sample filename.
    if len(sys.argv) > 1:
        filepath = sys.argv[1]
    else:
        filepath = "FiveBus_PQ.txt"
        print(f"No file provided via command line. Defaulting to '{filepath}'.")
        print("Usage for custom path: python your_script.py <path_to_FiveBus_PQ.txt>\n")
        
    try:
        # Run the solve and retain its returned voltage magnitudes and angles.
        voltages, angles = newton_raphson_power_flow(filepath)
    except FileNotFoundError:
        # Give a readable message for the common case of an invalid input path.
        print(f"Error: Could not find the file '{filepath}'. Please check the path and try again.")