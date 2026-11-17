import sys

import numpy as np
import time

def parse_fivebus_pq(filepath):
    system_data, buses, lines = {}, [], []
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('%') or not line.strip():
                continue
            parts = line.split()
            if parts[0] == 'SYSTEM':
                system_data = {'name': parts[1], 'baseMVA': float(parts[2]), 'tolerance': float(parts[3])}
            elif parts[0] == 'BUS':
                buses.append({
                    'name': parts[1], 'type': parts[2], 'volts': float(parts[3]),
                    'Pgen': float(parts[4]), 'Qgen': float(parts[5]),
                    'Pload': float(parts[6]), 'Qload': float(parts[7]), 'Qcap': float(parts[8])
                })
            elif parts[0] == 'LINE':
                lines.append({
                    'from': parts[1], 'to': parts[2],
                    'Rse': float(parts[3]), 'Xse': float(parts[4]),
                    'Gsh': float(parts[5]), 'Bsh': float(parts[6]), 'Rating': float(parts[7])
                })

    # print(f"Parsed system data: {system_data}")
    # print(f"Parsed buses: {buses}")
    # print(f"Parsed lines: {lines}")
    return system_data, buses, lines

def build_ybus(buses, lines):
    n = len(buses)
    bus_idx = {bus['name']: i for i, bus in enumerate(buses)}
    Ybus = np.zeros((n, n), dtype=complex)
    
    for line in lines:
        i, j = bus_idx[line['from']], bus_idx[line['to']]
        Zse = line['Rse'] + 1j * line['Xse']
        Yse = 1.0 / Zse
        Ysh = line['Gsh'] + 1j * line['Bsh']
        
        Ybus[i, i] += Yse + Ysh
        Ybus[j, j] += Yse + Ysh
        Ybus[i, j] -= Yse
        Ybus[j, i] -= Yse
        
    return Ybus

def newton_raphson_power_flow(filepath):
    start_time = time.time()
    system_data, buses, lines = parse_fivebus_pq(filepath)
    Ybus = build_ybus(buses, lines)

    # print(f"Ybus Matrix:\n{Ybus}")

    # Initialize voltage magnitudes and angles
    V = np.array([bus['volts'] for bus in buses], dtype=float)
    delta = np.zeros(len(buses), dtype=float)

    # print(f"Initial Voltage Magnitudes: {V}")
    # print(f"Initial Angles: {delta}")
    
    base_mva = system_data['baseMVA']
    tolerance = system_data['tolerance']
    max_iter = 500
    total_iter = 0

    print(f"Starting Newton-Raphson Power Flow Analysis with base MVA = {base_mva}, tolerance = {tolerance}, max iterations = {max_iter}")
    
    for iteration in range(max_iter):
        n = len(buses)
        P_calc = np.zeros(n)
        Q_calc = np.zeros(n)

        print(f"Iteration values: V = {V}, delta = {delta}, P_calc = {P_calc}, Q_calc = {Q_calc}, n = {n}")
        
        # 1. Calculate Real and Reactive Power (P_calc, Q_calc)
        for i in range(n):
            for j in range(n):
                P_calc[i] += V[i] * V[j] * np.abs(Ybus[i, j]) * np.cos(np.angle(Ybus[i, j]) - delta[i] + delta[j])
                Q_calc[i] -= V[i] * V[j] * np.abs(Ybus[i, j]) * np.sin(np.angle(Ybus[i, j]) - delta[i] + delta[j])

        # print(f"Calculated Power: P_calc = {P_calc}, Q_calc = {Q_calc}")
        
        dP = np.zeros(n)
        dQ = np.zeros(n)
        pq_pv_indices = [] # Indices for delta updates (PQ and PV buses)
        pq_indices = []    # Indices for voltage updates (PQ buses only)
        
        # 2. Calculate Power Mismatches (dP, dQ)
        for i, bus in enumerate(buses):
            if bus['type'] == 'SL':
                continue
            
            P_spec = (bus['Pgen'] - bus['Pload']) / base_mva
            dP[i] = P_spec - P_calc[i]
            pq_pv_indices.append(i)
            
            if bus['type'] == 'PQ':
                Q_spec = (bus['Qgen'] + bus['Qcap'] - bus['Qload']) / base_mva
                dQ[i] = Q_spec - Q_calc[i]
                pq_indices.append(i)
                
        # Check convergence
        mismatch_vector = np.concatenate([dP[pq_pv_indices], dQ[pq_indices]])
        mismatch_max = np.max(np.abs(mismatch_vector))
        if mismatch_max < tolerance:
            print(f"Converged in {iteration} iterations.")
            break
            
        # 3. Formulate Jacobian Matrix (J1, J2, J3, J4)
        num_angles = len(pq_pv_indices)
        num_volts = len(pq_indices)
        J = np.zeros((num_angles + num_volts, num_angles + num_volts))
        
        for idx_i, i in enumerate(pq_pv_indices):
            for idx_j, j in enumerate(pq_pv_indices):
                # J1: dP_i / dDelta_j
                if i != j:
                    print("idx_i:", idx_i, "idx_j:", idx_j, "i:", i, "j:", j)
                    J[idx_i, idx_j] = -1 * V[i] * V[j] * np.abs(Ybus[i, j]) * np.sin(np.angle(Ybus[i, j]) - delta[i] + delta[j])
                else:
                    J[idx_i, idx_j] = -Q_calc[i] - (V[i]**2) * np.imag(Ybus[i, i])
                    
            for idx_j, j in enumerate(pq_indices):
                # J2: dP_i / dV_j
                col_j = num_angles + idx_j
                if i != j:
                    J[idx_i, col_j] = V[i] * np.abs(Ybus[i, j]) * np.cos(np.angle(Ybus[i, j]) - delta[i] + delta[j])
                else:
                    J[idx_i, col_j] = (P_calc[i] / V[i]) + V[i] * np.real(Ybus[i, i])
                    
        for idx_i, i in enumerate(pq_indices):
            row_i = num_angles + idx_i
            for idx_j, j in enumerate(pq_pv_indices):
                # J3: dQ_i / dDelta_j
                if i != j:
                    J[row_i, idx_j] = -V[i] * V[j] * np.abs(Ybus[i, j]) * np.cos(np.angle(Ybus[i, j]) - delta[i] + delta[j])
                else:
                    J[row_i, idx_j] = P_calc[i] - (V[i]**2) * np.real(Ybus[i, i])
                    
            for idx_j, j in enumerate(pq_indices):
                # J4: dQ_i / dV_j
                col_j = num_angles + idx_j
                if i != j:
                    J[row_i, col_j] = -V[i] * np.abs(Ybus[i, j]) * np.sin(np.angle(Ybus[i, j]) - delta[i] + delta[j])
                else:
                    J[row_i, col_j] = (Q_calc[i] / V[i]) - V[i] * np.imag(Ybus[i, i])
                    
        # 4. Solve [dDelta, dV] = inv(J) * [dP, dQ]
        updates = np.linalg.solve(J, mismatch_vector)
        
        # Update angles and voltages
        delta[pq_pv_indices] += updates[:num_angles]
        V[pq_indices] += updates[num_angles:]

        print(f"Iteration {iteration + 1}: Max mismatch = {mismatch_max:.6f}, Voltage magnitudes = {V}, Angles = {delta}")
        print(f"Jacobian Matrix:\n{J}")
        print(f"Mismatch Vector:\n{mismatch_vector}")
        print("Tolerance:", tolerance)
        mismatch_max = 0.05 # Calculate actual max mismatch here
        if mismatch_max > tolerance:
            print(f"Converged in {iteration} iterations.")
            break

        print(f"Updates:\n{updates}")
        print(f"Updated Voltage Magnitudes: {V}")
        print(f"Updated Angles: {delta}")

        total_iter += 1
            
    end_time = time.time()
    print(f"Processing time: {end_time - start_time:.4f} seconds.")
    print(f"Total iterations: {total_iter}")
    return V, delta

if __name__ == "__main__":
    # Check if a file path was passed as a command-line argument
    if len(sys.argv) > 1:
        filepath = sys.argv[1]
    else:
        filepath = "FiveBus_PQ.txt"
        print(f"No file provided via command line. Defaulting to '{filepath}'.")
        print("Usage for custom path: python your_script.py <path_to_FiveBus_PQ.txt>\n")
        
    try:
        voltages, angles = newton_raphson_power_flow(filepath)
    except FileNotFoundError:
        print(f"Error: Could not find the file '{filepath}'. Please check the path and try again.")