# import numpy as np
import time

def parse_fivebus_pq(filepath):
    system_data, buses, lines = {}, [], []
    with open(filepath, 'r') as f:
        print("Parsing FiveBus_PQ data from file...")
        for line in f:
            if line.startswith('%') or not line.strip():
                continue
            print("Processing line:", line.strip())
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
    print("Inside newton_raphson_power_flow function")
    start_time = time.time()
    print(f"Start time: {start_time}")
    system_data, buses, lines = parse_fivebus_pq(filepath)
    Ybus = build_ybus(buses, lines)

    print(f"Ybus Matrix:\n{Ybus}")
    
    # Initialize voltage magnitudes and angles
    V = np.array([bus['volts'] for bus in buses], dtype=float)
    delta = np.zeros(len(buses), dtype=float)
    
    base_mva = system_data['baseMVA']
    tolerance = system_data['tolerance']
    max_iter = 10
    
    print(f"Starting N-R Power Flow for System: {system_data['name']}")
    
    for iteration in range(max_iter):
        # Placeholder for power mismatch calculations (dP, dQ)
        # Placeholder for Jacobian matrix formation (J1, J2, J3, J4)
        # Placeholder for solving: [dDelta, dV] = inv(J) * [dP, dQ]
        
        # Simulated convergence check
        mismatch_max = 0.0 # Calculate actual max mismatch here
        if mismatch_max < tolerance:
            print(f"Converged in {iteration} iterations.")
            break
            
    end_time = time.time()
    print(f"Processing time: {end_time - start_time:.4f} seconds.")
    return V, delta

if __name__ == "__main__":
    # Ensure a file named FiveBus_PQ.txt is in the same directory
    try:
        voltages, angles = newton_raphson_power_flow("FiveBus_PQ.txt")
    except FileNotFoundError:
        print("Please save the Chapman FiveBus_PQ data into a text file named 'FiveBus_PQ.txt' to execute.")