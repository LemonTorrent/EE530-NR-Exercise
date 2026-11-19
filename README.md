# Newton-Raphson Power Flow

This program solves a power-flow problem using the Newton-Raphson method. It
reads the system, bus, and transmission-line data from a text file and prints
calculation progress and convergence information in the terminal.

## Requirements

- Python 3
- NumPy and Matplotlib (listed in `requirements.txt`)

## Create and activate a virtual environment

From the project directory, create a virtual environment:

```powershell
python -m venv .venv
```

Activate it in PowerShell:

```powershell
.\.venv\Scripts\Activate.ps1
```

If PowerShell prevents activation because of its execution policy, you can
activate the environment in Command Prompt instead:

```bat
.venv\Scripts\activate.bat
```

With the environment active, install the project dependencies:

```powershell
python -m pip install -r requirements.txt
```

Run the script using the same active environment, following the command below.
When finished, deactivate the environment with:

```powershell
deactivate
```

## Run with an input file

Pass the input file path as the first argument to `main.py`. From the project
directory, run the included example with:

```powershell
python main.py input_files\initial-fivebus_PQ.txt
```

You can also provide an absolute path. Quote paths that contain spaces:

```powershell
python main.py "C:\path\to\my power system.txt"
```

The path is resolved from the current working directory when it is relative.
The program prints the bus-admittance matrix, iteration details, convergence
status, processing time, and total iterations. It does not write results to an
output file.

## Compare processing time and iterations

Run `main.py` once for each tolerance you want to compare. Set the tolerance
on the `SYSTEM` record in the input file, then run the solver and note the
tolerance and the two summary values printed at the end of the command-line
output. For example:

```text
Processing time: 0.0054 seconds.
Total iterations: 2
```

Record the processing time in seconds and the total iteration count for that
run. Repeat with each tolerance, changing only the tolerance on the `SYSTEM`
record so the results are comparable.

Add the collected values to the corresponding lists in
`figures\plot-graphs.py`:

```python
tolerances = [0.1, 0.01, 0.001]
times = [0.0054, 0.0062, 0.0073]
iterations = [2, 2, 3]
```

Each position across the three lists must describe the same run: for example,
the time and iteration count at index `0` belong to tolerance `0.1`. Keep all
three lists the same length. The script sorts the values by tolerance and
generates two plots: `figures\nr_time_performance.png` (tolerance versus
processing time) and `figures\nr_iterations_performance.png` (tolerance
versus total iterations).

From the project directory, generate the plots with:

```powershell
python figures\plot-graphs.py
```

## Input file format

Input files are whitespace-delimited text files. Lines starting with `%` and
blank lines are ignored. Records use these formats:

```text
SYSTEM <name> <baseMVA> <tolerance>
BUS <name> <type> <volts> <Pgen> <Qgen> <Pload> <Qload> <Qcap>
LINE <from> <to> <Rse> <Xse> <Gsh> <Bsh> <rating>
```

Bus types used by the solver are `SL` (slack), `PV`, and `PQ`. Bus names in
`LINE` records must match names in `BUS` records. The included
`input_files\initial-fivebus_PQ.txt` provides a complete example.

If no input path is supplied, the script tries to open `FiveBus_PQ.txt` in the
current directory. To use the included sample, provide its path explicitly as
shown above. If the supplied file cannot be found, the program prints an error
message.
