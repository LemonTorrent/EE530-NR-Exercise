# Newton-Raphson Power Flow

This program solves a power-flow problem using the Newton-Raphson method. It
reads the system, bus, and transmission-line data from a text file and prints
calculation progress and convergence information in the terminal.

## Requirements

- Python 3
- NumPy (listed in `requirements.txt`)

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

With the environment active, install the project dependency:

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
