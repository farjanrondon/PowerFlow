# Power Flow Solver (Python)

A Python program that solves the load flow of an electrical power system with Gauss-Seidel, Newton-Raphson, fast-decoupled and DC methods. The case is read from an Excel workbook and the results are written back to Excel.

More of my work is on my portfolio: [farjanrondon.github.io](https://farjanrondon.github.io).

## Features

- **Four load flow methods:** Gauss-Seidel, Newton-Raphson, fast-decoupled and DC. You choose which ones to run.
- **Network model:** the Y-bus is built from π-equivalent models of lines and of transformers with off-nominal taps, plus shunt elements.
- **Load model:** voltage-dependent ZIP loads (constant impedance, current and power).
- **Results for each method:** bus voltages and angles, calculated and generated power per bus, line power flows in both directions, line losses, the system power balance and convergence information.
- **Case studies:** each run saves a copy of the input workbook, under the output name you choose, with the results added. You can run several cases and compare them side by side.

## Requirements

- Python 3
- pandas, NumPy and openpyxl

```bash
pip install pandas numpy openpyxl
```

## How to run

1. Describe your system in `data_io.xlsx` (the input file read by `Main.py`).
2. Choose the methods and settings in the `CONFIG` sheet.
3. Run:

```bash
python Main.py
```

## Input workbook

| Sheet | Contents |
|---|---|
| `CONFIG` | Which methods to run (Y/N for GS, NR, FD and DC), convergence tolerance, maximum iterations and output file name. |
| `BUS` | Bus number and type (`SLACK`, `PV` or `PQ`), voltage magnitude (pu) and angle (degrees), generation and load P and Q (pu), ZIP fractions. |
| `LINES` | From and to bus, R, X and shunt susceptance B (pu). |
| `TRX` | From and to bus, R, X (pu) and tap. Leave it empty if the system has no transformers. |
| `SHUNT_ELEMENTS` | Bus and R, X (pu) of each shunt element. Leave it empty if there are none. |

The sheets are read by column position, with the first row as a header, so keep the column order of the template.

## Project structure

```
Main.py                    Entry point: reads the case, builds the Y-bus, runs the selected methods
readdata/ReadData.py       Reads the Excel sheets and builds line, transformer and shunt admittances
Ybus/Ybus.py               Nodal admittance matrix
OutputCopy/OutputCopy.py   Creates the output workbook as a copy of the input
MethodCalculus/Method.py   Gauss-Seidel, Newton-Raphson, fast-decoupled and DC solvers
MethodCalculus/S_Calculus.py   Line flows, bus powers, generation and power balance
MethodCalculus/WriteData.py    Writes the results to Excel
```

## Author

**Farjan Rondón**, Electrical Engineer, Universidad Simón Bolívar (Venezuela)
[farjanrondon.github.io](https://farjanrondon.github.io) · farjan.santos@gmail.com
