# CIDEP

Two modern ports of the original Chemically Induced Dynamic Electron
Polarisation simulator. Original source remains in [FORTRAN77](FORTRAN77).
See [MIGRATION.md](MIGRATION.md) for numerical decisions and compatibility limits.

## Acknowlegements

CIDEP was written at the University of Oxford Physical Chemistry Department
by the group led by Professor Keith McLauchlan FRS, Emeritus Professor.

Other faculty and students contributed to making this program run on a
PDP-11/23, Norsk Data mini-computer and later MS-DOS PCs.

This port has been done by Paul D.Smith, D.Phil, a former student of Professor
McLauchlan, using AI to generate implementations in Jupyter Notebooks and
Microsoft C#.  The original FORTRAN-77 code on which this it is based, which is
itself a port to a 16-bit MS-DOS PC, is present in this repo.

No attempt has been made to validate these implementations beyond a casual look
to see if they "appear correct".

## Build and Run

See [BUILD.md](BUILD.md) for prerequisites, environment setup, rebuilding the
notebook and C# applications, publishing, and verification. After setup, run
these commands from the repository root in PowerShell.

### Windows Application

```powershell
dotnet run --project CSharp/Cidep.App
```

File menus open/save JSON and import legacy binary spectra. The parameter panel
contains editable nuclear-group tables, field and physics settings, mechanism
weights, and relaxation/exchange parameters. Calculate commits the inputs and
runs off the UI thread. Select a component above the plot. Click the plot to
place the cursor; left/right arrows then move it one simulated sample at a time.
Delete a selected nuclear-group row with Delete; the empty row adds a group.
An empty Q exponent uses the original square-root-minus-linear dependence.

### Jupyter Notebook

```powershell
.\.venv\Scripts\python.exe -m jupyter lab --ServerApp.root_dir=Jupyter
```

Open [Jupyter/CIDEP.ipynb](Jupyter/CIDEP.ipynb) and run its cell. In VS Code,
select `.venv/Scripts/python.exe` as the notebook kernel. The workbench provides
the same scientific parameters through notebook widgets. Open / Import accepts
JSON or legacy files; the JSON file input is a save path relative to the notebook
working directory. Click the plot before using left/right arrows.

The plot uses ipympl's resizable canvas. Narrow notebook panes scroll
horizontally rather than scale coordinates; desktop JupyterLab is recommended.
The native Windows application is not a mobile application.

### Fresh Environment

Use Python 3.11 or newer and the .NET 10 SDK on Windows.

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r Jupyter/requirements.txt
dotnet run --project CSharp/Cidep.App
```

The regular Jupyter launch command above retains Jupyter's token authentication.
The temporary verification server is loopback-only with authentication disabled;
do not expose that server through a network interface or proxy.

## Verify

```powershell
.\.venv\Scripts\python.exe -m unittest discover -s Jupyter -v
dotnet test CSharp/Cidep.Tests --nologo
dotnet build CSharp/CIDEP.slnx --nologo
dotnet run --project CSharp/Cidep.App -- --smoke-test
.\.venv\Scripts\python.exe Jupyter/verify_ports.py --check-csharp artifacts/csharp-spectrum.json
```

The native smoke test exits automatically and writes screenshots, a result
summary, and a JSON spectrum into ignored `artifacts/`. Cross-port fixtures are
committed in [verification/fixtures.json](verification/fixtures.json). Regenerate
them deliberately with `python Jupyter/verify_ports.py` only after reviewing a
scientific change; regeneration is not proof of correctness.

## Compatibility Status

Implemented: first/second-order sticks, ST-zero RPM, ST-minus RPM with optional
net contribution, TM, weighted mixing, transient relaxation, slow degenerate
exchange, JSON interchange, and candidate legacy binary import.

Analytic tests and agreement between independent ports have passed. Historical
parity is not yet established: no original executable, reference spectra, or
real legacy binary samples were supplied. Binary import currently accepts only
strict little-endian int16/IEEE-real32 stream layouts inferred from the source.