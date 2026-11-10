# Build and Run CIDEP

Run the following commands from the repository root in PowerShell. Generated
environments, build output, notebook checkpoints, and verification artifacts
are ignored by Git and are not included in the repository.

## Prerequisites

- Python 3.11 or newer, available as `python`.
- The .NET 10 SDK, available as `dotnet` (not just the .NET runtime).
- Windows to build and run the WPF desktop application.
- Internet access for the initial Python and NuGet package restores.

The Python notebook and .NET core tests can also run on other platforms.
Use `.venv/bin/python` instead of `.venv/Scripts/python.exe` outside Windows.
No legacy FORTRAN compiler is needed for either modern port.

## Jupyter Notebook

Create the environment and install the dependencies:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r Jupyter/requirements.txt
.\.venv\Scripts\python.exe -m ipykernel install --user --name cidep --display-name "CIDEP (.venv)"
.\.venv\Scripts\python.exe -m jupyter lab --ServerApp.root_dir=Jupyter
```

Open `CIDEP.ipynb`, select the **CIDEP (.venv)** kernel, and run all cells.
The notebook loads the Python implementation and creates the interactive UI;
there is no separate compilation step. In VS Code, open
[Jupyter/CIDEP.ipynb](Jupyter/CIDEP.ipynb) and select the same kernel.
Keep Jupyter's default token authentication enabled.

To regenerate a separate executed notebook without changing the source:

```powershell
New-Item -ItemType Directory -Force artifacts | Out-Null
.\.venv\Scripts\python.exe -m jupyter nbconvert --to notebook --execute Jupyter/CIDEP.ipynb --ExecutePreprocessor.kernel_name=cidep --ExecutePreprocessor.timeout=180 --output CIDEP.executed.ipynb --output-dir artifacts
```

Open the source notebook in a live Jupyter session to use the widgets; a saved
executed notebook does not replace a running kernel.

## C# Applications

Restore dependencies and build the core library, WPF application, and tests:

```powershell
dotnet restore CSharp/CIDEP.slnx
dotnet build CSharp/CIDEP.slnx --configuration Release --no-restore
dotnet run --project CSharp/Cidep.App --configuration Release --no-build
```

Build results are under each project's `bin/Release` and intermediates under
`obj`. To create a standalone Windows x64 distribution including the runtime:

```powershell
dotnet publish CSharp/Cidep.App/Cidep.App.csproj --configuration Release --runtime win-x64 --self-contained true --output artifacts/publish/win-x64
```

Run `artifacts/publish/win-x64/Cidep.App.exe`. Distribute the entire publish
directory, not just the executable.

## Verification

```powershell
.\.venv\Scripts\python.exe -m unittest discover -s Jupyter -v
dotnet test CSharp/Cidep.Tests/Cidep.Tests.csproj --configuration Release --nologo
dotnet run --project CSharp/Cidep.App --configuration Release -- --smoke-test
.\.venv\Scripts\python.exe Jupyter/verify_ports.py --check-csharp artifacts/csharp-spectrum.json
```

The Windows smoke test exits automatically and regenerates screenshots, a
summary, and the spectrum consumed by the final cross-port check in `artifacts`.
Keep [verification/fixtures.json](verification/fixtures.json): these are reviewed
regression-test inputs, not disposable build output. Regenerate them with
`python Jupyter/verify_ports.py` only after deliberately reviewing a scientific
change; see [MIGRATION.md](MIGRATION.md) for compatibility limits.

## Cleaning

With Jupyter and the desktop application stopped, remove `artifacts`, project
`bin` and `obj` directories, Python `__pycache__` directories, and notebook
`.ipynb_checkpoints` directories to discard generated output. The `.venv`
environment can also be deleted and recreated with the commands above. A local
SDK under `.tools` is optional; the commands here use the installed `dotnet` SDK.
The files under [FORTRAN77](FORTRAN77) are preserved historical source and build
scripts, not outputs of the modern builds.