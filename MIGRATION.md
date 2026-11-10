# CIDEP Migration

## Status and Permission

Implementation authorized on 2026-11-10, including project files, dependency
downloads, tooling installation, builds, and tests. Original files under
FORTRAN77 remain unchanged. This document distinguishes implemented behavior
from planned work; incomplete features must not be interpreted as validated.

Both numerical engines, the interactive notebook, and the WPF application are
implemented. Scientific validation against historical output remains open.

## Shared Design

Preserve the FORTRAN scientific equations, field conventions (Gauss), and
frequency units (GHz). Separate nuclear-state generation, polarization,
discretization, convolution, persistence, and presentation. Replace COMMON
blocks with explicit configuration and result objects. Use double precision;
expect rounding differences from the original single-precision calculations.

Validate state degeneracies against exact multiplets, polarization against
small independently calculated examples, and both ports against the same
fixtures. Without a working original executable or reference datasets,
cross-port agreement alone does not establish historical numerical parity.

## Jupyter Design

Python modules under Jupyter provide the numerical engine and persistence.
NumPy/SciPy supply array, polynomial-root, and convolution operations. A notebook uses widgets
and interactive Matplotlib plots for inputs and mouse/arrow cursor placement.
The notebook uses browser controls, not native Windows menus. Keep numerical
work independent of the notebook so it can be tested without a display.
Independent lineshape offsets use a thread pool. Notebook calculations run
outside the event-loop thread. ipympl provides mouse and keyboard events.
An automatic-canvas resize experiment was removed after browser validation
failed; the standard resizable canvas uses horizontal scrolling in narrow panes.

## C# Design

CSharp contains a separately implemented numerical core and WPF application.
Use native Windows menus, input controls, and a plot with mouse and keyboard
cursor placement. Keep the core independent of WPF for automated testing.
Parallelize independent work only, avoiding nondeterministic shared sums.
The implementation targets .NET 10 and uses MathNet.Numerics polynomial roots
and OxyPlot.Wpf plotting. Parallel.For computes independent polarization rows,
lineshape offsets, and output convolution samples. Task.Run keeps WPF responsive.

## Persistence Design

Use a versioned JSON format with explicit units, configuration, mechanism
components, and sampled results. Do not write the legacy binary format.
The FORTRAN build uses /4I2, suggesting 16-bit integers; REAL defaults and
FORM='binary' layout must be checked before declaring import compatibility.
Different writers store stick and shaped spectra with different headers.
Reject malformed or ambiguous files rather than guessing silently.

The implemented version-1 JSON envelope has `version`, `units`, `source`,
`metadata`, nullable `config` and `shape`, and `data`. Data includes a uniformly
spaced increasing field axis, mixed/equilibrium samples, optional mechanism and
shaped arrays, and a scalar total_degeneracy when known. Import preserves header
values in metadata. Missing legacy configurations are not invented. Shaped
legacy imports have zero equilibrium placeholders because that component was
not stored separately; their original total is retained in mixed and shaped.

Supported binary layouts (no record markers, no trailing bytes):

- Stick / mixed: int16 count; real32 sweep, centre, frequency; count polarized
  real32 samples; count equilibrium real32 samples. Size: 14 + 8 * count bytes.
- Shaped: int16 count; real32 sweep, centre, frequency, microwave, delay,
  interval; int16 samples; count real32 values. Size: 28 + 4 * count bytes.

Both importers require 2..2048 samples, finite values, and positive field and
frequency metadata. Other compiler real encodings and record-marked files are
not supported. Only synthetic files have been verified so far.

## Decisions and Open Verification

- Keep FORTRAN77 in its existing location without modifying it.
- Preserve suspected scientific defects as documented decisions pending review;
  do not silently repair mathematics while claiming exact compatibility.
- Guard invalid inputs, zero normalization ranges, and excessive state counts.
- Initial discriminating check: one group of two spin-1/2 nuclei must produce
  projections -1, 0, 1 with degeneracies 1, 2, 1.
- Pending: representative historical output and actual legacy binary files.

## Scientific Mapping and Departures

- STICK.FOR first/second/ami map to group and radical states. Repeated polynomial
  convolution computes exact multiplet degeneracies without enumerating every
  individual spin function. Total-J multiplicities are differences of adjacent
  projection degeneracies. This handles higher nuclear spins without zero-weight
  states from the original hard-coded first multiplicity.
- POL_FLD.FOR pol/option/field2 supply Q dependence, the Adrian scale constants,
  resonance conversion 714.4844244, and the positive-energy/low-field convention.
  The original square-root-minus-linear dependence and configurable positive Q
  exponent are available. Pair scales are S=1.138871796, T=-1.138871796,
  F=-0.281368326.
- Second-order corrections and delta-g use the configured centre field before
  optional resonance recentering, as in the original call order. Second-order
  normal sticks activate only for radicals with a coupling above 10 G. ST-minus
  always uses total-J states and second-order shifts.
- The saved legacy spectrum represents radical 1. Both ports explicitly select
  observed radical 1 or 2, rather than silently persisting the temporary combined
  display of both radicals. Combined-radical concentration fitting is not added.
- STM.FOR supplies multiplet and counter-radical net terms. The per-group
  degeneracy weighting in sumstm is retained, including its non-product form
  for multiple groups. The original non-identical-radical loop repeats radical 2
  indefinitely; finite per-radical evaluation implements the stated purpose,
  not that defective control flow. ST-minus is not rescaled by S/F pair choice,
  matching its triplet-pair provenance in the original menu.
- TM is a signed multiple of equilibrium degeneracies. Components can be mixed
  with signed, including zero, coefficients. Optional normalization scales each
  component to a range of 100 without baseline subtraction. Zero ranges are
  left unchanged instead of dividing by zero. The initial first-radical component
  now follows the same normalization rule as subsequent mixtures; the legacy
  first-component cstkt assignment bypassed that rule.
- CALC_ETC.FOR const/calc/conste/calce supply the cubic/quartic coefficients,
  residues, steady-state term, and geometric sample sum. Generic library roots
  replace bespoke solvers. Near-coincident roots produce an explicit error,
  rather than using uninitialized residues. Small complex exponential increments
  use stable expm1 calculations. Time is microseconds; offset conversion is
  17.60844567 rad MHz/G.
- The even lineshape retains the duplicated zero-offset centre and legacy
  half-kernel trim convention, not a silently substituted Lorentzian/Gaussian.
  Output has the input sample count. Slow exchange groups occupied bins by
  equilibrium degeneracy and uses the full simulated degeneracy even if part of
  the spectrum lies outside the sweep. It requires a generated configuration;
  arbitrary imported spectra do not contain the total nuclear degeneracy.
- The legacy group routine's uninitialized line count, integer overflow limits,
  uninitialized COMMON values, printer output, and DOS graphics behavior are not
  reproduced. Two modern engines use double precision and explicit state.
- Inputs allow up to 4096 resolved states per radical and 16384 spectrum points;
  original fixed buffer limits were 1024 states and 2048 points. Group count,
  spin, sampling, finite values, and excessive state products are validated.

## Verification Record

- Python: exact spin-1/2 multiplet, total-J state conservation through spin 3/2,
  independent two-line RPM calculation, second-order shift, zero-coupling case,
  mixture identity, cubic/quartic finite responses, zero-delay limits, JSON
  roundtrip, malformed input, synthetic binary layouts, and precise cursor tests.
- C#: the corresponding analytic and persistence checks, plus seven shared
  fixtures covering defaults, mixing, distinct radicals, observed radical 2,
  nuclear spin 1, relaxation, and slow exchange. Absolute tolerance 1e-7 plus
  relative tolerance 1e-8; fixtures are generated with the Python engine and
  are regression references, not original FORTRAN output.
- WPF: clean build, default and shaped calculations, cursor key events,
  1200x820 and 800x620 screenshots with colored-pixel checks, and JSON export.
- Notebook: nbconvert execution with widget backend, browser calculation,
  nonblank canvas pixel check, mouse placement, and exact one-sample arrow step.
- Historical scientific parity and compiler-specific binary compatibility
  remain unverified. These require real spectra/files or a usable original
  executable and its unavailable DOS graphics dependencies.