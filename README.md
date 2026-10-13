# cemaSourceTermsDNS (branch `reactingDNS`)

OpenFOAM-8 post-processing utility that splits the species (`YEqn`) and
energy (`EEqn`) transport equations of the `reactingDNS` solver into
**chemical** and **non-chemical** source terms for one saved time step. The
output is intended as input for Chemical Explosive Mode Analysis (CEMA).

This branch is the `reactingDNS` counterpart of `cemaSourceTerms` (branch
`main`, for `reactingFoam`). `reactingDNS` does not use OpenFOAM's
thermophysical transport model. It computes mixture-averaged transport
properties itself (`solvers/reactingDNS/transport`), and its diffusion terms
differ from `divj` / `divq`. This utility includes the same transport headers
and rebuilds the left-hand side of `Y-hEqn_DD.H` (or `Y-hEqn.H` when
`differentialDiffusion off`) without the `ddt` and convection terms.

For each active species `Yi` and for the enthalpy `he`, the utility evaluates
the terms explicitly from the saved fields:

| Term | Expression | Output field | Units |
| --- | --- | --- | --- |
| Convection | ∇·(φ Yi) | `conv_<Yi>` | kg/m³/s |
| Diffusion | sum of the four terms below | `diff_<Yi>` | kg/m³/s |
| – Fick | −∇·((ρDᵢ + μt/Sct)∇Yi) (Le=1: −∇·((α + μt/Sct)∇Yi)) | `diffFick_<Yi>` | kg/m³/s |
| – correction velocity | ∇·(φ_c Yi), φ_c = (Σₖ ρVₖ)·S_f | `diffCorr_<Yi>` | kg/m³/s |
| – mean molecular weight | −∇·(ρDᵢ Yi/W ∇W) | `diffMW_<Yi>` | kg/m³/s |
| – Soret (H, H2 only) | −∇·(−ρDᵢ Θᵢ Yi/(T Xᵢ) ∇T) | `diffSoret_<Yi>` | kg/m³/s |
| Non-chemical | −conv − diff | `nonChem_<Yi>` | kg/m³/s |
| Chemical *(optional)* | ω̇ᵢ (from `reaction->R(Yi)`) | `chem_<Yi>` | kg/m³/s |
| Convection | ∇·(φ he) | `conv_enthalpy` | W/m³ |
| Diffusion | sum of the three terms below | `diff_enthalpy` | W/m³ |
| – conduction | −∇·((α + μt/Prt)∇h) | `diffCond_enthalpy` | W/m³ |
| – species enthalpy (α) | Σᵢ ∇·(α hs,ᵢ ∇Yi) (`diffsh`) | `diffsh_enthalpy` | W/m³ |
| – species enthalpy flux | −∇·(Σᵢ hs,ᵢ ρDᵢ…) (`−fvc::div(vk)`) | `diffFlux_enthalpy` | W/m³ |
| Pressure work | ∂p/∂t (`fvc::ddt(p)`; for `e`: −∇·(p u)) | `dpdt_enthalpy` | W/m³ |
| Non-chemical | −conv − diff + dpdt | `nonChem_enthalpy` | W/m³ |
| Chemical *(optional)* | heat release `reaction->Qdot()` | `chem_enthalpy` | W/m³ |

The correction-velocity, molecular-weight, Soret and species-enthalpy terms
are zero when `differentialDiffusion off`. The Soret term is non-zero only
when `thermalDiffusion on`. The component fields (`diffFick_` etc.) are
written only with `-writeDiffTerms`.

`rho` and `mixAlpha_post` (the recomputed `mixAlpha`, for comparison with the
solver output `mixAlpha`) are also written. All fields go into the selected
time folder.

The kinetic energy terms (`ddt(ρK)`, `div(φ,K)`) are not included in
`nonChem_enthalpy`.

`dpdt_enthalpy` is evaluated from the saved `p` and `p_0` with
`fvc::ddt(p)`, in the same way as `pEqn.H` (no pressure iteration is needed).
It is zero when `dpdt no` is set in `thermophysicalProperties`. The solver
computes `dpdt` at the end of each time step and uses it in the `EEqn` of the
next step. So the value written for time n is the one the solver uses in step
n+1, and it differs by O(Δt) from the value used in step n. On the first step
after a restart, the solver uses `dpdt = 0`.

## Requirements

- OpenFOAM-8 (the commands below use the `openfoam/openfoam8-paraview56:8`
  Docker image)
- A `reactingDNS` case that uses `psiReactionThermo`
- `constant/chemistryProperties` (`differentialDiffusion`, `thermalDiffusion`,
  `Prt`, `Sct`) and the transport data read by the solver: `trandat`,
  `binaryDiff`, `thermoDiff`, `speciesLambda`, `speciesMu`
- The `reactingDNS` source tree at `../../solvers/reactingDNS` (the transport
  headers are included from there at build time)
- In the selected time folder: `T`, `p`, `U`, the species mass fractions, and
  any turbulence fields the case needs. `phi` is read if it is there and
  computed from `rho` and `U` if it is not.
- `p_0` (the old-time pressure, written by the solver) in the selected time
  folder. Without it, `dpdt_enthalpy` is zero.
- The case's `system/fvSchemes` must define `div(phi,Yi_h)`. The convection
  terms are evaluated with that scheme.

## Environment (Docker)

```bash
docker run -it --rm \
  --user $(id -u):$(id -g) \
  -v /NAS/29/shioyoke/foam8/run:/home/openfoam/run \
  -v /NAS/29/shioyoke/foam8/applications:/home/openfoam/applications \
  -v /NAS/29/shioyoke/foam8/platforms:/home/openfoam/platforms \
  -v /NAS/29/shioyoke/KU_MILD:/home/openfoam/KU_MILD \
  -e DISPLAY=$DISPLAY \
  -v /tmp/.X11-unix:/tmp/.X11-unix \
  openfoam/openfoam8-paraview56:8 \
  bash
```

## Build

```bash
cd applications/utilities/cemaSourceTerms/
wmake
```

The executable is installed to `$FOAM_USER_APPBIN/cemaSourceTermsDNS`, so it
can coexist with `cemaSourceTerms` built from `main`.

## Usage

Run it from the case directory:

```bash
cemaSourceTermsDNS <time-folder> [-calculateChem] [-writeDiffTerms]
```

- `<time-folder>`: name of an existing time directory, for example `0.05`.
  The utility stops with an error if that directory does not exist.
- `-calculateChem`: also compute the chemical source terms (`chem_<Yi>`,
  `chem_enthalpy`). Without this flag, only the non-chemical terms are written.
- `-writeDiffTerms`: also write the individual diffusion terms
  (`diffFick_<Yi>`, `diffCorr_<Yi>`, `diffMW_<Yi>`, `diffSoret_<Yi>`,
  `diffCond_enthalpy`, `diffsh_enthalpy`, `diffFlux_enthalpy`).

Example using the tutorial case:

```bash
cd run/1D_H2_premixed_inletoutlet_MA/
cemaSourceTermsDNS 0.005                 # non-chemical terms only
cemaSourceTermsDNS 0.005 -calculateChem  # non-chemical + chemical terms
```

### Notes on `-calculateChem`

- The chemistry is integrated over the simulation time step, which is read
  from `deltaT` in `<time-folder>/uniform/time`. That file must exist, so
  use a time folder that the solver wrote.
- `chem_<Yi>` is taken from the explicit source of `reaction->R(Yi)`
  (`-source/V`). This is only exact for combustion models that add the
  reaction rate as a pure explicit source (`Su += RR`), such as `laminar`,
  `PaSR` and `EDC`. Models that add implicit `fvm::Sp()` terms
  (`infinitelyFastChemistry`, `diffusion`, `FSD`, semi-implicit treatments)
  need a different extraction. See the comment in
  [cemaSourceTerms.C](cemaSourceTerms.C).

### Differences from the solver

- The solver evaluates the transport properties, the correction velocity and
  `vk` at intermediate PIMPLE iterations, and it updates `Yi` species by
  species. This utility evaluates everything from the final saved fields, so
  the terms differ from the solver's by a small amount.
- The inert species has no transport equation in the solver
  (`Y_inert = 1 − ΣYi`). Its terms are evaluated with the same operator as
  the other species.
- In the solver's discretisation, `Σᵢ diff_<Yi>` is not exactly zero. The
  correction flux uses the `div(phi,Yi_h)` convection scheme, so mass is not
  conserved discretely.

## Files

| File | Purpose |
| --- | --- |
| [cemaSourceTerms.C](cemaSourceTerms.C) | Main program: argument handling, transport properties and chemical terms |
| [YEqnTerms.H](YEqnTerms.H) | Species convection / diffusion terms (reactingDNS form) |
| [EEqnTerms.H](EEqnTerms.H) | Enthalpy convection / diffusion terms (reactingDNS form) |
| [createFields.H](createFields.H) | Builds thermo, U, phi, turbulence, thermophysical transport and combustion models |
| [createFieldRefs.H](createFieldRefs.H) | Convenience references (`psi`, `T`, inert specie index) |
| [Make/](Make/) | `wmake` build files |
