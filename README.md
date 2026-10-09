# cemaSourceTerms

OpenFOAM-8 post-processing utility that splits the species (`YEqn`) and
energy (`EEqn`) transport equations of `reactingFoam` into **chemical** and
**non-chemical** source terms for one saved time step. The output is
intended as input for Chemical Explosive Mode Analysis (CEMA).

For each active species `Yi` and for the enthalpy `he`, the utility evaluates
the terms explicitly from the saved fields:

| Term | Expression | Output field | Units |
| --- | --- | --- | --- |
| Convection | ∇·(φ Yi) | `conv_<Yi>` | kg/m³/s |
| Diffusion | ∇·jᵢ (from `thermophysicalTransport->divj(Yi)`) | `diff_<Yi>` | kg/m³/s |
| Non-chemical | −conv − diff | `nonChem_<Yi>` | kg/m³/s |
| Chemical *(optional)* | ω̇ᵢ (from `reaction->R(Yi)`) | `chem_<Yi>` | kg/m³/s |
| Convection | ∇·(φ he) | `conv_enthalpy` | W/m³ |
| Diffusion | ∇·q (from `thermophysicalTransport->divq(he)`) | `diff_enthalpy` | W/m³ |
| Non-chemical | −conv − diff | `nonChem_enthalpy` | W/m³ |
| Chemical *(optional)* | heat release `reaction->Qdot()` | `chem_enthalpy` | W/m³ |

`rho` and `alphaEff` are also written. All fields go into the selected time
folder.

## Requirements

- OpenFOAM-8 (the commands below use the `openfoam/openfoam8-paraview56:8`
  Docker image)
- A `reactingFoam` case that uses `psiReactionThermo`
- In the selected time folder: `T`, `p`, `U`, the species mass fractions, and
  any turbulence fields the case needs. `phi` is read if it is there and
  computed from `rho` and `U` if it is not.
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

The executable is installed to `$FOAM_USER_APPBIN/cemaSourceTerms`.

## Usage

Run it from the case directory:

```bash
cemaSourceTerms <time-folder> [-calculateChem]
```

- `<time-folder>`: name of an existing time directory, for example `0.05`.
  The utility stops with an error if that directory does not exist.
- `-calculateChem`: also compute the chemical source terms (`chem_<Yi>`,
  `chem_enthalpy`). Without this flag, only the non-chemical terms are written.

Example using the tutorial case:

```bash
cd run/counterFlowFlame2D/
cemaSourceTerms 0.05                 # non-chemical terms only
cemaSourceTerms 0.05 -calculateChem  # non-chemical + chemical terms
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

## Files

| File | Purpose |
| --- | --- |
| [cemaSourceTerms.C](cemaSourceTerms.C) | Main program: argument handling and source-term evaluation |
| [createFields.H](createFields.H) | Builds thermo, U, phi, turbulence, thermophysical transport and combustion models |
| [createFieldRefs.H](createFieldRefs.H) | Convenience references (`psi`, `T`, inert specie index) |
| [Make/](Make/) | `wmake` build files |
