# cemaSourceTerms
Utility for obtaining chemical and non-chemical source terms in reactingFoam.  

Docker command for current setting is shown bellow.

```
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

Command too build the utilitie is shown bellow

```
cd applications/utilities/cemaSourceTerms/
wmake
```

To test utilitie in tutorial case, type bellow commands.

```
cd run/counterFlowFlame2D/
cemaSourceTerms 
```