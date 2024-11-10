

#### Python:

Before i forget: ChimeraX runs with python3.11 and the on this arch machine this python is activated per-shell via `pyenv shell 3.11.9`


#### CMD:

putting this whole thign together:
What needs to be called from python? 

- jump into the correct dir
- feed models by `.mmcif` file (in the load script)

Chimerax headless mode (pseudocode, figure out how to): `chimerax   { open ribrepr.py } --nogui  "chimerax_make_movie.cxc" $rcsb_id`