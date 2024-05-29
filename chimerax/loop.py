import os
from chimerax.core.commands import run

for i in ["7K00", "4UG0"]:
    run(session,"ribmovie {}".format(i))