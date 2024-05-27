from chimerax.core.session import Session
from chimerax.dist_monitor import _DistMonitorBundleAPI
from chimerax.atomic import initialize_atomic
from chimerax.core.fetch import fetch_file
from chimerax.core.commands import run, runscript

session = Session('cx standalone')
_DistMonitorBundleAPI.initialize(session)
initialize_atomic(session)

#? -------------------------------------------


