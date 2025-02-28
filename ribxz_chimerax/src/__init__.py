# src/ribxz_chimerax/__init__.py

from chimerax.core.toolshed import BundleAPI

class _RibosomeXBundleAPI(BundleAPI):
    api_version = 1
    
    @staticmethod
    def register_command(command_name, logger):
        # Register the command with ChimeraX
        from . import cmd
        cmd.register_command(command_name, logger)
    
    @staticmethod
    def start_tool(session, tool_name):
        # Not implementing a GUI tool in this example
        pass

# This is the critical line - it MUST match what ChimeraX expects
bundle_api = _RibosomeXBundleAPI()