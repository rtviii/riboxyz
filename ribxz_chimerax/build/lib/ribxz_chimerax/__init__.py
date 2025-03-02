from chimerax.core.toolshed import BundleAPI

class _RibosomeXBundleAPI(BundleAPI):
    api_version = 1
    
    @staticmethod
    def register_command(bi, ci, logger):
        # Register commands with ChimeraX
        # The ChimeraX API passes: bundle_info, command_info, and logger
        command_name = ci.name if hasattr(ci, 'name') else None
        if command_name == "ribxz":
            from . import cmd_registry
            cmd_registry.register_ribxz_command(logger)
        else:
            logger.error(f"Unknown command name: {command_name}")
    
    @staticmethod
    def start_tool(session, tool_name):
        # Not implementing a GUI tool in this example
        pass

# This is the critical line - it MUST match what ChimeraX expects
bundle_api = _RibosomeXBundleAPI()