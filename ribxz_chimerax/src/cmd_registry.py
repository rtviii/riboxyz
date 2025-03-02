"""
Command registry for ribxz ChimeraX bundle.
This module handles the registration of all commands under the 'ribxz' namespace.
"""

def register_ribxz_command(logger):
    """Register the 'ribxz' command with ChimeraX with all subcommands."""
    from chimerax.core.commands import CmdDesc, register, create_alias
    
    # Import all command modules
    from . import cmd_loci
    from . import cmd_polymers
    
    # Main ribxz command description
    desc = CmdDesc(
        synopsis="Ribosome visualization and analysis tools"
    )
    
    # Register the main ribxz command
    register("ribxz", desc, _ribxz_cmd, logger=logger)
    
    # Register all subcommands
    cmd_loci.register_loci_command(logger)
    from . import cmd_polymers
    cmd_polymers.register_command("ribxz polymers", logger)
    
    # Create aliases for backward compatibility if needed
    create_alias("rx", "ribxz", logger=logger)

def _ribxz_cmd(session):
    """The main 'ribxz' command shows usage help."""
    from chimerax.core.commands import Command
    lines = [
        "Ribosome visualization and analysis tools",
        "",
        "Subcommands:",
        "  ribxz loci <rcsb_id> <loci_name>   - Render specific loci (tunnel, PTC, constriction)",
        "  ribxz polymers <structure> [rcsb_id]  - Apply ribosome representation with custom coloring",
        "",
        "Examples:",
        "  ribxz loci 4UG0 ptc                - Visualize peptidyl transferase center for 4UG0",
        "  ribxz polymers #1 4V6X                - Apply ribosome coloring to model #1 using 4V6X metadata"
    ]
    session.logger.info("\n".join(lines))