from chimerax.core.session import Session
def connect_atoms(session, atoms, to_atoms = None, distance = 2.2):
    session:Session
    session.
    session.logger.status(f'Made {len(bonds)} bonds between'
                          f' {len(atoms)} and {len(to_atoms)} atoms',
                          log = True)
    return bonds

def register_command(logger):
    from chimerax.core.commands import register, CmdDesc, FloatArg
    from chimerax.core.session import Session
    
    desc = CmdDesc(
                   synopsis='Connect close atoms')
    register('connect', desc, connect_atoms, logger=logger)

register_command(session.logger)