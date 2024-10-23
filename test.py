import click
import sys

def is_pipe():
    """Check if the script is being piped to"""
    return not sys.stdin.isatty()

@click.command()
@click.option('--uppercase', '-u', is_flag=True, help='Convert text to uppercase')
@click.option('--prefix', '-p', default='', help='Add prefix to each line')
def process_text(uppercase, prefix):
    """Process text from stdin, applying specified transformations."""
    if not is_pipe():
        ctx = click.get_current_context()
        click.echo(ctx.get_help())
        ctx.exit()

    for line in sys.stdin:
        line = line.strip()

        if uppercase:
            line = line.upper()
        if prefix:
            line = f"{prefix}{line}"

        print(line)

if __name__ == '__main__':
    process_text()
