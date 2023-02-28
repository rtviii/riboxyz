import argparse

parser = argparse.ArgumentParser(description='Simulation presets')
parser .add_argument('s', '--struct', type=int, required=True,
                     help='')
parser .add_argument('db', '--initial_number', type=int,
                     help='Starting number of individuals')

args = parser .parse_args()
