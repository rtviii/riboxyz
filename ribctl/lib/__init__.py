import json

print(f'Invoking __init__.py for {__name__}')

# open both subunit maps:
LSU_map = {k: v for k, v in json.load(open('ribctl/data/subunit_map_LSU.json', 'r')).items()}
SSU_map = {k: v for k, v in json.load(open('ribctl/data/subunit_map_SSU.json', 'r')).items()}