from ribctl.etl.ribosome_assets import RibosomeAssets


ra = RibosomeAssets('8cd1').profile().model_dump_json()
print(ra)