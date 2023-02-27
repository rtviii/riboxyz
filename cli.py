from ribctl.lib.types.types_ribosome import RibosomeAssets


new = RibosomeAssets("4ug0")

new._verify_cif(True)
new._verify_json_profile(True)
new._verify_cif_modified(True)
new._verify_ligads_and_ligandlike_polys(True)
new._verify_chains_dir()

      
        