import json
from pprint import pprint
from ribctl.etl.ribosome_assets import RibosomeAssets
from ribctl.lib.ribosome_types.types_ribosome import MitochondrialRNAClass


def get_stats():
    struct_stat = {
        "taxid":[]
    }

    for struct in RibosomeAssets.list_all_structs():

    # struct = "2FTC"
        print(struct)
        ra                 = RibosomeAssets(struct)
        profile            = ra.profile()
        d                  = profile.json()
        d                  = json.loads(d)
        d['mitochondrial'] = False
        for rna in profile.rnas:
            if len( rna.nomenclature )>0 :
               if ( rna.nomenclature[0].value in [k.value for k in list(MitochondrialRNAClass)] ):
                    print("{} : Mitochodrial rna {}".format(struct, rna.nomenclature), )
                    d['mitochondrial'] = True
                    break
        ra.write_own_json_profile(d, overwrite=True)
        taxid = profile.src_organism_ids[0]




get_stats()
