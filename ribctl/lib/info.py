from ribctl.etl.ribosome_assets import RibosomeAssets
from ribctl.lib.ribosome_types.types_ribosome import MitochondrialRNAClass


def get_stats():
    struct_stat = {
        "taxid":[]
    }

    for struct in RibosomeAssets.list_all_structs():
        # print(struct)
        print(struct)
        ra = RibosomeAssets(struct)
        profile = ra.profile()
        for rna in profile.rnas:
            if len( rna.nomenclature )>0 :
               if ( rna.nomenclature[0].value in [k.value for k in list(MitochondrialRNAClass)] ):
                    print("{} : Mitochodrial rna {}".format(struct, rna.nomenclature), )
                    d = profile.dict()
                    d['mitochondrial'] = True
                    print(d)
                    ra.write_own_json_profile(d, overwrite=True)
                    continue
        exit()
        taxid = profile.src_organism_ids[0]




get_stats()
