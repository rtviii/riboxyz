
def verify_structure_profile_schema(rcsb_id:str):
    s = RibosomeAssets(rcsb_id).profile().model_dump_json()
    try:
        RibosomeStructure.model_validate_json(s)
        return True
    except Exception as e:
        print(e)
        return False

def verify_profile_exists(rcsb_id:str):
    return os.path.exists(RibosomeAssets(rcsb_id)._json_profile_filepath())