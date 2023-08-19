
ribosome = ETLPipeline(query_rcsb_api(rcsb_single_structure_graphql("3J7Z"))).process_structure()
print(ribosome)