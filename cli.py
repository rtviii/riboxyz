from pprint import pprint
from ribctl.lib.types.types_ribosome import RibosomeAssets, RibosomeStructure
from ribctl.neo4j.ingress import init_driver


driver = init_driver()
pprint(driver)



with driver.session(database="system") as s:
    print(s.run("show default database").single())
    print(s.run("show default database").single())


      