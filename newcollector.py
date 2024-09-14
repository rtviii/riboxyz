import asyncio
from ribctl.etl.collector  import ETLCollector



asyncio.run(ETLCollector('3j9m').process_structure(overwrite=True))

# ETLCollector('5AFI').process_structure()