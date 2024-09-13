import asyncio
from ribctl.etl.collector  import ETLCollector



asyncio.run(ETLCollector('7k00').process_structure())

# ETLCollector('5AFI').process_structure()