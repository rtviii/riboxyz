import asyncio
from ribctl.etl.etl_collector import ETLCollector


asyncio.run(ETLCollector('5afi').process_structure())