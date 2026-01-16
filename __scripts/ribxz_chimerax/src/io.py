# ribxz_chimerax/io.py

import aiohttp
import json
import asyncio

async def fetch_ribosome_data(pdb_id):
    """Fetch ribosome metadata from your server for a given PDB ID."""
    # Replace with your actual server URL
    server_url = "https://api.ribosome.xyz/api/structures/profile"
    
    async with aiohttp.ClientSession() as session:
        try:
            # Adjust this to match your API endpoint structure
            params = {'rcsb_id': pdb_id}
            async with session.get(server_url, params=params) as response:
                if response.status == 200:
                    data = await response.json()
                    return data
                else:
                    return {"error": f"Failed to fetch data: HTTP {response.status}"}
        except Exception as e:
            return {"error": f"Network error: {str(e)}"}
