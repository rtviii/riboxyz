import json
from urllib.request import Request
from ninja import Router, Body
import os
from pprint import pprint
from django.http import JsonResponse, HttpResponseServerError
from ninja import Path, Router, Schema
from django.http import FileResponse
from ninja.responses import Response
from ribctl import ASSETS, ASSETS_PATH, RIBETL_DATA
from ribctl.asset_manager.asset_manager import RibosomeAssetManager
from ribctl.asset_manager.asset_types import AssetType
from ribctl.lib.schema.types_ribosome import (
    ConstrictionSite,
    PTCInfo,
)
import tempfile
import subprocess
import os.path as osp
from enum import Enum
import trimesh
import numpy as np
from pathlib import Path

router_loci = Router()
TAG_LOCI = "Biologically Relevant Loci & Landmarks"

class MeshFormat(str, Enum):
    PLY = "ply"
    OBJ = "obj"
    STL = "stl"
    VTK = "vtk"

@router_loci.get("/tunnel_geometry", tags=[TAG_LOCI])
def get_shape(request, rcsb_id: str, is_ascii: bool = False, format: MeshFormat = MeshFormat.PLY):
    rcsb_id = rcsb_id.upper()
    
    # Default PLY file path
    file_path = AssetType.NPET_MESH_ASCII.get_path(rcsb_id)
    
    if not os.path.exists(file_path):
        return Response({"error": "Shape file not found"}, status=404)
    
    # If PLY format requested, return original file
    if format == MeshFormat.PLY:
        filename = (
            "{}_poisson_recon.ply".format(rcsb_id)
            if not is_ascii
            else "{}_poisson_recon_ascii.ply".format(rcsb_id)
        )
        try:
            file = open(file_path, "rb")
            return FileResponse(
                file, content_type="application/octet-stream", filename=filename
            )
        except IOError:
            return Response({"error": "Error reading the shape file"}, status=500)
    
    # For other formats, convert the PLY file
    try:
        # Load the PLY mesh using trimesh
        mesh = trimesh.load(file_path)
        
        # Create a temp file for the converted mesh
        with tempfile.NamedTemporaryFile(suffix=f".{format}", delete=False) as temp_file:
            temp_path = temp_file.name
        
        # Export to the requested format
        if format == MeshFormat.OBJ:
            mesh.export(temp_path, file_type='obj')
            content_type = "model/obj"
        elif format == MeshFormat.STL:
            mesh.export(temp_path, file_type='stl')
            content_type = "model/stl"
        elif format == MeshFormat.VTK:
            # Manual VTK export since trimesh doesn't support it directly
            try:
                # Write a legacy VTK file manually
                with open(temp_path, 'w') as f:
                    vertices = mesh.vertices
                    faces = mesh.faces
                    
                    # Write VTK header
                    f.write("# vtk DataFile Version 2.0\n")
                    f.write(f"Converted from PLY: {Path(file_path).name}\n")
                    f.write("ASCII\n")
                    f.write("DATASET POLYDATA\n")
                    
                    # Write vertices
                    f.write(f"POINTS {len(vertices)} float\n")
                    for v in vertices:
                        f.write(f"{v[0]} {v[1]} {v[2]}\n")
                    
                    # Write faces
                    f.write(f"POLYGONS {len(faces)} {len(faces) * 4}\n")
                    for face in faces:
                        f.write(f"3 {face[0]} {face[1]} {face[2]}\n")
                    
                    # If the mesh has vertex colors, add them
                    if hasattr(mesh, 'visual') and hasattr(mesh.visual, 'vertex_colors'):
                        colors = mesh.visual.vertex_colors
                        f.write(f"POINT_DATA {len(vertices)}\n")
                        f.write("COLOR_SCALARS vertex_colors 3\n")
                        for color in colors[:, :3]:  # Use RGB, ignoring alpha
                            normalized = [c/255.0 for c in color]
                            f.write(f"{normalized[0]} {normalized[1]} {normalized[2]}\n")
                
                content_type = "model/vtk"
            except Exception as e:
                return Response({"error": f"Error creating VTK file: {str(e)}"}, status=500)
        
        # Define the output filename
        output_filename = f"{rcsb_id}_poisson_recon.{format}"
        
        # Return the converted file
        converted_file = open(temp_path, "rb")
        response = FileResponse(
            converted_file, 
            content_type=content_type, 
            filename=output_filename,
            as_attachment=True
        )
        
        # Set cleanup function to delete temp file after sending
        response._resource_closers.append(lambda: os.unlink(temp_path) if os.path.exists(temp_path) else None)
        
        return response
        
    except Exception as e:
        return Response({"error": f"Error converting mesh: {str(e)}"}, status=500)


@router_loci.get("/cylinder_residues", tags=[TAG_LOCI], include_in_schema=False)
def cylinder_residues(request):
    try:
        with open("/Users/rtviii/dev/RIBETL_DATA/4UG0/artifacts/cylinder_residues.json", "rb") as f:
            map = json.load(f)
        return Response(map)
    except IOError:
        print("Couldn not read file")


@router_loci.get("/half_cylinder_residues", tags=[TAG_LOCI], include_in_schema=False)
def half_cylinder_residues(request):
    try:
        with open(
            "/home/rtviii/dev/npet-cg-sim/half_cylinder_residues.json", "rb"
        ) as f:
            map = json.load(f)
        return Response(map)
    except IOError:
        print("Couldn not read file")


@router_loci.get("/helices/{rcsb_id}", tags=[TAG_LOCI])
def get_helices(request, rcsb_id: str):
    rcsb_id = rcsb_id.upper()
    file_path = os.path.join("/home/rtviii/dev/riboxyz/7K00_rrna_helices.json")
    try:
        with open(file_path, "rb") as f:
            result = json.load(f)
        return Response(result)
    except IOError as e:
        return Response({"error": f"Error reading the helices file {e}"}, status=500)


@router_loci.get(
    "/ptc",
    response=PTCInfo,
    tags=[TAG_LOCI],
)
def structure_ptc(request, rcsb_id: str):
    params = dict(request.GET)
    rcsb_id = str.upper(params["rcsb_id"][0])
    ptc = RibosomeAssetManager().load_model(rcsb_id, AssetType.PTC)
    path = AssetType.PTC.get_path(rcsb_id)
    if not ptc:
        return JsonResponse({"error": "No PTC found for {}".format(rcsb_id)})
    return JsonResponse(ptc.model_dump())

@router_loci.get(
    "/constriction_site",
    response=ConstrictionSite,
    tags=[TAG_LOCI],
)
def constriction_site(request, rcsb_id: str,):
    params = dict(request.GET)
    rcsb_id = str.upper(params["rcsb_id"][0])
    try:
        constriction = RibosomeAssetManager().load_model(
            rcsb_id, AssetType.CONSTRICTION_SITE
        )
        if not constriction:
            return JsonResponse(
                {"error": "No constriction site found for {}".format(rcsb_id)}
            )
        return JsonResponse(constriction.model_dump())
    except Exception as e:
        return HttpResponseServerError(e)
