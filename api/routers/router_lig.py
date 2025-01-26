import json
import os
from pprint import pprint
from typing import List, Optional, Tuple, TypeAlias
from django.http import JsonResponse, HttpResponseServerError
from ninja import Field, Router
from pydantic import BaseModel
from ribctl import ASSETS_PATH
from ribctl.ribosome_ops import RibosomeOps
from ribctl.lib.libbsite import (
    bsite_ligand,
    bsite_transpose,
    lig_get_chemical_categories,
)
from ribctl.lib.schema.types_binding_site import BindingSite, LigandTransposition
from ribctl.lib.seq_project_many_to_one import compact_class, prepare_mapping_sources

router_lig = Router()
TAG = "Ligands, Antibitics & Small Molecules"


@router_lig.get("/binding_pocket", tags=[TAG], response=BindingSite)
def lig_nbhd(request, source_structure: str, chemical_id: str, radius: int = 5):
    """All members of the given RNA class: small and large subunit, cytosolic and mitochondrial RNA; tRNA."""
    path = RibosomeOps(source_structure).assets.paths.binding_site(chemical_id)
    if not os.path.exists(path):
        try:
            bsite = bsite_ligand(chemical_id, source_structure, radius)
            return JsonResponse(bsite.model_dump(), safe=False)
        except Exception as e:
            return HttpResponseServerError(e)

    with open(path, "r") as f:
        return JsonResponse(json.load(f), safe=False)


@router_lig.get("/transpose", tags=[TAG], response=LigandTransposition)
def lig_transpose(
    request,
    source_structure: str,
    target_structure: str,
    chemical_id: str,
    radius: float,
):
    """All members of the given RNA class: small and large subunit, cytosolic and mitochondrial RNA; tRNA."""
    bsite_path = RibosomeOps(source_structure).assets.paths.binding_site(chemical_id)
    prediction_path = RibosomeOps(
        target_structure
    ).assets.paths.binding_site_prediction(chemical_id, source_structure)

    if os.path.exists(prediction_path):
        with open(prediction_path, "r") as f:
            return JsonResponse(json.load(f), safe=False)

    bsite      = bsite_ligand(chemical_id, source_structure, radius)
    prediction = bsite_transpose(source_structure, target_structure, bsite)

    with open(prediction_path, "w") as f:
        json.dump(prediction.model_dump(), f, indent=4)
        pprint("Saved to file {}".format(prediction_path))

    with open(bsite_path, "w") as f:
        json.dump(bsite.model_dump(), f, indent=4)
        pprint("Saved to file {}".format(bsite_path))

    return JsonResponse(prediction.model_dump(), safe=False)


@router_lig.get("/bsite_composite", tags=[TAG], response=dict)
def bsite_composite(
    request,
):
    with open("bsite_composite.json", "r") as f:
        compacted_registry = json.load(f)
        pprint("Saved to file bsite_composite.json")
    return JsonResponse(compacted_registry, safe=False)



BindingSite: TypeAlias = List[Tuple[str, int]]

class LigandInput(BaseModel):
    """Input model for a single ligand"""
    chemicalId: str
    chemicalName: str
    purported_7K00_binding_site: Optional[List[List[str | int]]] = Field(default_factory=list)
    ribosome_ligand_categories: Optional[List[str]] = Field(default_factory=list)

class ProcessedLigand(BaseModel):
    """Model for processed ligand data"""
    chemicalId: str
    chemicalName: str
    purported_7K00_binding_site: List[List[str | int]] = Field(default_factory=list)

class CategoryData(BaseModel):
    """Model for data within each category"""
    items: List[ProcessedLigand] = Field(default_factory=list)
    composite_bsite: List[List[str | int]] = Field(default_factory=list)

class ProcessedLigands(BaseModel):
    """
    Root model for the processed ligands response
    The keys are category names and values are CategoryData objects
    """
    class Config:
        extra = "allow"  # Allows dynamic category names as keys
    
    @classmethod
    def from_dict(cls, data: dict) -> "ProcessedLigands":
        """Create a ProcessedLigands instance from a dictionary"""
        return cls(**{
            category: CategoryData(**category_data)
            for category, category_data in data.items()
        })



@router_lig.get("/demo_7k00", tags=[TAG], response=ProcessedLigands)
def demo_7k00(request):
    file = os.path.join(ASSETS_PATH, "ligands", "composite_bsites.json")

    if not os.path.exists(file):
        return HttpResponseServerError("File not found")
    with open(file, "r") as f:
        compacted_registry = json.load(f)
    return JsonResponse(compacted_registry, safe=False)