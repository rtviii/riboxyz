The overarchign idea is basically creating a coarse-grained simulation environment from the ribosome exit tunnel.
We do that by first extracting the empty space from the ribosome exit tunnel as voxels, then running a poisson surface reconstruction on it and then doing the same over the actual atom coordinates of the whole ribosome itself to tget the exterior surface. Then we combine the two. That's the method in the nutshell.

Currently it's been kind of hand-engineered and ad-hoc over the last 2 years and there are a lot of params that need to be tuned manually. Furthermore there are a lot of intermediate artifacts that get produced during the pipeline and at this point it's kind of slowing down the development a lot to inspect each to pinpoint what needs to change between runs of the pipeline. 

So my goal right now is to basically create a bit of a harness software for it that would both help me orchestrate the pipeline live but also allow me to inspect the intermediate artifacts more easily. I guess a small gui is what i mean where each stage and it's parameters can be run and inspected sequentially before advancing to the next stage and if something fails at one stage -- the user can see it, adjust the parameters and rerun this stage in isolation.

Let e show you the "orchestrator" i have going on at the moment. I'm mostly dissatisfied by it because it's slow, janky and doesn't really connect well to the artifacts that are already produced (i.e if i initialize the it with some strcutrure 7K00 -- which has already had all of its artifacts produced -- i want to be able to just load them in and inspect them rather than having to re-run the whole pipeline from scratch).


Also pyqt is a piece of shit and im wondering if we can switch to something web-based like nicegui. My concern here of course is that we use some pretty involved software there like pyvista which is finicky as fuck as is and putting it in the browser context seems like a taller task yet. What do u think?

Does that make sense? Can you help me with that?


Don't jump straight into the code, just talk to me about the tradeoffs. Let's chart a good course for improving this junk. Because eventually i want this to be the modus operandi for other methods where lots of intermediate artifacts are involved. 

I think the general architecture should be a viewer/orchestrator that is able to verify the completenes off predefined-path artifacts PER stage, knows how to visualize each and is able to rerun that stage with new parameters (ui-editable) on demand.



-----------

Ok cool. You are absolutely right in your assessment of my issues and your proposals are spot on. Let me try to show you some relevant parts of the code around my pipeline (some stage definitions, paths of the artifcats, the general computations) and we can try to refactor and re-engineer this "orchestrator" that i just attached according to the decoupled and modular architecture that you just described so well. 


First of all, this extraction pipeline is a part of my bigger database layout (but that doesn't matter much because the idea with filesystem-based asset managedement is the same, that's just to tell you that there might be other files there in the same folders that are not related to this pipeline). The primary data is organized as followsin in the `RIBETL_DATA` folder:
```
│   │   ├── 4P70_refined_clusters.png
│   │   └── 4P70_surface_points.png
│   └── TUNNELS
├── 4TUA
│   ├── 4TUA_PAR_STRUCTURE.cif
│   ├── 4TUA.cif
│   ├── 4TUA.png
│   ├── artifacts
│   │   ├── 4TUA_clusters.png
│   │   ├── 4TUA_mesh.png
│   │   ├── 4TUA_point_cloud.png
│   │   ├── 4TUA_refined_clusters.png
│   │   └── 4TUA_surface_points.png
│   └── TUNNELS
├── 4TUB
│   ├── 4TUB_PAR_STRUCTURE.cif
│   ├── 4TUB.cif
│   ├── 4TUB.png
│   ├── artifacts
│   │   ├── 4TUB_clusters.png
│   │   ├── 4TUB_mesh.png
│   │   ├── 4TUB_point_cloud.png
│   │   ├── 4TUB_refined_clusters.png
│   │   └── 4TUB_surface_points.png
│   └── TUNNELS
├── 4TUC
│   ├── 4TUC.cif
│   ├── 4TUC.png
│   ├── artifacts
│   │   ├── 4TUC_clusters.png
│   │   ├── 4TUC_point_cloud.png
ᢹ saeta.rtviii[ dev/RIBETL_DATA ]  tree 7K00  -L 5
7K00
├── 7K00_ALPHA_SHAPE_ascii.ply
├── 7K00_ALPHA_SHAPE.ply
├── 7K00_CONSTRICTION_SITE.json
├── 7K00_LIG_SPD.json
├── 7K00_NPET_MESH_ascii.ply
├── 7K00_NPET_MESH.ply
├── 7K00_PAR_STRUCTURE.cif
├── 7K00_PTC_COORDINATES.json
├── 7K00_PTC.json
├── 7K00_SPD_STRUCTURE.cif
├── 7K00_SPM_STRUCTURE.cif
├── 7K00.cif
├── 7K00.json
├── 7K00.png
├── artifacts
│   ├── 7K00_back_projected.npy
│   ├── 7K00_cluster_0.npy
│   ├── 7K00_cluster_1.npy
│   ├── 7K00_cluster_2.npy
│   ├── 7K00_cluster_3.npy
│   ├── 7K00_cluster_4.npy
│   ├── 7K00_clusters.png
│   ├── 7K00_constriction.npy
│   ├── 7K00_filtered_points.npy
│   ├── 7K00_interior_points.npy
│   ├── 7K00_largest_cluster.npy
│   ├── 7K00_mesh.png
│   ├── 7K00_normal_estimated_pcd.ply
│   ├── 7K00_point_cloud.png
│   ├── 7K00_ptc.npy
│   ├── 7K00_refined_cluster.npy
│   ├── 7K00_refined_clusters.png
│   ├── 7K00_surface_points.npy
│   ├── 7K00_surface_points.png
│   └── 7K00_transformed_points.npy
├── classification_report_7K00.json
├── test.py
└── TUNNELS
    └── tunnel_7K00.csv

3 directories, 288 files
```

And now here is the project code:
```
ᢹ saeta.rtviii[ dev/riboxyz ]  tree -I 'node_modules|venv|*__pycache__*|profiles|cache|debug_output|*json|*npy|*ply|*fasta|*csv|assets_*' -L 3
.
├── __scripts
│   ├── gettrna.py
│   ├── ligclass_distill_composite.py
│   ├── ligclass_distill.py
│   ├── ligclass.py
│   ├── maps_acquire.py
│   ├── masif_plugin.py
│   ├── merge_lig_info.py
│   ├── mesh_grid.py
│   ├── move_classification.sh
│   ├── move_tunnels.sh
│   ├── narstructs_tally
│   ├── npy_to_pdb.py
│   ├── process_prediction_7k00.py
│   ├── pymol_visualtion.py
│   ├── rcsb_cumulative_entries_barplot.py
│   ├── reclassify.py
│   ├── ribxz_chimerax
│   │   ├── build
│   │   ├── bundle_info.xml
│   │   ├── dist
│   │   ├── ribxz_chimerax.egg-info
│   │   └── src
│   ├── uniprot_seeds_query_record.py
│   └── williamson_assembly.py
├── api
│   ├── logs
│   ├── manage.py
│   ├── rbxz_bend
│   ├── reqs.txt
│   ├── ribxz_api
│   │   ├── asgi.py
│   │   ├── driver.py
│   │   ├── settings.py
│   │   ├── urls.py
│   │   └── wsgi.py
│   ├── routers
│   │   ├── router_lig.py
│   │   ├── router_loci.py
│   │   ├── router_mmcif.py
│   │   ├── router_polymers.py
│   │   └── router_struct.py
│   ├── staticfiles
│   │   ├── admin
│   │   ├── ninja
│   │   └── rest_framework
│   └── templates
│       └── swagger.html
├── chain_extractor.py
├── docker-compose.yml
├── Dockerfile-django
├── docs.md
├── kingdom_analysis.py
├── missing_meshes.py
├── neo4j_ribosome
│   ├── __archive
│   │   ├── cypher_ops
│   │   └── riboxyz_seed_data
│   ├── __init__.py
│   ├── cypher
│   │   └── list_filtered_structs.cypher
│   ├── db_driver.py
│   ├── db_lib_builder.py
│   ├── db_lib_reader.py
│   ├── node_ligand.py
│   ├── node_phylogeny.py
│   ├── node_polymer.py
│   ├── node_protein.py
│   ├── node_rna.py
│   └── node_structure.py
├── neo4j.conf
├── neo4j.old.conf
├── notes
│   ├── binding_affinities.md
│   ├── docs
│   │   ├── architecture_environment_configs.md
│   │   └── update.md
│   ├── factors.md
│   ├── functional_sites.md
│   ├── general_directions.md
│   ├── hmm-based-classification.md
│   └── pymol.md
├── npet_bulk_viewer.py
├── npet_inspect.py
├── npet_orchestrator.py
├── q.md
├── ribctl
│   ├── __init__.py
│   ├── asset_manager
│   │   ├── asset_manager.py
│   │   ├── asset_registry.py
│   │   ├── asset_types.py
│   │   ├── doc.md
│   │   └── parallel_acquisition.py
│   ├── assets
│   │   ├── 5T5H
│   │   ├── 7K00
│   │   ├── 85PD
│   │   ├── 8P5D
│   │   ├── 8P60
│   │   ├── classification_reports
│   │   ├── taxa.sqlite
│   │   ├── taxa.sqlite.traverse.pkl
│   │   └── taxdump.tar.gz
│   ├── etl
│   │   ├── __init__.py
│   │   ├── etl_collector.py
│   │   └── gql_querystrings.py
│   ├── global_ops.py
│   ├── lib
│   │   ├── __libseq.py
│   │   ├── chimerax
│   │   ├── enumunion.py
│   │   ├── info.py
│   │   ├── landmarks
│   │   ├── libbsite.py
│   │   ├── libhmm.py
│   │   ├── libmsa.py
│   │   ├── libseq.py
│   │   ├── libtax.py
│   │   ├── npet
│   │   ├── nsearch_gemmi.py
│   │   ├── ribosome_types
│   │   ├── schema
│   │   ├── seq_project_many_to_one.py
│   │   ├── thumbnail.py
│   │   ├── types
│   │   └── utils.py
│   ├── logger_config.py
│   ├── logs
│   │   ├── etl.log
│   │   └── loggers.py
│   ├── ribd.py
│   └── ribosome_ops.py
├── ribxz_logo_black.png
├── taxdump.tar.gz
└── visualize_mesh.py

42 directories, 96 files
```


and the pipeline code itself (npet lib):
```
ᢹ saeta.rtviii[ dev/riboxyz ]  tree -I 'node_modules|venv|*__pycache__*|profiles|cache|debug_output|*json|*npy|*ply|*fasta|*csv|assets_*' -L 5 ribctl/lib/npet/
ribctl/lib/npet/
├── _dbscan_pairs
├── alphalib.py
├── gpu_ops.py
├── kdtree_approach.py
├── npet_pipeline.py
├── pipeline
│   ├── base_stage.py
│   ├── clustering_stage.py
│   ├── constriction_identification_stage.py
│   ├── entity_filtering_stage.py
│   ├── exterior_mesh_stage.py
│   ├── landmark_identification_stage.py
│   ├── logs
│   ├── mesh_reconstruction_stage.py
│   ├── normal_estimation_stage.py
│   ├── point_cloud_processing_stage.py
│   ├── ptc_identification_stage.py
│   ├── refinement_stage.py
│   ├── setup_stage.py
│   ├── surface_extraction_stage.py
│   └── validation_stage.py
├── tunnel_asset_manager.py
└── visualization
    ├── bulk_artifact_viewer.py
    ├── dashboard_with_kingdom.html
    ├── dashboard.html
    ├── html_tally.py
    ├── landmark_ptcloud_collector.py
    ├── pipeline_status_tracker.py
    ├── various_visualization.py
    └── view_pipeline_artifacts.py

4 directories, 28 files
```


npet_pipeline:
```
from pathlib import Path
from typing import Dict, List, Optional, Any

import numpy as np

from ribctl.asset_manager.asset_types import AssetType
from ribctl.lib.npet.kdtree_approach import (
    filter_residues_parallel,
    landmark_constriction_site,
    landmark_ptc,
    ribosome_entities,
)
from ribctl.lib.npet.pipeline.setup_stage import SetupStage
from ribctl.lib.npet.pipeline.ptc_identification_stage import PTCIdentificationStage
from ribctl.lib.npet.pipeline.constriction_identification_stage import (
    ConstrictionIdentificationStage,
)
from ribctl.lib.npet.pipeline.exterior_mesh_stage import AlphaShapeStage
from ribctl.lib.npet.pipeline.landmark_identification_stage import (
    LandmarkIdentificationStage,
)
from ribctl.lib.npet.pipeline.entity_filtering_stage import EntityFilteringStage
from ribctl.lib.npet.pipeline.point_cloud_processing_stage import (
    PointCloudProcessingStage,
)
from ribctl.lib.npet.pipeline.clustering_stage import ClusteringStage
from ribctl.lib.npet.pipeline.refinement_stage import RefinementStage
from ribctl.lib.npet.pipeline.surface_extraction_stage import SurfaceExtractionStage
from ribctl.lib.npet.pipeline.normal_estimation_stage import NormalEstimationStage
from ribctl.lib.npet.pipeline.mesh_reconstruction_stage import MeshReconstructionStage
from ribctl.lib.npet.pipeline.validation_stage import ValidationStage
from ribctl.lib.npet.visualization.pipeline_status_tracker import (
    NPETProcessingTracker,
)

def generate_final_cluster_pdb(
    refined_cluster: np.ndarray,
    output_path: str,
    rcsb_id: str
) -> str:
    """
    Generate a PDB file for the final refined cluster used for mesh generation.
    
    Args:
        refined_cluster: The final cluster points array
        output_path: Where to write the PDB file
        rcsb_id: PDB ID for header
        
    Returns:
        str: Path to the generated PDB file
    """
    
    with open(output_path, 'w') as pdb_file:
        # Write PDB header
        pdb_file.write(f"HEADER    FINAL TUNNEL CLUSTER                    {rcsb_id}\n")
        pdb_file.write(f"TITLE     FINAL REFINED CLUSTER FOR {rcsb_id}\n")
        pdb_file.write("REMARK    Final cluster used for mesh generation\n")
        
        atom_num = 1
        connect_records = []
        
        print(f"Writing final cluster with {len(refined_cluster)} points")
        
        # Write each point as a CA atom in chain A, residue 1
        for point_idx, (x, y, z) in enumerate(refined_cluster):
            pdb_file.write(
                f"ATOM  {atom_num:5d}  CA  ALA A   1    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C  \n"
            )
            
            # Create bonds between consecutive atoms
            if point_idx > 0:
                connect_records.append(f"CONECT{atom_num-1:5d}{atom_num:5d}\n")
            
            atom_num += 1
        
        # Write chain terminator
        pdb_file.write(f"TER   {atom_num:5d}      ALA A   1\n")
        
        # Write all CONECT records
        for connect in connect_records:
            pdb_file.write(connect)
            
        # Write end record
        pdb_file.write("END\n")
    
    print(f"Generated final cluster PDB: {output_path}")
    return output_path

def generate_cylinder_residues_json(
    rcsb_id: str,
    cif_path: str,
    ptc_coords: np.ndarray,
    constriction_coords: np.ndarray,
    output_path: str,
    radius: float = 35.0,
    skip_nascent_chain: List[str] = [],
) -> dict:
    """
    Generate cylinder residues JSON for molstar visualization.

    Args:
        rcsb_id: PDB ID
        cif_path: Path to mmCIF file
        ptc_coords: PTC coordinates from pipeline context
        constriction_coords: Constriction site coordinates from pipeline context
        output_path: Where to write the JSON file
        radius: Cylinder radius in Angstroms
        skip_nascent_chain: Chain IDs to skip

    Returns:
        dict: The cylinder residues data structure
    """
    import json

    # Define cylinder between PTC and constriction site
    base_point = ptc_coords
    axis_point = constriction_coords

    # Calculate height as distance between landmarks plus some buffer
    height = np.linalg.norm(axis_point - base_point) + 10.0

    print(
        f"Cylinder: base={base_point}, axis={axis_point}, radius={radius}, height={height}"
    )

    # Load structure and get residues
    residues = ribosome_entities(
        rcsb_id=rcsb_id,
        cifpath=cif_path,
        level="R",  # Residue level
        skip_nascent_chain=skip_nascent_chain,
    )

    print(f"Loaded {len(residues)} residues from structure")

    # Filter residues that fall within cylinder
    cylinder_residues = filter_residues_parallel(
        residues=residues,
        base_point=base_point,
        axis_point=axis_point,
        radius=radius,
        height=120,
    )

    print(f"Found {len(cylinder_residues)} residues within cylinder")

    # Group by chain and extract residue numbers
    cylinder_data = {}
    for residue in cylinder_residues:
        chain_id = residue.get_parent().id  # auth_asym_id
        residue_num = residue.id[1]  # auth_seq_id

        if chain_id not in cylinder_data:
            cylinder_data[chain_id] = []
        cylinder_data[chain_id].append(residue_num)

    # Sort residue numbers within each chain
    for chain_id in cylinder_data:
        cylinder_data[chain_id].sort()

    print(
        f"Cylinder residues by chain: {dict((k, len(v)) for k, v in cylinder_data.items())}"
    )

    # Write to JSON file
    with open(output_path, "w") as f:
        json.dump(cylinder_data, f, indent=2)

    print(f"Wrote cylinder residues to {output_path}")
    return cylinder_data



def create_npet_mesh(
    rcsb_id: str, log_dir: Optional[Path] = None, force: bool = False
) -> NPETProcessingTracker:
    """
    Creates NPET mesh for a given RCSB ID with comprehensive tracking and logging.

    This function orchestrates the execution of all pipeline stages, handling setup,
    error management, and tracking of the overall process.

    Args:
        rcsb_id: The RCSB PDB identifier
        log_dir: Directory to store logs (default: logs)
        force: Whether to force regeneration of existing assets

    Returns:
        NPETProcessingTracker: Processing tracker with detailed status information
    """
    print(f"Creating NPET mesh for {rcsb_id}")
    meshpath = AssetType.NPET_MESH.get_path(rcsb_id)
    if meshpath.exists():
        print(f"NPET mesh already exists for {rcsb_id}")

    rcsb_id = rcsb_id.upper()
    tracker = NPETProcessingTracker(rcsb_id, log_dir)
    
    try:
        # Initialize pipeline context to share data between stages
        context: Dict[str, Any] = {}
        
        # Define stage parameters
        alpha_shape_params = {
            "d3d_alpha": 200,
            "d3d_tol": 10,
            "d3d_offset": 3,
            "kdtree_radius": 40,
            "max_nn": 60,
            "tangent_planes_k": 20,

            "PR_depth": 6,
            "PR_ptweight": 4,
        }
        
        point_cloud_params = {
            "radius": 35,
            "height": 120,
            "voxel_size": 1,
            "atom_size": 2,
        }
        
        clustering_params = {
            "epsilon": 5.5,
            "min_samples": 600,
        }
        
        refinement_params = {
            "epsilon": 3.5,
            "min_samples": 175,
        }
        
        surface_params = {
            "alpha": 2,
            "tolerance": 1,
            "offset": 2,
        }
        
        normal_params = {
            "kdtree_radius": 10,
            "kdtree_max_nn": 15,
            "correction_tangent_planes_n": 10,
        }
        
        mesh_params = {
            "depth": 6,
            "ptweight": 3,
        }
        
        # Run the Setup stage
        setup_stage = SetupStage(rcsb_id, tracker, None, force)
        setup_result = setup_stage.run(context)
        context.update(setup_result)
        
        # Create artifacts directory
        artifacts_dir = Path(context["ro"].assets.paths.dir) / "artifacts"
        artifacts_dir.mkdir(exist_ok=True, parents=True)
        
        # Run the PTC identification stage
        ptc_stage = PTCIdentificationStage(rcsb_id, tracker, artifacts_dir, force)
        ptc_result = ptc_stage.run(context)
        context.update(ptc_result)
        
        # Run the constriction identification stage
        constriction_stage = ConstrictionIdentificationStage(rcsb_id, tracker, artifacts_dir, force)
        constriction_result = constriction_stage.run(context)
        context.update(constriction_result)
        
        # Run the alpha shape stage
        alpha_stage = AlphaShapeStage(rcsb_id, tracker, artifacts_dir, True, alpha_shape_params)
        alpha_result = alpha_stage.run(context)
        context.update(alpha_result)
        
        # Run the landmark identification stage (verifies PTC and constriction)
        landmark_stage = LandmarkIdentificationStage(rcsb_id, tracker, artifacts_dir, force)
        landmark_result = landmark_stage.run(context)
        context.update(landmark_result)
        
        # Run the entity filtering stage
        entity_filtering_stage = EntityFilteringStage(
            rcsb_id, tracker, artifacts_dir, force,
            radius=point_cloud_params["radius"],
            height=point_cloud_params["height"]
        )
        entity_result = entity_filtering_stage.run(context)
        context.update(entity_result)
        
        # Run the point cloud processing stage
        point_cloud_stage = PointCloudProcessingStage(
            rcsb_id, tracker, artifacts_dir, force,
            radius=point_cloud_params["radius"],
            height=point_cloud_params["height"],
            voxel_size=point_cloud_params["voxel_size"],
            atom_size=point_cloud_params["atom_size"]
        )
        point_cloud_result = point_cloud_stage.run(context)
        context.update(point_cloud_result)
        
        # Run the clustering stage
        clustering_stage = ClusteringStage(
            rcsb_id, tracker, artifacts_dir, force,
            epsilon=clustering_params["epsilon"],
            min_samples=clustering_params["min_samples"]
        )
        clustering_result = clustering_stage.run(context)
        context.update(clustering_result)
        
        # Run the refinement stage
        refinement_stage = RefinementStage(
            rcsb_id, tracker, artifacts_dir, force,
            epsilon=refinement_params["epsilon"],
            min_samples=refinement_params["min_samples"]
        )
        refinement_result = refinement_stage.run(context)
        context.update(refinement_result)
        
        # Run the surface extraction stage
        surface_stage = SurfaceExtractionStage(
            rcsb_id, tracker, artifacts_dir, force,
            alpha=surface_params["alpha"],
            tolerance=surface_params["tolerance"],
            offset=surface_params["offset"]
        )
        surface_result = surface_stage.run(context)
        context.update(surface_result)
        
        # Run the normal estimation stage
        normal_stage = NormalEstimationStage(
            rcsb_id, tracker, artifacts_dir, force,
            kdtree_radius=normal_params["kdtree_radius"],
            kdtree_max_nn=normal_params["kdtree_max_nn"],
            correction_tangent_planes_n=normal_params["correction_tangent_planes_n"]
        )
        normal_result = normal_stage.run(context)
        context.update(normal_result)
        
        # Run the mesh reconstruction stage
        mesh_stage = MeshReconstructionStage(
            rcsb_id, tracker, artifacts_dir, force,
            depth=mesh_params["depth"],
            ptweight=mesh_params["ptweight"]
        )
        mesh_result = mesh_stage.run(context)
        context.update(mesh_result)
        
        # Run the validation stage
        validation_stage = ValidationStage(rcsb_id, tracker, artifacts_dir, force)
        validation_result = validation_stage.run(context)
        context.update(validation_result)
        
        # Complete pipeline and set final status
        tracker.complete_processing(True, watertight=context.get("watertight", False))
        
    except Exception as e:
        # Handle any uncaught exceptions
        if tracker.current_stage:
            tracker.end_stage(tracker.current_stage, False, error=e)
        tracker.complete_processing(False)
    
    return tracker
```


alphalib.py:
```
import os
from pathlib import Path
import sys
sys.path.append("/home/rtviii/dev/riboxyz")
from ribctl.asset_manager.asset_types import AssetType
from ribctl import RIBXZ_TEMP_FILES
from ribctl.lib.npet.kdtree_approach import apply_poisson_reconstruction

from Bio.PDB.MMCIFParser import MMCIFParser
from ribctl.lib.npet.visualization.various_visualization import visualize_mesh, visualize_pointcloud
from ribctl.ribosome_ops import RibosomeOps

import numpy as np
import pyvista as pv
import open3d as o3d

def cif_to_point_cloud(cif_path: str, chains: list[str] | None = None,  do_atoms:bool=False):
    parser = MMCIFParser()
    structure = parser.get_structure("structure", cif_path)
    coordinates = []

    first_model = structure[0]
    if do_atoms:
        for chain in first_model:
            if chains is not None and chain.id not in chains:
                continue
            for residue in chain:
                for atom in residue:
                    coordinates.append(atom.get_coord())    
    else:
        for chain in first_model:
            if chains is not None and chain.id not in chains:
                continue
            for residue in chain:
                coordinates.append(residue.center_of_mass())

    if not coordinates:
        raise ValueError(f"No coordinates found in {cif_path}")

    return np.array(coordinates)

def validate_mesh_pyvista(mesh_or_path, stage="unknown"):
    """
    Validates if a mesh is watertight.
    
    Parameters:
        mesh_or_path: Either a PyVista mesh object or a path to a mesh file
        stage: Optional identifier for logging
        
    Returns:
        bool: True if the mesh is watertight, False otherwise
    """
    import pyvista as pv
    import os
    from pathlib import Path

    # If mesh is a path, load it
    if isinstance(mesh_or_path, (str, Path)):
        if not os.path.exists(mesh_or_path):
            print(f"WARNING: Mesh file does not exist: {mesh_or_path}")
            return False
        
        try:
            mesh = pv.read(str(mesh_or_path))
        except Exception as e:
            print(f"WARNING: Failed to load mesh at {mesh_or_path}: {e}")
            return False
    else:
        mesh = mesh_or_path

    if mesh is None:
        print(f"WARNING: Null mesh at stage {stage}")
        return False

    print(f"\nMesh properties at stage: {stage}")

    try:
        # Check watertightness by looking for boundary edges
        # A mesh is watertight if it has no boundary edges
        edges = mesh.extract_feature_edges(
            boundary_edges=True,
            feature_edges=False,
            manifold_edges=False,
            non_manifold_edges=False,
        )
        is_watertight = edges.n_cells == 0

        print(f"- Is watertight: {is_watertight}")
        print("DONE\n\n ")
        return is_watertight
    except Exception as e:
        print(f"WARNING: Error checking watertightness: {e}")
        return False

def quick_surface_points(
    pointcloud: np.ndarray, alpha: float, tolerance: float, offset: float
) -> np.ndarray:
    cloud = pv.PolyData(pointcloud)
    # Using larger tolerance and smaller offset for faster computation
    grid    = cloud.delaunay_3d(alpha=alpha, tol=tolerance, offset=offset, progress_bar=True )
    surface = grid.extract_surface().cast_to_pointset()
    return surface.points

def fast_normal_estimation(
    surface_pts: np.ndarray,
    kdtree_radius,
    max_nn,
    tangent_planes_k, 
) -> o3d.geometry.PointCloud:
    """
    Returns:
        o3d.geometry.PointCloud: Point cloud with estimated normals
    """
    pcd        = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(surface_pts)

    # Use hybrid search with reduced parameters
    search_param = o3d.geometry.KDTreeSearchParamHybrid( radius=kdtree_radius, max_nn=max_nn )

    pcd.estimate_normals(search_param=search_param)
    pcd.orient_normals_consistent_tangent_plane(k=tangent_planes_k)

    return pcd

def alpha_contour_via_poisson_recon(rcsb_id:str, verbose:bool=False):
    ptcloudpath = os.path.join(RIBXZ_TEMP_FILES, '{}_ptcloud.npx'.format(rcsb_id))
    rops                  = RibosomeOps(rcsb_id)
    cifpath               = rops.assets.paths.cif


    print("Cifpath:", cifpath)
    if not os.path.exists(ptcloudpath):
        print("Extracting point cloud from CIF file")
        first_assembly_chains = rops.first_assembly_auth_asym_ids()
        ptcloud               = cif_to_point_cloud(cifpath, first_assembly_chains, do_atoms=True)
        np.save(ptcloudpath, ptcloud)
    else:
        print("Loaded.")
        ptcloud = np.load(ptcloudpath)


    output_normals_pcd = os.path.join(RIBXZ_TEMP_FILES, "{}_normal_estimated_pcd.ply".format(rcsb_id))
    output_mesh        = AssetType.ALPHA_SHAPE.get_path(rcsb_id)

    d3d_alpha  = 75    # Increase from 35 - be more aggressive
    d3d_tol    = 4     # Increase tolerance to smooth out small details
    d3d_offset = 3     # Slightly larger offset

    kdtree_radius   = 30    # Larger radius to catch more global structure
    max_nn          = 40    # More neighbors for more robust estimation
    tanget_planes_k = 15    # Actually decrease this to avoid over-smoothing

    PR_depth    = 6        # Reduced from 6
    PR_ptweight = 4

    print("Beginning Delaunay 3d reconstruction")
    surface_pts = quick_surface_points(ptcloud, d3d_alpha, d3d_tol, d3d_offset)

    if verbose:
        visualize_pointcloud(surface_pts, rcsb_id)

    normal_estimated_pcd = fast_normal_estimation(surface_pts, kdtree_radius, max_nn, tanget_planes_k)

    if verbose:
        o3d.visualization.draw_geometries([normal_estimated_pcd], point_show_normal=True)

    o3d.io.write_point_cloud(output_normals_pcd, normal_estimated_pcd)

    apply_poisson_reconstruction(
        output_normals_pcd,
        output_mesh,
        recon_depth=PR_depth,
        recon_pt_weight=PR_ptweight,
    )

    mesh = pv.read(output_mesh)
    labeled = mesh.connectivity(largest=True)  # This keeps only the largest component
    
    # Save the filtered mesh
    labeled.save(output_mesh)
    
    watertight = validate_mesh_pyvista(labeled)
    if verbose:
        visualize_mesh(output_mesh, rcsb_id)
    if not watertight:
        print("XXXX Watertightness check failed, removing", output_mesh , " XXXX")
        os.remove(output_mesh)
```


And of course NPET pipeline itself relies on a bunch of other stages that are defined in the pipeline folder that parse the ribosome, manipualte pointclouds, masks, clustering, correspond etc. loosely grouped into this stages. Let's try to reuse these stages as they are right now just for the sake of makign progress, but i'm also aware that they are super bloated and can each be basically one function as opposed a fucking file...

In any case, here are some:
```
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Dict, Optional

from ribctl.lib.npet.visualization.pipeline_status_tracker import NPETProcessingTracker, ProcessingStage



class NPETPipelineStage(ABC):
    """
    Base class for all NPET pipeline stages.
    
    Each stage handles a specific part of the mesh creation process,
    with standardized interfaces for input/output and error handling.
    """
    
    def __init__(self, 
                 rcsb_id: str, 
                 tracker: NPETProcessingTracker,
                 artifacts_dir: Path,
                 force: bool = False):
        """
        Initialize a pipeline stage.
        
        Args:
            rcsb_id: The RCSB PDB identifier
            tracker: Processing tracker for logging progress
            artifacts_dir: Directory to store output artifacts
            force: Whether to force regeneration even if artifacts exist
        """
        self.rcsb_id = rcsb_id.upper()
        self.tracker = tracker
        self.artifacts_dir = artifacts_dir
        self.force = force
    
    @property
    @abstractmethod
    def stage(self) -> ProcessingStage:
        """The processing stage this class handles."""
        pass
    
    @property
    @abstractmethod
    def stage_params(self) -> Dict[str, Any]:
        """Parameters used by this stage for tracking and reproducibility."""
        pass
    
    @abstractmethod
    def process(self, context: Dict[str, Any]) -> Dict[str, Any]:
        """
        Execute the processing for this stage.
        
        Args:
            context: Input data from previous stages
            
        Returns:
            Dict with outputs for next stages
        """
        pass
    
    def run(self, context: Dict[str, Any]) -> Dict[str, Any]:
        """
        Run the stage with progress tracking and error handling.
        
        Args:
            context: Input data from previous stages
            
        Returns:
            Dict with outputs for next stages
        """
        should_process = self.tracker.begin_stage(self.stage, self.stage_params)
        
        if not should_process and not self.force:
            print(f"Skipping {self.stage} (artifacts exist and parameters unchanged)")
            return {}
            
        try:
            results = self.process(context)
            self.tracker.end_stage(self.stage, True)
            return results
        except Exception as e:
            self.tracker.end_stage(self.stage, False, error=e)
            # Re-raise to stop pipeline
            raise
```

entity_filtering_stage.py:
```
import numpy as np
from pathlib import Path
from typing import Any, Dict, List

from ribctl.lib.npet.kdtree_approach import (
    ribosome_entities,
    filter_residues_parallel,
)
from ribctl.lib.npet.pipeline.base_stage import NPETPipelineStage
from ribctl.lib.npet.visualization.pipeline_status_tracker import NPETProcessingTracker, ProcessingStage
from ribctl.lib.npet.visualization.various_visualization import visualize_filtered_residues


class EntityFilteringStage(NPETPipelineStage):
    """
    Stage for filtering ribosome entities to only those within the tunnel region.
    """
    
    def __init__(self, 
                 rcsb_id: str, 
                 tracker: NPETProcessingTracker,
                 artifacts_dir: Path,
                 force: bool = False,
                 radius: float = 35,
                 height: float = 120):
        """
        Initialize the entity filtering stage.
        
        Args:
            rcsb_id: The RCSB PDB identifier
            tracker: Processing tracker
            artifacts_dir: Directory for artifacts
            force: Whether to force regeneration
            radius: Radius of the tunnel cylinder
            height: Height of the tunnel cylinder
        """
        super().__init__(rcsb_id, tracker, artifacts_dir, force)
        self.radius = radius
        self.height = height
    
    @property
    def stage(self) -> ProcessingStage:
        return ProcessingStage.ENTITY_FILTERING
    
    @property
    def stage_params(self) -> Dict[str, Any]:
        return {
            "radius": self.radius,
            "height": self.height,
        }
    
    def process(self, context: Dict[str, Any]) -> Dict[str, Any]:
        """
        Filter ribosome entities to only those relevant to the tunnel.
        """
        cifpath = context["cifpath"]
        ptc_pt = context["ptc_pt"]
        constriction_pt = context["constriction_pt"]
        profile = context["profile"]
        
        # Determine tunnel debris chains to exclude
        tunnel_debris = {
            "3J7Z": ["a", "7"],
            "5GAK": ["z"],
            "5NWY": ["s"],
            "7A5G": ["Y2"],
            "9F1D": ["BK"],
        }
        
        # For mitochondrial ribosomes, also exclude mL45 chain
        if profile.mitochondrial:
            try:
                chain = context["ro"].get_poly_by_polyclass("mL45")
                tunnel_debris[self.rcsb_id] = [chain.auth_asym_id]
            except Exception as chain_error:
                print("Mitochondrial mL45 chain not found:", str(chain_error))
        
        # Get all ribosome entities
        residues = ribosome_entities(
            self.rcsb_id,
            cifpath,
            "R",
            tunnel_debris[self.rcsb_id] if self.rcsb_id in tunnel_debris else [],
        )
        
        # Filter to only those within tunnel region
        filtered_residues = filter_residues_parallel(
            residues, ptc_pt, constriction_pt, self.radius, self.height
        )
        
        # Extract atom coordinates
        filtered_points = np.array([
            atom.get_coord()
            for residue in filtered_residues
            for atom in residue.child_list
        ])
        
        # Save filtered points
        filtered_points_path = self.artifacts_dir / f"{self.rcsb_id}_filtered_points.npy"
        np.save(filtered_points_path, filtered_points)
        self.tracker.add_artifact(self.stage, filtered_points_path)
        
        # Save visualization if possible
        try:
            viz_path = self.artifacts_dir / f"{self.rcsb_id}_filtered_residues.png"
            visualize_filtered_residues(
                filtered_residues, residues, ptc_pt, constriction_pt, self.radius, self.height,
                output_path=str(viz_path)
            )
            self.tracker.add_artifact(self.stage, viz_path)
        except Exception as viz_error:
            print(f"Warning: Could not save visualization: {str(viz_error)}")
        
        return {
            "filtered_residues": filtered_residues,
            "filtered_points": filtered_points
        }
```

ptc_identification_stage.py:
```
import numpy as np
from pathlib import Path
from typing import Any, Dict

from ribctl.asset_manager.asset_types import AssetType
from ribctl.lib.landmarks.ptc_via_trna import PTC_location
from ribctl.lib.npet.pipeline.base_stage import NPETPipelineStage
from ribctl.lib.npet.visualization.pipeline_status_tracker import ProcessingStage


class PTCIdentificationStage(NPETPipelineStage):
    """
    Stage for identifying the PTC (Peptidyl Transferase Center) location.
    """
    
    @property
    def stage(self) -> ProcessingStage:
        return ProcessingStage.PTC_IDENTIFICATION
    
    @property
    def stage_params(self) -> Dict[str, Any]:
        # No configurable parameters for this stage
        return {}
    
    def process(self, context: Dict[str, Any]) -> Dict[str, Any]:
        """
        Identify the PTC location based on the ribosome structure.
        """
        try:
            # Get PTC location
            ptc_info = PTC_location(self.rcsb_id)
            ptc_pt = np.array(ptc_info.location)
            
            # Save PTC info as JSON
            ptc_json_path = AssetType.PTC.get_path(self.rcsb_id)
            with open(ptc_json_path, 'w') as f:
                f.write(ptc_info.model_dump_json())
            
            # Track the artifact
            self.tracker.add_artifact(self.stage, ptc_json_path)
            
            return {
                "ptc_info": ptc_info,
                "ptc_pt": ptc_pt
            }
        except Exception as e:
            raise RuntimeError(f"Failed to identify PTC: {str(e)}") from e
```



No matter how these stages are ultimately organized or encoded (whether as separate classes etc. ) conceptually i think the separation is correct because i'd like to be able to basically associate each particular "stage" with an artifact and a way to visualize that artifact, which is lacking right now (or at least works like shit in the currnet orchestrator). For example the mesh producing code at the end should obviosuly display the mesh s.t. it can be interacted with. But also something like "PTC idenfication" stage should have a visualizaer that loads some auxilary data (ex. riboosme ptcloud) and then a claer mkrker fro the location of the ptc that has just been calcualted. Don't worry about this too much for now, but engineere it to be simple and extensible.

Furthermore let's stick to the pyqt and repair what we can about this threading issue (mind you im on macos so it gets finicky with sticking vtk visualizer into the pyqt thread). Furthrmore lets reduce the bloat where we can: no need to be overprotective and create 50-line catch blocks and stuff, no need for elaborrate warnings and sophisticated typechecking and loggign beyond what's there already. If the thing fails on one run -- it's no problem, it's only me who's using it right now nad i'd prefer that the code is lean much more than that the software to be pretty. Everythign should be minimal and functional. No need for crazy code comments or docstrings beyond what is absolutely necessary to understand the code etc.


OH and here the is the assets thingy that specific the tunnel artifact paths that we can rpobably reuse:
tunnel_asset_manager.py:
```
import os

from ribctl import ASSETS_PATH, RIBETL_DATA, RIBXZ_TEMP_FILES
from ribctl.ribosome_ops import RibosomeOps


class TunnelMeshAssetsManager:

    rcsb_id: str

    def __init__(self, rcsb_id: str) -> None:
        self.rcsb_id = rcsb_id.upper()
        RO = RibosomeOps(rcsb_id)
        self.structpath = RO.assets.paths.cif
        pass

    @property
    def cif_struct(self):
        return self.structpath

    @property
    def tunnel_pcd_normal_estimated(self):
        return os.path.join( RIBXZ_TEMP_FILES, "{}_tunnel_pcd_normal_estimated.ply".format(self.rcsb_id) )

    @property
    def tunnel_half_mesh(self):
        return os.path.join(
            self.structpath, "{}_tunnel_half_poisson_recon.ply".format(self.rcsb_id)
        )

    @property
    def ashape_half_mesh(self):
        return os.path.join(
            self.structpath, "{}_half_ashape_watertight.ply".format(self.rcsb_id)
        )

    # These are just to visualize lamps trajectories/pointclouds as pdb files.

    @property
    def lammps_traj_tunnel(self):
        return os.path.join(
            self.structpath, "{}_tunnel.lammpstraj".format(self.rcsb_id)
        )

    @property
    def lammps_traj_ashape(self):
        return os.path.join(
            self.structpath, "{}_ashape.lammpstraj".format(self.rcsb_id)
        )

    @property
    def lammps_traj_tunnel_as_pdb(self):
        return os.path.join(
            self.structpath, "{}_tunnel.lammpstraj.pdb".format(self.rcsb_id)
        )

    @property
    def lammps_traj_ashape_as_pdb(self):
        return os.path.join(
            self.structpath, "{}_ashape.lammpstraj.pdb".format(self.rcsb_id)
        )

    @property
    def dbscan_clusters_PDB_noise(self):
        return os.path.join(
            self.structpath, "{}_dbscan_clusters.noise.pdb".format(self.rcsb_id)
        )

    @property
    def dbscan_clusters_PDB_refined(self):
        return os.path.join(
            self.structpath, "{}_dbscan_clusters.refined.pdb".format(self.rcsb_id)
        )

    @property
    def dbscan_clusters_PDB_largest(self):
        return os.path.join(
            self.structpath, "{}_dbscan_clusters.largest.pdb".format(self.rcsb_id)
        )

    @property
    def dbscan_clusters_mmcif_noise(self):
        return os.path.join(
            self.structpath, "{}_dbscan_clusters.noise.cif".format(self.rcsb_id)
        )

    @property
    def dbscan_clusters_mmcif_refined(self):
        return os.path.join(
            self.structpath, "{}_dbscan_clusters.refined.cif".format(self.rcsb_id)
        )

    @property
    def dbscan_clusters_mmcif_largest(self):
        return os.path.join(
            self.structpath, "{}_dbscan_clusters.largest.cif".format(self.rcsb_id)
        )

    @property
    def dbscan_clusters_xyz_noise(self):
        return os.path.join(
            self.structpath, "{}_dbscan_clusters.noise.xyz".format(self.rcsb_id)
        )

    @property
    def dbscan_clusters_xyz_refined(self):
        return os.path.join(
            self.structpath, "{}_dbscan_clusters.refined.xyz".format(self.rcsb_id)
        )
```


here are some computational methods that are used all over the fucking place by these stages, i gave up keeping track of what's what.. lots of redundant shit and i'll need to clean it up at some point..:
```
from concurrent.futures import ProcessPoolExecutor
import multiprocessing
import os
from typing import Callable, List, Literal, Optional, Tuple, TypeVar
import open3d as o3d
import numpy as np
from scipy.spatial import cKDTree
import sys
from scipy.spatial import cKDTree
import numpy as np
from Bio.PDB.MMCIFParser import FastMMCIFParser
from typing import List, Tuple

from ribctl.lib.schema.types_ribosome import ResidueSummary
data_dir = os.getenv("DATA_DIR")
sys.dont_write_bytecode = True
from Bio.PDB.Entity import Entity
from Bio.PDB.MMCIFParser import FastMMCIFParser
import numpy as np
from Bio.PDB.Chain import Chain
import pyvista as pv
from pathlib import Path
from pprint import pprint
import subprocess
from matplotlib import pyplot as plt
import open3d as o3d
import pyvista as pv
import numpy as np
import plyfile
import warnings
from ribctl import POISSON_RECON_BIN
import numpy as np
from sklearn.cluster import DBSCAN as skDBSCAN
import requests

warnings.filterwarnings("ignore")
import os


def landmark_constriction_site(rcsb_id: str) -> np.ndarray:
    """
    Fetches the constriction site location for a given RCSB PDB ID.

    Args:
        rcsb_id (str): The RCSB PDB identifier (e.g., "4UG0")

    Returns:
        np.ndarray: Array containing the x, y, z coordinates of the constriction site
    """
    url = f"http://localhost:8000/loci/constriction_site?rcsb_id={rcsb_id}"
    headers = {"Accept": "application/json"}

    response = requests.get(url, headers=headers)
    response.raise_for_status()  # Raises an exception for 4XX/5XX status codes

    data = response.json()
    return np.array(data["location"])


def landmark_ptc(rcsb_id: str) -> np.ndarray:
    """
    Fetches the peptidyl transferase center (PTC) location for a given RCSB PDB ID.

    Args:
        rcsb_id (str): The RCSB PDB identifier (e.g., "4UG0")

    Returns:
        np.ndarray: Array containing the x, y, z coordinates of the PTC site
    """
    url = f"http://localhost:8000/loci/ptc?rcsb_id={rcsb_id}"
    headers = {"Accept": "application/json"}

    response = requests.get(url, headers=headers)
    response.raise_for_status()

    data = response.json()
    return np.array(data["location"])

def DBSCAN_capture(
    ptcloud: np.ndarray,
    eps,
    min_samples,
    metric: str = "euclidean",
):

    u_EPSILON     = eps
    u_MIN_SAMPLES = min_samples
    u_METRIC      = metric

    print(
        "Running DBSCAN on {} points. eps={}, min_samples={}, distance_metric={}".format(
            len(ptcloud), u_EPSILON, u_MIN_SAMPLES, u_METRIC
        )
    )
    db = skDBSCAN(eps=eps, min_samples=min_samples, metric=metric).fit(ptcloud)
    labels = db.labels_

    CLUSTERS_CONTAINER = {}
    for point, label in zip(ptcloud, labels):
        if label not in CLUSTERS_CONTAINER:
            CLUSTERS_CONTAINER[label] = []
        CLUSTERS_CONTAINER[label].append(point)

    CLUSTERS_CONTAINER = dict(sorted(CLUSTERS_CONTAINER.items()))
    return db, CLUSTERS_CONTAINER


def DBSCAN_pick_largest_cluster(
    clusters_container: dict[int, list], pick_manually: bool = False
) -> tuple[np.ndarray, int]:
    DBSCAN_CLUSTER_ID = 0
    # if pick_manually:
    #     print("-------------------------------")
    #     print("Running Manual Cluster Selection")
    #     picked_id =  int(input("Enter Cluster ID to proceed the reconstruction with\n (options:[{}]):".format(list( clusters_container.keys() ))))
    #     print("Choise cluster # {}".format(picked_id))
    #     if picked_id == -2:
    #         # if picked -2 ==> return largest
    #         for k, v in clusters_container.items():
    #             if int(k) == -1:
    #                 continue
    #             elif len(v) > len(clusters_container[DBSCAN_CLUSTER_ID]):
    #                 DBSCAN_CLUSTER_ID = int(k)
    #         return np.array(clusters_container[DBSCAN_CLUSTER_ID])

    #     return np.array(clusters_container[picked_id])

    for k, v in clusters_container.items():
        # print(f"Cluster {k} has {len(v)} points.")
        if int(k) == -1:
            continue
        elif len(v) > len(clusters_container[DBSCAN_CLUSTER_ID]):
            DBSCAN_CLUSTER_ID = int(k)
    return np.array(clusters_container[DBSCAN_CLUSTER_ID]), DBSCAN_CLUSTER_ID


def apply_poisson_reconstruction(
    surf_estimated_ptcloud_path: str,
    output_path: Path,
    recon_depth: int = 6,
    recon_pt_weight: int = 3,
):
    # The documentation can be found at https://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version16.04/ in "PoissonRecon" binary
    print( "Rolling Poisson Reconstruction: {} -> {}".format( surf_estimated_ptcloud_path, output_path ) )
    command = [
        POISSON_RECON_BIN,
        "--in",
        surf_estimated_ptcloud_path,
        "--out",
        output_path,
        "--depth",
        str(recon_depth),
        "--pointWeight",
        str(recon_pt_weight),
    ]
    process = subprocess.run(command, capture_output=True, text=True)
    if process.returncode == 0:
        data = plyfile.PlyData.read(output_path)
        data.text = True
        ascii_duplicate = output_path.as_posix().split(".")[0] + "_ascii.ply"
        data.write(ascii_duplicate)
        print(">>Wrote {} and the _ascii version.".format(output_path))
    else:
        print(">>Error:", process.stderr)


def ptcloud_convex_hull_points(
    pointcloud: np.ndarray, ALPHA: float, TOLERANCE: float, OFFSET: float
) -> np.ndarray:
    assert pointcloud is not None
    cloud = pv.PolyData(pointcloud)
    grid = cloud.delaunay_3d(
        alpha=ALPHA, tol=TOLERANCE, offset=OFFSET, progress_bar=True
    )
    convex_hull = grid.extract_surface().cast_to_pointset()
    return convex_hull.points


def estimate_normals(
    convex_hull_surface_pts: np.ndarray,
    kdtree_radius=None,
    kdtree_max_nn=None,
    correction_tangent_planes_n=None,
):
    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(convex_hull_surface_pts)
    pcd.estimate_normals(
        search_param=o3d.geometry.KDTreeSearchParamHybrid(
            radius=kdtree_radius, max_nn=kdtree_max_nn
        )
    )
    pcd.orient_normals_consistent_tangent_plane(k=correction_tangent_planes_n)
    # o3d.visualization.draw_geometries([pcd], point_show_normal=True)
    return pcd


T = TypeVar("T")


def generate_voxel_centers(radius: float, height: float, voxel_size: float) -> tuple:
    nx = ny = int(2 * radius / voxel_size) + 1
    nz = int(height / voxel_size) + 1
    x = np.linspace(-radius, radius, nx)
    y = np.linspace(-radius, radius, ny)
    z = np.linspace(0, height, nz)

    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    voxel_centers = np.column_stack((X.ravel(), Y.ravel(), Z.ravel()))
    return voxel_centers, (X.shape, x, y, z)


def create_point_cloud_mask(
    points: np.ndarray,
    radius: float,
    height: float,
    voxel_size: float = 1.0,
    radius_around_point: float = 2.0,
):

    voxel_centers, (grid_shape, x, y, z) = generate_voxel_centers(
        radius, height, voxel_size
    )
    tree = cKDTree(points)
    indices = tree.query_ball_point(voxel_centers, radius_around_point)

    point_cloud_mask = np.zeros(len(voxel_centers), dtype=bool)
    point_cloud_mask[[i for i, idx in enumerate(indices) if idx]] = True

    point_cloud_mask = point_cloud_mask.reshape(grid_shape)

    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    cylinder_mask = np.sqrt(X**2 + Y**2) <= radius
    hollow_cylinder = ~cylinder_mask

    final_mask = hollow_cylinder | point_cloud_mask
    return final_mask, (x, y, z)


def get_transformation_to_C0(
    base_point: np.ndarray, axis_point: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute transformation matrices to move arbitrary cylinder to C0 configuration.
    Returns translation vector and rotation matrix.
    """
    # Get cylinder axis vector
    axis = axis_point - base_point
    axis_length = np.linalg.norm(axis)
    axis_unit = axis / axis_length

    # Get rotation that aligns axis_unit with [0, 0, 1]
    z_axis = np.array([0, 0, 1])

    # Use Rodrigues rotation formula to find rotation matrix
    # that rotates axis_unit to z_axis
    if np.allclose(axis_unit, z_axis):
        R = np.eye(3)
    elif np.allclose(axis_unit, -z_axis):
        R = np.diag([1, 1, -1])  # 180-degree rotation around x-axis
    else:
        v = np.cross(axis_unit, z_axis)
        s = np.linalg.norm(v)
        c = np.dot(axis_unit, z_axis)
        v_skew = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        R = np.eye(3) + v_skew + (v_skew @ v_skew) * (1 - c) / (s * s)

    return -base_point, R


def transform_points_to_C0(
    points: np.ndarray, base_point: np.ndarray, axis_point: np.ndarray
) -> np.ndarray:
    translation, rotation = get_transformation_to_C0(base_point, axis_point)
    points_translated = points + translation
    points_transformed = points_translated @ rotation.T

    return points_transformed


def transform_points_from_C0(
    points: np.ndarray, base_point: np.ndarray, axis_point: np.ndarray
) -> np.ndarray:

    translation, rotation = get_transformation_to_C0(base_point, axis_point)
    points_unrotated = points @ rotation
    points_untranslated = points_unrotated - translation

    return points_untranslated


def clip_pointcloud_with_mesh(points: np.ndarray, mesh_path: str) -> np.ndarray:
    """
    Clips a point cloud to keep only points that lie inside a mesh using point-by-point checking.

    Parameters:
        points (np.ndarray): Nx3 array of points to clip
        mesh_path (str): Path to the mesh file (PLY format)

    Returns:
        np.ndarray: Points that lie inside the mesh
    """
    # Load the mesh
    mesh = pv.read(mesh_path)

    # Convert to PolyData if needed
    if not isinstance(mesh, pv.PolyData):
        mesh = mesh.extract_surface()

    # Ensure mesh is triangulated
    if not mesh.is_all_triangles:
        mesh = mesh.triangulate()
    else:
        print("Mesh is triangulated")

    # Initialize mask array
    mask = np.zeros(len(points), dtype=bool)
    print("Got mask", mask.shape)

    # Check each point
    for i, point in enumerate(points):
        # Create a single-point PolyData
        point_data = pv.PolyData(point.reshape(1, 3))
        # Check if point is inside
        selection = mesh.select_enclosed_points(point_data, check_surface=False)
        mask[i] = bool(selection.point_data["SelectedPoints"][0])

        # Progress indicator
        if i % 1000 == 0:
            print(f"Processed {i}/{len(points)} points")

    # Return filtered points
    clipped_points = points[mask]
    print(f"Kept {len(clipped_points)}/{len(points)} points")

    return clipped_points


def verify_mesh_quality(mesh) -> dict:
    """
    Verifies the quality of the input mesh and returns diagnostics.
    """
    stats = {
        "n_points": mesh.n_points,
        "n_faces": mesh.n_faces,
        "is_manifold": mesh.is_manifold,
        "bounds": mesh.bounds,
        "open_edges": mesh.n_open_edges,
    }

    try:
        stats["volume"] = mesh.volume
    except:
        stats["volume"] = None
        print("Warning: Could not compute mesh volume")

    return stats


def clip_pcd_via_ashape(
    pcd: np.ndarray, mesh: pv.PolyData
) -> Tuple[np.ndarray, np.ndarray]:
    points_poly = pv.PolyData(pcd)
    select = points_poly.select_enclosed_points(mesh)
    mask = select["SelectedPoints"]
    ashape_interior = pcd[mask == 1]
    ashape_exterior = pcd[mask == 0]
    return ashape_interior, ashape_exterior


def ribosome_entities(
    rcsb_id: str,
    cifpath: str,
    level=Literal["R"] | Literal["A"],
    skip_nascent_chain: List[str] = [],
) -> list[Entity]:
    structure = FastMMCIFParser(QUIET=True).get_structure(rcsb_id, cifpath)
    residues = []
    for chain in structure.child_list[0]:
        if chain.id in skip_nascent_chain:
            print("Skipping nascent chain ", chain.id)
            continue
        residues.extend(chain)

    for residue in residues:
        if not ResidueSummary.is_canonical(residue.get_resname()):
            residues.remove(residue)

    if level == "R":
        return residues
    elif level == "A":
        atoms = []
        [atoms.extend(a) for a in residues]
        return atoms
    else:
        raise


def is_point_in_cylinder(
    point: np.ndarray,
    base_point: np.ndarray,
    axis_point: np.ndarray,
    radius: float,
    height: float,
) -> bool:

    point = np.asarray(point)
    base_point = np.asarray(base_point)
    axis_point = np.asarray(axis_point)

    # Calculate cylinder axis direction vector
    axis = axis_point - base_point
    axis_length = np.linalg.norm(axis)
    axis_unit = axis / axis_length

    # Calculate vector from base to point
    point_vector = point - base_point

    # Project vector onto cylinder axis
    projection = np.dot(point_vector, axis_unit)

    # Calculate perpendicular vector from axis to point
    projection_point = base_point + projection * axis_unit
    radial_vector = point - projection_point

    # Calculate radial distance
    radial_distance = np.linalg.norm(radial_vector)

    # Check if point is inside cylinder
    return (radial_distance <= radius) and (-10 <= projection <= height)


def make_cylinder_predicate(
    base_point: np.ndarray,
    axis_point: np.ndarray,
    radius: float,
    height: float,
    position_getter: Callable[[T], np.ndarray],
) -> Callable[[T], bool]:
    def predicate(obj: T) -> bool:
        position = position_getter(obj)
        return is_point_in_cylinder(position, base_point, axis_point, radius, height)

    return predicate


def get_residue_position(residue):
    return residue.center_of_mass()


def _worker_process_chunk(chunk_data):
    """
    Worker function that processes a chunk of residues.
    Takes a tuple of (residue_positions, base_point, axis_point, radius, height, indices)
    Returns indices of residues that are inside the cylinder.
    """
    positions, base_point, axis_point, radius, height, indices = chunk_data

    results = []
    for i, pos in enumerate(positions):
        if is_point_in_cylinder(pos, base_point, axis_point, radius, height):
            results.append(indices[i])

    return results


def filter_residues_parallel(
    residues: List[T],
    base_point: np.ndarray,
    axis_point: np.ndarray,
    radius: float,
    height: float,
    chunk_size: Optional[int] = None,
    max_workers: Optional[int] = None,
) -> List[T]:
    """
    Filter residues in parallel using ProcessPoolExecutor.

    Parameters:
    -----------
    residues : List[T]
        List of residue objects to filter
    base_point : np.ndarray
        Center point of cylinder base
    axis_point : np.ndarray
        Point defining cylinder axis direction
    radius : float
        Radius of cylinder
    height : float
        Height of cylinder
    chunk_size : Optional[int]
        Size of chunks to process in parallel. If None, calculated automatically
    max_workers : Optional[int]
        Maximum number of worker processes. If None, uses CPU count

    Returns:
    --------
    List[T]
        Filtered list of residues whose positions lie within the cylinder
    """
    # Set defaults for parallel processing parameters
    if max_workers is None:
        max_workers = multiprocessing.cpu_count()

    if chunk_size is None:
        # Aim for ~4 chunks per worker for better load balancing
        chunk_size = max(1, len(residues) // (max_workers * 4))

    # Pre-compute all positions and create index mapping
    positions = np.array([get_residue_position(r) for r in residues])

    indices = list(range(len(residues)))
    index_chunks = [
        indices[i : i + chunk_size] for i in range(0, len(indices), chunk_size)
    ]

    # Create data chunks for processing
    chunks_data = [
        (positions[idx], base_point, axis_point, radius, height, idx)
        for idx in index_chunks
    ]

    # Process chunks in parallel
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        results = list(executor.map(_worker_process_chunk, chunks_data))

    # Flatten results and get corresponding residues
    filtered_indices = [idx for chunk_result in results for idx in chunk_result]
    return [residues[i] for i in filtered_indices]


def clip_tunnel_by_chain_proximity(
    tunnel_points: np.ndarray,
    rcsb_id: str,
    cif_path: str,
    chain_id: str = "Y2",
    n_start_residues: int = 7,
    start_proximity_threshold: float = 10.0,  # Tighter threshold for start residues
    rest_proximity_threshold: float = 15.0,  # Wider threshold for rest of chain
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Clips a tunnel point cloud using different proximity thresholds for the start
    and rest of the chain.

    Args:
        tunnel_points: np.ndarray
            Nx3 array of tunnel point cloud coordinates
        rcsb_id: str
            PDB/MMCIF ID of the structure
        cif_path: str
            Path to the MMCIF file
        chain_id: str
            Chain ID of the nascent peptide/reference chain
        n_start_residues: int
            Number of N-terminal residues to apply tighter threshold to
        start_proximity_threshold: float
            Maximum distance (Å) for points near start residues
        rest_proximity_threshold: float
            Maximum distance (Å) for points near rest of chain

    Returns:
        Tuple[np.ndarray, np.ndarray]:
            - Filtered point cloud (points within thresholds)
            - Removed points (points outside thresholds)
    """
    # Parse structure and extract specified chain
    parser = FastMMCIFParser(QUIET=True)
    structure = parser.get_structure(rcsb_id, cif_path)

    # Get the first model and find the specified chain
    model = structure[0]
    try:
        chain = model[chain_id]
    except KeyError:
        raise ValueError(f"Chain {chain_id} not found in structure {rcsb_id}")

    # Calculate centers of mass for each residue in the chain
    residue_centers = []
    residue_indices = []  # Track residue numbers
    for residue in chain:
        coords = np.array([atom.get_coord() for atom in residue])
        center = coords.mean(axis=0)
        residue_centers.append(center)
        residue_indices.append(residue.id[1])  # Get residue number

    residue_centers = np.array(residue_centers)
    residue_indices = np.array(residue_indices)

    # Split centers into start and rest based on residue numbers
    start_centers = residue_centers[residue_indices <= n_start_residues]
    rest_centers = residue_centers[residue_indices > n_start_residues]

    # Build KD-trees for both parts
    start_tree = cKDTree(start_centers)
    rest_tree = cKDTree(rest_centers) if len(rest_centers) > 0 else None

    # Find distances from each tunnel point to nearest residue center
    start_distances, _ = start_tree.query(tunnel_points)
    if rest_tree is not None:
        rest_distances, _ = rest_tree.query(tunnel_points)
    else:
        rest_distances = np.full_like(start_distances, np.inf)

    # Point is kept if it's within threshold of either part
    within_start = start_distances <= start_proximity_threshold
    within_rest = rest_distances <= rest_proximity_threshold
    keep_mask = within_start | within_rest

    # Split points into kept and removed
    kept_points = tunnel_points[keep_mask]
    removed_points = tunnel_points[~keep_mask]

    print(f"Kept {len(kept_points)} points:")
    print(
        f"- {np.sum(within_start)} points within {start_proximity_threshold}Å of first {n_start_residues} residues"
    )
    print(
        f"- {np.sum(within_rest)} points within {rest_proximity_threshold}Å of remaining residues"
    )
    print(f"Removed {len(removed_points)} points")

    return kept_points, removed_points

```