Alrighty, i have an NPET extraction pipeline in my riboxyz codebase that i want to first extract (limit its integration with the actual riboxyz application to only a few key points/coordinates) and then heavily refactor and optimize (bubble up the main computational parameters, implement fine resolution and more flexbile parametrization of operations where appropriate, gpu operations where possible).

For that thought i first want to reorganize this fucking mess. It grew somewhat adhoc.. well, very adhoc when i was building it and currently relies on a mix of a bunch of different disjointed packages and a flow of data that's not robust to a few things that happen pretty often (ex. if the exterior mesh of the ribosome is not watertight then we can't use it to clip the interior mesh, which is pretty important). There is not a shred of formal reasoning about parameters we pick for various computational tools during this work (alpha shape, poisson reconstruction, dbscan etc.) so to succeed even half the time everything has to be hand-picked "just so" and my current paramset is a result of heavy trial-and-error.

Also i'm not super duper satisfied with the way the pipeline produces a thousand different artifacts, or rather the way they are organized. I frequently need to visualize one or a combination of them but these visualizations have been hitherto interspersed into the pipeline code at different stages to basically "have a glance" at this or that stage and go ahead.. I think this is also bad design and unsustainable -- we should have a registry of artifacts PER a given structure (per parameters with which they were produces) and separate librray of visualizations that knows how to access this registry and visualize each thing correctly. I made some efforts to this end, but it's still pretty garbage.

So that's basically my conundrum right now and i need your help to make this a publishable and usable piece of software as opposed to a bloated and convoluted buncha scripts...

For later -- i'm saying this for context, but basically when we clean this shit up a little what i mean to do is 
- to make things like PTC and Constriction site (coordinates vital to our construction) be ingestable either from the local files or from the riboxyz api
- add a configurable atom radius representation -- right now all atoms are uniformly represented with 2A radius sphere
- add an option to create the grid at other resolution than 1 angstrom wide (i'm not even sure if that's what we use currently...). But basically my hunch is that if the voxel grid is made finer we can actually be more robust to various extraneous clusters and have better separation, have better resolutin in the end. Actually this would be a huge and central improvement but will also be incredibly costly for any finer grid so we should reason very carefully and intentianally about this. This would also of course influence our choice of dbscan parameters and i'd like to have at least some notion of formal/back-of-the-envelope correlation between grid size and the dbscan/poisson reconstruction parameters that i could pass on to the users.
- with all of the above in mind -- i'd like to enable this thing to run on gpu and also parallellize whatever can be parallellized... I'm happy to try to move away from the kdtree to octree or whatever, use other langugaes than python even depending on the tradeoffs..


Last two points i'd be very excited to hear your thoughts on. Cheers.

Ok here is the repo and _some_ of the code...


```
(venv) ᢹ saeta.rtviii[ dev/riboxyz ]  tree -L 6 -I 'node_modules|venv|__pycache__|profiles|cache|debug_output|*.npy|*.ply|*.fasta|*.csv|assets_*|staticfiles|api|assets|*.png|TUBETL_DATA|*.pkl|*hmm|*fasta|npet|*.mdx|*.ts.map|*.d.ts|nightingale' .  e
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
│   │   │   ├── bdist.macosx-10.9-universal2
│   │   │   └── lib
│   │   │       └── ribxz_chimerax
│   │   │           ├── __init__.py
│   │   │           ├── cmd_loci.py
│   │   │           ├── cmd_polymers.py
│   │   │           ├── cmd_registry.py
│   │   │           └── io.py
│   │   ├── bundle_info.xml
│   │   ├── dist
│   │   │   └── ribxz_chimerax-0.1-py3-none-any.whl
│   │   ├── ribxz_chimerax.egg-info
│   │   │   ├── dependency_links.txt
│   │   │   ├── PKG-INFO
│   │   │   ├── requires.txt
│   │   │   ├── SOURCES.txt
│   │   │   └── top_level.txt
│   │   └── src
│   │       ├── __init__.py
│   │       ├── cmd_loci.py
│   │       ├── cmd_polymers.py
│   │       ├── cmd_registry.py
│   │       └── io.py
│   ├── stats.json
│   ├── uniprot_seeds_query_record.py
│   └── williamson_assembly.py
├── dirtocontext.py
├── docker-compose.yml
├── Dockerfile-django
├── docs.md
├── muscle3.8.31_src.tar.gz
├── ncbi_taxonomy.sqlite
├── neo4j_ribosome
│   ├── __archive
│   │   ├── cypher_ops
│   │   │   ├── cypher_exec
│   │   │   ├── neo4j_commit_structure.sh
│   │   │   └── neo4j_seed_db_ontology.sh
│   │   └── riboxyz_seed_data
│   │       ├── ban-pfam-map-lsu.json
│   │       ├── ban-pfam-map-ssu.json
│   │       ├── interpro-base.json
│   │       ├── interpro-go-1.json
│   │       ├── interpro-go-2.json
│   │       ├── interpro-go-3.json
│   │       ├── interpro-go-4.json
│   │       ├── package-lock.json
│   │       ├── package.json
│   │       ├── pfam-interpro-1.json
│   │       ├── pfam-interpro-2.json
│   │       ├── pfam-interpro-3.json
│   │       ├── pfam-interpro-4.json
│   │       └── pfam-to-interpro-map.json
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
├── npet_orchestrator.py
├── pipeline_manager.py
├── PLAN_refactor_npet_pipeline.md
├── q.md
├── ribctl
│   ├── __init__.py
│   ├── asset_manager
│   │   ├── asset_manager.py
│   │   ├── asset_registry.py
│   │   ├── asset_types.py
│   │   ├── doc.md
│   │   └── parallel_acquisition.py
│   ├── classifyre.json
│   ├── etl
│   │   ├── __init__.py
│   │   ├── etl_collector.py
│   │   └── gql_querystrings.py
│   ├── global_ops.py
│   ├── lib
│   │   ├── __libseq.py
│   │   ├── chimerax
│   │   │   ├── _cmd_ribrepr.py
│   │   │   ├── cmd_chainsplitter.py
│   │   │   ├── cmd_ligvis.py
│   │   │   ├── cmd_ribetl.py
│   │   │   ├── cmd_ribmovie.py
│   │   │   ├── cmd_ribrepr.py
│   │   │   ├── cmds_all.py
│   │   │   ├── ffmpeg_convert.sh
│   │   │   ├── ffmpeg_firstframe.sh
│   │   │   ├── gen_movies.py
│   │   │   ├── loop_ligvis.py
│   │   │   ├── loop_movies.py
│   │   │   ├── loop_split_chains.py
│   │   │   ├── notes.md
│   │   │   ├── produce_gif.py
│   │   │   └── thumbnails_from_ribetl.sh
│   │   ├── enumunion.py
│   │   ├── info.py
│   │   ├── landmarks
│   │   │   ├── __ptc_via_doris.py
│   │   │   ├── constriction_site.py
│   │   │   ├── notes.md
│   │   │   ├── ptc_via_trna.py
│   │   │   └── rrna_helices
│   │   │       ├── convert.py
│   │   │       ├── ecoli_7K00.json
│   │   │       ├── rrna_helices.py
│   │   │       ├── thermus_1VY4.json
│   │   │       └── yeast_4V88.json
│   │   ├── libbsite.py
│   │   ├── libhmm.py
│   │   ├── libmsa.py
│   │   ├── libseq.py
│   │   ├── libtax.py
│   │   ├── nsearch_gemmi.py
│   │   ├── ribosome_types
│   │   ├── schema
│   │   │   ├── __init__.py
│   │   │   ├── primitives.py
│   │   │   ├── types_binding_site.py
│   │   │   └── types_ribosome.py
│   │   ├── seq_project_many_to_one.py
│   │   ├── thumbnail.py
│   │   ├── types
│   │   │   └── polymer
│   │   │       ├── __init__.py
│   │   │       ├── base.py
│   │   │       ├── hierarchies.py
│   │   │       └── types.py
│   │   └── utils.py
│   ├── logger_config.py
│   ├── logs
│   │   ├── etl.log
│   │   └── loggers.py
│   ├── ribd.py
│   └── ribosome_ops.py
└── taxdump.tar.gz
e  [error opening dir]

29 directories, 146 files

```

ribctl/lib/npet/kdtree_approach.py
```py
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
    print(
        "Rolling Poisson Reconstruction: {} -> {}".format(
            surf_estimated_ptcloud_path, output_path
        )
    )
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

ribctl/lib/npet/npet_pipeline.py
```py
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
from ribctl.lib.npet.pipeline_status_tracker import (
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

ribctl/lib/npet/landmark_ptcloud_collector.py
```py
#!/usr/bin/env python
import os
import sys
import glob
import json
import argparse
import numpy as np
from pathlib import Path
import multiprocessing as mp
import scipy.spatial.distance as distance
from functools import partial

from Bio.PDB.Chain import Chain
from ribctl.lib.types.polymer.base import CytosolicProteinClass, MitochondrialProteinClass
from ribctl.lib.utils import find_closest_pair_two_sets, midpoint
from ribctl.ribosome_ops import RibosomeOps
from ribctl.lib.landmarks.ptc_via_trna import PTC_location


def get_constriction(rcsb_id: str) -> np.ndarray:
    """Calculate constriction site based on uL4 and uL22 protein chains"""
    ro = RibosomeOps(rcsb_id)
    is_mitochondrial = ro.profile.mitochondrial
    if is_mitochondrial:
        uL4 = ro.get_poly_by_polyclass(MitochondrialProteinClass.uL4m)
        uL22 = ro.get_poly_by_polyclass(MitochondrialProteinClass.uL22m)
    else:
        uL4 = ro.get_poly_by_polyclass(CytosolicProteinClass.uL4)
        uL22 = ro.get_poly_by_polyclass(CytosolicProteinClass.uL22)
    if uL4 is None or uL22 is None:
        raise ValueError(f"Could not find uL4 or uL22 in {rcsb_id}")
    structure = ro.assets.biopython_structure()
    
    uL4_c: Chain = structure[0][uL4.auth_asym_id]
    uL22_c: Chain = structure[0][uL22.auth_asym_id]
    uL4_coords = [r.center_of_mass() for r in uL4_c.child_list]
    uL22_coords = [r_.center_of_mass() for r_ in uL22_c.child_list]
    return midpoint(*find_closest_pair_two_sets(uL4_coords, uL22_coords))


def get_closest_residues(residues, target_location, n=15):
    """Get the n residues with centers closest to the target location"""
    residue_centers = []
    
    for residue in residues:
        try:
            # Calculate center of residue
            coords = np.array([atom.coord for atom in residue.get_atoms()])
            if len(coords) > 0:
                center = coords.mean(axis=0)
                residue_centers.append(center)
        except Exception:
            continue
    
    if not residue_centers:
        return np.array([])
    
    # Calculate distances to target location
    residue_centers = np.array(residue_centers)
    distances = distance.cdist([target_location], residue_centers)[0]
    
    # Get indices of n closest residues
    closest_indices = np.argsort(distances)[:n]
    
    # Return the coordinates of the centers of the n closest residues
    return residue_centers[closest_indices]


def process_structure(rcsb_id, output_dir=None, skip_existing=True):
    """Process a single structure and extract landmark data with zero transformations"""
    # Check if output file already exists (if output_dir and skip_existing are provided)
    if output_dir and skip_existing:
        output_file = os.path.join(output_dir, f"{rcsb_id}_landmarks_pts.json")
        if os.path.exists(output_file):
            print(f"Skipping {rcsb_id}, file already exists")
            return None
    
    print(f"Processing {rcsb_id}...")
    
    landmarks = {
        "ptc": None,
        "constriction": None,
        "uL4": None,
        "uL22": None
    }
    
    try:
        # Initialize RibosomeOps
        ro = RibosomeOps(rcsb_id)
        
        # Compute PTC location using the function from ptc_via_trna
        try:
            ptc_info = PTC_location(rcsb_id)
            landmarks["ptc"] = ptc_info.location
            print(f"  Found PTC at {landmarks['ptc']}")
        except Exception as e:
            print(f"  Error finding PTC: {e}")
            return rcsb_id, landmarks
        
        # Compute constriction site location using the provided function
        try:
            constriction_location = get_constriction(rcsb_id)
            landmarks["constriction"] = constriction_location.tolist()
            print(f"  Found constriction site at {landmarks['constriction']}")
        except Exception as e:
            print(f"  Error finding constriction site: {e}")
            return rcsb_id, landmarks
        
        # Find uL4 and uL22 protein chains
        profile = ro.profile
        ul4_chain = None
        ul22_chain = None
        

        if ro.profile.mitochondrial:
            ul4_chain = ro.get_poly_by_polyclass(MitochondrialProteinClass.uL4m).auth_asym_id
            ul22_chain = ro.get_poly_by_polyclass(MitochondrialProteinClass.uL22m).auth_asym_id
        else:
            ul4_chain = ro.get_poly_by_polyclass(CytosolicProteinClass.uL4).auth_asym_id
            ul22_chain = ro.get_poly_by_polyclass(CytosolicProteinClass.uL22).auth_asym_id

        # Process protein chains if found
        model = ro.assets.biopython_structure()[0]
        constriction_location = np.array(landmarks["constriction"])
        
        # Process uL4 chain
        if ul4_chain and ul4_chain in model.child_dict:
            ul4_residues = list(model[ul4_chain].get_residues())
            ul4_points = get_closest_residues(ul4_residues, constriction_location, 15)
            
            if len(ul4_points) > 0:
                landmarks["uL4"] = ul4_points.tolist()
            else:
                print(f"  Warning: No valid uL4 residues found for {rcsb_id}")
        else:
            print(f"  Warning: No uL4 chain found for {rcsb_id}")
        
        # Process uL22 chain
        if ul22_chain and ul22_chain in model.child_dict:
            ul22_residues = list(model[ul22_chain].get_residues())
            ul22_points = get_closest_residues(ul22_residues, constriction_location, 15)
            
            if len(ul22_points) > 0:
                landmarks["uL22"] = ul22_points.tolist()
            else:
                print(f"  Warning: No valid uL22 residues found for {rcsb_id}")
        else:
            print(f"  Warning: No uL22 chain found for {rcsb_id}")
    
    except Exception as e:
        print(f"  Error processing {rcsb_id}: {str(e)}")
    
    # Save to file if output_dir is provided
    if output_dir:
        output_file = os.path.join(output_dir, f"{rcsb_id}_landmarks_pts.json")
        with open(output_file, 'w') as f:
            json.dump(landmarks, f, indent=2)
        print(f"Saved landmark data for {rcsb_id}")
    
    return rcsb_id, landmarks


def process_worker(rcsb_id, output_dir, skip_existing):
    """Worker function for multiprocessing"""
    try:
        result = process_structure(rcsb_id, output_dir, skip_existing)
        if result:
            return result[0]  # Return rcsb_id if processed successfully
        return None  # Return None if skipped
    except Exception as e:
        print(f"Error in worker processing {rcsb_id}: {str(e)}")
        return None


def main():
    parser = argparse.ArgumentParser(description='Extract landmark data from ribosome structures')
    parser.add_argument('--data-dir', type=str, help='Directory containing structure data')
    parser.add_argument('--output-dir', type=str, help='Directory to save output files')
    parser.add_argument('--ids', type=str, nargs='+', help='Specific structure IDs to process')
    parser.add_argument('--processes', type=int, default=mp.cpu_count(), 
                      help='Number of processes to use (default: number of CPU cores)')
    parser.add_argument('--force', action='store_true', 
                      help='Force processing even if output file exists')
    args = parser.parse_args()
    
    # Set default data directory if not provided
    data_dir = args.data_dir if args.data_dir else os.environ.get("RIBETL_DATA", ".")
    
    # Set default output directory if not provided
    output_dir = args.output_dir if args.output_dir else os.path.join(data_dir, "landmarks")
    os.makedirs(output_dir, exist_ok=True)
    
    # Get structure IDs
    if args.ids:
        structure_ids = args.ids
    else:
        # Find all structure directories that have NPET meshes
        mesh_files = glob.glob(os.path.join(data_dir, "*/*_NPET_MESH.ply"))
        structure_ids = [os.path.basename(f).split('_')[0] for f in mesh_files]
        structure_ids = list(set(structure_ids))  # Remove duplicates
    
    print(f"Found {len(structure_ids)} structures to process")
    
    # Adjust number of processes if needed
    num_processes = min(args.processes, len(structure_ids))
    if num_processes < 1:
        num_processes = 1
    
    # Skip check based on force flag
    skip_existing = not args.force
    
    print(f"Using {num_processes} processes")
    
    # Process structures in parallel
    if num_processes > 1:
        # Create a partial function with fixed arguments
        worker_func = partial(process_worker, output_dir=output_dir, skip_existing=skip_existing)
        
        # Create pool and map worker function to structure IDs
        with mp.Pool(processes=num_processes) as pool:
            processed_ids = list(pool.map(worker_func, structure_ids))
        
        # Count successfully processed structures (non-None results)
        successful = [pid for pid in processed_ids if pid is not None]
        print(f"Successfully processed {len(successful)} structures")
    else:
        # Process sequentially if only one process requested
        for rcsb_id in structure_ids:
            process_structure(rcsb_id, output_dir, skip_existing)
    
    print("Done!")


if __name__ == "__main__":
    main()

```

ribctl/lib/npet/pipeline_status_tracker.py
```py
import copy
import json
import os
import time
import traceback
from enum import Enum
from pathlib import Path
from typing import Dict, Any, Optional, List

class ProcessingStage(str, Enum):
    """Enum defining the stages of the NPET mesh pipeline"""
    SETUP = "setup"

    PTC_IDENTIFICATION = "ptc_identification"           # New stage
    CONSTRICTION_IDENTIFICATION = "constriction_identification"  # New stage
    ALPHA_SHAPE = "alpha_shape"
    LANDMARK_IDENTIFICATION = "landmark_identification"
    ENTITY_FILTERING = "entity_filtering"
    POINT_CLOUD_PROCESSING = "point_cloud_processing"
    CLUSTERING = "clustering"
    REFINEMENT = "refinement"
    SURFACE_EXTRACTION = "surface_extraction"
    NORMAL_ESTIMATION = "normal_estimation"
    MESH_RECONSTRUCTION = "mesh_reconstruction"
    VALIDATION = "validation"
    COMPLETE = "complete"

class ProcessingStatus(str, Enum):
    """Status of a processing stage"""
    PENDING     = "pending"
    IN_PROGRESS = "in_progress"
    SUCCESS = "success"
    FAILURE = "failure"
    SKIPPED = "skipped"

class NPETProcessingTracker:
    """Tracks the processing status of NPET mesh generation for a structure"""
    
    def __init__(self, rcsb_id: str, output_dir: Optional[Path] = None):
        self.rcsb_id = rcsb_id.upper()
        self.start_time = time.time()
        self.end_time: Optional[float] = None
        self.current_stage: ProcessingStage = ProcessingStage.SETUP
        self.output_dir = output_dir or Path(f"logs")
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize stage tracking
        self.stages: Dict[ProcessingStage, Dict[str, Any]] = {
            stage: {
                "status": ProcessingStatus.PENDING,
                "start_time": None,
                "end_time": None,
                "duration": None,
                "error": None,
                "artifacts": [],
                "parameters": {}  # New field to store parameters
            } for stage in ProcessingStage
        }
        
        # Initialize overall status
        self.status: Dict[str, Any] = {
            "rcsb_id": self.rcsb_id,
            "overall_status": ProcessingStatus.PENDING,
            "duration": None,
            "stages": self.stages,
            "summary": {
                "success": False,
                "failed_stage": None,
                "error_summary": None,
                "watertight": False,
                "artifacts_generated": []
            }
        }
        
        # Save initial status
        self._save_status()
    
    def begin_stage(self, stage: ProcessingStage, parameters: Optional[Dict[str, Any]] = None) -> bool:
        """
        Mark the beginning of a processing stage.
        
        Args:
            stage: The processing stage
            parameters: Optional dictionary of parameters for this stage
            
        Returns:
            bool: True if the stage should be processed, False if it can be skipped
        """
        self.current_stage = stage
        now = time.time()
        
        # Print stage start
        print(f"\n⏳ Starting stage: {stage}")
        
        # Check if we can skip processing based on parameters
        can_skip = False
        if parameters is not None:
            # Get previous parameters and artifacts
            prev_params = self.stages[stage].get("parameters", {})
            prev_artifacts = self.stages[stage].get("artifacts", [])
            
            # Check if parameters match and artifacts exist
            if prev_params and prev_artifacts and prev_params == parameters:
                # Check that artifacts actually exist on disk
                artifacts_exist = all(os.path.exists(artifact) for artifact in prev_artifacts)
                if artifacts_exist:
                    can_skip = True
                    print(f"⏩ Stage {stage} skipped (parameters unchanged, artifacts exist)")
                    
                    # Update status but keep existing data
                    self.stages[stage].update({
                        "status": ProcessingStatus.SUCCESS,
                        "start_time": now,
                        "end_time": now,
                        "duration": 0.0,  # Zero duration for skipped stages
                    })
                    self._save_status()
                    return False  # Skip processing
        
        # Proceed with processing
        self.stages[stage].update({
            "status": ProcessingStatus.IN_PROGRESS,
            "start_time": now,
            "parameters": parameters or {}  # Store new parameters
        })
        self._save_status()
        return True  # Process the stage
    
    def end_stage(self, stage: ProcessingStage, success: bool, 
                artifacts: Optional[List[Path]] = None,
                error: Optional[Exception] = None) -> None:
        """Mark the end of a processing stage"""
        now = time.time()
        start_time = self.stages[stage]["start_time"] or now
        
        status = ProcessingStatus.SUCCESS if success else ProcessingStatus.FAILURE
        
        error_info = None
        if error:
            error_info = {
                "type": type(error).__name__,
                "message": str(error),
                "traceback": traceback.format_exc()
            }
        
        # Create a dict with the updates, but don't include artifacts unless provided
        update_dict = {
            "status": status,
            "end_time": now,
            "duration": now - start_time,
            "error": error_info,
        }
        
        # Only update the artifacts list if explicitly provided
        if artifacts is not None:
            update_dict["artifacts"] = [str(a) for a in artifacts]
        
        # Apply the updates
        self.stages[stage].update(update_dict)
        
        if not success:
            # Update summary with failure info
            self.status["summary"]["success"] = False
            self.status["summary"]["failed_stage"] = stage
            self.status["summary"]["error_summary"] = str(error) if error else "Unknown error"
            
            # Mark remaining stages as skipped
            all_stages = list(ProcessingStage)
            current_index = all_stages.index(stage)
            for skip_stage in all_stages[current_index+1:]:
                if skip_stage != ProcessingStage.COMPLETE:
                    self.stages[skip_stage]["status"] = ProcessingStatus.SKIPPED
                    
            # Print failure message with emoji
            print(f"❌ Stage FAILED: {stage} - {str(error) if error else 'Unknown error'}")
        else:
            # Print success message with emoji
            duration = now - start_time
            print(f"✅ Stage completed: {stage} ({duration:.2f}s)")
        
        # Print artifact count
        artifact_count = len(self.stages[stage].get("artifacts", []))
        if artifact_count > 0:
            print(f"   {artifact_count} artifacts generated")
        
        self._save_status()
    
    def complete_processing(self, success: bool, watertight: bool = False) -> None:
        """Mark the overall processing as complete"""
        now = time.time()
        self.end_time = now
        total_duration = now - self.start_time
        
        # Get all generated artifacts with debugging
        artifacts = []
        for stage_name, stage_info in self.stages.items():
            stage_artifacts = stage_info.get("artifacts", [])
            if stage_artifacts:
                print(f"Found {len(stage_artifacts)} artifacts in stage {stage_name}: {stage_artifacts}")
                artifacts.extend(stage_artifacts)
        
        print(f"Total artifacts collected: {len(artifacts)}")
        
        self.status.update({
            "overall_status": ProcessingStatus.SUCCESS if success else ProcessingStatus.FAILURE,
            "duration": total_duration,
            "summary": {
                "success": success,
                "watertight": watertight,
                "artifacts_generated": artifacts
            }
        })
        
        # Update the complete stage
        self.stages[ProcessingStage.COMPLETE].update({
            "status": ProcessingStatus.SUCCESS if success else ProcessingStatus.FAILURE,
            "duration": 0
        })
        
        # Don't print summary here - let the caller handle it
        # self.print_summary()
        
        self._save_status()
    
    def add_artifact(self, stage: ProcessingStage, artifact_path: Path) -> None:
        """Add an artifact to the specified stage"""
        artifact_str = str(artifact_path)
        
        if artifact_str not in self.stages[stage].get("artifacts", []):
            self.stages[stage].setdefault("artifacts", []).append(artifact_str)
            self._save_status()
    
    def print_summary(self) -> None:
        """Print a summary of all stages with their status"""
        print("\n" + "=" * 60)
        print(f"PROCESSING SUMMARY FOR {self.rcsb_id}")
        print("=" * 60)
        
        # Calculate the widest stage name for nice formatting
        stage_width = max(len(str(stage)) for stage in ProcessingStage) + 2
        
        # Print status for each stage
        for stage in ProcessingStage:
            if stage != ProcessingStage.COMPLETE:
                status = self.stages[stage]["status"]
                
                # Select appropriate emoji and color based on status
                emoji = {
                    ProcessingStatus.SUCCESS: "✅",
                    ProcessingStatus.FAILURE: "❌",
                    ProcessingStatus.SKIPPED: "⏩",
                    ProcessingStatus.IN_PROGRESS: "⏳",
                    ProcessingStatus.PENDING: "⏱️",
                }.get(status, "❓")
                
                # Get duration if available
                duration_str = ""
                if self.stages[stage]["duration"] is not None:
                    duration = self.stages[stage]["duration"]
                    if duration < 0.001:  # Skipped stages
                        duration_str = "(skipped)"
                    else:
                        duration_str = f"({duration:.2f}s)"
                
                # Get artifacts if available
                artifact_count = len(self.stages[stage].get("artifacts", []))
                artifact_str = f"{artifact_count} artifacts" if artifact_count > 0 else ""
                
                # Get error if available
                error_str = ""
                if self.stages[stage]["error"]:
                    error_str = f"ERROR: {self.stages[stage]['error']['message']}"
                
                # Format the line
                stage_name = f"{stage}".ljust(stage_width)
                status_name = f"{status.name}".ljust(10)
                duration_str = duration_str.ljust(15)
                artifact_str = artifact_str.ljust(15)
                
                print(f"{emoji} {stage_name} {status_name} {duration_str} {artifact_str} {error_str}")
        
        print("-" * 60)
        
        # Print overall status
        total_duration = self.status.get("duration", 0)
        duration_str = f"Total time: {total_duration:.2f}s" if total_duration else ""
        
        success = self.status["summary"]["success"]
        if success:
            print(f"✅ PIPELINE COMPLETED SUCCESSFULLY  {duration_str}")
            if self.status["summary"]["watertight"]:
                mesh_files = [a for a in self.status["summary"]["artifacts_generated"] 
                            if a.endswith('.ply') and not a.endswith('_normal_estimated_pcd.ply')]
                if mesh_files:
                    print(f"   Watertight mesh: {os.path.basename(mesh_files[0])}")
        else:
            failed_stage = self.status["summary"]["failed_stage"]
            error = self.status["summary"]["error_summary"]
            print(f"❌ PIPELINE FAILED at stage {failed_stage}  {duration_str}")
            if error:
                print(f"   Error: {error}")
        
        print("=" * 60)
    
    def _save_status(self) -> None:
        """Save the current status to a JSON file"""
        log_file = self.output_dir / f"{self.rcsb_id}_processing_log.json"
        
        # Create a simplified version of the status for JSON output
        output_status = copy.deepcopy(self.status)
        
        # Simplify the time representation - remove start_time and end_time
        for stage_name, stage_info in output_status["stages"].items():
            if "start_time" in stage_info:
                del stage_info["start_time"]
            if "end_time" in stage_info:
                del stage_info["end_time"]
        
        with open(log_file, 'w') as f:
            json.dump(output_status, f, indent=2, default=str)
```

ribctl/lib/npet/tunnel_asset_manager.py
```py
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

ribctl/lib/npet/pipeline/base_stage.py
```py
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Dict, Optional

from ribctl.lib.npet.pipeline_status_tracker import NPETProcessingTracker, ProcessingStage



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

ribctl/lib/npet/pipeline/clustering_stage.py
```py
import numpy as np
from pathlib import Path
from typing import Any, Dict

from ribctl.lib.npet.kdtree_approach import (
    DBSCAN_capture,
    DBSCAN_pick_largest_cluster,
)
from ribctl.lib.npet.pipeline.base_stage import NPETPipelineStage
from ribctl.lib.npet.pipeline_status_tracker import (
    NPETProcessingTracker,
    ProcessingStage,
)
from ribctl.lib.npet.various_visualization import (
    visualize_DBSCAN_CLUSTERS_particular_eps_minnbrs,
)


def generate_clusters_pdb(
    clusters_container: dict, output_path: str, rcsb_id: str
) -> str:
    """
    Generate a PDB file with each DBSCAN cluster as a separate chain.
    """

    def get_chain_id(cluster_idx: int) -> str:
        """Generate unique chain IDs: A, B, C... Z, AA, AB, AC..."""
        if cluster_idx < 26:
            return chr(65 + cluster_idx)  # A-Z
        else:
            first = chr(65 + (cluster_idx - 26) // 26)
            second = chr(65 + (cluster_idx - 26) % 26)
            return first + second

    with open(output_path, "w") as pdb_file:
        # Write PDB header
        pdb_file.write(f"HEADER    DBSCAN CLUSTERS                         {rcsb_id}\n")
        pdb_file.write(f"TITLE     DBSCAN CLUSTERS FOR {rcsb_id}\n")
        pdb_file.write("REMARK    Generated from DBSCAN clustering\n")

        atom_num = 1
        valid_cluster_idx = 0
        connect_records = []  # Store CONECT records

        # Process each cluster (skip noise cluster -1)
        for cluster_id, points in clusters_container.items():
            if cluster_id == -1:  # Skip noise
                continue

            chain_id = get_chain_id(valid_cluster_idx)
            points_array = np.array(points)

            print(
                f"Writing cluster {cluster_id} as chain {chain_id} with {len(points_array)} points"
            )

            cluster_start_atom = atom_num

            # Write each point as a CA atom in the SAME residue (residue 1)
            for point_idx, (x, y, z) in enumerate(points_array):
                pdb_file.write(
                    f"ATOM  {atom_num:5d}  CA  ALA {chain_id}   1    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C  \n"
                )

                # Create bonds between consecutive atoms in the cluster
                if point_idx > 0:
                    connect_records.append(f"CONECT{atom_num-1:5d}{atom_num:5d}\n")

                atom_num += 1

            # Write chain terminator
            pdb_file.write(f"TER   {atom_num:5d}      ALA {chain_id}   1\n")
            atom_num += 1
            valid_cluster_idx += 1

        # Write all CONECT records at the end
        for connect in connect_records:
            pdb_file.write(connect)

        # Write end record
        pdb_file.write("END\n")

    print(f"Generated clusters PDB: {output_path}")
    print(f"Total clusters written: {valid_cluster_idx}")
    return output_path


class ClusteringStage(NPETPipelineStage):
    """
    Stage for applying DBSCAN clustering to identify the tunnel points.
    """

    def __init__(
        self,
        rcsb_id: str,
        tracker: NPETProcessingTracker,
        artifacts_dir: Path,
        force: bool = False,
        epsilon: float = 5.5,
        min_samples: int = 600,
    ):
        """
        Initialize the clustering stage.

        Args:
            rcsb_id: The RCSB PDB identifier
            tracker: Processing tracker
            artifacts_dir: Directory for artifacts
            force: Whether to force regeneration
            epsilon: DBSCAN epsilon parameter (neighborhood radius)
            min_samples: DBSCAN min_samples parameter (min points in neighborhood)
        """
        super().__init__(rcsb_id, tracker, artifacts_dir, force)
        self.epsilon = epsilon
        self.min_samples = min_samples

    @property
    def stage(self) -> ProcessingStage:
        return ProcessingStage.CLUSTERING

    @property
    def stage_params(self) -> Dict[str, Any]:
        return {
            "epsilon": self.epsilon,
            "min_samples": self.min_samples,
        }

    def process(self, context: Dict[str, Any]) -> Dict[str, Any]:
        """
        Apply DBSCAN clustering to identify the tunnel points.

        Steps:
        1. Apply DBSCAN to the interior points
        2. Extract the largest cluster (representing the tunnel)
        3. Save the cluster data and visualization
        """
        empty_in_world_coords = context["empty_in_world_coords"]
        ptc_pt = context["ptc_pt"]
        constriction_pt = context["constriction_pt"]

        # Apply DBSCAN clustering
        db, clusters_container = DBSCAN_capture(
            empty_in_world_coords, self.epsilon, self.min_samples
        )

        # Extract the largest cluster
        largest_cluster, largest_cluster_id = DBSCAN_pick_largest_cluster(
            clusters_container
        )

        # Save largest cluster
        largest_cluster_path = (
            self.artifacts_dir / f"{self.rcsb_id}_largest_cluster.npy"
        )
        np.save(largest_cluster_path, largest_cluster)
        self.tracker.add_artifact(self.stage, largest_cluster_path)

        # Save visualization if possible
        try:
            cluster_viz_path = self.artifacts_dir / f"{self.rcsb_id}_clusters.png"
            visualize_DBSCAN_CLUSTERS_particular_eps_minnbrs(
                clusters_container,
                self.epsilon,
                self.min_samples,
                ptc_pt,
                constriction_pt,
                largest_cluster,
                context.get("radius", 35),  # Default radius if not in context
                context.get("height", 120),  # Default height if not in context
                output_path=str(cluster_viz_path),
            )
            self.tracker.add_artifact(self.stage, cluster_viz_path)
        except Exception as viz_error:
            print(f"Warning: Could not save cluster visualization: {str(viz_error)}")

        # Save individual clusters
        try:
            # Save each cluster as a separate file
            for cluster_id, cluster_points in clusters_container.items():
                if cluster_id != -1:  # Skip noise
                    cluster_path = (
                        self.artifacts_dir / f"{self.rcsb_id}_cluster_{cluster_id}.npy"
                    )
                    np.save(cluster_path, np.array(cluster_points))
                    if cluster_id == largest_cluster_id:
                        self.tracker.add_artifact(self.stage, cluster_path)
        except Exception as cluster_error:
            print(f"Warning: Could not save individual clusters: {str(cluster_error)}")

        return {
            "largest_cluster": largest_cluster,
            "clusters_container": clusters_container,
            "largest_cluster_id": largest_cluster_id,
        }

```

ribctl/lib/npet/pipeline/constriction_identification_stage.py
```py
import json
import numpy as np
from pathlib import Path
from typing import Any, Dict

from ribctl.asset_manager.asset_types import AssetType
from ribctl.lib.landmarks.constriction_site import get_constriction
from ribctl.lib.npet.pipeline.base_stage import NPETPipelineStage
from ribctl.lib.npet.pipeline_status_tracker import ProcessingStage
from ribctl.lib.schema.types_ribosome import ConstrictionSite


class ConstrictionIdentificationStage(NPETPipelineStage):
    """
    Stage for identifying the constriction site in the ribosome tunnel.
    """
    
    @property
    def stage(self) -> ProcessingStage:
        return ProcessingStage.CONSTRICTION_IDENTIFICATION
    
    @property
    def stage_params(self) -> Dict[str, Any]:
        # No configurable parameters for this stage
        return {}
    
    def process(self, context: Dict[str, Any]) -> Dict[str, Any]:
        """
        Identify the constriction site within the ribosome tunnel.
        """
        try:
            # Get constriction point
            constriction_pt = get_constriction(self.rcsb_id)
            
            # Save constriction point as JSON
            constriction_path = AssetType.CONSTRICTION_SITE.get_path(self.rcsb_id)
            with open(constriction_path, 'w') as f:
                json.dump(ConstrictionSite(location=constriction_pt.tolist()).model_dump(), f)
            
            # Track the artifact
            self.tracker.add_artifact(self.stage, constriction_path)
            
            return {
                "constriction_pt": constriction_pt
            }
        except Exception as e:
            raise RuntimeError(f"Failed to identify constriction site: {str(e)}") from e
```

ribctl/lib/npet/pipeline/entity_filtering_stage.py
```py
import numpy as np
from pathlib import Path
from typing import Any, Dict, List

from ribctl.lib.npet.kdtree_approach import (
    ribosome_entities,
    filter_residues_parallel,
)
from ribctl.lib.npet.pipeline.base_stage import NPETPipelineStage
from ribctl.lib.npet.pipeline_status_tracker import NPETProcessingTracker, ProcessingStage
from ribctl.lib.npet.various_visualization import visualize_filtered_residues


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

ribctl/lib/npet/pipeline/exterior_mesh_stage.py
```py
import os
import numpy as np
import pyvista as pv
import open3d as o3d
from pathlib import Path
from typing import Any, Dict

from ribctl.asset_manager.asset_types import AssetType
from ribctl.lib.npet.alphalib import cif_to_point_cloud, fast_normal_estimation, quick_surface_points, validate_mesh_pyvista
from ribctl.lib.npet.kdtree_approach import apply_poisson_reconstruction
from ribctl.lib.npet.pipeline.base_stage import NPETPipelineStage
from ribctl.lib.npet.pipeline_status_tracker import NPETProcessingTracker, ProcessingStage
from ribctl.lib.npet.various_visualization import visualize_pointcloud, visualize_mesh


class AlphaShapeStage(NPETPipelineStage):
    """
    Stage for generating the alpha shape representation of the ribosome structure.
    """
    
    def __init__(self, 
                 rcsb_id: str, 
                 tracker: NPETProcessingTracker,
                 artifacts_dir: Path,
                 force: bool = False,
                 alpha_params: Dict[str, Any] = None):
        """
        Initialize the alpha shape stage with specific parameters.
        """
        super().__init__(rcsb_id, tracker, artifacts_dir, force)
        if alpha_params is None:
            raise ValueError("Alpha shape parameters must be provided")
        self.params = alpha_params     
    @property
    def stage(self) -> ProcessingStage:
        return ProcessingStage.ALPHA_SHAPE
    
    @property
    def stage_params(self) -> Dict[str, Any]:
        return self.params

    def process(self, context: Dict[str, Any]) -> Dict[str, Any]:
        """
        Generate the alpha shape representation of the ribosome structure.
        """
        # Get paths from context
        cifpath = context["cifpath"]
        ashapepath = context["ashapepath"]
        
        # Check if alpha shape already exists
        regenerate_alpha = self.force or not os.path.exists(ashapepath)
        
        if not regenerate_alpha:
            print(f"Using existing alpha shape mesh: {ashapepath}")
            self.tracker.add_artifact(self.stage, Path(ashapepath))
        else:
            print(f"Generating alpha shape mesh for {self.rcsb_id}")
            
            # Generate point cloud
            ptcloudpath = self.artifacts_dir / f"{self.rcsb_id}_structure_ptcloud.npy"
            
            if not os.path.exists(ptcloudpath) or self.force:
                print("Extracting point cloud from CIF file")
                first_assembly_chains = context["ro"].first_assembly_auth_asym_ids()
                ptcloud = cif_to_point_cloud(cifpath, first_assembly_chains, do_atoms=True)
                np.save(ptcloudpath, ptcloud)
            else:
                print(f"Loading existing point cloud from {ptcloudpath}")
                ptcloud = np.load(ptcloudpath)
            
            self.tracker.add_artifact(self.stage, ptcloudpath)
            
            # Surface extraction
            print("Beginning Delaunay 3D reconstruction")
            surface_pts = quick_surface_points(
                ptcloud, 
                self.params["d3d_alpha"], 
                self.params["d3d_tol"], 
                self.params["d3d_offset"]
            )
            
            # Save surface points
            surface_pts_path = self.artifacts_dir / f"{self.rcsb_id}_alpha_surface_points.npy"
            np.save(surface_pts_path, surface_pts)
            self.tracker.add_artifact(self.stage, surface_pts_path)
            
            # Visualize surface points
            try:
                surface_viz_path = self.artifacts_dir / f"{self.rcsb_id}_alpha_surface.png"
                visualize_pointcloud(surface_pts, self.rcsb_id, output_path=str(surface_viz_path))
                self.tracker.add_artifact(self.stage, surface_viz_path)
            except Exception as viz_error:
                print(f"Warning: Could not save visualization: {str(viz_error)}")
            
            # Normal estimation
            normal_estimated_pcd = fast_normal_estimation(
                surface_pts, 
                self.params["kdtree_radius"], 
                self.params["max_nn"], 
                self.params["tangent_planes_k"]
            )

            # --- CRITICAL FIX 1: Robust Normal Orientation ---
            # Checkerboard patterns in your viz show flipped normals.
            # We force all normals to point OUTWARD from the center.
            center = normal_estimated_pcd.get_center()
            # First, point them all at the center (inward)
            normal_estimated_pcd.orient_normals_towards_camera_location(camera_location=center)
            # Then flip them to face out
            normal_estimated_pcd.normals = o3d.utility.Vector3dVector(-np.asarray(normal_estimated_pcd.normals))
            
            # Save normal-estimated point cloud
            alpha_normals_path = self.artifacts_dir / f"{self.rcsb_id}_alpha_normals.ply"
            o3d.io.write_point_cloud(str(alpha_normals_path), normal_estimated_pcd)
            self.tracker.add_artifact(self.stage, alpha_normals_path)
            
            # Apply Poisson reconstruction
            apply_poisson_reconstruction(
                str(alpha_normals_path),
                ashapepath,
                recon_depth=self.params["PR_depth"],
                recon_pt_weight=self.params["PR_ptweight"],
            )
            
            # --- CRITICAL FIX 2: Explicit Hole Filling ---
            # Load the generated mesh and physically seal the gaps.
            mesh = pv.read(ashapepath)
            
            # We use a high threshold (2000) to bridge that specific edge gap.
            mesh = mesh.fill_holes(2000) 
            
            # Extract largest component and finalize geometry
            labeled = mesh.connectivity(largest=True)
            labeled = labeled.triangulate()
            labeled.save(ashapepath)
            
            # Check watertightness
            watertight = validate_mesh_pyvista(labeled)
            
            if not watertight:
                print("Warning: Alpha shape mesh is still not watertight, attempting aggressive repair...")
                # Optional: extra smoothing/repair pass if needed
                # labeled = labeled.smooth(n_iter=10).fill_holes(5000)
                # labeled.save(ashapepath)
            
            # Save final mesh visualization
            try:
                mesh_viz_path = self.artifacts_dir / f"{self.rcsb_id}_alpha_mesh.png"
                visualize_mesh(ashapepath, self.rcsb_id, output_path=str(mesh_viz_path))
                self.tracker.add_artifact(self.stage, mesh_viz_path)
            except Exception as viz_error:
                print(f"Warning: Could not save visualization: {str(viz_error)}")
            
            # Also add ASCII version as artifact if it exists
            ascii_path = Path(str(ashapepath).split(".")[0] + "_ascii.ply")
            if ascii_path.exists():
                self.tracker.add_artifact(self.stage, ascii_path)
        
        # Add the final alpha shape as artifact
        self.tracker.add_artifact(self.stage, Path(ashapepath))
        
        # Verify alpha shape exists for the remaining pipeline
        if not os.path.exists(ashapepath):
            raise FileNotFoundError(f"Alpha shape file {ashapepath} not found after generation")
            
        return {
            "ashapepath": ashapepath
        }
```

ribctl/lib/npet/pipeline/landmark_identification_stage.py
```py
import json
import numpy as np
from typing import Any, Dict

from ribctl.asset_manager.asset_types import AssetType
from ribctl.lib.npet.pipeline.base_stage import NPETPipelineStage
from ribctl.lib.npet.pipeline_status_tracker import ProcessingStage
from ribctl.lib.schema.types_ribosome import PTCInfo, ConstrictionSite


class LandmarkIdentificationStage(NPETPipelineStage):
    """
    Stage for loading and verifying landmark points (PTC and constriction site).
    
    This stage serves as a verification point to ensure both the PTC and
    constriction site are available for subsequent stages.
    """
    
    @property
    def stage(self) -> ProcessingStage:
        return ProcessingStage.LANDMARK_IDENTIFICATION
    
    @property
    def stage_params(self) -> Dict[str, Any]:
        # This stage doesn't have configurable parameters
        return {}
    
    def process(self, context: Dict[str, Any]) -> Dict[str, Any]:
        """
        Load and verify landmark points from previous stages or files.
        
        The stage either uses the PTC and constriction points from previous
        stages or loads them from disk if they weren't already computed.
        """
        try:
            # Check if we already have the landmarks from previous stages
            ptc_pt = context.get("ptc_pt")
            constriction_pt = context.get("constriction_pt")
            
            # If not available in context, try to load from files
            if ptc_pt is None or constriction_pt is None:
                # Load PTC and constriction points from previous stages
                ptc_path = AssetType.PTC.get_path(self.rcsb_id)
                constriction_path = AssetType.CONSTRICTION_SITE.get_path(self.rcsb_id)
                
                # If they exist, use them directly
                if ptc_path.exists() and constriction_path.exists():
                    with open(ptc_path, 'r') as f:
                        ptc_info = PTCInfo.model_validate(json.load(f))
                    with open(constriction_path, 'r') as f:
                        constriction_site = ConstrictionSite.model_validate(json.load(f))
                    ptc_pt = np.array(ptc_info.location)
                    constriction_pt = np.array(constriction_site.location)
                else:
                    raise FileNotFoundError("PTC or constriction site not found")
            
            # Add existing files as artifacts if they weren't already tracked
            ptc_path = AssetType.PTC.get_path(self.rcsb_id)
            if ptc_path.exists():
                self.tracker.add_artifact(self.stage, ptc_path)
                
            constriction_path = AssetType.CONSTRICTION_SITE.get_path(self.rcsb_id)
            if constriction_path.exists():
                self.tracker.add_artifact(self.stage, constriction_path)
    
            return {
                "ptc_pt": ptc_pt,
                "constriction_pt": constriction_pt
            }
            
        except Exception as e:
            raise RuntimeError(f"Failed to identify landmarks: {str(e)}") from e
```

ribctl/lib/npet/pipeline/mesh_reconstruction_stage.py
```py
from pathlib import Path
from typing import Any, Dict
import os

from ribctl.lib.npet.kdtree_approach import apply_poisson_reconstruction
from ribctl.lib.npet.pipeline.base_stage import NPETPipelineStage
from ribctl.lib.npet.pipeline_status_tracker import NPETProcessingTracker, ProcessingStage


class MeshReconstructionStage(NPETPipelineStage):
    """
    Stage for reconstructing the mesh from the point cloud with normals.
    
    This stage applies Poisson surface reconstruction to create a 
    3D mesh representing the tunnel surface.
    """
    
    def __init__(self, 
                 rcsb_id: str, 
                 tracker: NPETProcessingTracker,
                 artifacts_dir: Path,
                 force: bool = False,
                 depth: int = 6,
                 ptweight: int = 3):
        """
        Initialize the mesh reconstruction stage.
        
        Args:
            rcsb_id: The RCSB PDB identifier
            tracker: Processing tracker
            artifacts_dir: Directory for artifacts
            force: Whether to force regeneration
            depth: Depth parameter for Poisson reconstruction
            ptweight: Point weight for Poisson reconstruction
        """
        super().__init__(rcsb_id, tracker, artifacts_dir, force)
        self.depth = depth
        self.ptweight = ptweight
    
    @property
    def stage(self) -> ProcessingStage:
        return ProcessingStage.MESH_RECONSTRUCTION
    
    @property
    def stage_params(self) -> Dict[str, Any]:
        return {
            "depth": self.depth,
            "ptweight": self.ptweight,
        }
    
    def process(self, context: Dict[str, Any]) -> Dict[str, Any]:
        """
        Reconstruct the mesh from the point cloud with normals.
        
        Steps:
        1. Apply Poisson reconstruction to create the mesh
        2. Save the mesh and related files
        """
        normals_pcd_path = context["normals_pcd_path"]
        meshpath = context["meshpath"]
        
        # Apply Poisson reconstruction
        apply_poisson_reconstruction(
            str(normals_pcd_path),
            meshpath,
            recon_depth=self.depth,
            recon_pt_weight=self.ptweight,
        )
        
        # Add mesh as artifact
        self.tracker.add_artifact(self.stage, Path(meshpath))
        
        # Also save ASCII version as artifact if it exists
        ascii_path = Path(str(meshpath).split(".")[0] + "_ascii.ply")
        if ascii_path.exists():
            self.tracker.add_artifact(self.stage, ascii_path)
        
        return {
            "meshpath": meshpath
        }
```

ribctl/lib/npet/pipeline/normal_estimation_stage.py
```py
import open3d as o3d
from pathlib import Path
from typing import Any, Dict

from ribctl.lib.npet.kdtree_approach import estimate_normals
from ribctl.lib.npet.pipeline.base_stage import NPETPipelineStage
from ribctl.lib.npet.pipeline_status_tracker import NPETProcessingTracker, ProcessingStage


class NormalEstimationStage(NPETPipelineStage):
    """
    Stage for estimating normals for the surface points.
    
    This stage estimates surface normals for the point cloud,
    which are required for the Poisson surface reconstruction.
    """
    
    def __init__(self, 
                 rcsb_id: str, 
                 tracker: NPETProcessingTracker,
                 artifacts_dir: Path,
                 force: bool = False,
                 kdtree_radius: float = 10,
                 kdtree_max_nn: int = 15,
                 correction_tangent_planes_n: int = 10):
        """
        Initialize the normal estimation stage.
        
        Args:
            rcsb_id: The RCSB PDB identifier
            tracker: Processing tracker
            artifacts_dir: Directory for artifacts
            force: Whether to force regeneration
            kdtree_radius: Radius for kdtree search during normal estimation
            kdtree_max_nn: Maximum number of neighbors for normal estimation
            correction_tangent_planes_n: Number of neighbors for tangent plane fitting
        """
        super().__init__(rcsb_id, tracker, artifacts_dir, force)
        self.kdtree_radius = kdtree_radius
        self.kdtree_max_nn = kdtree_max_nn
        self.correction_tangent_planes_n = correction_tangent_planes_n
    
    @property
    def stage(self) -> ProcessingStage:
        return ProcessingStage.NORMAL_ESTIMATION
    
    @property
    def stage_params(self) -> Dict[str, Any]:
        return {
            "kdtree_radius": self.kdtree_radius,
            "kdtree_max_nn": self.kdtree_max_nn,
            "correction_tangent_planes_n": self.correction_tangent_planes_n,
        }
    
    def process(self, context: Dict[str, Any]) -> Dict[str, Any]:
        """
        Estimate normals for the surface points.
        
        Steps:
        1. Apply normal estimation to the surface points
        2. Save the point cloud with normals
        """
        surface_pts = context["surface_pts"]
        
        # Estimate normals
        normal_estimated_pcd = estimate_normals(
            surface_pts,
            kdtree_radius=self.kdtree_radius,
            kdtree_max_nn=self.kdtree_max_nn,
            correction_tangent_planes_n=self.correction_tangent_planes_n,
        )
        
        # Save the point cloud with normals
        normals_pcd_path = self.artifacts_dir / f"{self.rcsb_id}_normal_estimated_pcd.ply"
        o3d.io.write_point_cloud(str(normals_pcd_path), normal_estimated_pcd)
        self.tracker.add_artifact(self.stage, normals_pcd_path)
        
        return {
            "normal_estimated_pcd": normal_estimated_pcd,
            "normals_pcd_path": normals_pcd_path
        }
```

ribctl/lib/npet/pipeline/point_cloud_processing_stage.py
```py
import numpy as np
import pyvista as pv
from pathlib import Path
from typing import Any, Dict

from ribctl.lib.npet.kdtree_approach import (
    transform_points_to_C0,
    create_point_cloud_mask,
    transform_points_from_C0,
)
from ribctl.lib.npet.pipeline.base_stage import NPETPipelineStage
from ribctl.lib.npet.pipeline_status_tracker import (
    NPETProcessingTracker,
    ProcessingStage,
)
from ribctl.lib.npet.various_visualization import visualize_pointcloud


class PointCloudProcessingStage(NPETPipelineStage):
    """
    Stage for transforming points to C0 space, creating the tunnel mask,
    and generating points representing the empty space inside the ribosome.
    """

    def __init__(
        self,
        rcsb_id: str,
        tracker: NPETProcessingTracker,
        artifacts_dir: Path,
        force: bool = False,
        radius: float = 35,
        height: float = 120,
        voxel_size: float = 1,
        atom_size: float = 2,
    ):
        """
        Initialize the point cloud processing stage.

        Args:
            rcsb_id: The RCSB PDB identifier
            tracker: Processing tracker
            artifacts_dir: Directory for artifacts
            force: Whether to force regeneration
            radius: Radius of the tunnel cylinder
            height: Height of the tunnel cylinder
            voxel_size: Size of voxels for discretization
            atom_size: Size of atoms for creating the mask
        """
        super().__init__(rcsb_id, tracker, artifacts_dir, force)
        self.radius = radius
        self.height = height
        self.voxel_size = voxel_size
        self.atom_size = atom_size

    @property
    def stage(self) -> ProcessingStage:
        return ProcessingStage.POINT_CLOUD_PROCESSING

    @property
    def stage_params(self) -> Dict[str, Any]:
        return {
            "radius": self.radius,
            "height": self.height,
            "voxel_size": self.voxel_size,
            "atom_size": self.atom_size,
        }

    def process(self, context: Dict[str, Any]) -> Dict[str, Any]:
        """
        Process the filtered points to create a representation of the tunnel.

        Steps:
        1. Transform filtered points to canonical space (C0)
        2. Create a mask representing occupied space
        3. Find points in empty space within the tunnel region
        4. Transform back to world coordinates
        5. Select only points inside the alpha shape mesh
        """
        filtered_points = context["filtered_points"]
        ptc_pt           = context["ptc_pt"]
        constriction_pt = context["constriction_pt"]
        ashapepath = context["ashapepath"]

        # Transform points to canonical space (C0)
        transformed_points = transform_points_to_C0(
            filtered_points, ptc_pt, constriction_pt
        )

        # Save transformed points
        transformed_points_path = (
            self.artifacts_dir / f"{self.rcsb_id}_transformed_points.npy"
        )
        np.save(transformed_points_path, transformed_points)
        self.tracker.add_artifact(self.stage, transformed_points_path)

        # Create mask representing occupied space
        mask, (x, y, z) = create_point_cloud_mask(
            transformed_points,
            radius=self.radius,
            height=self.height,
            voxel_size=self.voxel_size,
            radius_around_point=self.atom_size,
        )

        # Find coordinates of empty space
        points = np.where(~mask)
        empty_coordinates = np.column_stack((x[points[0]], y[points[1]], z[points[2]]))

        # Transform back to world coordinates
        back_projected = transform_points_from_C0(
            empty_coordinates, ptc_pt, constriction_pt
        )

        # Save back projected points
        back_projected_path = self.artifacts_dir / f"{self.rcsb_id}_back_projected.npy"
        np.save(back_projected_path, back_projected)
        self.tracker.add_artifact(self.stage, back_projected_path)

        ashape_watertight_mesh = pv.read(ashapepath)
        select = pv.PolyData(back_projected).select_enclosed_points(
            ashape_watertight_mesh
        )
        mask = select["SelectedPoints"]
        interior = back_projected[mask == 1]
        empty_in_world_coords = np.array(interior)

        # Save interior points
        interior_points_path = (
            self.artifacts_dir / f"{self.rcsb_id}_interior_points.npy"
        )
        np.save(interior_points_path, empty_in_world_coords)
        self.tracker.add_artifact(self.stage, interior_points_path)

        # Save visualization if possible
        try:
            pc_viz_path = self.artifacts_dir / f"{self.rcsb_id}_point_cloud.png"
            visualize_pointcloud(
                empty_in_world_coords, self.rcsb_id, output_path=str(pc_viz_path)
            )
            self.tracker.add_artifact(self.stage, pc_viz_path)
        except Exception as viz_error:
            print(
                f"Warning: Could not save point cloud visualization: {str(viz_error)}"
            )

        return {"empty_in_world_coords": empty_in_world_coords}

```

ribctl/lib/npet/pipeline/refinement_stage.py
```py
import numpy as np
from pathlib import Path
from typing import Any, Dict

from ribctl.lib.npet.kdtree_approach import (
    DBSCAN_capture,
    DBSCAN_pick_largest_cluster,
)
from ribctl.lib.npet.pipeline.base_stage import NPETPipelineStage
from ribctl.lib.npet.pipeline_status_tracker import NPETProcessingTracker, ProcessingStage
from ribctl.lib.npet.various_visualization import visualize_DBSCAN_CLUSTERS_particular_eps_minnbrs


class RefinementStage(NPETPipelineStage):
    """
    Stage for refining the cluster to improve tunnel representation.
    
    This stage applies a second round of DBSCAN with stricter parameters
    to the largest cluster from the previous stage, in order to further 
    refine the tunnel representation.
    """
    
    def __init__(self, 
                 rcsb_id: str, 
                 tracker: NPETProcessingTracker,
                 artifacts_dir: Path,
                 force: bool = False,
                 epsilon: float = 3.5,
                 min_samples: int = 175):
        """
        Initialize the refinement stage.
        
        Args:
            rcsb_id: The RCSB PDB identifier
            tracker: Processing tracker
            artifacts_dir: Directory for artifacts
            force: Whether to force regeneration
            epsilon: DBSCAN epsilon parameter for refinement
            min_samples: DBSCAN min_samples parameter for refinement
        """
        super().__init__(rcsb_id, tracker, artifacts_dir, force)
        self.epsilon = epsilon
        self.min_samples = min_samples
    
    @property
    def stage(self) -> ProcessingStage:
        return ProcessingStage.REFINEMENT
    
    @property
    def stage_params(self) -> Dict[str, Any]:
        return {
            "epsilon": self.epsilon,
            "min_samples": self.min_samples,
        }
    
    def process(self, context: Dict[str, Any]) -> Dict[str, Any]:
        """
        Apply a second round of DBSCAN to refine the tunnel representation.
        
        Steps:
        1. Apply DBSCAN with stricter parameters to the largest cluster
        2. Extract the largest refined cluster
        3. Save the refined cluster data and visualization
        """
        largest_cluster = context["largest_cluster"]
        ptc_pt = context["ptc_pt"]
        constriction_pt = context["constriction_pt"]
        
        # Apply second DBSCAN for refinement
        db_2, refined_clusters_container = DBSCAN_capture(
            largest_cluster, self.epsilon, self.min_samples
        )
        
        # Extract the largest refined cluster
        refined_cluster, refined_cluster_id = DBSCAN_pick_largest_cluster(
            refined_clusters_container
        )
        
        # Save refined cluster
        refined_cluster_path = self.artifacts_dir / f"{self.rcsb_id}_refined_cluster.npy"
        np.save(refined_cluster_path, refined_cluster)
        self.tracker.add_artifact(self.stage, refined_cluster_path)
        
        # Save visualization if possible
        try:
            refined_viz_path = self.artifacts_dir / f"{self.rcsb_id}_refined_clusters.png"
            visualize_DBSCAN_CLUSTERS_particular_eps_minnbrs(
                refined_clusters_container,
                self.epsilon,
                self.min_samples,
                ptc_pt,
                constriction_pt,
                refined_cluster,
                context.get("radius", 35),  # Default radius if not in context
                context.get("height", 120),  # Default height if not in context
                output_path=str(refined_viz_path)
            )
            self.tracker.add_artifact(self.stage, refined_viz_path)
        except Exception as viz_error:
            print(f"Warning: Could not save refined cluster visualization: {str(viz_error)}")
        
        # Save individual refined clusters
        try:
            for cluster_id, cluster_points in refined_clusters_container.items():
                if cluster_id != -1:  # Skip noise
                    cluster_path = self.artifacts_dir / f"{self.rcsb_id}_refined_cluster_{cluster_id}.npy"
                    np.save(cluster_path, np.array(cluster_points))
                    if cluster_id == refined_cluster_id:
                        self.tracker.add_artifact(self.stage, cluster_path)
        except Exception as cluster_error:
            print(f"Warning: Could not save individual refined clusters: {str(cluster_error)}")
        
        return {
            "refined_cluster": refined_cluster,
            "refined_clusters_container": refined_clusters_container,
            "refined_cluster_id": refined_cluster_id
        }
```

ribctl/lib/npet/pipeline/ptc_identification_stage.py
```py
import numpy as np
from pathlib import Path
from typing import Any, Dict

from ribctl.asset_manager.asset_types import AssetType
from ribctl.lib.landmarks.ptc_via_trna import PTC_location
from ribctl.lib.npet.pipeline.base_stage import NPETPipelineStage
from ribctl.lib.npet.pipeline_status_tracker import ProcessingStage


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
            ptc_info = PTC_location(self.rcsb_id)
            ptc_pt = np.array(ptc_info.location)
            
            ptc_json_path = AssetType.PTC.get_path(self.rcsb_id)
            with open(ptc_json_path, 'w') as f:
                f.write(ptc_info.model_dump_json())
            
            self.tracker.add_artifact(self.stage, ptc_json_path)
            
            return {
                "ptc_info": ptc_info,
                "ptc_pt": ptc_pt
            }

        except Exception as e:
            raise RuntimeError(f"Failed to identify PTC: {str(e)}") from e
```

ribctl/lib/npet/pipeline/setup_stage.py
```py
from pathlib import Path
from typing import Any, Dict

from ribctl.asset_manager.asset_types import AssetType
from ribctl.lib.npet.pipeline.base_stage import NPETPipelineStage
from ribctl.lib.npet.pipeline_status_tracker import ProcessingStage
from ribctl.lib.npet.tunnel_asset_manager import TunnelMeshAssetsManager
from ribctl.ribosome_ops import RibosomeOps


class SetupStage(NPETPipelineStage):
    """
    Setup stage - initializes resources and validates inputs.
    """
    
    @property
    def stage(self) -> ProcessingStage:
        return ProcessingStage.SETUP
    
    @property
    def stage_params(self) -> Dict[str, Any]:
        # This stage doesn't have configurable parameters
        return {}
    
    def process(self, context: Dict[str, Any]) -> Dict[str, Any]:
        """
        Initialize the pipeline by setting up paths and resources.
        """
        # Initialize asset managers
        assets = TunnelMeshAssetsManager(self.rcsb_id)
        ro = RibosomeOps(self.rcsb_id)
        
        # Set up paths
        cifpath = AssetType.MMCIF.get_path(self.rcsb_id)
        ashapepath = AssetType.ALPHA_SHAPE.get_path(self.rcsb_id)
        meshpath = AssetType.NPET_MESH.get_path(self.rcsb_id)
        
        # Track the CIF file as an artifact
        self.tracker.add_artifact(self.stage, Path(cifpath))
        
        return {
            "assets": assets,
            "ro": ro,
            "profile": ro.profile,
            "cifpath": cifpath,
            "ashapepath": ashapepath,
            "meshpath": meshpath
        }
```

ribctl/lib/npet/pipeline/surface_extraction_stage.py
```py
import numpy as np
from pathlib import Path
from typing import Any, Dict

from ribctl.lib.npet.kdtree_approach import ptcloud_convex_hull_points
from ribctl.lib.npet.pipeline.base_stage import NPETPipelineStage
from ribctl.lib.npet.pipeline_status_tracker import NPETProcessingTracker, ProcessingStage
from ribctl.lib.npet.various_visualization import visualize_pointcloud


class SurfaceExtractionStage(NPETPipelineStage):
    """
    Stage for extracting the surface points from the refined cluster.
    
    This stage extracts points on the surface of the refined cluster
    using a convex hull approach, which will be used for mesh generation.
    """
    
    def __init__(self, 
                 rcsb_id: str, 
                 tracker: NPETProcessingTracker,
                 artifacts_dir: Path,
                 force: bool = False,
                 alpha: float = 2,
                 tolerance: float = 1,
                 offset: float = 2):
        """
        Initialize the surface extraction stage.
        
        Args:
            rcsb_id: The RCSB PDB identifier
            tracker: Processing tracker
            artifacts_dir: Directory for artifacts
            force: Whether to force regeneration
            alpha: Alpha value for surface extraction
            tolerance: Tolerance for surface extraction
            offset: Offset for surface extraction
        """
        super().__init__(rcsb_id, tracker, artifacts_dir, force)
        self.alpha = alpha
        self.tolerance = tolerance
        self.offset = offset
    
    @property
    def stage(self) -> ProcessingStage:
        return ProcessingStage.SURFACE_EXTRACTION
    
    @property
    def stage_params(self) -> Dict[str, Any]:
        return {
            "alpha": self.alpha,
            "tolerance": self.tolerance,
            "offset": self.offset,
        }
    
    def process(self, context: Dict[str, Any]) -> Dict[str, Any]:
        """
        Extract surface points from the refined cluster.
        
        Steps:
        1. Apply convex hull algorithm to extract surface points
        2. Save the surface points and visualization
        """
        refined_cluster = context["refined_cluster"]
        
        # Extract surface points
        surface_pts = ptcloud_convex_hull_points(
            refined_cluster, self.alpha, self.tolerance, self.offset
        )
        
        # Save surface points
        surface_pts_path = self.artifacts_dir / f"{self.rcsb_id}_surface_points.npy"
        np.save(surface_pts_path, surface_pts)
        self.tracker.add_artifact(self.stage, surface_pts_path)
        
        # Save visualization if possible
        try:
            surface_viz_path = self.artifacts_dir / f"{self.rcsb_id}_surface_points.png"
            visualize_pointcloud(surface_pts, self.rcsb_id, output_path=str(surface_viz_path))
            self.tracker.add_artifact(self.stage, surface_viz_path)
        except Exception as viz_error:
            print(f"Warning: Could not save surface point visualization: {str(viz_error)}")
        
        return {
            "surface_pts": surface_pts
        }
```

ribctl/lib/npet/pipeline/validation_stage.py
```py
import os
from pathlib import Path
from typing import Any, Dict
from ribctl.lib.npet.alphalib import validate_mesh_pyvista
from ribctl.lib.npet.pipeline.base_stage import NPETPipelineStage
from ribctl.lib.npet.pipeline_status_tracker import ProcessingStage
from ribctl.lib.npet.various_visualization import visualize_mesh


class ValidationStage(NPETPipelineStage):
    """
    Stage for validating the final mesh.
    
    This stage checks if the generated mesh is watertight and meets
    quality standards. It also generates visualizations of the final mesh.
    """
    
    @property
    def stage(self) -> ProcessingStage:
        return ProcessingStage.VALIDATION
    
    @property
    def stage_params(self) -> Dict[str, Any]:
        return {}
    
    def process(self, context: Dict[str, Any]) -> Dict[str, Any]:
        """
        Validate the final mesh and ensure it meets quality standards.
        
        Steps:
        1. Generate visualization of the final mesh
        2. Check if the mesh is watertight
        3. If not watertight, remove the mesh files and raise an error
        """
        meshpath = context["meshpath"]
        
        # Save mesh visualization if possible
        try:
            mesh_viz_path = self.artifacts_dir / f"{self.rcsb_id}_mesh.png"
            visualize_mesh(meshpath, self.rcsb_id, output_path=str(mesh_viz_path))
            self.tracker.add_artifact(self.stage, mesh_viz_path)
        except Exception as viz_error:
            print(f"Warning: Could not save mesh visualization: {str(viz_error)}")
        
        # Check watertightness
        watertight = validate_mesh_pyvista(meshpath)
        
        if not watertight:
            error_msg = "Mesh is not watertight"
            print("XXXX Watertightness check failed, removing", meshpath, " XXXX")
            os.remove(meshpath)
            
            # Also remove ASCII version if it exists
            ascii_path = str(meshpath).split(".")[0] + "_ascii.ply"
            if os.path.exists(ascii_path):
                os.remove(ascii_path)
                
            raise ValueError(error_msg)
        
        return {
            "watertight": watertight
        }
```

