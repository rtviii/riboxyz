from pprint import pprint
import typing
from matplotlib import pyplot as plt
import pyvista as pv
import json
import numpy as np
import numpy as np
from mesh_generation.paths import *
from mesh_generation.voxelize import (
    expand_atomcenters_to_spheres_threadpool,
    normalize_atom_coordinates,
)
from ribctl import EXIT_TUNNEL_WORK, POISSON_RECON_BIN, RIBETL_DATA
from ribctl.etl.ribosome_assets import RibosomeAssets
from ribctl.lib.libmsa import Taxid


diagram_tunnels = {
    "bacteria": [
        "4W29", # GOOD
        "6WD4", # GOOD
        "7UNV",
        "5DM6",

        "5O60",
        "8BUU",
        "6HMA",
        "7RYH", #GOOD

        "7MSZ",
        "7P7T",
        "5MYJ", # GOOD, weird structure
        "7JIL",
        # 3j9y -- good
    ],
    "eukaryota": [
        "6P5N",
        "7QGG",
        "4UG0", # GOOD
        "4U3M",

        "7OYB",
        "7OLC",
        "8EUI",
        "5XXB",

        "4V91",
        "6XU8",
        "4V7E",
        "8P5D",

        "8BTR",
        "3JBO",
        "7CPU",
        "7Q08",
        "6AZ3", #GOOD

        "5T5H",
        "5XY3",
        "7QEP",
    ],
    "archaea": ["4V6U",  # GOOD
                "4V9F",] # GOOD
}

FONT                  = 'courier'
CHAIN_PT_SIZE         = 8
PTC_PT_SIZE           = 20
CHAIN_LANDMARK_COLORS = ["purple","orange", "cornflowerblue", "cornsilk", "crimson", "darkblue", "darkcyan", "darkgoldenrod", "darkgray", "darkgreen", "darkkhaki", "darkmagenta", "darkolivegreen", "darkorange", "darkorchid", "darkred",
"rebeccapurple",
"rosybrown",
"royalblue",
"saddlebrown",
"salmon",
"sandybrown",
"seagreen"]

available_tunnels = {
    "archaea": ["1W2B", "4V6U", "4V9F"],
    "bacteria": [
        "1NKW",
        "1NWX",
        "1NWY",
        "1SM1",
        "1VVJ",
        "1VY4",
        "1VY5",
        "1VY6",
        "1VY7",
        "1XBP",
        "2ZJP",
        "2ZJQ",
        "2ZJR",
        "3CF5",
        "3DLL",
        "3J3V",
        "3J3W",
        "3J7Z",
        "3J9W",
        "3J9Y",
        "3J9Z",
        "3JA1",
        "3JBU",
        "3JCD",
        "3JCJ",
        "3PIO",
        "3PIP",
        "4IO9",
        "4IOA",
        "4IOC",
        "4L47",
        "4L71",
        "4LEL",
        "4LFZ",
        "4LNT",
        "4LSK",
        "4LT8",
        "4P6F",
        "4P70",
        "4TUA",
        "4TUB",
        "4TUC",
        "4TUD",
        "4TUE",
        "4U1U",
        "4U1V",
        "4U20",
        "4U24",
        "4U25",
        "4U26",
        "4U27",
        "4UY8",
        "4V4H",
        "4V4I",
        "4V4J",
        "4V4Q",
        "4V50",
        "4V51",
        "4V52",
        "4V53",
        "4V54",
        "4V56",
        "4V57",
        "4V5A",
        "4V5B",
        "4V5C",
        "4V5D",
        "4V5E",
        "4V5F",
        "4V5G",
        "4V5J",
        "4V5K",
        "4V5L",
        "4V5P",
        "4V5Q",
        "4V5R",
        "4V5S",
        "4V63",
        "4V64",
        "4V67",
        "4V6A",
        "4V6C",
        "4V6D",
        "4V6E",
        "4V6F",
        "4V6G",
        "4V7L",
        "4V7M",
        "4V7P",
        "4V7S",
        "4V7T",
        "4V7U",
        "4V7V",
        "4V7W",
        "4V7X",
        "4V7Y",
        "4V7Z",
        "4V83",
        "4V84",
        "4V85",
        "4V89",
        "4V8A",
        "4V8G",
        "4V8H",
        "4V8I",
        "4V8J",
        "4V8N",
        "4V8O",
        "4V8Q",
        "4V8U",
        "4V8X",
        "4V90",
        "4V95",
        "4V97",
        "4V9A",
        "4V9B",
        "4V9C",
        "4V9D",
        "4V9H",
        "4V9I",
        "4V9J",
        "4V9K",
        "4V9L",
        "4V9N",
        "4V9O",
        "4V9P",
        "4V9Q",
        "4V9R",
        "4V9S",
        "4W29",
        "4W2E",
        "4W2F",
        "4W2G",
        "4W2H",
        "4W2I",
        "4W4G",
        "4WF1",
        "4WOI",
        "4WPO",
        "4WQ1",
        "4WQF",
        "4WQR",
        "4WQU",
        "4WQY",
        "4WR6",
        "4WRA",
        "4WRO",
        "4WSD",
        "4WSM",
        "4WT1",
        "4WT8",
        "4WU1",
        "4WWW",
        "4WZD",
        "4WZO",
        "4XEJ",
        "4Y4O",
        "4Y4P",
        "4YBB",
        "4YPB",
        "4YZV",
        "4Z3S",
        "4Z8C",
        "4ZER",
        "4ZSN",
        "5AFI",
        "5CZP",
        "5DM6",
        "5DM7",
        "5DOX",
        "5DOY",
        "5E7K",
        "5E81",
        "5EL4",
        "5EL5",
        "5EL6",
        "5EL7",
        "5F8K",
        "5FDU",
        "5FDV",
        "5GAD",
        "5GAE",
        "5GAG",
        "5GAH",
        "5H5U",
        "5HAU",
        "5HCP",
        "5HCQ",
        "5HCR",
        "5HD1",
        "5IB8",
        "5IBB",
        "5IMQ",
        "5IT8",
        "5J30",
        "5J3C",
        "5J4B",
        "5J4C",
        "5J4D",
        "5J5B",
        "5J7L",
        "5J88",
        "5J8A",
        "5J8B",
        "5J91",
        "5JC9",
        "5JTE",
        "5JU8",
        "5KCR",
        "5KCS",
        "5KPS",
        "5KPW",
        "5KPX",
        "5L3P",
        "5LZA",
        "5LZD",
        "5LZE",
        "5MDV",
        "5MDW",
        "5MDY",
        "5MDZ",
        "5MGP",
        "5MYJ",
        "5NDJ",
        "5NDK",
        "5NGM",
        "5NP6",
        "5NWY",
        "5O2R",
        "5O60",
        "5OT7",
        "5T7V",
        "5U4I",
        "5U9F",
        "5U9G",
        "5UQ7",
        "5UQ8",
        "5UYK",
        "5UYL",
        "5UYM",
        "5UYP",
        "5UYQ",
        "5V7Q",
        "5V8I",
        "5V93",
        "5VP2",
        "5VPO",
        "5VPP",
        "5W4K",
        "5WDT",
        "5WE4",
        "5WE6",
        "5WF0",
        "5WFK",
        "5WFS",
        "5WIS",
        "5WIT",
        "5XYM",
        "5ZEB",
        "5ZEP",
        "5ZET",
        "5ZLU",
        "6B4V",
        "6BOH",
        "6BOK",
        "6BU8",
        "6BUW",
        "6BY1",
        "6BZ6",
        "6BZ7",
        "6BZ8",
        "6C4I",
        "6C5L",
        "6CAE",
        "6CFJ",
        "6CFK",
        "6CFL",
        "6CZR",
        "6DNC",
        "6DZI",
        "6DZP",
        "6ENF",
        "6ENJ",
        "6ENU",
        "6FKR",
        "6GC8",
        "6GSJ",
        "6GSK",
        "6GSL",
        "6GWT",
        "6GXM",
        "6GXN",
        "6GXO",
        "6GZQ",
        "6HA1",
        "6HA8",
        "6HMA",
        "6HTQ",
        "6I0Y",
        "6I7V",
        "6N1D",
        "6N9E",
        "6N9F",
        "6ND5",
        "6ND6",
        "6NDK",
        "6NSH",
        "6NTA",
        "6NUO",
        "6NWY",
        "6O3M",
        "6O97",
        "6O9J",
        "6OF1",
        "6OF6",
        "6OFX",
        "6OG7",
        "6OGF",
        "6OGI",
        "6OJ2",
        "6OM6",
        "6OPE",
        "6ORD",
        "6ORE",
        "6ORL",
        "6OSK",
        "6OSQ",
        "6OT3",
        "6OTR",
        "6OUO",
        "6OXA",
        "6OXI",
        "6PJ6",
        "6Q95",
        "6Q97",
        "6Q9A",
        "6QNQ",
        "6QNR",
        "6S0K",
        "6S0X",
        "6S0Z",
        "6S12",
        "6S13",
        "6SZS",
        "6TBV",
        "6TC3",
        "6U48",
        "6UCQ",
        "6UO1",
        "6V39",
        "6V3A",
        "6V3B",
        "6V3D",
        "6VU3",
        "6VWL",
        "6VWM",
        "6VWN",
        "6VYQ",
        "6VYR",
        "6VYS",
        "6WD0",
        "6WD1",
        "6WD2",
        "6WD3",
        "6WD4",
        "6WD5",
        "6WD6",
        "6WD7",
        "6WD8",
        "6WD9",
        "6WDA",
        "6WDD",
        "6WDE",
        "6WDF",
        "6WDG",
        "6WDJ",
        "6WDK",
        "6WDL",
        "6WDM",
        "6WRU",
        "6X7F",
        "6X7K",
        "6XDQ",
        "6XHV",
        "6XHW",
        "6XHX",
        "6XHY",
        "6XQD",
        "6XQE",
        "6XZ7",
        "6XZA",
        "6XZB",
        "6Y69",
        "6YSR",
        "6YSS",
        "6YST",
        "6YSU",
        "7AQC",
        "7AQD",
        "7ASO",
        "7ASP",
        "7AZO",
        "7AZS",
        "7B5K",
        "7BL2",
        "7BL3",
        "7BL4",
        "7BV8",
        "7D6Z",
        "7F0D",
        "7JQL",
        "7JQM",
        "7JSS",
        "7JSW",
        "7JSZ",
        "7JT1",
        "7JT2",
        "7JT3",
        "7K00",
        "7K50",
        "7K51",
        "7K52",
        "7K53",
        "7K54",
        "7K55",
        "7KGB",
        "7LH5",
        "7LV0",
        "7LVK",
        "7M4V",
        "7M4W",
        "7M4X",
        "7M4Y",
        "7M4Z",
        "7MD7",
        "7MSC",
        "7MSH",
        "7MSM",
        "7MSZ",
        "7MT2",
        "7MT3",
        "7MT7",
        "7N1P",
        "7N2C",
        "7N2U",
        "7N2V",
        "7N30",
        "7N31",
        "7NSO",
        "7NSP",
        "7NSQ",
        "7NWT",
        "7NWW",
        "7O5B",
        "7OIF",
        "7OIG",
        "7OII",
        "7OIZ",
        "7OJ0",
        "7OPE",
        "7OT5",
        "7P3K",
        "7P7T",
        "7PJS",
        "7PJV",
        "7PJY",
        "7QG8",
        "7QGN",
        "7QGU",
        "7QH4",
        "7QV1",
        "7QV3",
        "7RQ8",
        "7RQ9",
        "7RQA",
        "7RQB",
        "7RQC",
        "7RQD",
        "7RQE",
        "7RYF",
        "7RYG",
        "7RYH",
        "7S0S",
        "7S1G",
        "7S1H",
        "7S1I",
        "7S1J",
        "7S1K",
        "7SFR",
        "7SS9",
        "7SSD",
        "7SSL",
        "7SSN",
        "7SSO",
        "7SSW",
        "7ST2",
        "7ST6",
        "7ST7",
        "7TOS",
        "7U2H",
        "7U2I",
        "7U2J",
        "7UG7",
        "7UNR",
        "7UNU",
        "7UNV",
        "7UNW",
        "7XAM",
        "7Y41",
        "8BUU",
    ],
    "eukaryota": [
        "3J6X",
        "3J6Y",
        "3J77",
        "3J78",
        "3J79",
        "3J7O",
        "3J7P",
        "3J7Q",
        "3J7R",
        "3J92",
        "3JAG",
        "3JAH",
        "3JAI",
        "3JAJ",
        "3JAN",
        "3JBN",
        "3JBO",
        "3JBP",
        "3JCS",
        "3JCT",
        "4D5Y",
        "4D67",
        "4U3M",
        "4U3N",
        "4U3U",
        "4U4N",
        "4U4O",
        "4U4Q",
        "4U4R",
        "4U4U",
        "4U4Y",
        "4U4Z",
        "4U50",
        "4U51",
        "4U52",
        "4U53",
        "4U55",
        "4U56",
        "4U6F",
        "4UG0",
        "4UJC",
        "4UJD",
        "4UJE",
        "4V3P",
        "4V6I",
        "4V6W",
        "4V6X",
        "4V7E",
        "4V7H",
        "4V7R",
        "4V88",
        "4V8M",
        "4V8P",
        "4V8T",
        "4V8Y",
        "4V8Z",
        "4V91",
        "5AJ0",
        "5APN",
        "5APO",
        "5DAT",
        "5DC3",
        "5DGE",
        "5DGF",
        "5DGV",
        "5FCI",
        "5FCJ",
        "5GAK",
        "5H4P",
        "5I4L",
        "5IT7",
        "5JCS",
        "5JUO",
        "5JUP",
        "5JUS",
        "5JUT",
        "5JUU",
        "5LKS",
        "5LYB",
        "5LZS",
        "5LZT",
        "5LZU",
        "5LZV",
        "5LZW",
        "5LZX",
        "5LZY",
        "5LZZ",
        "5M1J",
        "5MC6",
        "5MEI",
        "5NDG",
        "5NDV",
        "5NDW",
        "5OBM",
        "5ON6",
        "5T2A",
        "5T2C",
        "5T5H",
        "5T62",
        "5T6R",
        "5TBW",
        "5TGA",
        "5TGM",
        "5UMD",
        "5XXB",
        "5XY3",
        "6AZ3",
        "6D90",
        "6D9J",
        "6GQ1",
        "6GQB",
        "6GQV",
        "6GZ3",
        "6GZ4",
        "6GZ5",
        "6HCF",
        "6HCJ",
        "6HCM",
        "6HCQ",
        "6HHQ",
        "6IP5",
        "6IP6",
        "6IP8",
        "6M62",
        "6N8J",
        "6N8K",
        "6N8L",
        "6N8M",
        "6N8N",
        "6N8O",
        "6OIG",
        "6OLE",
        "6OLF",
        "6OLG",
        "6OLI",
        "6OLZ",
        "6OM0",
        "6OM7",
        "6P5I",
        "6P5J",
        "6P5K",
        "6P5N",
        "6QIK",
        "6QT0",
        "6QTZ",
        "6R5Q",
        "6R6G",
        "6R6P",
        "6R7Q",
        "6R84",
        "6R86",
        "6R87",
        "6RI5",
        "6RM3",
        "6RZZ",
        "6S05",
        "6SGC",
        "6SWA",
        "6T59",
        "6UZ7",
        "6W6L",
        "6WOO",
        "6XIQ",
        "6XIR",
        "6XU6",
        "6XU7",
        "6XU8",
        "6Y0G",
        "6Y2L",
        "6Y57",
        "6YLX",
        "6Z6J",
        "6Z6K",
        "6Z6L",
        "6Z6M",
        "6Z6N",
        "6ZU5",
        "7A01",
        "7AZY",
        "7BHP",
        "7CPU",
        "7CPV",
        "7LS1",
        "7LS2",
        "7MDZ",
        "7NFX",
        "7NRC",
        "7NRD",
        "7NWG",
        "7NWH",
        "7NWI",
        "7OF1",
        "7OH3",
        "7OHQ",
        "7OLC",
        "7OLD",
        "7OSA",
        "7OYA",
        "7OYB",
        "7OYC",
        "7PWG",
        "7PWO",
        "7PZY",
        "7Q08",
        "7Q0F",
        "7Q0P",
        "7Q0R",
        "7QCA",
        "7QEP",
        "7QGG",
        "7QWQ",
        "7QWR",
        "7QWS",
        "7RR5",
        "7TM3",
        "7TOP",
        "7TOQ",
        "7TOR",
        "7TUT",
        "7UG6",
        "7Z34",
        "7ZPQ",
        "7ZRS",
        "7ZS5",
        "7ZUW",
        "7ZUX",
        "7ZW0",
        "8BIP",
        "8BJQ",
        "8BQD",
        "8BQX",
        "8BR8",
        "8BRM",
        "8BSI",
        "8BSJ",
        "8BTD",
        "8BTR",
        "8EUG",
        "8EUI",
        "8HFR",
        "8IR3",
        "8P5D",
        "8P60",
        "8Q5I",
    ],
}

dbscan_pairs = [
    (2,33),
    (2.3, 57),
    (3, 123),
    (3.2, 145),
    (3.5, 175),
    (4.2, 280),
    (5, 490),
    (5.5,600)
]


def move_cords_to_normalized_cord_frame(
    grid_dimensions    : np.ndarray,
    translation_vectors: np.ndarray,
    original_cords     : np.ndarray,
)                  : 
    """this is a helper function for plotting to move additional atom coordinates into the cord frame of the mesh (vox grid indices)"""
    normalized_original_cords = ( original_cords - translation_vectors[0] + translation_vectors[1] )
    voxel_size = 1
    normalized_original_cords_quantized = np.round(
        normalized_original_cords / voxel_size
    ).astype(int)
    vox_grid = np.zeros(grid_dimensions)
    vox_grid[
        normalized_original_cords_quantized[:, 0],
        normalized_original_cords_quantized[:, 1],
        normalized_original_cords_quantized[:, 2],
    ] = 1
    __xyz_v_positive_ix = np.asarray( np.where(vox_grid == 1) )
    return __xyz_v_positive_ix.T

def retrieve_ptc_and_chain_atoms(rcsb_id):
        with open( tunnel_atom_encoding_path(rcsb_id), "r", ) as infile:
            bbox_atoms: list[dict] = json.load(infile)
            _atom_centers       = np.array(list(map(lambda x: x["coord"], bbox_atoms)))
            _vdw_radii          = np.array(list(map(lambda x: x["vdw_radius"], bbox_atoms)))

            normalized_sphere_cords, translation_vectors = normalize_atom_coordinates(_atom_centers)
            voxel_size = 1
            sphere_cords_quantized = np.round( np.array(normalized_sphere_cords / voxel_size) ).astype(int)
            max_values      = np.max(sphere_cords_quantized, axis=0)
            grid_dimensions = max_values + 1

        with open( ptc_data_path(rcsb_id), "r", ) as infile:
            ptc_data = json.load(infile)

        atom_coordinates_by_chain: dict[str, list] = {}
        for atom in bbox_atoms:
            if len(atom["chain_nomenclature"]) < 1:
                # print( "atom ", atom, "has no chain nomenclature", atom["chain_nomenclature"] )
                continue
            if atom["chain_nomenclature"][0] not in atom_coordinates_by_chain:
                atom_coordinates_by_chain[atom["chain_nomenclature"][0]] = []
            atom_coordinates_by_chain[atom["chain_nomenclature"][0]].extend([atom["coord"]])

        ptc_midpoint = np.array(ptc_data["midpoint_coordinates"])

        return ptc_midpoint, atom_coordinates_by_chain, grid_dimensions, translation_vectors


#! For figure only
def DBSCAN_CLUSTERS_particular_eps_minnbrs( dbscan_cluster_dict: dict[int, list], eps, min_nbrs):
    plotter               = pv.Plotter()
    # plotter               = pv.Plotter(shape=(1, 2))

    plotter.subplot(0,0)
    #? Visualize all clusters
    for k, v in dbscan_cluster_dict.items():
        print("Cluster {} has {} points.".format(k, len(v)))
    clusters_palette = dict(zip(range(-1, 60), plt.cm.terrain(np.linspace(0, 1, 60))))
    for k, v in clusters_palette.items():
        clusters_palette[k] = [*v[:3], 0.5]

    combined_cluster_colors = []
    combined_cluster_points = []

    for dbscan_label, coordinates in dbscan_cluster_dict.items():
        combined_cluster_points.extend(coordinates)
        combined_cluster_colors.extend( [clusters_palette[( dbscan_label * 5 )%len(clusters_palette)]   if dbscan_label != -1 else [0, 0, 0, 0.1]] * len(coordinates) )

    ptcloud_all_clusters         = pv.PolyData(combined_cluster_points)
    ptcloud_all_clusters["rgba"] = combined_cluster_colors

    plotter.add_mesh(ptcloud_all_clusters, scalars="rgba", rgb=True, show_scalar_bar=False)

    plotter.add_text('eps: {} \nmin_nbrs: {}'.format(eps, min_nbrs), position='upper_left', font_size=20, shadow=True, font=FONT, color='black')
    # plotter.add_text('Volume: {}'.format(round(mesh_.volume, 3)), position='lower_left', font_size=8, shadow=True, font=FONT, color='black')

    # # ? Visualize selected cluster
    # plotter.subplot(0,1)
    # rgbas_cluster = [[15, 10, 221, 1] for datapoint in selected_cluster]
    # rgbas_positive = np.array([[205, 209, 228, 0.2] for _ in positive_space])
    # combined = np.concatenate([selected_cluster, positive_space])
    # rgbas_combined = np.concatenate([rgbas_cluster, rgbas_positive])

    # point_cloud = pv.PolyData(combined)
    # point_cloud["rgba"] = rgbas_combined
    # plotter.add_mesh(point_cloud, scalars="rgba", rgb=True, show_scalar_bar=False)


    plotter.open_gif("dbscan.gif")
    viewup = [0, 0, 1]
    orbit  = plotter.generate_orbital_path( factor=2.0, n_points=48, shift=0.0, viewup=viewup )
    plotter.orbit_on_path( orbit, write_frames=True, viewup=viewup, step=0.02 )

    plotter.show()

def DBSCAN_CLUSTERS_visualize_largest(positive_space: np.ndarray, dbscan_cluster_dict: dict[int, list], selected_cluster: np.ndarray):
    # plotter               = pv.Plotter()
    plotter               = pv.Plotter(shape=(1, 2))

    plotter.subplot(0,0)
    #? Visualize all clusters
    for k, v in dbscan_cluster_dict.items():
        print("Cluster {} has {} points.".format(k, len(v)))
    clusters_palette = dict(zip(range(-1, 60), plt.cm.terrain(np.linspace(0, 1, 60))))
    for k, v in clusters_palette.items():
        clusters_palette[k] = [*v[:3], 0.5]
    combined_cluster_colors = []
    combined_cluster_points = []
    for dbscan_label, coordinates in dbscan_cluster_dict.items():
        combined_cluster_points.extend(coordinates)
        combined_cluster_colors.extend( [clusters_palette[( dbscan_label * 5 )%len(clusters_palette)]   if dbscan_label != -1 else [0, 0, 0, 0.1]] * len(coordinates) )

    ptcloud_all_clusters         = pv.PolyData(combined_cluster_points)
    ptcloud_all_clusters["rgba"] = combined_cluster_colors

    plotter.add_mesh(ptcloud_all_clusters, scalars="rgba", rgb=True, show_scalar_bar=False)

    # ? Visualize selected cluster
    plotter.subplot(0,1)
    rgbas_cluster = [[15, 10, 221, 1] for datapoint in selected_cluster]
    rgbas_positive = np.array([[205, 209, 228, 0.2] for _ in positive_space])
    combined = np.concatenate([selected_cluster, positive_space])
    rgbas_combined = np.concatenate([rgbas_cluster, rgbas_positive])
    point_cloud         = pv.PolyData(combined)
    point_cloud["rgba"] = rgbas_combined
    plotter.add_mesh(point_cloud, scalars="rgba", rgb=True, show_scalar_bar=False)



    plotter.show()

def plot_multiple_surfaces(rcsb_id:str):

    rcsb_id   = rcsb_id.upper()
    src_taxid = RibosomeAssets(rcsb_id).get_taxids()[0][0]
    taxname   = list( Taxid.get_name(str(src_taxid)).items() )[0][1]
    plotter               = pv.Plotter(shape=(2, 4))


    ptc_midpoint, atom_coordinates_by_chain, grid_dimensions, mean_abs_vectors = retrieve_ptc_and_chain_atoms(rcsb_id)

    #! * FOR EACH SUBPLOT *
    for (i,j) in [(0,0),(0,1),(0,2),(0,3),(1,0),(1,1),(1,2),(1,3)]:
        plotter.subplot(i,j)

        eps, min_nbrs = dbscan_pairs[i*4+j]
        # ? Add mesh to the plotter
        mesh_  = pv.read(custom_cluster_recon_path(rcsb_id, eps, min_nbrs))
        plotter.add_mesh(mesh_, opacity=0.5)


        pprint(list(atom_coordinates_by_chain.keys()))

        for i, ( chain_name, coords ) in enumerate(atom_coordinates_by_chain.items()):
            # print("Plotting " + chain_name, "with index", i ,)
            # ? Adding coordinates to the plotter for each chain( coordinates and color )
            plotter.add_points(
                move_cords_to_normalized_cord_frame(grid_dimensions, mean_abs_vectors, np.array(coords)),
                  point_size               = 6 if chain_name in ["eL39","uL4","uL22","uL23"] else 2 if "rRNA" in chain_name else 4 ,
                  color                    =  'gray' if "rRNA" in chain_name else "cyan" if chain_name == "eL39" else "lightgreen" if chain_name == "uL4" else "gold" if chain_name =="uL22" else 'magenta' if chain_name =='uL23' else CHAIN_LANDMARK_COLORS[i],
                  opacity                  = 0.1 if chain_name not in ["eL39","uL4","uL22",'uL23'] else 1 ,
                  render_points_as_spheres = True ,
            )

        for i, (label, color) in enumerate([( 'eL39','cyan' ),( 'uL4','lightgreen' ),( 'uL22','gold' ), ('uL23','magenta')]):
            offset   = i * 50  # Adjust the offset as needed
            position = (20, 200 - offset, 0)
            plotter.add_text( label, position=position, font_size=20, font=FONT,color=color, shadow=True )


        plotter.add_points( move_cords_to_normalized_cord_frame( grid_dimensions, mean_abs_vectors, np.array([ptc_midpoint]) ), point_size=PTC_PT_SIZE, color="red", render_points_as_spheres=True, )
        plotter.add_text('RCSB_ID:{}'.format(rcsb_id), position='upper_right', font_size=14, shadow=True, font=FONT, color='black')
        plotter.add_text('eps: {} \nmin_nbrs: {}'.format(eps, min_nbrs), position='upper_left', font_size=8, shadow=True, font=FONT, color='black')
        plotter.add_text('Volume: {}'.format(round(mesh_.volume, 3)), position='lower_left', font_size=8, shadow=True, font=FONT, color='black')
        plotter.add_text('{}'.format(taxname), position='lower_right', font_size=8, shadow=True, font=FONT, color='black') 
        

    plotter.show()

def plot_multiple_by_kingdom(kingdom:typing.Literal['bacteria','archaea','eukaryota'], eps:float, min_nbrs:int):

    plotter               = pv.Plotter(shape=(2, 4))


    for (i,j) in [(0,0),(0,1),(0,2),(0,3),(1,0),(1,1),(1,2),(1,3)]:
        try:
            plotter.subplot(i,j)

            rcsb_id   = diagram_tunnels[kingdom][i*4+j]
            src_taxid = RibosomeAssets(rcsb_id).get_taxids()[0][0]
            taxname   = list( Taxid.get_name(str(src_taxid)).items() )[0][1]

            ptc_midpoint, atom_chains_dict, grid_dimensions, translation_vectors = retrieve_ptc_and_chain_atoms(rcsb_id)


            # ? Add mesh to the plotter
            mesh_  = pv.read(custom_cluster_recon_path(rcsb_id, eps, min_nbrs))
            plotter.add_mesh(mesh_, opacity=0.5)


            for i, ( chain_name, coords ) in enumerate(atom_chains_dict.items()):
                # print("Plotting " + chain_name, "with index", i ,)
                # ? Adding coordinates to the plotter for each chain( coordinates and color )
                plotter.add_points(
                    move_cords_to_normalized_cord_frame(grid_dimensions, translation_vectors, np.array(coords)),
                      point_size               = 8 if chain_name in ["eL39","uL4","uL22", "uL23"] else 2 if "rRNA" in chain_name else 4 ,
                      color                    =  'gray' if "rRNA" in chain_name else 'pink' if chain_name =="uL23" else "cyan" if chain_name == "eL39" else "lightgreen" if chain_name == "uL4" else "gold" if chain_name =="uL22" else CHAIN_LANDMARK_COLORS[i],
                      opacity                  = 0.1 if chain_name not in ["eL39","uL4","uL22", "uL23"] else 1 ,
                      render_points_as_spheres = True ,
                )

            # ? Adding PTC coordinatesfor each chain( coordinates and color )
            plotter.add_points( move_cords_to_normalized_cord_frame( grid_dimensions, translation_vectors, np.array([ptc_midpoint]) ), point_size=PTC_PT_SIZE, color="red", render_points_as_spheres=True, )

            #? Add text labels to the plotter
            plotter.add_text('RCSB_ID:{}'.format(rcsb_id), position='upper_right', font_size=14, shadow=True, font=FONT, color='black')
            plotter.add_text('eps: {} \nmin_nbrs: {}'.format(eps, min_nbrs), position='upper_left', font_size=8, shadow=True, font=FONT, color='black')
            plotter.add_text('Volume: {}'.format(round(mesh_.volume, 3)), position='lower_left', font_size=8, shadow=True, font=FONT, color='black')
            plotter.add_text('{}'.format(taxname), position='lower_right', font_size=8, shadow=True, font=FONT, color='black') 

        except Exception as e:
            print(e)
            continue
            

    plotter.show()

def plot_with_landmarks( rcsb_id: str, eps, min_nbrs,poisson_recon_custom_path:str|None=None, ):
    """
    @translation_vectors is a np.ndarray of shape (2,3) where
        - the first row is the means of the coordinate set
        - the second row is the deviations of the normalized coordinate set
        (to be used to reverse the normalization process or to travel to this coordinate frame)

    """

    src_taxid = RibosomeAssets(rcsb_id).get_taxids()[0][0]
    taxname   = list( Taxid.get_name(str(src_taxid)).items() )[0][1]


    ptc_midpoint,atom_coordinates_by_chain, grid_dimensions, mean_abs_vectors= retrieve_ptc_and_chain_atoms(rcsb_id)

    if poisson_recon_custom_path == None:
        poisson_recon = poisson_recon_path(rcsb_id)
    else:
        poisson_recon = poisson_recon_custom_path

    print("Opened poisson recon file at \033[32m{}\033[0m".format(poisson_recon))
    mesh_   = pv.read(poisson_recon)
    plotter = pv.Plotter()
    plotter.add_mesh(mesh_, opacity=1)


    for i, ( chain_name, coords ) in enumerate(atom_coordinates_by_chain.items()):
        # print("Plotting " + chain_name, "with index", i ,)
        # ? Adding coordinates to the plotter for each chain( coordinates and color )
        plotter.add_points(
            # move_cords_to_normalized_cord_frame(grid_dimensions, mean_abs_vectors, np.array(coords)),
                np.array(coords),
              point_size               = 8 if chain_name in ["eL39","uL4","uL22", "uL23"] else 2 if "rRNA" in chain_name else 4 ,
              color                    =  'gray' if "rRNA" in chain_name else "cyan" if chain_name == "eL39" else 'pink' if chain_name=='uL23' else "lightgreen" if chain_name == "uL4" else "gold" if chain_name =="uL22" else CHAIN_LANDMARK_COLORS[i],
              opacity                  = 0.1 if chain_name not in ["eL39","uL4","uL22", 'uL23'] else 1 ,
              render_points_as_spheres = True ,
        )
        #!ALL CHAINS
        # plotter.add_points(
        #     move_cords_to_normalized_cord_frame(grid_dimensions, mean_abs_vectors, np.array(coords)),
        #       point_size               = 7,
        #       color                    =  'gray' if "rRNA" in chain_name else "cyan" if chain_name == "eL39" else "lightgreen" if chain_name == "uL4" else "gold" if chain_name =="uL22" else CHAIN_LANDMARK_COLORS[i],
        #       opacity                  = 0.9,
        #       render_points_as_spheres = True ,
        # )
        # #! RNA ONLY
        # if "rRNA" in chain_name:
        #     plotter.add_points(
        #         move_cords_to_normalized_cord_frame(grid_dimensions, mean_abs_vectors, np.array(coords)),
        #           point_size               = 8 ,
        #           color                    =  'blue' ,
        #           opacity                  = 1,
        #           render_points_as_spheres = True ,
        #     )

    for i, (label, color) in enumerate([( 'eL39','cyan' ),( 'uL4','lightgreen' ),( 'uL22','gold' )]):
        offset   = i * 50  # Adjust the offset as needed
        position = (20, 200 - offset, 0)
        plotter.add_text( label, position=position, font_size=20, font=FONT,color=color, shadow=True )

    plotter.add_points( 
        # move_cords_to_normalized_cord_frame( grid_dimensions, mean_abs_vectors, np.array([ptc_midpoint]) ),
        np.array([ptc_midpoint]),
                        point_size=PTC_PT_SIZE, color="red", render_points_as_spheres=True )

    #!--- Labels ----
    plotter.add_text('RCSB_ID:{}'.format(rcsb_id), position='upper_right', font_size=14, shadow=True, font=FONT, color='black')
    plotter.add_text('eps: {} \nmin_nbrs: {}'.format(eps, min_nbrs), position='upper_left', font_size=8, shadow=True, font=FONT, color='black')
    plotter.add_text('Volume: {}'.format(round(mesh_.volume, 3)), position='lower_left', font_size=8, shadow=True, font=FONT, color='black')
    plotter.add_text('{}'.format(taxname), position='lower_right', font_size=8, shadow=True, font=FONT, color='black') 

    # plotter.open_gif("just_chains.gif")
    # viewup = [1, 0, 0]
    # orbit  = plotter.generate_orbital_path( factor=4.0, n_points=72, shift=2.0, viewup=viewup )
    # plotter.orbit_on_path( orbit, write_frames=True, viewup=viewup, step=0.02 )

    plotter.show(auto_close=False)