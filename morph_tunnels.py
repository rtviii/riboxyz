import numpy as np
import pyvista as pv
from scipy.spatial.distance import cdist
from itertools import permutations
import imageio
import os
import traceback
import ot


from ribctl import RIBETL_DATA
tunnelset = ["6OFX", "6OG7", "6OGF", "6OGG", "6OGI"]
mesh_files = list(
    map(lambda x: os.path.join(RIBETL_DATA, x, "{}_NPET_MESH_ASCII.ply".format(x)), tunnelset)
)


# --- Morphing Parameters ---
NUM_SAMPLES = 10000
NUM_INTERPOLATION_STEPS = 40
GIF_FILENAME = "morphing_ot.gif"
GIF_FPS = 25
OT_REGULARIZATION = 0.005


def sample_surface_manually(mesh, n_samples):
    vertices = mesh.points
    try:
        faces = mesh.faces.reshape(-1, 4)[:, 1:]
    except ValueError:
        faces = mesh.regular_faces
    triangles = vertices[faces]
    v0, v1, v2 = triangles[:, 0], triangles[:, 1], triangles[:, 2]
    areas = 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0), axis=1)
    total_area = np.sum(areas)
    if total_area == 0:
        indices = np.random.choice(len(vertices), size=n_samples, replace=True)
        return vertices[indices]
    probabilities = areas / total_area
    chosen_triangle_indices = np.random.choice(
        len(faces), size=n_samples, p=probabilities
    )
    r1, r2 = np.random.rand(n_samples, 1), np.random.rand(n_samples, 1)
    u, v, w = 1.0 - np.sqrt(r1), np.sqrt(r1) * (1.0 - r2), np.sqrt(r1) * r2
    chosen_triangles = triangles[chosen_triangle_indices]
    points = (
        u * chosen_triangles[:, 0]
        + v * chosen_triangles[:, 1]
        + w * chosen_triangles[:, 2]
    )
    return points


def load_and_sample_meshes(files, num_samples):
    print(f"Loading and sampling {num_samples} points from each mesh...")
    point_clouds = []
    for file in files:
        if not os.path.exists(file):
            raise FileNotFoundError(f"Error: '{file}' not found.")
        mesh = pv.read(file)
        if isinstance(mesh, pv.MultiBlock):
            mesh = mesh[0]
        print(f"  -> Sampling surface for {os.path.basename(file)}...")
        point_clouds.append(sample_surface_manually(mesh, num_samples))
    return point_clouds


def calculate_distance_matrix(point_clouds):
    """Calculates the Wasserstein distance between each pair of point clouds."""
    print("Calculating pairwise Wasserstein distances to find optimal path...")
    num_clouds = len(point_clouds)
    dist_matrix = np.zeros((num_clouds, num_clouds))

    # Get the number of points and define a uniform distribution
    num_points = point_clouds[0].shape[0]
    a = np.ones(num_points) / num_points
    b = np.ones(num_points) / num_points

    for i in range(num_clouds):
        for j in range(i, num_clouds):
            if i == j:
                continue

            # Manually compute the cost matrix (squared Euclidean distance)
            M = cdist(point_clouds[i], point_clouds[j], "sqeuclidean")

            # Provide the distributions (a, b) and the cost matrix (M)
            dist = ot.emd2(a, b, M)

            dist_matrix[i, j] = dist_matrix[j, i] = dist

    return dist_matrix


def find_optimal_path(dist_matrix):
    print("Finding the most logical sequence for morphing...")
    num_nodes = len(dist_matrix)
    if num_nodes < 2:
        return list(range(num_nodes))
    min_path, min_dist = None, float("inf")
    for p in permutations(range(1, num_nodes)):
        path = [0] + list(p)
        current_dist = sum(
            dist_matrix[path[i], path[i + 1]] for i in range(num_nodes - 1)
        )
        if current_dist < min_dist:
            min_dist, min_path = current_dist, path
    return min_path


def create_morph_animation_ot(point_clouds, path, steps, reg, output_filename, fps):
    """Generates a GIF using Optimal Transport for interpolation."""
    print(f"\nGenerating animation with Optimal Transport (regularization={reg})...")
    plotter = pv.Plotter(off_screen=True, window_size=[800, 800])
    plotter.background_color = "white"
    plotter.add_points(
        point_clouds[path[0]], color="teal", point_size=5, render_points_as_spheres=True
    )
    plotter.view_isometric()
    plotter.camera.zoom(1.4)
    frames = []

    num_points = point_clouds[0].shape[0]
    a = np.ones(num_points) / num_points
    b = np.ones(num_points) / num_points

    for i in range(len(path) - 1):
        start_cloud = point_clouds[path[i]]
        end_cloud = point_clouds[path[i + 1]]

        print(
            f"  Calculating OT Plan: {os.path.basename(mesh_files[path[i]])} -> {os.path.basename(mesh_files[path[i+1]])}"
        )

        M = cdist(start_cloud, end_cloud, "sqeuclidean")
        M /= M.max()

        gamma = ot.sinkhorn(a, b, M, reg)

        print(f"  Morphing frames...")
        for t in np.linspace(0, 1, steps):
            mapped_end_cloud = (gamma / a[:, None]) @ end_cloud
            interpolated = (1 - t) * start_cloud + t * mapped_end_cloud

            plotter.update_coordinates(interpolated, render=False)
            plotter.render()
            frames.append(plotter.screenshot(None, return_img=True))

    frames.extend([frames[-1]] * (fps // 2))
    print(f"Saving GIF '{output_filename}' with {len(frames)} frames at {fps} FPS...")
    imageio.mimsave(output_filename, frames, fps=fps)
    print("✨ Done!")


# --- Main execution ---
if __name__ == "__main__":
    try:
        clouds = load_and_sample_meshes(mesh_files, NUM_SAMPLES)
        distance_matrix = calculate_distance_matrix(clouds)
        optimal_path = find_optimal_path(distance_matrix)

        ordered_files = [os.path.basename(mesh_files[i]) for i in optimal_path]
        print(f"\nOptimal morphing order determined:\n{' -> '.join(ordered_files)}")

        create_morph_animation_ot(
            clouds,
            optimal_path,
            NUM_INTERPOLATION_STEPS,
            OT_REGULARIZATION,
            GIF_FILENAME,
            GIF_FPS,
        )

    except FileNotFoundError as e:
        print(e)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        traceback.print_exc()
