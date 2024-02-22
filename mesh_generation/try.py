import concurrent.futures
from functools import partial

import numpy as np

from mesh_generation.voxelize import get_sphere_indices_voxelized

# def your_task(container, task_index, additional_arg1, additional_arg2):
#     point_list.append((task_index, additional_arg1, additional_arg2))


# Function to be executed by each worker
def sphere_task(container_sink:list, atom_center_coordinate:np.ndarray, vdw_R=2):
    result = get_sphere_indices_voxelized(atom_center_coordinate, 2)
    container_sink.append(result)
    return result



def expand_atomcenters_to_spheres_threadpool(sink_container:list, sphere_sources:list[tuple[np.ndarray, int]]):
  with concurrent.futures.ThreadPoolExecutor(max_workers=20) as executor:

        futures = []
        for (center_coordinate, vdw_R) in sphere_sources:
            partial_task = partial(sphere_task, sink_container, center_coordinate, vdw_R)
            future = executor.submit(partial_task)
            futures.append(future)

        concurrent.futures.wait(futures)

  return sink_container

sink = []
k = expand_atomcenters_to_spheres_threadpool(sink,sphere_sources)
