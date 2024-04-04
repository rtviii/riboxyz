# import numpy as np


# def compute_naive_map_to_map_distance_matrix(self_distance_matrices, loss_fun):
#     n_maps = 100
#     results = np.zeros((n_maps,n_maps))

#     for idx_1 in range(n_maps):
#         for idx_2 in range(n_maps):
#             if idx_1 < idx_2: 
#                 self_distance_matrix_1 = self_distance_matrices[idx_1]
#                 self_distance_matrix_2 = self_distance_matrices[idx_2]
#                 dist = map_to_map_distance(self_distance_matrix_1, self_distance_matrix_2, loss_fun=loss_fun) # returns a scalar
#                 results[idx_1, idx_2] = dist
#                 results[idx_2, idx_1] = dist
#     return results


import numpy as np

def compare_dicts(dict1, dict2, compare_func):
    """Compare two dictionaries using the provided comparison function."""
    return compare_func(dict1, dict2)

def task(task_id, dict1, dict2, compare_func):
    """Example task function."""
    print(f"Task {task_id} started")
    result = compare_dicts(dict1, dict2, compare_func)
    print(f"Task {task_id} completed with result {result}")

# Example comparison function
def custom_compare(dict1, dict2):
    """Custom comparison function."""
    return dict1 == dict2

# Define two lists of integers
n_maps = 100
list1 = list(np.linspace(0, 1, n_maps))
list2 = list(np.linspace(0, 1, n_maps))

# Generate the Cartesian product using list comprehension
cartesian_product = [(x, y) for x, y in product(list1, list2)]

# Sample dictionaries
dict1 = {'a': 1, 'b': 2, 'c': 3}
dict2 = {'a': 1, 'b': 2, 'c': 3}

# Create a thread pool executor with a specified number of threads
with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
    # Submit the tasks to the executor
    future_tasks = {executor.submit(task, index, dict1, dict2, custom_compare): index for index in range(num_tasks)}
    
    # Wait for all tasks to complete
    for future in concurrent.futures.as_completed(future_tasks):
        # Get the task id associated with the completed task
        task_id = future_tasks[future]
        try:
            result = future.result()
        except Exception as e:
            print(f"Task {task_id} generated an exception: {e}")
        else:
            print(f"Task {task_id} returned result: {result}")

# ---------------

import concurrent.futures
from itertools import product
from functools import partial

nmaps = len(self_distance_matrices)
sink_container = np.zeros(nmaps, nmaps)

def distance_task(global_results_matrix:list, self_distance_matrices:list, idx_1:int, idx_2:int, loss_function):
    distance = loss_function(self_distance_matrices[idx_1], self_distance_matrices[idx_2])
    global_results_matrix[idx_1, idx_2] = distance
    return distance

def distances_threadpool(sink_container:list, ...):
  map_indices = list(product(range(nmaps), range(nmaps))) # Cartesian product of map indices
  map_indices = list(filter( lambda x,y : x >= y, map_indices)) # Filter out pairs of indices where the first index is less than the second index

  with concurrent.futures.ThreadPoolExecutor(max_workers=CORES-2) as executor:
        futures = []
        for (idx1, idx2) in map_indices:
            partial_task = partial(distance_task_task, global_results_matrix, self_distance_matrices, idx_1, idx_2, loss_function)
            future = executor.submit(partial_task)
            futures.append(future)

        concurrent.futures.wait(futures)

  return sink_container