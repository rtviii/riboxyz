
# Old
def __np_workflow(C):
    rescaled_coordinates, dim = normalize_atom_coordinates(C)
    x, y, z = np.indices((dim, dim, dim))
    xc = midpoints(x)
    yc = midpoints(y)
    zc = midpoints(z)
    filled = xc + yc + zc < -1
    filled = visualize_source_coordinates(filled, rescaled_coordinates)
    return filled