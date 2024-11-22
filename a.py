import numpy as np

def compare_meshgrid_indexing():
    # Create simple coordinate arrays
    x = np.array([1, 2, 3])  # 3 points
    y = np.array([10, 20])   # 2 points
    
    print("Input arrays:")
    print("x:", x)
    print("y:", y)
    print("\n" + "="*50 + "\n")
    
    # Using 'ij' indexing (matrix indexing)
    Xi, Yi = np.meshgrid(x, y, indexing='ij')
    print("'ij' indexing:")
    print("\nXi shape:", Xi.shape)  # (3, 2)
    print("Xi array:")
    print(Xi)
    print("\nYi array:")
    print(Yi)
    print("\nPoints (as row, col pairs):")
    for i in range(Xi.shape[0]):
        for j in range(Xi.shape[1]):
            print(f"({i},{j}) -> ({Xi[i,j]}, {Yi[i,j]})")
    
    print("\n" + "="*50 + "\n")
    
    # Using 'xy' indexing (cartesian indexing)
    Xx, Yx = np.meshgrid(x, y, indexing='xy')
    print("'xy' indexing:")
    print("\nXx shape:", Xx.shape)  # (2, 3)
    print("Xx array:")
    print(Xx)
    print("\nYx array:")
    print(Yx)
    print("\nPoints (as row, col pairs):")
    for i in range(Xx.shape[0]):
        for j in range(Xx.shape[1]):
            print(f"({i},{j}) -> ({Xx[i,j]}, {Yx[i,j]})")

compare_meshgrid_indexing()