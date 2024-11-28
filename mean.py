import mrcfile
import numpy as np

def analyze_mrc(file_path):
    """
    Open and analyze an MRC file to calculate basic statistics.
    
    Parameters:
    file_path (str): Path to the MRC file
    
    Returns:
    dict: Dictionary containing statistical measurements
    """
    try:
        # Open the MRC file in read-only mode
        with mrcfile.open(file_path, mode='r', permissive=True) as mrc:
            # Get the data array
            data = mrc.data
            
            # Calculate statistics
            stats = {
                'mean': float(np.mean(data)),
                'min': float(np.min(data)),
                'max': float(np.max(data)),
                'std': float(np.std(data)),
                'shape': data.shape
            }
            
            return stats
            
    except Exception as e:
        print(f"Error processing MRC file: {str(e)}")
        return None

def print_stats(stats):
    """
    Print the statistics in a formatted way
    """
    if stats is not None:
        print("\nMRC File Statistics:")
        print("-" * 20)
        print(f"Mean value: {stats['mean']:.4f}")
        print(f"Minimum value: {stats['min']:.4f}")
        print(f"Maximum value: {stats['max']:.4f}")
        print(f"Standard deviation: {stats['std']:.4f}")
        print(f"Map dimensions: {stats['shape']}")
    else:
        print("No statistics available")

# Example usage
if __name__ == "__main__":
    file_path = "EMD-6058.map"  # Replace with your MRC file path
    stats = analyze_mrc(file_path)
    print_stats(stats)