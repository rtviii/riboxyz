# Make scripts executable
chmod +x __scripts/npet2_viz_enhanced.py
chmod +x __scripts/npet2_slice_viewer.py

# Basic viz: show shell, ROI, and empty points
python __scripts/npet2_viz_enhanced.py --rcsb 7K00 --shell --roi --points1

# Show occupancy (occupied vs empty voxels) for level 0
python __scripts/npet2_viz_enhanced.py --rcsb 7K00 --occupancy0 --empty_mask0 --grid_bounds0 --downsample 2

# Show full pipeline: grid bounds, ROI, clusters, refined result
python __scripts/npet2_viz_enhanced.py --rcsb 7K00 \
    --shell --roi --grid_bounds0 --grid_bounds1 \
    --stage50_coarse --refined --downsample 3

# Slice viewer for occupancy grid
python __scripts/npet2_slice_viewer.py --rcsb 7K00 --grid occupancy_grid_level_0

# Slice viewer for empty mask
python __scripts/npet2_slice_viewer.py --rcsb 7K00 --grid empty_mask_level_1