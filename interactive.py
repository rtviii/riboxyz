from ribctl.lib.libtax import get_ncbi
from ribctl.ribosome_ops import RibosomeOps
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import umap
import os
import sys
from dash import Dash, dcc, html, Input, Output, State, callback_context, no_update
import dash_bootstrap_components as dbc
import dash_vtk
from dash_vtk.utils import to_mesh_state
import base64
from io import BytesIO
import vtk

# File paths
ids_file = "wenjun_umap/ids.txt"
dist_file = "wenjun_umap/gw_distances.npy"
meshes_dir = "wenjun_data_meshes"

# Initialize NCBITaxa database
ncbi = get_ncbi()

def get_taxonomy_for_id(ribosome_id):
    try:
        ribosome_metadata = RibosomeOps(ribosome_id).profile
        tax_id = ribosome_metadata.src_organism_ids[0]
        
        # Get lineage info using NCBITaxa
        lineage = ncbi.get_lineage(tax_id)
        names = ncbi.get_taxid_translator(lineage)
        ranks = ncbi.get_rank(lineage)
        
        # Create taxonomy dictionary
        taxonomy = {ranks[taxid]: names[taxid] for taxid in lineage}
        taxonomy['taxid'] = tax_id
        
        # Get specific taxonomic levels
        superkingdom_id = next((id for id, rank in ranks.items() if rank == 'superkingdom'), None)
        taxonomy['superkingdom'] = names.get(superkingdom_id, 'Unknown')
        
        phylum_id = next((id for id, rank in ranks.items() if rank == 'phylum'), None)
        taxonomy['phylum'] = names.get(phylum_id, 'Unknown')
        
        species_id = next((id for id, rank in ranks.items() if rank == 'species'), None)
        taxonomy['species'] = names.get(species_id, 'Unknown')
        
        return taxonomy
    except Exception as e:
        print(f"Error getting taxonomy for {ribosome_id}: {e}")
        return {'taxid': 0, 'superkingdom': 'Unknown', 'phylum': 'Unknown', 'species': 'Unknown'}

# Load the data
K_GW_full = np.load(dist_file)
K_GW_full = K_GW_full/K_GW_full.max()

# Load the list of IDs
with open(ids_file, 'r') as f:
    ribosome_ids = [line.strip() for line in f.readlines()]

# Ensure we have the right number of IDs
assert len(ribosome_ids) == K_GW_full.shape[0], "Number of IDs doesn't match matrix dimensions"

# Fit UMAP with 2 components
fit = umap.UMAP(metric="precomputed", n_components=2, random_state=52)
mapper = fit.fit_transform(K_GW_full)

# Create a DataFrame with the UMAP coordinates, IDs, and taxonomy
df = pd.DataFrame({
    'UMAP_1': mapper[:, 0],
    'UMAP_2': mapper[:, 1],
    'ID': ribosome_ids
})

# Add taxonomy data to the DataFrame
print("Loading taxonomy data...")
taxonomy_data = []
for ribosome_id in ribosome_ids:
    taxonomy_data.append(get_taxonomy_for_id(ribosome_id))

# Extract specific taxonomy levels to the DataFrame
df['Taxid'] = [t['taxid'] for t in taxonomy_data]
df['Superkingdom'] = [t['superkingdom'] for t in taxonomy_data]
df['Phylum'] = [t['phylum'] for t in taxonomy_data]
df['Species'] = [t['species'] for t in taxonomy_data]

# Define taxonomy levels for coloring
taxonomy_levels = ['Superkingdom', 'Phylum', 'Species']

# Function to read a PLY file and convert it to a format for dash_vtk
def read_ply_to_vtk_mesh(mesh_id):
    mesh_path = os.path.join(meshes_dir, f'{mesh_id}_NPET_MESH.ply')
    
    if not os.path.exists(mesh_path):
        print(f"Mesh file not found: {mesh_path}")
        return None
    
    try:
        # Read the PLY file using VTK
        reader = vtk.vtkPLYReader()
        reader.SetFileName(mesh_path)
        reader.Update()
        
        # Convert to mesh state for dash_vtk
        return to_mesh_state(reader.GetOutput())
    except Exception as e:
        print(f"Error reading PLY file for {mesh_id}: {e}")
        return None

# Cache for mesh data to avoid reloading
mesh_cache = {}

# Create the Dash app
app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

# Define the VTK view component
def get_vtk_view(mesh_id=None):
    if mesh_id is None:
        # Return empty view if no mesh is selected
        return html.Div("Select a point to view the mesh")
    
    # Check if we have the mesh in cache
    if mesh_id in mesh_cache:
        mesh_state = mesh_cache[mesh_id]
    else:
        mesh_state = read_ply_to_vtk_mesh(mesh_id)
        if mesh_state:
            mesh_cache[mesh_id] = mesh_state
    
    if mesh_state is None:
        return html.Div(f"Mesh not found for {mesh_id}")
    
    # Unique key forces a complete re-render, which resets the camera
    unique_key = f"vtk-view-{mesh_id}-{np.random.randint(10000)}"
    
    return dash_vtk.View(
        id=unique_key,
        children=[
            dash_vtk.GeometryRepresentation([
                dash_vtk.Mesh(state=mesh_state)
            ], 
            property={"color": [0.5, 0.8, 1.0], "representation": 2}),
        ], 
        background=[0.9, 0.9, 0.9],
        cameraPosition=[0, 0, 1],
    )

# Define the layout
app.layout = dbc.Container([
    dbc.Row([
        dbc.Col([
            html.H4("Visualization Settings"),
            html.Label("Color By:"),
            dcc.Dropdown(
                id='color-by',
                options=[{'label': level, 'value': level} for level in taxonomy_levels],
                value='Phylum',
                clearable=False
            ),
            html.Div(id='selected-point-info', className="mt-3"),
            
            html.Div([
                html.H5("3D Mesh Visualization", className="mt-4"),
                html.Div(id='vtk-container', style={'height': '400px', 'width': '100%'}),
                html.Div([
                    dbc.Button("Reset Camera", id="reset-camera-btn", color="secondary", size="sm", className="mt-2"),
                    dbc.Button("Toggle Wireframe", id="toggle-wireframe-btn", color="secondary", size="sm", className="mt-2 ms-2"),
                ], className="d-flex"),
            ], className="mt-3"),
            
        ], width=4),
        
        dbc.Col([
            dcc.Graph(
                id='umap-scatter',
                config={'displayModeBar': True},
                style={'height': '80vh'}
            ),
        ], width=8),
    ]),
    
    # Hidden div to store the currently selected ID
    dcc.Store(id='selected-id-store'),
    
    # Load indicator
    dbc.Spinner(html.Div(id="loading-output"), type="grow", fullscreen=True, 
                fullscreen_style={"backgroundColor": "rgba(0, 0, 0, 0.3)"}),
], fluid=True)

# Callback to update the scatter plot
@app.callback(
    Output('umap-scatter', 'figure'),
    [Input('color-by', 'value')]
)
def update_graph(color_by):
    # Create the scatter plot colored by the selected taxonomy level
    fig = px.scatter(
        df, 
        x='UMAP_1', 
        y='UMAP_2', 
        color=color_by,
        hover_data=['ID', 'Superkingdom', 'Phylum', 'Species'],
        labels={'UMAP_1': 'UMAP 1', 'UMAP_2': 'UMAP 2'},
        title=f'UMAP Projection Colored by {color_by}',
        custom_data=['ID', 'Superkingdom', 'Phylum', 'Species', 'Taxid']  # Include all data we need for callbacks
    )
    
    fig.update_traces(marker=dict(size=10))
    fig.update_layout(clickmode='event+select')
    
    return fig

# Global variables to track visualization state
wireframe_mode = False

# Callback for click events
@app.callback(
    [Output('selected-point-info', 'children'),
     Output('vtk-container', 'children'),
     Output('selected-id-store', 'data'),
     Output('loading-output', 'children')],
    [Input('umap-scatter', 'clickData'),
     Input('reset-camera-btn', 'n_clicks'),
     Input('toggle-wireframe-btn', 'n_clicks')],
    [State('selected-id-store', 'data')]
)
def handle_interaction(clickData, reset_btn, toggle_btn, current_id):
    global wireframe_mode
    
    # Determine which input triggered the callback
    ctx = callback_context
    if not ctx.triggered:
        trigger_id = 'no-trigger'
    else:
        trigger_id = ctx.triggered[0]['prop_id'].split('.')[0]
    
    # Handle reset camera button
    if trigger_id == 'reset-camera-btn' and current_id:
        # Force a re-render of the VTK view with the same mesh to reset camera
        vtk_view = get_vtk_view(current_id)
        return no_update, vtk_view, current_id, ''
    
    # Handle toggle wireframe button
    if trigger_id == 'toggle-wireframe-btn' and current_id:
        wireframe_mode = not wireframe_mode
        
        # Re-render the VTK view with updated representation
        mesh_state = mesh_cache.get(current_id, read_ply_to_vtk_mesh(current_id))
        if mesh_state is None:
            return no_update, html.Div(f"Mesh not found for {current_id}"), current_id, ''
        
        representation = 1 if wireframe_mode else 2  # 1=wireframe, 2=surface
        
        unique_key = f"vtk-view-{current_id}-{np.random.randint(10000)}"
        vtk_view = dash_vtk.View(
            id=unique_key,
            children=[
                dash_vtk.GeometryRepresentation([
                    dash_vtk.Mesh(state=mesh_state)
                ], 
                property={"color": [0.5, 0.8, 1.0], "representation": representation}),
            ], 
            background=[0.9, 0.9, 0.9],
            cameraPosition=[0, 0, 1],
        )
        
        return no_update, vtk_view, current_id, ''
    
    # Handle click on scatter plot
    if trigger_id == 'umap-scatter' and clickData:
        # Extract data from the clicked point
        point_data = clickData['points'][0]['customdata']
        ribosome_id = point_data[0]
        superkingdom = point_data[1]
        phylum = point_data[2]
        species = point_data[3]
        taxid = point_data[4]
        
        # Only update if a different point is clicked
        if ribosome_id == current_id:
            return no_update, no_update, current_id, ''
        
        # Create info display
        info_div = [
            html.H5(f"Selected Point: {ribosome_id}"),
            html.P(f"Superkingdom: {superkingdom}"),
            html.P(f"Phylum: {phylum}"),
            html.P(f"Species: {species}"),
            html.P(f"Taxid: {taxid}"),
        ]
        
        # Reset wireframe mode for new selection
        wireframe_mode = False
        
        # Create VTK visualization
        vtk_view = get_vtk_view(ribosome_id)
        
        return info_div, vtk_view, ribosome_id, ''
    
    # Default case - no point selected yet
    if not current_id:
        return (
            html.Div([html.H5("Click on a point to see details")]), 
            html.Div("Click on a point to view the mesh"), 
            None, 
            ''
        )
    
    # No change needed
    return no_update, no_update, current_id, ''

# Run the app
if __name__ == '__main__':
    app.run(debug=True)