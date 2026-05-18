import dash
from dash import dcc, html, Input, Output
import plotly.graph_objects as go
import networkx as nx
import numpy as np
import pandas as pd
import re
from collections import defaultdict
import pickle
import os

# =============== UTILITY FUNCTIONS ====================

def psamm_to_cobra_id(s: str) -> str:
    """Apply character replacements for package compatibility."""
    renames = {
        ".": "_DOT_", "(": "_LPAREN_", ")": "_RPAREN_", "-": "_DASH_",
        "[": "_LSQBKT", "]": "_RSQBKT", ",": "_COMMA_", ":": "_COLON_",
        ">": "_GT_", "<": "_LT", "/": "_FLASH", "\\": "_BSLASH",
        "+": "_PLUS_", "=": "_EQ_", " ": "_SPACE_", "'": "_SQUOT_",
        '"': "_DQUOT_",
    }
    result = s
    for char, replacement in renames.items():
        result = result.replace(char, replacement)
    return result

def parse_formula(formula: str) -> dict:
    """Parse chemical formula string to element counts."""
    pattern = r'([A-Z][a-z]?)(\d*)'
    matches = re.findall(pattern, formula)
    
    result = defaultdict(int)
    for element, count in matches:
        count = int(count) if count else 1
        result[element] += count
    
    return dict(result)

def interpolate_color(value, vmin, vmax, vcenter=0):
    """
    Interpolate between blue-grey-red for a diverging colormap.
    Returns rgba string.
    """
    # Normalize value to [-1, 1] range
    if value < vcenter:
        if vmin == vcenter:
            norm_value = 0
        else:
            norm_value = -1 * (vcenter - value) / (vcenter - vmin)
    else:
        if vmax == vcenter:
            norm_value = 0
        else:
            norm_value = (value - vcenter) / (vmax - vcenter)
    
    # Clamp to [-1, 1]
    norm_value = max(-1, min(1, norm_value))
    
    # Blue to grey to red interpolation
    if norm_value < 0:
        # Blue (33, 102, 172) to grey (211, 211, 211)
        t = abs(norm_value)
        r = int(211 * (1 - t) + 33 * t)
        g = int(211 * (1 - t) + 102 * t)
        b = int(211 * (1 - t) + 172 * t)
    else:
        # Grey (211, 211, 211) to red (178, 24, 43)
        t = norm_value
        r = int(211 * (1 - t) + 178 * t)
        g = int(211 * (1 - t) + 24 * t)
        b = int(211 * (1 - t) + 43 * t)
    
    return f'rgba({r}, {g}, {b}, 1)'


def load_labels_from_file(filepath):
    """
    Load labels from a TSV/CSV file.
    Expected columns: text, x, y, font_size
    Optional columns: font_color, font_weight, angle
    """
    # Try to determine separator
    with open(filepath, 'r') as f:
        first_line = f.readline()
        sep = '\t' if '\t' in first_line else ','
    
    df = pd.read_csv(filepath, sep=sep)
    labels = []
    for _, row in df.iterrows():
        # Replace \n with <br> for HTML line breaks
        text = str(row['text']).replace('\\n', '<br>')
        
        label = {
            'text': text,
            'x': row['x'],
            'y': row['y'],
            'font_size': row.get('font_size', 14),
            'font_color': row.get('font_color', '#2c3e50'),
            'font_weight': row.get('font_weight', 'normal'),
            'angle': row.get('angle', 0)  # Default to 0 if not specified
        }
        labels.append(label)
    return labels

# =============== PARAMETERS ====================
NODE_SIZE_HAS_DATA = 105
NODE_SIZE_NO_DATA = 25
NODE_EDGEWIDTH = 1
NODE_EDGE_COLOR = "black"
NODE_NOTSIG_COLOR = "lightgrey"
NODE_NODATA_COLOR = "white"

EDGE_WIDTH_MIN = 2
EDGE_WIDTH_SCALING = 10
EDGE_COLOR_BOUNDS = (-2, 2)
NODE_COLOR_BOUNDS = (-5, 5)
DEFAULT_EDGE_COLOR = "lightgrey"

COLOR_METABOLITE_DATA = True

NODE_SIZE_SCALE_HAS_DATA = 0.2
NODE_SIZE_SCALE_NO_DATA = 0.5
EDGE_WIDTH_SCALE = 1.5

PADDING = 0.1

# Reactions to ignore for width calculation
IGNORE_WIDTH = ["SerK", "R00582"]

# Strain mappings
STRAIN_DISPLAY_TO_INTERNAL = {
    'BDO-ALS': 'MW698',
    'BDO': 'MW697',
    'Parent-ALS': 'MW268',
    'Parent-COM': 'COM1c'
}

STRAIN_INTERNAL_TO_DISPLAY = {v: k for k, v in STRAIN_DISPLAY_TO_INTERNAL.items()}

# Strains with metabolomics data
STRAINS_WITH_METABOLOMICS = ['MW698', 'COM1c']  # BDO-ALS and Parent-COM

# Normalization options
NORMALIZE_OPTIONS = {
    "Biomass": "Biomass_Pfu",
    "Carbon uptake": "cellobioseABC",
    "None": None
}

# =============== LOAD STATIC DATA ====================

def load_network_data():
    """Load and process all network data."""
    
    # Load graph
    G = nx.read_graphml("/app/data/GRAPH.graphml")
    
    # Load node positions
    node_table = pd.read_csv("/app/data/NODE_POSITIONS.csv", index_col=0)
    node_table.index = node_table.index.astype("str")
    
    # Set node positions (flip y-axis for plotting)
    pos = {}
    for node in G.nodes():
        if node in node_table.index:
            G.nodes[node]['x'] = node_table.loc[node, 'x']
            G.nodes[node]['y'] = node_table.loc[node, 'y']
            pos[node] = (node_table.loc[node, 'x'], -node_table.loc[node, 'y'])
    
    # Load flux data
    with open('/app/data/samples_dict.pkl', 'rb') as f:
        samples = pickle.load(f)
    
    # Load labels
    labels = load_labels_from_file('/app/data/labels.tsv')
    
    return G, pos, samples, labels

# Load static data once at startup
G_ORIGINAL, pos, samples, LABELS = load_network_data()

# Calculate plot bounds once
x_values = [pos[node][0] for node in G_ORIGINAL.nodes()]
y_values = [pos[node][1] for node in G_ORIGINAL.nodes()]
x_min, x_max = min(x_values), max(x_values)
y_min, y_max = min(y_values), max(y_values)

x_range = x_max - x_min
y_range = y_max - y_min
x_min -= PADDING * x_range
x_max += PADDING * x_range
y_min -= PADDING * y_range
y_max += PADDING * y_range

# =============== DATA PROCESSING FUNCTIONS ====================

def calculate_flux_changes(samples, group1, group2, normalize_by):
    """Calculate flux fold changes between groups."""
    
    # Calculate normalized mean fluxes
    if normalize_by:
        group1_norm = samples[group1].div(samples[group1][normalize_by], axis=0).median()
        group2_norm = samples[group2].div(samples[group2][normalize_by], axis=0).median()
    else:
        group1_norm = samples[group1].median()
        group2_norm = samples[group2].median()
    
    # Calculate log2 fold change with small epsilon to avoid log(0)
    epsilon = 1e-9
    flux_log2fc = np.log2(
        (np.abs(group1_norm) + epsilon) / (np.abs(group2_norm) + epsilon)
    ).replace([np.inf, -np.inf], np.nan).dropna()
    
    return flux_log2fc, samples[group1]

def process_edges(G, flux_log2fc, group1_samples):
    """Process edge styling, colors, and widths."""
    
    # Collect edge properties
    edge_styles = []
    edge_values = []
    edge_has_value = []
    edge_widths_raw = []
    edge_width_has_value = []
    
    for u, v, edge_data in G.edges(data=True):
        # Style
        if edge_data.get('edge_type') == 'link':
            edge_styles.append('dotted')
        else:
            edge_styles.append('solid')
        
        # Values for coloring
        if 'reaction' in edge_data:
            rxn_id = psamm_to_cobra_id(edge_data['reaction'])
            if rxn_id in flux_log2fc.index:
                edge_values.append(flux_log2fc[rxn_id])
                edge_has_value.append(True)
            else:
                edge_values.append(0)
                edge_has_value.append(False)
        else:
            edge_values.append(0)
            edge_has_value.append(False)
        
        # Widths based on carbon flux
        if ('reaction' in edge_data) and ('element_transfer' in edge_data):
            rxn_id = psamm_to_cobra_id(edge_data['reaction'])
            molecules = parse_formula(edge_data['element_transfer'])
            
            if rxn_id not in IGNORE_WIDTH and rxn_id in group1_samples.columns:
                carbon_flux = molecules.get("C", 0) * np.abs(group1_samples[rxn_id].median())
                edge_widths_raw.append(carbon_flux)
                edge_width_has_value.append(True)
            else:
                edge_widths_raw.append(0)
                edge_width_has_value.append(False)
        else:
            edge_widths_raw.append(0)
            edge_width_has_value.append(False)
    
    # Color mapping
    vmin, vmax = EDGE_COLOR_BOUNDS
    edge_colors = []
    for i, has_val in enumerate(edge_has_value):
        if has_val:
            edge_colors.append(interpolate_color(edge_values[i], vmin, vmax, vcenter=0))
        else:
            edge_colors.append(DEFAULT_EDGE_COLOR)
    
    # Width normalization
    edge_widths_array = np.array(edge_widths_raw)
    if any(edge_width_has_value):
        valid_widths = edge_widths_array[edge_width_has_value]
        width_min = valid_widths.min()
        width_max = valid_widths.max()
        width_range = width_max - width_min if width_max != width_min else 1
        
        edge_widths = []
        for i, w in enumerate(edge_widths_raw):
            if edge_width_has_value[i]:
                normalized = (w - width_min) / width_range
                edge_widths.append(EDGE_WIDTH_MIN + normalized * EDGE_WIDTH_SCALING)
            else:
                edge_widths.append(EDGE_WIDTH_MIN)
    else:
        edge_widths = [EDGE_WIDTH_MIN] * len(edge_widths_raw)
    
    return {
        'styles': edge_styles,
        'values': edge_values,
        'has_value': edge_has_value,
        'colors': edge_colors,
        'widths': edge_widths,
        'widths_raw': edge_widths_raw,
        'width_has_value': edge_width_has_value
    }

def process_nodes(G, group1_internal, group2_internal):
    """Process node colors and sizes."""
    
    # Check if we have metabolomics data for both strains
    has_metabolomics = (group1_internal in STRAINS_WITH_METABOLOMICS and 
                       group2_internal in STRAINS_WITH_METABOLOMICS)
    
    # Determine if we need to invert the fold changes
    # Original data is BDO-ALS vs Parent-COM (MW698 vs COM1c)
    # If comparison is reversed, invert the fold changes
    invert_fc = False
    same_group = False
    
    if has_metabolomics:
        if group1_internal == group2_internal:
            same_group = True
        elif group1_internal == 'COM1c' and group2_internal == 'MW698':
            # Parent-COM vs BDO-ALS (reversed from original)
            invert_fc = True
    
    def flip_significance(sig):
        """Flip Up/Down significance when inverting fold changes."""
        if sig == 'Up':
            return 'Down'
        elif sig == 'Down':
            return 'Up'
        else:
            return sig  # 'NS' stays as 'NS'
    
    node_values = []
    node_has_value = []
    node_sizes = []
    node_sig = []
    node_pvals = []
    
    for node in G.nodes():
        node_data = G.nodes[node]
        
        # Only use metabolomics data if both strains have it
        if has_metabolomics and 'logFC' in node_data and node_data['logFC'] != "None":
            if same_group:
                # Same group comparison: FC = 0, p-value = 1
                node_values.append(0)
                node_has_value.append(True)
                node_sizes.append(NODE_SIZE_HAS_DATA)
                node_sig.append('NS')
                node_pvals.append(1.0)
            else:
                # Normal or inverted comparison
                fc_value = node_data['logFC']
                sig_value = node_data.get('significance', 'NS')
                
                if invert_fc:
                    fc_value = -fc_value
                    sig_value = flip_significance(sig_value)
                
                node_values.append(fc_value)
                node_has_value.append(True)
                node_sizes.append(NODE_SIZE_HAS_DATA)
                node_sig.append(sig_value)
                node_pvals.append(node_data.get('adj.P.Val', None))
        else:
            node_values.append(0)
            node_has_value.append(False)
            node_sizes.append(NODE_SIZE_NO_DATA)
            node_sig.append(None)
            node_pvals.append(None)
    
    # Color mapping
    vmin, vmax = NODE_COLOR_BOUNDS
    node_colors = []
    
    for i, node in enumerate(G.nodes()):
        if node_has_value[i] and COLOR_METABOLITE_DATA:
            if node_sig[i] == 'NS':
                node_colors.append(NODE_NOTSIG_COLOR)
            else:
                node_colors.append(interpolate_color(node_values[i], vmin, vmax, vcenter=0))
        else:
            node_colors.append(NODE_NODATA_COLOR)
    
    return {
        'values': node_values,
        'has_value': node_has_value,
        'sizes': node_sizes,
        'sig': node_sig,
        'colors': node_colors,
        'has_metabolomics': has_metabolomics,
        'pvals': node_pvals
    }

# =============== CREATE FIGURE ====================

def create_figure(group1_display, group2_display, normalize_display):
    """Create the network figure based on selected parameters."""
    
    # Convert display names to internal IDs
    group1_internal = STRAIN_DISPLAY_TO_INTERNAL[group1_display]
    group2_internal = STRAIN_DISPLAY_TO_INTERNAL[group2_display]
    normalize_by = NORMALIZE_OPTIONS[normalize_display]
    
    # Make a copy of the graph to modify
    G = G_ORIGINAL.copy()
    
    # Calculate flux changes
    flux_log2fc, group1_samples = calculate_flux_changes(
        samples, group1_internal, group2_internal, normalize_by
    )
    
    # Process edges and nodes
    edge_data = process_edges(G, flux_log2fc, group1_samples)
    node_data = process_nodes(G, group1_internal, group2_internal)
    
    # Edge traces
    edge_traces = []
    for u, v, edge_attr in G.edges(data=True):
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        
        edge_idx = list(G.edges()).index((u, v))
        
        edge_trace = go.Scatter(
            x=[x0, x1, None],
            y=[y0, y1, None],
            mode='lines',
            line=dict(
                width=edge_data['widths'][edge_idx] * EDGE_WIDTH_SCALE,
                color=edge_data['colors'][edge_idx],
                dash='dot' if edge_data['styles'][edge_idx] == 'dotted' else 'solid'
            ),
            hoverinfo='skip',
            showlegend=False
        )
        edge_traces.append(edge_trace)
    
    # Edge hover points
    edge_hover_x = []
    edge_hover_y = []
    edge_hover_text = []
    
    for u, v, edge_attr in G.edges(data=True):
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        
        mid_x = (x0 + x1) / 2
        mid_y = (y0 + y1) / 2
        edge_hover_x.append(mid_x)
        edge_hover_y.append(mid_y)
        
        edge_idx = list(G.edges()).index((u, v))
        
        hover_text = ""
        if 'reaction' in edge_attr:
            hover_text += f"<b>Reaction:</b> {edge_attr['reaction']}<br>"
        
        # Add reaction_name if available
        if 'reaction_name' in edge_attr:
            hover_text += f"<b>Name:</b> {edge_attr['reaction_name']}<br>"
        
        if edge_data['has_value'][edge_idx]:
            hover_text += f"<b>Flux Log2FC:</b> {edge_data['values'][edge_idx]:.3f}<br>"
        else:
            hover_text += f"<b>Flux Log2FC:</b> No data<br>"
        
        if 'element_transfer' in edge_attr:
            hover_text += f"<b>Transfer:</b> {edge_attr['element_transfer']}<br>"
        
        if edge_data['width_has_value'][edge_idx]:
            hover_text += f"<b>Carbon Flux:</b> {edge_data['widths_raw'][edge_idx]:.2f}<br>"
        
        hover_text += f"<b>Type:</b> {edge_attr.get('edge_type', 'reaction')}"
        edge_hover_text.append(hover_text)
    
    edge_hover_trace = go.Scatter(
        x=edge_hover_x,
        y=edge_hover_y,
        mode='markers',
        marker=dict(size=10, color='rgba(0,0,0,0)'),
        hoverinfo='text',
        hovertext=edge_hover_text,
        showlegend=False
    )
    
    # Node trace
    node_x = []
    node_y = []
    node_colors_plot = []
    node_sizes_plot = []
    node_text = []
    
    for i, node in enumerate(G.nodes()):
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)
        node_colors_plot.append(node_data['colors'][i])
        
        if node_data['has_value'][i]:
            node_sizes_plot.append(node_data['sizes'][i] * NODE_SIZE_SCALE_HAS_DATA)
        else:
            node_sizes_plot.append(node_data['sizes'][i] * NODE_SIZE_SCALE_NO_DATA)
        
        node_attr = G.nodes[node]
        node_name = node_attr.get('name', node)
        
        hover_text = f"<b>{node_name}</b><br>"
        hover_text += f"<b>ID:</b> {node}<br>"
        
        if node_data['has_value'][i]:
            hover_text += f"<b>Log2FC:</b> {node_data['values'][i]:.3f}<br>"
            hover_text += f"<b>Significance:</b> {node_data['sig'][i]}<br>"
            
            # Add adj.P.Val if available
            if node_data['pvals'][i] is not None:
                hover_text += f"<b>Adj. P-value:</b> {node_data['pvals'][i]:.4e}<br>"
        else:
            if not node_data['has_metabolomics']:
                hover_text += "<i>No metabolomics data for selected strains</i><br>"
            else:
                hover_text += "<i>No metabolomics data</i><br>"
        
        if 'formula' in node_attr:
            hover_text += f"<b>Formula:</b> {node_attr['formula']}<br>"
        
        node_text.append(hover_text)
    
    node_trace = go.Scatter(
        x=node_x,
        y=node_y,
        mode='markers',
        hoverinfo='text',
        text=node_text,
        marker=dict(
            size=node_sizes_plot,
            color=node_colors_plot,
            line=dict(color=NODE_EDGE_COLOR, width=NODE_EDGEWIDTH)
        ),
        showlegend=False
    )
    
    # Create figure
    fig = go.Figure(data=edge_traces + [edge_hover_trace, node_trace])
    
    # Add annotations (labels)
    annotations = []
    for label in LABELS:
        annotation = dict(
            x=label['x'],
            y=label['y'],
            text=label['text'],
            showarrow=False,
            font=dict(
                size=label['font_size'],
                color=label.get('font_color', '#2c3e50'),
                family='Arial, sans-serif'
            ),
            textangle=label.get('angle', 0),  # Add rotation angle
            xref='x',
            yref='y',
            align='center',
            bgcolor='rgba(255, 255, 255, 0.7)',
            borderpad=4
        )
        
        if label.get('font_weight') == 'bold':
            annotation['font']['family'] = 'Arial Black, sans-serif'
        
        annotations.append(annotation)
    
    # Update layout
    fig.update_layout(
        title=dict(
            text=f"Interactive Network: {group1_display} vs {group2_display}",
            font=dict(size=20)
        ),
        annotations=annotations,
        showlegend=False,
        hovermode='closest',
        margin=dict(b=0, l=0, r=0, t=40),
        xaxis=dict(
            showgrid=False,
            zeroline=False,
            showticklabels=False,
            range=[x_min, x_max]
        ),
        yaxis=dict(
            showgrid=False,
            zeroline=False,
            showticklabels=False,
            range=[y_min, y_max]
        ),
        plot_bgcolor='white',
        width=1200,
        height=1100
    )
    
    return fig

# =============== DASH APP ====================
url_prefix = os.environ.get('SHINYPROXY_PUBLIC_PATH', '/')


app = dash.Dash(__name__,
                requests_pathname_prefix=url_prefix,
                routes_pathname_prefix=url_prefix,
                serve_locally=True
)
app.title = "Interactive Network Visualization"

app.layout = html.Div([
    html.H1(
        "Metabolic Network Analysis",
        style={'textAlign': 'center', 'color': '#2c3e50', 'padding': '20px'}
    ),
    
    # Control panel
    html.Div([
        html.Div([
            html.Label('Group 1 (numerator):', style={'fontWeight': 'bold', 'marginRight': '10px'}),
            dcc.Dropdown(
                id='group1-dropdown',
                options=[{'label': k, 'value': k} for k in STRAIN_DISPLAY_TO_INTERNAL.keys()],
                value='BDO-ALS',
                clearable=False,
                style={'width': '200px'}
            )
        ], style={'display': 'inline-block', 'marginRight': '20px'}),
        
        html.Div([
            html.Label('Group 2 (denominator):', style={'fontWeight': 'bold', 'marginRight': '10px'}),
            dcc.Dropdown(
                id='group2-dropdown',
                options=[{'label': k, 'value': k} for k in STRAIN_DISPLAY_TO_INTERNAL.keys()],
                value='Parent-COM',
                clearable=False,
                style={'width': '200px'}
            )
        ], style={'display': 'inline-block', 'marginRight': '20px'}),
        
        html.Div([
            html.Label('Normalize by:', style={'fontWeight': 'bold', 'marginRight': '10px'}),
            dcc.Dropdown(
                id='normalize-dropdown',
                options=[{'label': k, 'value': k} for k in NORMALIZE_OPTIONS.keys()],
                value='Biomass',
                clearable=False,
                style={'width': '200px'}
            )
        ], style={'display': 'inline-block'})
        
    ], style={
        'textAlign': 'center',
        'padding': '20px',
        'backgroundColor': '#f8f9fa',
        'borderRadius': '5px',
        'margin': '0 20px'
    }),
    
    dcc.Graph(
        id='network-graph',
        style={'height': '90vh'}
    )
])

@app.callback(
    Output('network-graph', 'figure'),
    [Input('group1-dropdown', 'value'),
     Input('group2-dropdown', 'value'),
     Input('normalize-dropdown', 'value')]
)
def update_figure(group1, group2, normalize):
    """Update the figure when dropdown values change."""
    return create_figure(group1, group2, normalize)

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=8050)
