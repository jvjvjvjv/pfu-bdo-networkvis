from app import create_figure

# Render the default view: BDO-ALS vs Parent-COM, no normalization
fig = create_figure('BDO-ALS', 'Parent-COM', 'None')

fig.write_html(
    'index.html',
    include_plotlyjs='cdn',
    full_html=True,
    config={
        'toImageButtonOptions': {
            'format': 'png',
            'filename': 'network_visualization',
            'width': 1200,
            'height': 1100,
            'scale': 3
        },
        'displaylogo': False,
    }
)
