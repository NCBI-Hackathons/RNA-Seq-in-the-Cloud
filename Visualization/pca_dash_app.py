import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import seaborn as sns
import plotly.graph_objs as go
from sklearn.decomposition import PCA
import numpy as np
# from sys import argv
import argparse

def get_args():
    '''
    Returns args from command line that are two data files:
        - count data
        - metadata file
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--count_data', required=True)
    parser.add_argument('-a', '--attributes', required=True)
    args = parser.parse_args()
    return args

args = get_args()

def data_process(data_file,metadata_file):
    data_df = pd.read_csv(data_file, sep='\t') #'~/Downloads/testRNAseqData.csv'
    data_df.index = data_df['geneID']
    metadata_df = pd.read_csv(metadata_file, sep='\t') #'~/Downloads/DESeq2_POE_data.csv'
    metadata_df.index = metadata_df['Run']
    col_not_counts = 8
    norm_df = data_df.iloc[:,col_not_counts:]
    return norm_df, metadata_df

data_df, metadata_df = data_process(args.count_data, args.attributes)

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__)# external_stylesheets=external_stylesheets)

def run_pca(data_df,m_df):
    '''
    Runs PCA on count data.

    Parameters:
        data_df: count data dataframe
        m_df: dataframe of metadata associated with an SRA Run
    Returns:
        m_pca_df: m_df with PCA components 1 - 3 added on
        var_explained: list of strings denoting variation explained used for
            labeling axes of PCA plot
    '''
    pca=PCA(n_components=3)
    pca.fit(data_df)
    m_pca_df = m_df.copy()
    var_explained = ['PC'+str((i+1))+'('+str(round(e*100, 1))+'% var. explained)' for i, e in enumerate(pca.explained_variance_ratio_)]
    m_pca_df["PC1"] = pca.components_[0]
    m_pca_df["PC2"] = pca.components_[1]
    m_pca_df["PC3"] = pca.components_[2]
    return m_pca_df, var_explained

def gen_metadata_colors(series):
    '''
    Returns appropriate assigning of colors based on user's dropdown selection.

    Parameters:
        series: series whose unique entries will determine the number
            of colors plotted
    Returns:
        colors: series with color hex code associated with a SRA RUN
    '''

    unique = series.unique()
    num_colors = len(unique)
    pal = sns.color_palette("Set2",num_colors)
    hex_codes = pal.as_hex()

    _d = dict(zip(unique, hex_codes))
    colors = series.map(_d)
    return colors

def separate_traces(metadata_df,metadata_OI):
    '''
    Separates traces needed by plotly into their various colors based on
    user's dropdown selection. This is a bit of a hack so that hoverboxes are
    colored correctly.

    Parameters:
        metadata_df: metadata dataframe
        metadata_OI: metadata "of interest" that the user has selected from
            the dropdown menu
    Returns:
        unique_df: dataframe of the unique type of metadata from the
            catagory the user has selected
        unique_c: list of hex codes used to color each SRA run
    '''
    colors = gen_metadata_colors(metadata_df[metadata_OI])
    unique_df =metadata_df[metadata_OI].unique()
    unique_c = colors.unique()
    if len(unique_c) < len(unique_df):
        diff = int(np.floor(len(unique_df) / len(unique_c)))
        unique_c = unique_c.tolist()*diff + unique_c[:len(unique_df)-len(unique_c)*diff].tolist()
    return unique_df, unique_c


def generate_traces(data_df,metadata_df,metadata_OI):
    '''
    Generates the traces needed for plotly to render the PCA plot.

    Parameters:
        data_df: count data dataframe
        metadata_df: metadata dataframe
        metadata_OI: metadata "of interest" that the user has selected from
            the dropdown menu
    Returns:
        data: list of plotly traces separated by their metadata "of interest"
        var_explained: list of strings denoting variation explained used for
            labeling axes of PCA plot
    '''
    metadata_df_pca, var_explained = run_pca(data_df,metadata_df)
    unique_df, unique_c = separate_traces(metadata_df_pca,metadata_OI)
    data = []
    for idx,_type in enumerate(unique_df):
        u_df = metadata_df_pca[metadata_df_pca[metadata_OI]==_type]
        trace = go.Scatter3d(x=u_df["PC1"],
                    y=u_df["PC2"],
                    z=u_df["PC3"],
                    mode='markers',
                    hoverinfo='text',
                    text=['<br>'.join(['{key}: {value}'.format(**locals()) for key, value in rowData.items()]) for index, rowData in metadata_df.loc[data_df.columns].rename_axis('Signature ID').reset_index().iterrows()],
                    marker=dict(color=unique_c[idx], size=15, opacity=0.8),
                    name=str(_type))
        #             marker={'size': 15, 'showscale': False, 'colorscale': 'Jet', 'color': metadata_dataframe.loc[signature_dataframe.columns]['pert_time']})
        data.append(trace)
    return data, var_explained

# grab traces
traces, var = generate_traces(data_df,metadata_df,'condition')

# dash html for rendering plot
app.layout = html.Div([
    # html.Label('Dropdown'),
    dcc.Dropdown(
        id ='attribute',
        options=[{'label': i, 'value': i} for i in metadata_df.columns
        ],
        value='condition'
    ),

    dcc.Graph(
        id='scatter3d' ,
        figure={
            'data': traces,# generate_traces(data_df,metadata_df,'condition')[0],

            'layout': go.Layout(
                        title='PCA Analysis | Scatter Plot<br><span style="font-style: italic;">Colored by '+'condition'+' </span>',
                        hovermode='closest',
                        width=900,
                        height = 900,
                        scene=dict(xaxis=dict(title=var[0]),
                            yaxis=dict(title=var[1]),
                            zaxis=dict(title=var[2])),
                        showlegend=True
            )
        }
    )
])



# handles dropdown selection
@app.callback(
    dash.dependencies.Output('scatter3d', 'figure'),
    [dash.dependencies.Input('attribute', 'value')])
def update_graph(v):
    n_traces,n_var = generate_traces(data_df,metadata_df,v)
    return {
                'data': n_traces,
                'layout': go.Layout(
                            title='PCA Analysis | Scatter Plot<br><span style="font-style: italic;">Colored by '+v+' </span>',
                            hovermode='closest',
                            width=900,
                            height = 900,
                            scene=dict(xaxis=dict(title=n_var[0]),
                                yaxis=dict(title=n_var[1]),
                                zaxis=dict(title=n_var[2])),
                            showlegend=True
                )
            }
        # )


if __name__ == '__main__':
    app.run_server(debug=False)
