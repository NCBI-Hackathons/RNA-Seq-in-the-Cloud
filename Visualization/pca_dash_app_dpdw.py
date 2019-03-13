import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import seaborn as sns
import plotly.graph_objs as go
from sklearn.decomposition import PCA
# from sys import argv
import argparse

def data_process(data_file,metadata_file,metadata_OI):
    data_df = pd.read_csv(data_file, sep='\t') #'~/Downloads/testRNAseqData.csv'
    # print(data_df.head())
    data_df.index = data_df['geneID']

    metadata_df = pd.read_csv(metadata_file, sep='\t') #'~/Downloads/DESeq2_POE_data.csv'
    # print(metadata_df.head())
    # metadata_df.index = metadata_df['Run']

    col_not_counts = 8
    norm_df = data_df.iloc[:,col_not_counts:]
    return norm_df, metadata_df,metadata_OI

def run_pca(data_df,metadata_df):
    pca=PCA(n_components=3)
    pca.fit(data_df)
    var_explained = ['PC'+str((i+1))+'('+str(round(e*100, 1))+'% var. explained)' for i, e in enumerate(pca.explained_variance_ratio_)]
    metadata_df["PC1"] = pca.components_[0]
    metadata_df["PC2"] = pca.components_[1]
    metadata_df["PC3"] = pca.components_[2]
    return metadata_df, var_explained

def gen_metadata_colors(series):
    '''
    Parameters:
        series: Pandas series that will be colored
    '''

    unique = series.unique()
    num_colors = len(unique)
    pal = sns.color_palette("Set2",num_colors)
    hex_codes = pal.as_hex()

    _d = dict(zip(unique, hex_codes))
    colors = series.map(_d)

    return colors

def separate_traces(metadata_df,metadata_OI):
    colors = gen_metadata_colors(metadata_df[metadata_OI])
    unique_df =metadata_df[metadata_OI].unique()
    unique_c = colors.unique()
    return colors, unique_df, unique_c

def generate_traces(data_df,metadata_df,metadata_OI):
    metadata_df, var_explained = run_pca(data_df,metadata_df)
    colors, unique_df, unique_c = separate_traces(metadata_df,metadata_OI)
    data = []
    for idx,_type in enumerate(unique_df):
        u_df = metadata_df[metadata_df[metadata_OI]==_type]
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

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-cd', '--count_data', required=True)
    parser.add_argument('-a', '--attributes', required=True)
    parser.add_argument('-m', '--metadata', required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = get_args()
    data_df, metadata_df, metadata_OI = data_process(args.count_data, args.attributes, args.metadata)

    external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

    app = dash.Dash(__name__)# external_stylesheets=external_stylesheets)

    app.layout = html.Div([
        # html.Label('Dropdown'),
        dcc.Dropdown(
            options=[
                {'label': 'Run', 'value': 'NYC'},
                {'label': 'Project', 'value': 'MTL'},
                {'label': 'Condition', 'value': 'SF'}
            ],
            value='MTL'
        ),

        dcc.Graph(
            id='scatter3d',
            figure={
                'data': generate_traces(data_df,metadata_df,metadata_OI)[0],

                'layout': go.Layout(
                            title='PCA Analysis | Scatter Plot<br><span style="font-style: italic;">Colored by '+metadata_OI+' </span>',
                            hovermode='closest',
                            width=900,
                            height = 900,
                            scene=dict(xaxis=dict(title=generate_traces(data_df,metadata_df,metadata_OI)[1][0]),
                                yaxis=dict(title=generate_traces(data_df,metadata_df,metadata_OI)[1][1]),
                                zaxis=dict(title=generate_traces(data_df,metadata_df,metadata_OI)[1][2])),
                            showlegend=True
                )
            }
        )
    ])
    app.run_server()#debug=True)
