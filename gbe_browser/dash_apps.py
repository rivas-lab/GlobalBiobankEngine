import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go

import os

import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
import scipy.spatial.distance as dist
import scipy.cluster.hierarchy as sch

def initialize():
    # Minimum z-score to display data for.
    min_z = 3
    
    fn = os.path.join('static/gcorr/opt_corr.tsv.gz')
    data = pd.read_table(fn, index_col=0)
    t = data.copy(deep=True)
    t.index = ['_'.join(list(reversed(x.split('_')))) for x in t.index]
    for columns in [['p1_code', 'p2_code'], ['tau1', 'tau2'], ['p1', 'p2'],
                    ['p1_num_cases', 'p2_num_cases']]:
        a = list(t[columns[0]])
        b = list(t[columns[1]])
        t.loc[:, columns[1]] = a
        t.loc[:, columns[0]] = b
    data = pd.concat([data, t])
    
    fn = os.path.join('static/gcorr/traits.tsv.gz')
    phenos = pd.read_table(fn, index_col=0)
    phenos = phenos.loc[sorted(list(set(data.p1_code) | set(data.p2_code)))]
    phenos['phenotype'] = (phenos['phenotype'].apply(lambda x: x.replace('_', ' ')).values
                           + ' (' + pd.Series(phenos.index) + ')').values
    data.loc[:, 'p1'] = phenos.loc[data['p1_code'], 'phenotype'].values
    data.loc[:, 'p2'] = phenos.loc[data['p2_code'], 'phenotype'].values
    pheno_categories = list(set(phenos['category']))
    name_to_code = pd.Series(
        dict(zip(list(data['p1']) + list(data['p2']),
                 list(data['p1_code']) + list(data['p2_code']))),
    )
    vc = pd.Series(list(data.p1) + list(data.p2)).value_counts()

    starting_phenos = list(vc.head(40).index)
    return(data, phenos, starting_phenos, pheno_categories, min_z, name_to_code)

DATA, PHENOS, STARTING_PHENOS, PHENO_CATEGORIES, MIN_Z, NAME_TO_CODE = initialize()

gcorr_layout = html.Div(children=[
    html.Div([
        dcc.Graph(
            id='gcorr-scatter',
            figure={'layout': {'title': 'Dash Data Visualization'}}
        ),
    ], style={'height':'900px', 'width':'1300px'}
    ),
    html.Div([
        html.Div([
            html.Div([
                html.P('SELECT phenotype code categories to include.'),
            ], style={'margin-left': '10px'}),
            dcc.Checklist(
                id='pheno-categories',
                options=[{'label': x, 'value': x} for x in PHENO_CATEGORIES],
                values=PHENO_CATEGORIES,
            )
        ], className='two columns'
        ),
        html.Div([
            html.Div([
                html.P('SELECT clustering method.'),
            ], style={'margin-left': '10px'}),
            dcc.RadioItems(
                id='cluster-method',
                options=[
                    {'label': 'single', 'value': 'single'},
                    {'label': 'complete', 'value': 'complete'},
                    {'label': 'average', 'value': 'average'},
                    {'label': 'ward', 'value': 'ward'},
                ],
                value='complete',
                labelStyle={'display': 'inline-block'}
            ),
        ], className='two columns'
        ),
        html.Div([
            html.Div([
                html.P('Minimum z-score.')
            ], style={'margin-left': '10px'}),
            dcc.Input(
                id='z-cutoff',
                placeholder='z >=', 
                type='number',
                value=MIN_Z,
            )
            ], className='two columns'
        ),
        html.Div([
            html.Div([
                html.P('Show phenotypes with no estimates.')
            ], style={'margin-left': '10px'}),
            dcc.Checklist(
                id='show-zero-estimates',
                options=[{'label': 'show', 'value': 'show'}],
                values=['show'],
            )
            ], className='two columns'
        ),
        html.Div([
            html.Div([
                html.P('Size scatter points based on:')
            ], style={'margin-left': '10px'}),
            dcc.RadioItems(
                id='size-var',
                options=[
                    {'label': 'z-score', 'value': 'z-score'},
                    {'label': 'membership', 'value': 'membership'},
                ],
                value='z-score',
            )
            ], className='two columns'
        ),
    ], className='row',
    ),
    html.Div([
        html.Div([
            html.Div([
                html.P('Correlation window.')
            ], style={'margin-left': '10px'}),
            dcc.RangeSlider(
                id='gcorr-range',
                count=1,
                min=-1,
                max=1,
                step=0.01,
                value=[-1, 1],
                marks=dict(zip(np.arange(-1, 1.25, 0.25), 
                               np.arange(-1, 1.25, 0.25))),
            ),
        ], className='four columns',
        ),
        html.Div([
            dcc.RadioItems(
                id='gcorr-radio',
                options=[{'label': 'include', 'value': 'include'},
                         {'label': 'exclude', 'value': 'exclude'}],
                value='include',
            )
        ], className='two columns', 
        ),
    ], className='row', style={'margin-top': '30px'},
    ),
    html.Div([
        html.Div([
            html.Div(id='gcorr-range-values', style={'margin-top': '20px'}),
        ], className='four columns',
        ),
    ], className='row', 
    ),
    html.Div([
        html.Div([
            html.Div([
                html.P('Membership window.')
            ], style={'margin-left': '10px'}),
            dcc.RangeSlider(
                id='pi2-range',
                count=1,
                min=0,
                max=0.2,
                step=0.005,
                value=[0, 0.2],
                marks=dict(zip(np.arange(0, 0.25, 0.05), 
                               [str(x) for x in np.arange(0, 0.25, 0.05)])),
            ),
        ], className='four columns',
        ),
        html.Div([
            dcc.RadioItems(
                id='pi2-radio',
                options=[{'label': 'include', 'value': 'include'},
                         {'label': 'exclude', 'value': 'exclude'}],
                value='include',
            )
        ], className='two columns', 
        ),
    ], className='row', style={'margin-top': '30px'},
    ),
    html.Div([
        html.Div([
            html.Div(id='pi2-range-values', style={'margin-top': '20px'}),
        ], className='four columns',
        ),
    ], className='row', 
    ),
    html.Div([
        html.Div([
            html.Div([
                html.P('SELECT phenotypes in the dropdown menu to add them to '
                       'the plot.')
            ], style={'margin-left': '10px'}),
            dcc.Dropdown(id='pheno-dropdown',
                         multi=True,
                         value=STARTING_PHENOS,
                        ),
        ], className='twelve columns'),
    ], className='row', style={'margin-top': '30px'},
    ),
    ])

def make_plot_df(
    selected_phenos, 
    cluster_method, 
    pheno_categories, 
    z_cutoff,
    gcorr_range, 
    gcorr_radio,
    pi2_range, 
    pi2_radio,
    show_zero_estimates,
    size_var,
):
    show_zero_estimates = show_zero_estimates == ['show']
    # Check whether z-score cutoff is a valid number and greater than the
    # minimum.
    try:
        z_cutoff = float(z_cutoff)
    except ValueError:
        z_cutoff = MIN_Z
    # Restrict selected phenotypes to those in the selected code categories if
    # needed. 
    if len(pheno_categories) < len(PHENO_CATEGORIES):
        c = PHENOS.loc[NAME_TO_CODE.loc[selected_phenos], 
                       'category'].isin(pheno_categories)
        selected_phenos = [selected_phenos[i] for i in
                           range(len(selected_phenos)) if c[i]]
    # Restrict to data from the selected phenotypes.
    tdf = DATA[DATA.p1.isin(selected_phenos) & DATA.p2.isin(selected_phenos)]
    # Filter using z-score cutoff.
    tdf = tdf[tdf['drawz'] >= z_cutoff]
    # Filter according to gcorr slider and radio.
    if gcorr_radio == 'include':
        tdf = tdf[(tdf['omegacor21'] <= gcorr_range[1]) & 
                  (tdf['omegacor21'] >= gcorr_range[0])]
    elif gcorr_radio == 'exclude':
        tdf = tdf[(tdf['omegacor21'] >= gcorr_range[1]) | 
                  (tdf['omegacor21'] <= gcorr_range[0])]
    # Filter according to pi2 slider and radio.
    if pi2_radio == 'include':
        tdf = tdf[(tdf['pi2'] <= pi2_range[1]) & 
                  (tdf['pi2'] >= pi2_range[0])]
    elif pi2_radio == 'exclude':
        tdf = tdf[(tdf['pi2'] >= pi2_range[1]) | 
                  (tdf['pi2'] <= pi2_range[0])]
    # tdf.shape[0] will always be even because each pair is represented twice in
    # the DATA table. If tdf.shape[0] == 2, then we only have one phenotype
    # which means that we won't have anything to display. Similarly,
    # tdf.shape[0] == 0 means we don't have anything to display.
    if tdf.shape[0] > 2:
        tdf_m = tdf.pivot_table(values='omegacor21', index='p1', columns='p2').fillna(0)
        if tdf_m.shape[0] < len(selected_phenos) and show_zero_estimates:
            missing = list(set(selected_phenos) - set(tdf_m.index))
            tdf_m = pd.concat([tdf_m, pd.DataFrame(0, index=missing,
                                                       columns=missing)]).fillna(0)
        distmat = dist.pdist(tdf_m, metric='euclidean')
        hclust = sch.linkage(distmat, method=cluster_method)
        dend = sch.dendrogram(hclust, no_plot=True)
        selected_phenos = pd.Series(tdf_m.index,
                                    index=range(tdf_m.shape[0]))[dend['leaves']]
    else:
        if not show_zero_estimates:
            selected_phenos = []
    
    sindex = pd.Series(1 + np.arange(len(selected_phenos)),
                       index=selected_phenos)
    tdf['xind'] = sindex[tdf.p1].values
    tdf['yind'] = sindex[tdf.p2].values
    # Set sizes for scatter points
    if size_var == 'membership':
        # Size based on membership
        # mins = np.log(DATA['pi2'].min())
        # maxs = np.log(DATA['pi2'].max())
        mins = DATA['pi2'].min()
        maxs = DATA['pi2'].max()
        minsize = 50.
        maxsize = 130
        m = (maxsize - minsize) / (maxs - mins)
        b = maxsize - m * maxs
        # tdf.loc[:, 'size'] = np.log(tdf['pi2']) * m + b
        tdf.loc[:, 'size'] = tdf['pi2'] * m + b
        tdf.loc[tdf['size'] > maxsize, 'size'] = maxsize
    elif size_var == 'z-score':
        # Default is z-score
        mins = np.log(3)
        maxs = np.log(4.5)
        minsize = 50.
        maxsize = 130
        m = (maxsize - minsize) / (maxs - mins)
        b = maxsize - m * maxs
        tdf.loc[:, 'size'] = np.log(tdf.drawz) * m + b
        tdf.loc[tdf['size'] > maxsize, 'size'] = maxsize
    return(tdf, sindex)

def gcorr_scatter(
    selected_phenos, 
    cluster_method, 
    pheno_categories, 
    z_cutoff,
    gcorr_range, 
    gcorr_radio, 
    pi2_range, 
    pi2_radio, 
    show_zero_estimates,
    size_var,
):
    if len(selected_phenos) > 0:
        tdf,sindex = make_plot_df(
            selected_phenos, 
            cluster_method,
            pheno_categories, 
            z_cutoff, 
            gcorr_range,
            gcorr_radio,
            pi2_range,
            pi2_radio,
            show_zero_estimates,
            size_var,
        )
    else:
        # If there are no selected phenos, we'll just show the starting data
        # until the user makes a selection.
        tdf,sindex = make_plot_df(
            STARTING_PHENOS, 
            cluster_method,
            pheno_categories, 
            z_cutoff, 
            gcorr_range,
            gcorr_radio,
            pi2_range,
            pi2_radio,
            show_zero_estimates,
            size_var,
        )

    text = ('gcorr: ' + tdf.omegacor21.round(3).astype(str) + '<BR>' + 
            'gcorr SE: ' + tdf.omegacor21_se.round(3).astype(str) + '<BR>' + 
            'z: ' + tdf.drawz.round(3).astype(str) + '<BR>' + 
            'pi2: ' + tdf.pi2.round(3).astype(str) + '<BR>' + 
            'tau1: ' + tdf.tau1.round(3).astype(str) +  '<BR>' + 
            'tau2: ' + tdf.tau2.round(3).astype(str))
    layout = go.Layout(
        height=800,
        width=1200,
        hovermode='closest',
        xaxis=dict(
            autorange='reversed',
            ticktext=sindex.index,
            tickvals=sindex,
            showticklabels=True,
            fixedrange=True,
            range=[0, sindex.shape[0] + 1.5],
            zeroline=False,
        ),
        yaxis=dict(
            ticktext=sindex.index,
            tickvals=sindex,
            showticklabels=True,
            fixedrange=True,
            range=[0, sindex.shape[0] + 1.5],
            zeroline=False,
        ),
        margin=dict(
            l=250,
            r=250,
            b=150,
        )
    )
    data = [go.Scatter(
        y = tdf['yind'],
        x = tdf['xind'],
        mode='markers',
        text=list(text),
        marker=dict(
            size=tdf['size'] / 10.,
            color = tdf['omegacor21'],
            colorbar=dict(title='Genetic correlation'),
            colorscale=[
                [0.0, 'rgb(165,0,38)'], 
                [0.1111111111111111, 'rgb(215,48,39)'],
                [0.2222222222222222, 'rgb(244,109,67)'],
                [0.3333333333333333, 'rgb(253,174,97)'],
                [0.4444444444444444, 'rgb(254,224,144)'],
                [0.5555555555555556, 'rgb(224,243,248)'],
                [0.6666666666666666, 'rgb(171,217,233)'],
                [0.7777777777777778, 'rgb(116,173,209)'],
                [0.8888888888888888, 'rgb(69,117,180)'], 
                [1.0, 'rgb(49,54,149)']
            ],
            cmin=-1,
            cmax=1,
            showscale=True,
        )
    )]
    fig = go.Figure(data=data, layout=layout)
    return(fig)
