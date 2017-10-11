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
    fn = os.path.join('static/gcorr/opt_corr.tsv.gz')
    data = pd.read_table(fn, index_col=0)
    t = data.copy(deep=True)
    t.index = ['_'.join(list(reversed(x.split('_')))) for x in t.index]
    for columns in [['p1_code', 'p2_code'], ['tau1', 'tau2'], ['p1', 'p2'],
                    ['p1_num_cases', 'p2_num_cases']]:
        a = list(t[columns[0]])
        b = list(t[columns[1]])
        t[columns[1]] = a
        t[columns[0]] = b
    data = pd.concat([data, t])
    data['p1'] = (data['p1'].apply(lambda x: x.replace('_', ' ')) + ' (' +
                  data['p1_code'] + ')')
    data['p2'] = (data['p2'].apply(lambda x: x.replace('_', ' ')) + ' (' + 
                  data['p2_code'] + ')')
    phenos = sorted(list(set(data.p1) | set(data.p2)))
    vc = pd.Series(list(data.p1) + list(data.p2)).value_counts()
    starting_phenos = list(vc.head(40).index)
    return(data, phenos, starting_phenos)

def make_plot_df(selected_phenos):
    tdf = DATA[DATA.p1.isin(selected_phenos) & DATA.p2.isin(selected_phenos)]
    if tdf.shape[0] > 2:
        tdf_m = tdf.pivot_table(values='omegacor21', index='p1', columns='p2').fillna(0)
        if tdf_m.shape[0] < len(selected_phenos):
            missing = list(set(selected_phenos) - set(tdf_m.index))
            tdf_m = pd.concat([tdf_m, pd.DataFrame(0, index=missing,
                                                   columns=missing)]).fillna(0)
        distmat = dist.pdist(tdf_m, metric='euclidean')
        # Allow the following methods 'single', 'complete', 'average', 'weighted'
        # and 'ward' because they are O(n^2)
        hclust = sch.linkage(distmat, method='complete')
        dend = sch.dendrogram(hclust, no_plot=True)
        order = pd.Series(tdf_m.index, index=range(tdf_m.shape[0]))[dend['leaves']]
    else:
        order = selected_phenos
    # Include all phenotypes in sindex regardless of whether they have any
    # estimates in tdf.
    sindex = pd.Series(1 + np.arange(len(selected_phenos)), index=order)
    tdf['xind'] = sindex[tdf.p1].values
    tdf['yind'] = sindex[tdf.p2].values
    minz = np.log(3)
    maxz = np.log(4.5)
    minsize = 50.
    maxsize = 130
    m = (maxsize - minsize) / (maxz - minz)
    b = maxsize - m * maxz
    tdf.loc[:, 'size'] = np.log(tdf.drawz) * m + b
    tdf.loc[tdf['size'] > maxsize, 'size'] = maxsize
    return(tdf, sindex)

def gcorr_scatter(selected_phenos):
    if len(selected_phenos) > 0:
        tdf,sindex = make_plot_df(selected_phenos)
    else:
        # If there are no selected phenos, we'll just show the starting data
        # until the user makes a selection.
        tdf,sindex = make_plot_df(STARTING_PHENOS)

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
            ticktext=sindex.index,
            tickvals=sindex,
            showticklabels=True,
            fixedrange=True,
            range=[0, sindex.shape[0] + 1.5],
        ),
        yaxis=dict(
            ticktext=sindex.index,
            tickvals=sindex,
            showticklabels=True,
            fixedrange=True,
            range=[0, sindex.shape[0] + 1.5],
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
            colorscale='RdBu',
            cmin=-1,
            cmax=1,
            showscale=True
        )
    )]
    fig = go.Figure(data=data, layout=layout)
    return(fig)

DATA, PHENOS, STARTING_PHENOS = initialize()
