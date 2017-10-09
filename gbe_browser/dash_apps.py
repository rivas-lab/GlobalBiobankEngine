import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go

import os

import numpy as np
import pandas as pd

def initialize():
    fn = os.path.join('static/gcorr/opt_corr.tsv')
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
    # starting_phenos = ['asthma', 'diabetes', 'rheumatoid_arthritis', 
    #                    'high_cholesterol']
    vc = pd.Series(list(data.p1) + list(data.p2)).value_counts()
    starting_phenos = list(vc.head(40).index)
    return(data, phenos, starting_phenos)

def make_plot_df(df):
    ind = sorted(list(set(df['p1']) | set(df['p2'])))
    tdf = df[df.apply(lambda x: x['p1'] in ind and x['p2'] in ind, axis=1)]
    sindex = pd.Series(range(len(ind)), index=ind)
    tdf['xind'] = sindex[tdf.p1].values
    tdf['yind'] = sindex[tdf.p2].values
    # tdf['size'] = -np.log10(tdf.drawp)
    tdf['size'] = np.log(tdf.drawz)
    tdf['size'] = tdf['size'] - tdf['size'].min()
    tdf['size'] = tdf['size'] / tdf['size'].max()
    tdf['size'] = tdf['size'] * 100 + 50
    return(tdf, sindex)

def gcorr_scatter(selected_phenos):
    tdf = DATA[DATA.p1.isin(selected_phenos) & DATA.p2.isin(selected_phenos)]
    tdf,sindex = make_plot_df(tdf)
    text = ('gcorr: ' + tdf.omegacor21.round(3).astype(str) + '<BR>' + 
            'gcorr SE: ' + tdf.omegacor21_se.round(3).astype(str) + '<BR>' + 
            'z: ' + tdf.drawz.round(3).astype(str) + '<BR>' + 
            'pi2: ' + tdf.pi2.round(3).astype(str) + '<BR>' + 
            'tau1: ' + tdf.tau1.round(3).astype(str) +  '<BR>' + 
            'tau2: ' + tdf.tau2.round(3).astype(str))
    layout = go.Layout(
        height=1000,
        width=1200,
        hovermode='closest',
        xaxis=dict(
            ticktext=sindex.index,
            tickvals=sindex,
            showticklabels=True,
            fixedrange=True,
        ),
        yaxis=dict(
            ticktext=sindex.index,
            tickvals=sindex,
            showticklabels=True,
            fixedrange=True,
        ),
        margin=dict(
            l=300,
            r=80,
            b=300,
            t=10,
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
            showscale=True
        )
    )]
    fig = go.Figure(data=data, layout=layout)
    return(fig)

DATA, PHENOS, STARTING_PHENOS = initialize()
