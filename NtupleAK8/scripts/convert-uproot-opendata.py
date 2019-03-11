
# coding: utf-8

import uproot
import pandas
import numpy as np
import pandas as pd
import h5py
import tables
import sys

filters = tables.Filters(complevel=7, complib='blosc')

infile = sys.argv[1]
upfile = uproot.open(infile)

tree = upfile['deepntuplizer/tree']

other_branches = ['event_no',  'npv', 'ntrueInt', 'rho', 'sample_isQCD', 'n_pfcands', 'npfcands', 'n_tracks', 'ntracks', 'n_sv', 'nsv']
fj_branches = []
pfcand_branches = []
track_branches = []
sv_branches = []

for branch_name in tree.iterkeys():
    if 'pfcand_' in branch_name:
        pfcand_branches.append(branch_name)
    elif 'track_' in branch_name or 'trackBTag_' in branch_name:
        track_branches.append(branch_name)
    elif 'sv_' in branch_name:
        sv_branches.append(branch_name)
    elif 'fj_' in branch_name:
        fj_branches.append(branch_name)
        
print(other_branches)
print(fj_branches)
print(pfcand_branches)
print(track_branches)
print(sv_branches)

def _write_carray(a, h5file, name, group_path='/', **kwargs):
    h5file.create_carray(group_path, name, obj=a, filters=filters, createparents=True, **kwargs)
    
def _transform(dataframe, max_particles=100, start=0, stop=-1):
    from collections import OrderedDict
    v = OrderedDict()
    
    df = dataframe.iloc[start:stop]
    max_jets = len(df.index.get_level_values(0).unique())
    
    for column in df.columns:
        c = np.zeros((max_jets, max_particles))
        for i, new_df in df[column].groupby(level=0):
            max_part = min(max_particles,len(new_df))
            c[i,:max_part] = new_df[:max_part].values
        v[column] = c
    return v

df_other = tree.pandas.df(branches=other_branches, entrystart=0, entrystop = None)
df_sv = tree.pandas.df(branches=sv_branches, entrystart=0, entrystop = None)
df_track = tree.pandas.df(branches=track_branches, entrystart=0, entrystop = None)
df_pfcand = tree.pandas.df(branches=pfcand_branches, entrystart=0, entrystop = None)
df_fj = tree.pandas.df(branches=fj_branches, entrystart=0, entrystop = None)


with tables.open_file(infile.replace('.root','.h5'), mode='w') as h5file:
    
    max_tracks = len(df_track.index.get_level_values(1).unique())
    max_sv = len(df_sv.index.get_level_values(1).unique())
    max_pfcand = len(df_pfcand.index.get_level_values(1).unique())
    print("max_tracks",max_tracks)
    print("max_sv",max_sv)
    print("max_pfcand",max_pfcand)
    
    v_track = _transform(df_track, max_particles = 60)
    for k in v_track.keys():
        _write_carray(v_track[k], h5file, name=k)
        
    v_sv = _transform(df_sv, max_particles = 5)
    for k in v_sv.keys():
        _write_carray(v_sv[k], h5file, name=k)
        
    v_pfcand = _transform(df_pfcand, max_particles = 100)
    for k in v_pfcand.keys():
        _write_carray(v_pfcand[k], h5file, name=k)
        
    for k in df_fj.columns:
        _write_carray(df_fj[k].values, h5file, name=k)
        
    for k in df_other.columns:
        _write_carray(df_other[k].values, h5file, name=k)


f = tables.open_file(infile.replace('.root','.h5'))
print(f)
f.close()

