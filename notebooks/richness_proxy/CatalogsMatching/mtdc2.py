import numpy as np
import pandas as pd

def mtrmdc2(member_data, cluster_data, truth_data):
    memberdf = member_data.to_pandas().rename(columns={'id_member': 'id'})
    truthdf = truth_data.to_pandas().rename(columns={'galaxy_id': 'id'})
    cluster_df=cluster_data.to_pandas()[['richness', 'cluster_id', 'redshift']]

    Mt_df = pd.merge(memberdf['id'], truthdf['id'], how='inner', on=['id'])
    mt_member = pd.merge(memberdf, Mt_df, how='inner', on=['id'])

    idc = mt_member[['id', 'cluster_id_member']].rename(columns={'cluster_id_member': 'cluster_id'})
    mt_cluster = pd.merge(idc, cluster_df, how='inner', on=['cluster_id'])

    mt_catalog = pd.merge(mt_cluster, truthdf, how='inner', on=['id'])

    clusters_id = mt_catalog[mt_catalog['is_central'] == True]['cluster_id'].unique()
    mt_catalog['freq'] = mt_catalog.groupby(['cluster_id', 'halo_id'])['halo_id'].transform('count')
    iscentral = mt_catalog[mt_catalog['is_central'] == True]

    match_dataframe = pd.DataFrame()
    for cl in clusters_id:
        gcut = iscentral.groupby('cluster_id').get_group(cl)
        match_dataframe = pd.concat([match_dataframe, gcut[gcut['freq'] == gcut['freq'].max()]], ignore_index=True)
    
    return match_dataframe