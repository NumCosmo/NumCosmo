import numpy as np
import pandas as pd

#-------------------------------------------------------------------------------------------------#
#MatchID: Class
#-------------------------------------------------------------------------------------------------#

class MatchID:
    def __init__(self, cluster_member_data, cluster_data, truth_data):
        self.member_data = cluster_member_data
        self.cluster_data = cluster_data
        self.truth_data = truth_data
        self.table_to_dataframe = self.table_to_dataframe()
        self.matching_by_id = self.matching_by_id()

#-------------------------------------------------------------------------------------------------#
# table_to_dataframe()
#
# Returns: 
# memberdf: Pandas dataframe of members data.
# truthdf: Pandas dataframe of halos data; 
# cluster_df: Pandas dataframe of clusters data.
#-------------------------------------------------------------------------------------------------# 

    def table_to_dataframe(self):
       
        memberdf = self.member_data.to_pandas().rename(columns={'id_member': 'id'})
        truthdf = self.truth_data.to_pandas().rename(columns={'galaxy_id': 'id'})
        cluster_df= self.cluster_data.to_pandas()[['richness', 'richness_err', 'cluster_id', 'cluster_redshift', 'cluster_redshift_err', 'cluster_ra', 'cluster_dec']]

        return memberdf, truthdf, cluster_df

#-------------------------------------------------------------------------------------------------#
# matching_by_id()
#
# Returns: 
# mt_catalog : Pandas dataframe with catalogs match.
#-------------------------------------------------------------------------------------------------# 
   
    def matching_by_id(self):

        memberdf, truthdf, cluster_df = self.table_to_dataframe
        
        #This selects commom members ID from cDC2 and cDC2 + RM data
        Mt_df = pd.merge(memberdf['id'], truthdf['id'], how='inner', on=['id']) 
        
        #This selects the matching members using Mt_df IDs
        mt_member = pd.merge(memberdf, Mt_df, how='inner', on=['id']) 
        
        #This selects the members and clusters IDs
        idc = mt_member[['id', 'cluster_id_member']].rename(columns={'cluster_id_member': 'cluster_id'})
        
        #This selects the clusters data by idc IDS 
        mt_cluster = pd.merge(idc, cluster_df, how='inner', on=['cluster_id']) 
        
        #This creates the matching catalog with clustes and halos
        mt_catalog = pd.merge(mt_cluster, truthdf, how='inner', on=['id'])

        return mt_catalog


#-------------------------------------------------------------------------------------------------#
# unique_id_match()
#
# Returns: 
# halos_bin_mz : Pandas dataframe with unique association between halos and clusters.
#-------------------------------------------------------------------------------------------------# 

    def unique_id_match(self):

        mt_catalog = self.matching_by_id
        
        #This selects the clusters IDs
        clusters_id = mt_catalog[mt_catalog['is_central'] == True]['cluster_id'].unique() 
       
        # This Adds a frequency column for each halo ID in mt_catalog
        mt_catalog['freq'] = mt_catalog.groupby(['cluster_id', 'halo_id'])['halo_id'].transform('count') 
        
        #This selects the central members
        iscentral = mt_catalog[mt_catalog['is_central'] == True] 
   
        # This loop creates the unique matching catalog using the halo ID frequency. 
        # It selects the halos by ra, dec, maximum frequency.
        
        idgroups = iscentral.groupby('cluster_id') # It groups data by cluster_id. 

        match_dataframe = pd.DataFrame()                                                                                                                                                                            
        for cl in clusters_id:                                                                                        
            
            gcut = idgroups.get_group(cl)    # For each halo in cluster_id group, it selects halos in a range of ra e dec and then 
                                             # select the halo with minimum proximity and maximum frequency of members.
            
            gcut['distance_2d'] = np.sqrt( (gcut['ra'] - gcut['cluster_ra']) ** 2 + (gcut['dec'] - gcut['cluster_dec']) ** 2  )

            print(len(gcut['distance_2d']))
                       
            dist_cut = gcut['distance_2d'] == gcut['distance_2d'].min()
            freq_cut = gcut['freq'] == gcut['freq'].max()
                        
            gcut = gcut[dist_cut]

            if len(gcut) > 1:
                gcut = gcut[freq_cut]
            else: 
                pass
            
            match_dataframe = pd.concat([match_dataframe, gcut], ignore_index=True)


        return match_dataframe