import numpy as np
import pandas as pd
from astropy.table import Table
import networkx as nx
import random

class IDMatch:
    
    def __init__(self, catalog1, catalog2, columns1, columns2):

        self.cat1 = catalog1
        self.cat2 = catalog2
        self.columns1 = columns1
        self.columns2 = columns2
        self.prepare_catalogs = self.prepare_catalogs()
          
    
    def prepare_catalogs(self):
        ''' It prepares the table to the matching. '''
        
        cat1 = self.cat1.to_pandas().rename(columns = self.columns1)    
        cat2 = self.cat2.to_pandas().rename(columns = self.columns2)

        #Member numbers:
        nmem_cat1 = pd.DataFrame(cat1['cat1_id'].value_counts())
        nmem_cat2 = pd.DataFrame(cat2['cat2_id'].value_counts())
        
        #It adds the 'members number' column to the catalogs:
        cat1_prepared = pd.merge(nmem_cat1, cat1, how='inner', on=['cat1_id']).rename(columns={'count': 'nmem_cat1'})
        cat2_prepared = pd.merge(nmem_cat2, cat2, how='inner', on=['cat2_id']).rename(columns={'count': 'nmem_cat2'})

        return cat1_prepared, cat2_prepared
    
    
    def get_shared_members_fraction_catalog(self, matched_catalog, option):
        ''' It calculates the share fraction between the objects of the two catalogs. '''
    
        shared_count = matched_catalog.groupby(['cat1_id', 'cat2_id'])['cat2_id'].transform('count')
        
        if option == 'pmem1':
            matched_catalog['pmem_sum_cat1'] = matched_catalog.groupby(['cat1_id'])['pmem'].transform('sum') 
            frac_catalog1 = matched_catalog.groupby(['cat1_id', 'cat2_id'])['pmem'].transform('sum') / matched_catalog['pmem_sum_cat1']
            frac_catalog2 = shared_count / matched_catalog['nmem_cat2']
    
        elif option == 'pmem2':
            matched_catalog['pmem_sum_cat2'] = matched_catalog.groupby(['cat2_id'])['pmem'].transform('sum') 
            frac_catalog1 = shared_count / matched_catalog['nmem_cat1']
            frac_catalog2 = matched_catalog.groupby(['cat1_id', 'cat2_id'])['pmem'].transform('sum') / matched_catalog['pmem_sum_cat2']
    
        elif option == 'pmem':
            matched_catalog['pmem_sum_cat1'] = matched_catalog.groupby(['cat1_id'])['pmem'].transform('sum') 
            matched_catalog['pmem_sum_cat2'] = matched_catalog.groupby(['cat2_id'])['pmem'].transform('sum') 
            frac_catalog1 = matched_catalog.groupby(['cat1_id', 'cat2_id'])['pmem'].transform('sum') / matched_catalog['pmem_sum_cat1']
            frac_catalog2 = matched_catalog.groupby(['cat1_id', 'cat2_id'])['pmem'].transform('sum') / matched_catalog['pmem_sum_cat2']
    
        else:
            frac_catalog1 = shared_count / matched_catalog['nmem_cat1']
            frac_catalog2 = shared_count / matched_catalog['nmem_cat2']
    
        
        return shared_count, frac_catalog1, frac_catalog2

    
    def unique(self, matched_catalog, n_components = 5, nodes_per_side = 10, edge_prob = 0.3):

        all_combinations = matched_catalog.copy().drop_duplicates(subset=['cat1_id', 'cat2_id'])      
        tuples = list(all_combinations[['cat1_id', 'cat2_id', 'cross_frac']].itertuples(index=False, name=None))
       
        G = nx.Graph()
        for id1, id2, w in tuples:
            G.add_edge(("cat1_id",id1), ("cat2_id",id2), weight=w)
        
        random.seed(42)
        total_weight = 0.0
        global_matching = set()
        
        for i, nodes in enumerate(nx.connected_components(G)):
            
            subG = G.subgraph(nodes)
            m = nx.max_weight_matching(subG, maxcardinality=False) 
            w = sum(subG[u][v]['weight'] for u, v in m)
            
            print("\nNodes in component:")
            
            for id1 in subG:
                if id1[0] == 'cat1_id':
                    for id2 in subG[id1]:
                        if id2[0] == 'cat2_id':
                            print(id1, id2, subG[id1][id2]["weight"])
            total_weight += w
            global_matching |= m
            
            print(f"Component {i}: {len(subG)} nodes, {len(m)} pairs, weight = {w:.3f}")

            
            for u, v in m:
                if u[0] == "cat1_id":
                    id1, id2 = u, v
                else:
                    id1, id2 = v, u
                    
                print(id1, id2)

        
        unique_combinations = []
        for row_data in global_matching:
            row_dict = dict(row_data) 
            unique_combinations.append(row_dict)
        
        
        # unique matched catalog
        return pd.merge(all_combinations, pd.DataFrame(unique_combinations), how='inner', on=['cat1_id', 'cat2_id']) 

    
    def matching_by_id(self, option = None):
        ''' Returns the multiples matches by ID '''

        catalogs = self.prepare_catalogs
        
        matched_catalog = pd.merge(catalogs[0], catalogs[1], how='inner', on=['id'])   
        
        # Shared members fraction
        shared_fraction = self.get_shared_members_fraction_catalog(matched_catalog, option)
    
        # It adds shared members fractions columns
        matched_catalog['shared_num'] = shared_fraction[0]
        matched_catalog['shared_frac_cat1'] = shared_fraction[1]
        matched_catalog['shared_frac_cat2'] = shared_fraction[2]
    
        frac1_array = np.array(matched_catalog['shared_frac_cat1'])
        frac2_array = np.array(matched_catalog['shared_frac_cat2'])
         
        matched_catalog['cross_frac'] = frac1_array * ( frac1_array + frac2_array ) / 2 

        unique = self.unique(matched_catalog)       
        
        return matched_catalog, unique

    


