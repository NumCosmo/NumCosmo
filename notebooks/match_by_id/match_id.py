def match_ID(cat_1, cat_2, option = 'right'):

    mt_catalog = pd.merge(cat_1, cat_2, how='inner', on=['ID_members'])
       

    match option:

                case "left":
                    
                        group_cat2 = mt_catalog.groupby(['ID_cat2']) # This groups the catalog by same "ID_cat2" objects.
                
                        cat2_id = mt_catalog['ID_cat2'].unique() # This get a "ID_cat2" list (without duplicates).
                    
                        match = []
                        
                        for id2 in cat2_id:   # For each id2 in "ID_cat2" list we create a Tuple with id2 and 
                                              #              the matched ids in "ID_cat1": (id2, [matched_ids]).
                    
                            mt_cat2 = group_cat2.get_group(id2)[['ID_cat1']]
                    
                            if len(mt_cat2) > 1:
                                match.append((id2, mt_cat2.squeeze().unique().tolist()))
                        
                            else:
                                match.append((id2, [mt_cat2.squeeze()]))
                    

                case "right":
                    
                       group_cat1 = mt_catalog.groupby(['ID_cat1']) # This groups the catalog by same "ID_cat1" objects.

                       cat1_id = mt_catalog['ID_cat1'].unique() # This get a "ID_cat1" list (without duplicates).
                
                       match = []
                        
                       for id1 in cat1_id:    # For each id1 in "ID_cat1" list we create a Tuple with id1 and 
                                              #              the matched ids in "ID_cat2": (id1, [matched_ids]).
                    
                            
                            mt_cat1 = group_cat1.get_group(id1)[['ID_cat2']]
                    
                            if len(mt_cat1) > 1:
                                match.append((id1, mt_cat1.squeeze().unique().tolist()))
                        
                            else:
                                match.append((id1, [mt_cat1.squeeze()]))

    return match
