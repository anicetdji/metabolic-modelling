import pandas as pd
import cobra
from math import fabs


def bilan_exchanges(dict_react_data, model, threshold=0, biomass_react='BIOMASS_Ec_iML1515_WT_75p37M', 
                    mets_c=False, name_formula=True, scen=None):
    """
    Returns pandas dataframe with all exchange metabolites for a given scenario.
    dict_react_data must be an Dict with keys scenarios(labels) and values reaction data SERIES.
    arg mets_c: if True metabolites containing Carbon are taken into account.
    """
    if scen is None:
        scen = dict_react_data.keys()
        
    df_exchanges = pd.DataFrame(columns=scen, index=['gr'])
    for s in dict_react_data:
        # Add growth rate
        df_exchanges.loc['gr', s] = float(format(dict_react_data[s].loc[biomass_react], '.2f'))
        # Add exchanges
        for r_exch in model.exchanges:
            m_exch = list(r_exch.metabolites.keys())[0]
            # Take into account Carbon exchanges only if mets_c is True.
            if m_exch.elements.keys().isdisjoint(['C']) or mets_c: 
                if fabs(dict_react_data[s].loc[r_exch.id]) < 10:
                    df_exchanges.loc[m_exch.id, s] = float(format(dict_react_data[s].loc[r_exch.id], '.1f'))
                else:
                    df_exchanges.loc[m_exch.id, s] = float(format(dict_react_data[s].loc[r_exch.id], '.0f'))
                    
        # Clean df_exchanges
        to_clean = list()
        for idx in df_exchanges.index:
            if sum([fabs(df_exchanges.loc[idx, s]) for s in scen]) < threshold and idx!='gr':
                to_clean.append(idx)
        df_exchanges = df_exchanges.drop(to_clean, axis=0)
        
        # Add name and formula
        if name_formula:
            mets_nogr = list(df_exchanges.index)
            mets_nogr.remove('gr')
            df_exchanges['name']=['Growth Rate']+[model.metabolites.get_by_id(m).name for m in mets_nogr]
            df_exchanges['formula']=['C41H64O16N11P']+[model.metabolites.get_by_id(m).formula for m in mets_nogr]
                                                
    return df_exchanges
    


def bilan_metabolite(met, dict_react_data, model, threshold=1, show_genes=True, gene_names=True, show_name=False, show_reaction=True, scen=None):
    """
    Returns pandas dataframe with all reactions for a given metabolite.
    Flux corresponds to creation (positive) or consumption (negative) of each metabolite 
    (its stoichiometric coefficient is taken into account).
    ** react_data : dictionary with labels as keys and values reaction data SERIES.
    If one wants a specific order it is better to use collections.OrderedDict.
    """
    
    if scen is None:
        scen = list(dict_react_data.keys())
    
    # Create Dataframe out
    df_out = pd.DataFrame(index=[r.id for r in model.metabolites.get_by_id(met).reactions], 
                             columns=scen)
    if show_genes:
        if gene_names:
            df_out['genes']=list(map(lambda x:', '.join([g.name for g in model.reactions.get_by_id(x).genes]), df_out.index))
        else:
            df_out['genes']=list(map(lambda x:', '.join([g.id for g in model.reactions.get_by_id(x).genes]), df_out.index))
            
    # Remplir DataFrame out
        
    if show_name:
        df_out['name'] = list(map(lambda x:model.reactions.get_by_id(x).name, df_out.index))
        
    if show_reaction:
        for r in df_out.index:
            if r.count('BIOMASS')==0:
                df_out.loc[r, 'reaction'] = model.reactions.get_by_id(r).reaction 
                #list(map(lambda x:model.reactions.get_by_id(x).reaction, df_out.index))
            else:
                df_out.loc[r, 'reaction'] = ''
    
    for react in df_out.index:
        stoich_coeff = model.reactions.get_by_id(react).get_coefficient(met)
        for s in scen:
            try:
                df_out.set_value(react, s, float(format(dict_react_data[s].loc[react]*stoich_coeff, '.1f')))
            except:
                df_out.set_value(react, s, 0)

    # Clean Dataframe out
    to_clean = list()
    for idx in df_out.index:
        if sum([fabs(df_out.loc[idx, s]) for s in scen]) < threshold:
            to_clean.append(idx)
    df_out = df_out.drop(to_clean, axis=0)
    
    return df_out


def reactions2genes(reacts, model):
    """
    Returns a dictionary with keys loci and names for the same list of genes (same order).
    """
    genes = list()
    for r in reacts:
        try:
            react = model.reactions.get_by_id(r)
            for g in react.genes:
                genes.append(g.id)
        except:
            print(' !!! Error: reaction' + r + 'not found in model')
            
    genes = list(set(genes))
    genes_names = list()
    for g in genes:
        g = model.genes.get_by_id(g)
        genes_names.append(g.name)
        
    return {'loci':genes, 'names':genes_names} 