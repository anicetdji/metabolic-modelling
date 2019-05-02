import cobra
from math import fabs
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import seaborn as sns


def recap_ferm_data(data_ferm, mets, scen=None, q_min=0.4, q_max=0.6):
    """
    Returns a dataframe with multiindex recapitulating fermentation data.
    All metabolites must be in positive terms, it returns consummed ones in negative terms.
    All dataframes must contain the same metabolites.
    Input: dictionary of dataframes with keys describing different scenarios analyzed.
    """
    if scen is None:
        scen = list(data_ferm.keys())
    hier_idx = pd.MultiIndex.from_product([mets, ['min', 'max', 'med', 'mean']], names=['mets', 'vals'])
    data_out = pd.DataFrame(index=hier_idx, columns=scen)
    
    for s in scen:
        for met in mets:
            if met=='glc' or met=='gly' or met=='o2' or met=='glc__D' or met=='glyc' or met=='sucr':
                data_out.loc[met, 'min'][s] = min(0,-1*data_ferm[s][met+'_sp'].quantile(q_max))
                data_out.loc[met, 'max'][s] = min(0,-1*data_ferm[s][met+'_sp'].quantile(q_min))
                data_out.loc[met, 'med'][s] = min(0,-1*data_ferm[s][met+'_sp'].quantile(0.5))
                data_out.loc[met, 'mean'][s] = min(0,-1*data_ferm[s][met+'_sp'].mean())
            else:
                data_out.loc[met, 'min'][s] = max(0,data_ferm[s][met+'_sp'].quantile(q_min))
                data_out.loc[met, 'max'][s] = max(0,data_ferm[s][met+'_sp'].quantile(q_max))
                data_out.loc[met, 'med'][s] = max(0,data_ferm[s][met+'_sp'].quantile(0.5))
                data_out.loc[met, 'mean'][s] = max(0,data_ferm[s][met+'_sp'].mean())
                
    return data_out


def bilanC(data_ferm=None, data_model=None, model=None, plot_title='Carbon Balance Fermentation', threshold=1, scen_model=[], figsize=(11,8)):
    """
    data_ferm and data_model must be dictionaries of SERIES 
    Plots and returns a figure containing Carbon Balance. If data_model provided, a cobra model must also be provided.
    """
    
    ###########################
    ## Fermentation Data ######
    ###########################
    
    data_plot_ferm = pd.DataFrame()
    data_plot_model = pd.DataFrame()
    scen = list()
    
    if not data_ferm is None:
    
        # Paramètres
        scen_ferm = list(data_ferm.columns)
        scen += scen_ferm

        # Input Carbon Sources
        input_c = list()
        idx = data_ferm.index
        if idx.contains('glc') or idx.contains('glc__D') or idx.contains('Glucose') or idx.contains('D-Glucose'):
            input_c.append('Glucose')
        if idx.contains('sucr') or idx.contains('Sucrose'):
            input_c.append('Sucrose')
        if idx.contains('gly') or idx.contains('glyc') or idx.contains('Glycerol'):
            input_c.append('Glycerol')


        # Parametres Fermentation
        names_ferm = {'glc':'Glucose', 'ac':'Acetate', 'acetone':'Acetone', 'o2':'O2', 'co2':'CO2', 'gly':'Glycerol', 'pyr':'Pyruvate', 'succ':'Succinate', 'lac':'Lactate', 'for':'Formate', 'MCA':'Prenate', 'HMG':'HMG', 'IBN':'Isobutene', 'sucr': 'Sucrose', 'iprop':'Isopropanol'}
        coeffs_c={'glc':6, 'ac':2, 'acetone':3, 'co2':1, 'bm':41, 'gly':3, 'pyr':3, 'lac':3, 'succ':4, 'for':1, 'MCA':5 , 'HMG':6 , 'IBN':4, 'sucr':12, 'o2':0, 'iprop':3}

        data_plot_ferm = pd.DataFrame(columns=scen)

        # Calculate biomass in mol
        bm_mw = 999.7 # gr/mol
        data_bm = list()
        data_bm_exch = list()
        data_plot_ferm.loc['Biomass'] = coeffs_c['bm']*data_ferm.loc['bm']

        # Calculate consumption for each metabolite
        for m in data_ferm.index:
            if m != 'bm':
                data_m = data_ferm.loc[m]*coeffs_c[m]      
                if sum([fabs(i) for i in data_m]) > 0:
                    data_plot_ferm.loc[names_ferm[m]] = data_m

        # Take all to 100
        sum_c_in = [fabs(sum(data_plot_ferm.loc[input_c, s])) for s in scen]
        data_plot_ferm = 100*data_plot_ferm/sum_c_in
    
    ###########################
    ###### Model Data #########
    ###########################        
    
    # Paramètres Modélisation
    if not data_model is None and model is None:
        print('Error for Model Data, Model not found')
        return None
    
    elif not data_model is None:
        if len(scen_model) == 0:
            scen_model = list(data_model.keys())
            
        scen += scen_model
            
        biomass_react = 'BIOMASS_Ec_iML1515_WT_75p37M'
        data_plot_model = pd.DataFrame(columns=scen_model)

        # Add Biomass reaction
        bmr = model.reactions.get_by_id(biomass_react)
        coeff_bmr = 0
        for m in bmr.metabolites:
            if list(m.elements.keys()).count('C') > 0:
                coeff_c = m.elements['C']
                coeff_react = bmr.get_coefficient(m.id)
                coeff_bmr += coeff_c*coeff_react
        data_plot_model.loc['Biomass'] = [fabs(coeff_bmr*data_model[s].loc[bmr.id]) for s in scen_model]

        # Other exchange reactions
        for r in model.exchanges:
            data_m = list()
            
            for m in r.metabolites:
                if list(m.elements.keys()).count('C') > 0:
                    coeff_m = r.get_coefficient(m.id)
                    coeff_c = m.elements['C']
                    for s in scen_model:
                        react_data_scen = data_model[s]
                        flux_c = -1*react_data_scen.loc[r.id]*coeff_c*coeff_m
                        data_m.append(flux_c)
                    break
                    
            if sum([fabs(i) for i in data_m]) > threshold:
                if m.name == 'D-Glucose':
                    data_plot_model.loc['Glucose'] = data_m
                else:
                    data_plot_model.loc[m.name] = data_m
            
        # Normalise to carbon intake (%)
        for s in scen_model:
            input_c = 0
            for m in data_plot_model.index:
                if data_plot_model.loc[m, s]<0:
                    input_c += fabs(data_plot_model.loc[m, s])
            data_plot_model[s] = 100*data_plot_model[s]/input_c

    ###########################
    ########   Plot   #########
    ###########################        
    
    #### A verifier !!!
    
    data_plot = pd.concat([data_plot_model, data_plot_ferm])
    data_table = data_plot.values.astype(np.float64).round(2).tolist()
    rows = data_plot.index
    data_plot[data_plot<0] = 0
    
    to_drop = list()
    for idx in data_plot.index:
        if sum([fabs(x) for x in data_plot.loc[idx]]) == 0:
            to_drop.append(idx)
    data_plot.drop(to_drop, axis=0, inplace=True)
    
    
    # Calculate number of metabolites to plot
    n_mets = len(data_plot.index)

    # Plot 
    fig = plt.figure(figsize=figsize)
    cols = list(cm.Spectral(np.linspace(0,1,n_mets)))
    ind = range(len(data_plot.columns))
    width = 0.7
    btm = [0 for _ in data_plot.columns]

    for m in data_plot.index:
        plt.bar(ind, data_plot.loc[m], width, color=cols.pop(), bottom=btm, label=m)
        btm += data_plot.loc[m]

    plt.ylabel('% of Carbon Intake')
    plt.title(plot_title)
    plt.xticks([])
    plt.legend(bbox_to_anchor=(1,1))

    
    plt.subplots_adjust(bottom=0.3)
    plt.table(cellText=data_table, rowLabels=rows, colLabels=scen, loc='bottom')
    return fig             
                  
        
        
        