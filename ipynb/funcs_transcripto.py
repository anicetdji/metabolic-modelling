import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from os import listdir
import numpy as np

def dataframe_diff_expression(time='T2', strain_base='WT'):
    """
    Returns a dataframe with differential expression with respect to strain_base.
    """
    df_out = None
    if time == 'T1':
        scen=['WT_T1', 'pkt_T1', 'goret_T1', 'Saccha_T1', 'gap_T1']
    elif time == 'T2':
        scen=['WT_T2', 'pkt_T2', 'goret_T2', 'Saccha_T2', 'gap_T2']
    
    scen.remove(strain_base+'_'+time)
    wd_data = '/Users/aalarcon/Documents/Data/IBN-Transcriptomique-7/HELIXIO/Isoforms-diff/'
    data_files = listdir(wd_data)
    
    for s in scen:
        name1 = 'Isoforms_Diff_Samples_' + s + '_vs_' + strain_base + '_' + time + '.xlsx'
        name2 = 'Isoforms_Diff_Samples_' + strain_base + '_' + time + '_vs_' + s + '.xlsx'
        if data_files.count(name1) == 1:
            tr_data = pd.read_excel(wd_data + name1, index_col=1)
            
        elif data_files.count(name2) == 1:
            tr_data = pd.read_excel(wd_data + name2, index_col=1)
        else:
            tr_data = None
            print('Error !!!!!!!!!!!!! no file found for', s2, ' vs ', s1)
        
        tr_data.columns = [c.replace(' ','').replace('[','').replace(']','').replace('FPKM','') 
                               for c in tr_data.columns]
            
        if df_out is None:
            df_out = pd.DataFrame(index=tr_data.index)
            df_out['Gene_name'] = tr_data['Gene_name']
        
        df_out[s] = 0

        for idx in tr_data.index:
            if tr_data.q_value.loc[idx] <= 0.05:
                if tr_data[s].loc[idx] > tr_data[strain_base+'_'+time].loc[idx]:
                    df_out.set_value(idx, s, 1)
                elif tr_data[s].loc[idx] < tr_data[strain_base+'_'+time].loc[idx]:
                    df_out.set_value(idx, s, -1)
    return df_out
           
    

def comp_genes_transcripto(loci, time='T2', plot=False, title_plot=None):
    """
    Compare transcriptomics of different genes given their locus.
    if arg plot=False, returns a dataframe, else returns nothing and plots the dataframe
    as a heatmap. Title of plot can be modified with title_plot argument.
    """
    if time == 'T1':
        scen=['WT_T1', 'pkt_T1', 'goret_T1', 'Saccha_T1', 'gap_T1']
    elif time == 'T2':
        scen=['WT_T2', 'pkt_T2', 'goret_T2', 'Saccha_T2', 'gap_T2']
    df_comp = pd.DataFrame(index=['time']+scen[:-2], columns=scen, data=np.nan)
    wd_data = '/Users/aalarcon/Documents/Data/IBN-Transcriptomique-7/HELIXIO/Isoforms-diff/'
    data_files = listdir(wd_data)
    for s1 in scen:
        for s2 in scen[scen.index(s1)+1:]:
            name1 = 'Isoforms_Diff_Samples_' + s2 + '_vs_' + s1 + '.xlsx'
            name2 = 'Isoforms_Diff_Samples_' + s1 + '_vs_' + s2 + '.xlsx'
            if data_files.count(name1) == 1:
                tr_data = pd.read_excel(wd_data + name1, index_col=1)
            elif data_files.count(name2) == 1:
                tr_data = pd.read_excel(wd_data + name2, index_col=1)
            else:
                tr_data = None
                print('Error !!!!!!!!!!!!! no file found for', s2, 'vs', s1)
            tr_data.columns = [c.replace(' ','').replace('[','').replace(']','').replace('FPKM','') 
                               for c in tr_data.columns]
            count = 0
            for locus in loci:
                if tr_data.q_value.loc[locus] <= 0.05:
                    if tr_data[s2].loc[locus] > tr_data[s1].loc[locus]:
                        count += 1
                    elif tr_data[s2].loc[locus] < tr_data[s1].loc[locus]:
                        count -= 1
            df_comp.set_value(s1, s2, count/len(loci))

    for s2 in scen:
            if s2.split('_')[1] == 'T1':
                s1 = s2.split('_')[0] + '_T2'
            else:
                s1 = s2.split('_')[0] + '_T1'
            name1 = 'Isoforms_Diff_Samples_' + s2 + '_vs_' + s1 + '.xlsx'
            name2 = 'Isoforms_Diff_Samples_' + s1 + '_vs_' + s2 + '.xlsx'
            if data_files.count(name1) == 1:
                tr_data = pd.read_excel(wd_data + name1, index_col=1)
            elif data_files.count(name2) == 1:
                tr_data = pd.read_excel(wd_data + name2, index_col=1)
            else:
                tr_data = None
                print('Error !!!!!!!!!!!!! no file found for', s2, 'vs', s1)
            tr_data.columns = [c.replace(' ','').replace('[','').replace(']','') for c in tr_data.columns]
            count = 0.0
            for locus in loci:
                if tr_data.q_value.loc[locus] <= 0.05:
                    if tr_data[s2].loc[locus] > tr_data[s1].loc[locus]:
                        count += 1
                    elif tr_data[s2].loc[locus] < tr_data[s1].loc[locus]:
                        count -= 1
            df_comp.set_value('time', s2, count/len(loci))
    if plot:
        fig = plt.figure(figsize=(5,5))
        sns.heatmap(df_comp, cmap='coolwarm', annot=True, linewidth=2, cbar=False, center=0)
        if title_plot is None:
            plt.title('Comparison Transcriptomics')
        else:
            plt.title(title_plot)
        return fig
    else:
        return df_comp