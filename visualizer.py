import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import f_oneway
import numpy as np
import os
from sklearn.metrics import r2_score
###
data_path = 'data/graphs'

def load_data() -> tuple[pd.DataFrame,pd.DataFrame]:
    
    main_df = pd.read_csv('data/csvs/LOY and LOX Master.csv')
    
    '''
    for _, row in main_df.iterrows():
        # if row['ID'] == float('nan'): continue
        #.strip()[0:str(row['ID']).index('.')]
        # age = str(row['Age'])[0:str(row['Age']).index('.')]
        try:
            id_to_age[str(int(row['ID']))] = row['Age']
            id_to_sex[str(int(row['ID']))] = row['Sex']
        except ValueError:
            pass
    '''
    # id_to_age,id_to_sex = __meta_per_id(main_df)        

    cell_type_yscr_df = pd.read_csv('data/csvs/Donor-CellType-Yscore 2000UMI (no X genes).csv')
    cell_type_yscr_df.drop(
        cell_type_yscr_df[cell_type_yscr_df['Cell Type'] == 'unassigned'].index, inplace = True)  # get rid of unassigned cells

    return main_df, cell_type_yscr_df

def __meta_per_id(df:pd.DataFrame) -> tuple[dict,dict]:
    '''Takes a df with columns ID, Age, and Sex, and returns a tuple
    of dictionaries aligning ID with Age and ID with Sex
    '''
    
    id_to_sex = {}
    id_to_age = {}
    for _, row in df.iterrows():
        try:
            id_to_age[str(int(row['ID']))] = row['Age']
            id_to_sex[str(int(row['ID']))] = row['Sex']
        except ValueError:
            pass
    return id_to_age,id_to_sex

def scr_against_misc(df:pd.DataFrame, save:bool, med:bool, chr:str, benchmark:int,only_age=False,dir=''):    
    '''Plot individuals by various against chr-score.

    Parameters
    ----------
    df - pandas.DataFrame
    
    save - bool, whether it saves to a file or displays
    
    med - bool, True uses median y, False uses average y

    chr - str, which chromosome looking at ('X' or 'Y')

    benchmark - int, floor number of UMIs for cells to include

    only_age - bool, default False, set True to only plot against age

    dir - str, include the directory path where to save files if saving files
    
    Notes
    -----
    Columns ID, date, and sex will be ignored for the axes, as well as any column
    that starts with x or y other than the specified chromosome_med(or avg).

    Plot points on the scatter graph will be colored according to sex.

    The df must include the specified chromosome_med(or avg) column, sex, and
    at least one column that will not be excluded.

    Warnings
    --------
    If one of the names of a column that will be plotted contains a character
    that cannot be in a file name (other than parantheses), and save=True,
    only_age=False, an error will be thrown.
    '''
    
    # split to two dataframes, one per sex
    female_df = df.loc[df['Sex'] == 'female']
    male_df = df.loc[df['Sex'] == 'male']
    
    norm = 'Med' if med else 'Avg'
    exclude = ['ID', 'Date', 'Sex']
    exclude.extend([x for x in df.columns if x.lower().startswith('y') or x.lower().startswith('x')])
    
    for col in df.columns:
        if only_age: col = 'Age'
        if col in exclude:
            continue
        if col == 'Age':
            x = col
            y = f'{chr}_{norm}_{benchmark}'.lower()
        else:
            x = f'{chr}_{norm}_{benchmark}'.lower()
            y = col
        
        plt.scatter(male_df[x], male_df[y], color = 'blue', label='Male')
        plt.scatter(female_df[x], female_df[y], color='red', label='Female')

        __add_best_fit(male_df[x],male_df[y])
        __add_best_fit(female_df[x],female_df[y])

        plt.legend(loc="center left")
        plt.xlabel(x)
        plt.ylabel(y)
        plt.title(f'Indiv. {norm}. {chr}-Score x {col}, Benchmark {benchmark}')
        if save:
            try:
                plt.savefig(f'{dir}/Indiv. {norm}. {chr}-Score x {col}, Benchmark {benchmark}.png')
                if only_age: break
            except:
                col = col[0:col.index('(')]
                plt.savefig(f'{dir}/Indiv. {norm}. {chr}-Score x {col}, Benchmark {benchmark}.png')
        else:
            plt.show()
            break

        plt.close()

def __add_best_fit(x,y):
        '''Add to a graph the line of best fit labelled with R-Squared value
        '''

        # Calculate the regression line
        x_range = np.array(range(int(min(x)),int(max(x))+1))
        slope, intercept = np.polyfit(x, y, 1)
        regression_line_m = slope * x_range + intercept
        
        # Calculate R-squared value
        y_pred = slope * np.array(x) + intercept
        r_squared = r2_score(y, y_pred)

        plt.plot(x_range, regression_line_m, color='black',
                 linestyle = 'dotted' if r_squared <.05 else 'solid')

        
        # Add R-squared value as text
        r_squared_label = f"R-squared: {r_squared:.4f}" # include up to 4 dec. places
        plt.text(25, slope * 25 + intercept, r_squared_label, fontsize=8, color='black',
                 bbox=dict(facecolor=(1,1,1,.85), edgecolor='black', boxstyle='round,pad=0.5'))

def cell_type_indiv(cell_df:pd.DataFrame, meta_df:pd.DataFrame,save:bool, chr:str, dir=''):
    '''Generate and save violin plots for every donor with with a separate violin per gene.
    
    The p-value displayed is taken using a one-way ANOVA test excluding the
    unknown cells.

    Parameters
    ----------
    cell_df - pd.DataFrame of individual cells
    
    meta_df - pd.DataFrame with information about the age and sex of each donor
    
    save - bool, if True will save all graphs, otherwise will display one then stop
    
    chr - str, 'X' or 'Y'
        
    dir - if save, then enter in the directory of where the file is to be saved

    Formatting
    -----
    cell_df must include columns 'Donor ID', 'X Score' and/or 'Y Score'
    
    meta_df must include columns 'ID' (referring to Donor ID), 'Age', and 'Sex'

    Notes
    -----
    If save, then a new folder will be created in the provided dir with the graphs.
    '''
    _,id_to_sex = __meta_per_id(meta_df)
    
    # df.sort_values("Donor ID", inplace=True)
    id_list = cell_df["Donor ID"].unique()
    
    # df.drop(df[df['Cell Type'] == 'unassigned'].index, inplace = True)  # get rid of unassigned cells
    if save: os.mkdir(f'{dir}/indiv_type_{chr.lower()}scr')
    
    for id in id_list:
        if id == 'unassigned': continue  # don't need an unassigned graph
        plt.figure(figsize=(15,6))
        subset = cell_df.loc[cell_df['Donor ID'] == id]  # subset df for rows for a given donor
        
        sns.violinplot(x ="Cell Type",
             y =f"{chr} Score",
             data = subset)
        
        plt.xlabel(xlabel="Cell Type",fontsize=14)
        plt.ylabel(ylabel=f"{chr} Score",fontsize=14)

        pvalue = ctype_pvalue(subset,chr)
        legend = plt.legend(loc=1,labels=[f'p-value = {pvalue}'], handletextpad=0, handlelength=0, markerscale=0,framealpha=1)
        legend.get_frame().set_facecolor('lightgray')
        for handle in legend.legend_handles:
            handle.set_linestyle('')

        sex = 'm' if id_to_sex[(str(id))] == 'male' else 'f'
        plt.title(f'Donor {id} ({sex}) {chr}-Score by Gene')
        if not save:
            plt.show()
            break
        plt.savefig(f'{dir}/indiv_type_{chr}scr/Donor {id} {chr}.png')
        break
        plt.close()
        
def ctype_pvalue(df:pd.DataFrame, chr:str) -> float:
    '''Generate p-value for differences between cell type using one-way ANOVA test.

    Helper method for cell_type_sex and cell_type_indiv
    '''
    ctypes = df['Cell Type'].unique()
    ctypes = np.delete(ctypes,np.where(ctypes=='Unknown'))
    # ctypes = np.delete(ctypes,np.where(ctypes=='Unknown'))
    ctypes = ctypes[~pd.isnull(ctypes)]
    # filter(lambda v: v==v, ctypes)
    datasets = [df.loc[df['Cell Type'] == ctype][f'{chr} Score']
                 for ctype in ctypes]
    _, pvalue = f_oneway(*datasets)
    return pvalue

def cell_type_sex(cell_df:pd.DataFrame, meta_df:pd.DataFrame,save:bool, chr:str, sex:str, dir=''):
    '''Generate a violin plot for each cell type with aggregate data from
    all cells from a specified sex.

    The p-value displayed is taken using a one-way ANOVA test excluding the
    unknown cells.

    Parameters
    ----------
    cell_df - pd.DataFrame of individual cells
    
    meta_df - pd.DataFrame with information about the age and sex of each donor
    
    save - bool, if True will save, otherwise will display
    
    chr - str, 'X' or 'Y'
    
    sex - str, 'male' or 'female'
    
    dir - if save, then enter in the directory of where the file is to be saved

    Formatting
    -----
    cell_df must include columns 'Donor ID', 'X Score' and/or 'Y Score'
    
    meta_df must include columns 'ID' (referring to Donor ID), 'Age', and 'Sex'
    '''
    _,id_to_sex = __meta_per_id(meta_df)

    plt.figure(figsize=(15,6))
    # subset = df.loc[id_to_sex[df['ID']] == 'male']
    id_to_sex['unassigned'] = ' ' # get rid of bug in subsetting
    
    subset = cell_df[cell_df['Donor ID'].apply(lambda x: id_to_sex[str(x)] == sex)] # subset of rows of a certain sex

    # confirm that no values are below 0
    # for gene in subset['Cell Type'].unique():
    #     print(str(gene)+': '+str(subset[subset['Cell Type'] == gene]['Y Score'].min()))

    sns.violinplot(x ="Cell Type",
             y =f"{chr} Score",
             data = subset)
    plt.title(f'Violin Plot {chr} {sex} b2000 Total by Type')
    
    plt.xlabel(xlabel="Cell Type",fontsize=14)
    plt.ylabel(ylabel=f"{chr} Score",fontsize=14)

    pvalue = ctype_pvalue(cell_df,chr)
    legend = plt.legend(loc=1,labels=[f'p-value = {pvalue}'], handletextpad=0, handlelength=0, markerscale=0,framealpha=1)
    legend.get_frame().set_facecolor('lightgray')
    for handle in legend.legend_handles:
        handle.set_linestyle('')
    
    if not save:
        plt.show()
    else:
        plt.savefig(f'{dir}/all {sex} {chr} b2000 violin plot by ctype.png')
    plt.close()

def y_zero_cell(cell_df:pd.DataFrame, meta_df:pd.DataFrame, save:bool, dir=''):
    '''Generates a scatter plot of all donors and the percentage of their cells
    that return a 0 Y-UMI count plotted against their age. Male and female donors
    are separated by color.

    Parameters
    ----------
    cell_df - pd.DataFrame of individual cells
    
    meta_df - pd.DataFrame with information about the age and sex of each donor
    
    save - bool, if True will save, otherwise will display
    
    dir - if save, then enter in the directory of where the file is to be saved

    '''
    id_to_age,id_to_sex = __meta_per_id(meta_df)

    id_list = cell_df["Donor ID"].unique()
    age_f = []
    age_m = []
    percent_zero_f = []
    percent_zero_m = []
    for id in id_list:
        if id == 'unassigned': continue  # don't need an unassigned graph
        subset = cell_df.loc[cell_df['Donor ID'] == id]
        total = len(subset['Donor ID'])
        zero = len(subset.loc[cell_df['Y Score'] == 0]['Donor ID'])
        percent = zero/total
        if id_to_sex[str(id)] == 'male':
            age_m.append(id_to_age[str(int(id))])
            percent_zero_m.append(percent)
        else:
            age_f.append(id_to_age[str(int(id))])
            percent_zero_f.append(percent)
    plt.scatter(age_m, percent_zero_m, color = 'blue', label='Male')
    plt.scatter(age_f, percent_zero_f, color='red', label='Female')
    
    __add_best_fit(age_f,percent_zero_f)
    __add_best_fit(age_m,percent_zero_m)
    
    plt.xlabel('Age')
    plt.ylabel('% Cells with 0 Y-Chr UMI')
    plt.title('Percent 0 Y-UMIs by Age')
    plt.legend(loc="center left")
    if save:
        plt.savefig(f'{dir}/y-chr percent 0 per person.png')
    else:
        plt.show()
    plt.close()

def main():
    main_df, type_df = load_data()
    # cell_type_indiv(type_df, main_df,True, 'Y', dir=data_path)
    cell_type_sex(type_df, main_df,True, 'Y', 'male',dir=data_path)
    # scr_against_misc(main_df, True, False, 'Y', 2000,only_age=True,dir=data_path)
    # y_zero_cell(type_df, main_df, True,dir=data_path)
    # print(compare_scores(old_scores_df, new_y_cells)[0:50])

if __name__ == '__main__':
    main()
