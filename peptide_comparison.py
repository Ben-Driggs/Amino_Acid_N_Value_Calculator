import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def main():
    mods = ['m', 'q', '1', '2', '3', '4', 's', 't', 'y', 'k', 'c', 'o']
    args = sys.argv[1:]
    
    fig = plt.figure(figsize=(10, 10))
    
    count = 1
    for d in args:
    
        d = d.split(',')
        group = d[0][:6]
        
        # set plot colors
        if group == "Diet_A":
            color = "blue"
            color2 = "cyan"
            # color3 = "navy"
            # color4 = "deepskyblue"
        elif group == "Diet_C":
            color = "red"
            color2 = "mediumorchid"
            # color3 = "darkred"
            # color4 = "darkorchid"
        elif group == "Diet_F":
            color = "green"
            color2 = "springgreen"
            # color3 = "darkgreen"
            # color4 = "mediumseagreen"
        else:
            color = "orange"
            color2 = "orangered"
            # color3 = "crimson"
            # color4 = "mediumvioletred"
        
        emp_df = pd.read_csv(d[0], sep='\t', low_memory=False)
        lit_df = pd.read_csv(d[1], sep='\t', low_memory=False)
        
        emp_df = emp_df[emp_df.loc[:, "time"] > 25.0]
        lit_df = lit_df[lit_df.loc[:, "time"] > 25.0]
        
        columns = ['Protein ID', 'Sequence', 'cf', 'mz', 'n_value', 'abundances', 'n_value_stddev']
        l_columns = ['Protein ID', 'Sequence', 'cf', 'mz', 'abundances', 'n_value']
        emp_df = emp_df[columns]
        lit_df = lit_df[l_columns]
        
        emp_df = emp_df.drop_duplicates()
        lit_df = lit_df.drop_duplicates()
        
        # remove any rows with modified peptides
        e_mask = emp_df['Sequence'].str.contains('|'.join(mods))
        l_mask = lit_df['Sequence'].str.contains('|'.join(mods))
        emp_df = emp_df[~e_mask]
        lit_df = lit_df[~l_mask]
        
        # filter by abundance
        # sum abundances and take top 50%
        emp_df.loc[:, 'sum_abundances'] = emp_df['abundances'].apply(lambda v: sum([float(value) for value in v[1:-1].split(', ')]))
        lit_df.loc[:, 'sum_abundances'] = lit_df['abundances'].apply(lambda v: sum([float(value) for value in v[1:-1].split(', ')]))
        
        emp_threshold = emp_df['sum_abundances'].quantile(0.25)
        lit_threshold = lit_df['sum_abundances'].quantile(0.25)
        
        emp_df = emp_df[emp_df['sum_abundances'] >= emp_threshold]
        lit_df = lit_df[lit_df['sum_abundances'] >= lit_threshold]
        
        # filter out noise by setting a limit on n_value standard deviation
        emp_df = emp_df[emp_df['n_value'] != "no valid time points"]
        emp_df['n_value_stddev'] = emp_df['n_value_stddev'].astype(float)
        emp_df['n_value'] = emp_df['n_value'].astype(float)
        # emp_df = emp_df[emp_df['n_value'] <= 100]
        filtered_emp_df = emp_df[emp_df.loc[:, 'n_value_stddev'] <= 2]
        
        emp_df['identifier'] = emp_df['Protein ID'] + '_' + emp_df['Sequence'] + '_' + emp_df['cf']
        filtered_emp_df.loc[:, 'identifier'] = filtered_emp_df.loc[:, 'Protein ID'] + '_' + filtered_emp_df.loc[:, 'Sequence'] + '_' + filtered_emp_df.loc[:, 'cf']
        lit_df['identifier'] = lit_df['Protein ID'] + '_' + lit_df['Sequence'] + '_' + lit_df['cf']
        
        df = pd.merge(emp_df, lit_df, on=['identifier'], how='left', suffixes=('_emp', '_lit'))
        filtered_df = pd.merge(filtered_emp_df, lit_df, on=['identifier'], how='left', suffixes=('_emp', '_lit'))
        
        df['n_value_emp'] = df['n_value_emp'].astype(float).round(2)
        filtered_df['n_value_emp'] = df['n_value_emp'].astype(float).round(2)
        
        x = df['n_value_lit']
        y = df['n_value_emp']
        f_x = filtered_df['n_value_lit']
        f_y = filtered_df['n_value_emp']
        
        plt.subplot(2, 2, count)
        plt.scatter(x, y, color=color, alpha=0.10, label='Unfiltered')
        coeffs1 = np.polyfit(x, y, 1)
        fit1 = np.poly1d(coeffs1)
        plt.plot(x, fit1(x), color=color, linewidth=2, label="Unfiltered Linear Fit")
        
        plt.scatter(f_x, f_y, color=color2, alpha=0.2, label='Filtered')
        coeffs2 = np.polyfit(f_x, f_y, 1)
        fit2 = np.poly1d(coeffs2)
        plt.plot(f_x, fit2(f_x), color=color2, linewidth=2, label="Filtered Linear Fit")
        
        plt.xlim(0, 100)
        plt.ylim(0, 100)
        plt.title(f"{group} - Literature vs Empirical Peptide N-Values")
        plt.xlabel("Peptide Literature N-Values")
        plt.ylabel("Peptide Empirical N-Values")
        count += 1
    
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
    plt.tight_layout()
    plt.show()
    sys.exit()


if __name__ == "__main__":
    main()
    
# filter out high standard deviations
# make tablular form of the AA n-values
# recalculate peptide n-values from empirical AA n-values
# filter top 50% abundant peptides
