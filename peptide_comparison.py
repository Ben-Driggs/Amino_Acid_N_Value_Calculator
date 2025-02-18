import sys
import pandas as pd
import matplotlib.pyplot as plt


def main():
    args = sys.argv[1:]
    
    fig = plt.figure(figsize=(10, 10))
    
    count = 1
    for d in args:
    
        d = d.split(',')
        group = d[0][:6]
        
        # set plot colors
        if group == "Diet_A":
            color = "blue"
        elif group == "Diet_C":
            color = "red"
        elif group == "Diet_F":
            color = "green"
        else:
            color = "orange"
        
        emp_df = pd.read_csv(d[0], sep='\t', low_memory=False)
        lit_df = pd.read_csv(d[1], sep='\t', low_memory=False)
        
        emp_df = emp_df[emp_df.loc[:, "time"] > 25.0]
        lit_df = lit_df[lit_df.loc[:, "time"] > 25.0]
        
        columns = ['Protein_ID', 'Sequence', 'cf', 'mz', 'n_value']
        emp_df = emp_df[columns]
        lit_df = lit_df[columns]
        
        emp_df = emp_df.drop_duplicates()
        lit_df = lit_df.drop_duplicates()
        
        emp_df['identifier'] = emp_df['Protein_ID'] + '_' + emp_df['Sequence'] + '_' + emp_df['cf']
        lit_df['identifier'] = lit_df['Protein_ID'] + '_' + lit_df['Sequence'] + '_' + lit_df['cf']
        
        df = pd.merge(emp_df, lit_df, on=['identifier'], how='left', suffixes=('_emp', '_lit'))
        
        df['n_value_emp'] = df['n_value_emp'].astype(float).round(2)
        
        x = df['n_value_lit']
        y = df['n_value_emp']
        
        plt.subplot(2, 2, count)
        plt.scatter(x, y, color=color)
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
