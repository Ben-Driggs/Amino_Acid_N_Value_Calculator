import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares


def graph_scatterplot(title, x_label, x_axis, y_label, y_axis, color, count):
    plt.subplot(4, 2, count)
    plt.scatter(x_axis, y_axis, color=color)
    plt.title(title)
    plt.xlim(0, 5)
    plt.ylim(-1, 10)
    plt.xlabel(x_label, labelpad=5)
    plt.ylabel(y_label, labelpad=5)
    return plt


def main():
    mods = ['m', 'q', '1', '2', '3', '4', 's', 't', 'y', 'k', 'c', 'o']
    args = sys.argv[1:]
    results = {}
    
    fig = plt.figure(figsize=(10, 10))
    std_distances = {}
    
    count = 1
    for d in args:
        
        d = d.split(',')
        group = d[0][:6]
        # group = "Diet_A"
        
        # set plot colors
        if group == "Diet_A":
            color = "grey"
            color2 = "blue"
            # color3 = "navy"
            # color4 = "deepskyblue"
        elif group == "Diet_C":
            color = "grey"
            color2 = "red"
            # color3 = "darkred"
            # color4 = "darkorchid"
        elif group == "Diet_F":
            color = "grey"
            color2 = "green"
            # color3 = "darkgreen"
            # color4 = "mediumseagreen"
        else:
            color = "grey"
            color2 = "orange"
            # color3 = "crimson"
            # color4 = "mediumvioletred"
        
        emp_df = pd.read_csv(d[0], sep='\t', low_memory=False)
        lit_df = pd.read_csv(d[1], sep='\t', low_memory=False)
        
        emp_df = emp_df[emp_df.loc[:, "time"] > 25.0]
        lit_df = lit_df[lit_df.loc[:, "time"] > 25.0]
        
        columns = ['Protein ID', 'Sequence', 'cf', 'mz', 'n_value', 'abundances', 'n_value_stddev']
        l_columns = ['Protein ID', 'Sequence', 'cf', 'mz', 'abundances', 'n_value']
        # columns = ['Protein ID', 'sequence', 'Elemental Composition', 'n_value', 'stddev']
        # l_columns = ['Protein ID', 'sequence', 'Elemental Composition', 'number of possible labeling sites', 'stddev']
        emp_df = emp_df[columns]
        lit_df = lit_df[l_columns]
        
        # merge dataframes
        emp_df['id'] = emp_df['Protein ID'] + emp_df['Sequence'] + str(emp_df['mz']) + emp_df['cf']
        lit_df['id'] = lit_df['Protein ID'] + lit_df['Sequence'] + str(lit_df['mz']) + lit_df['cf']
        emp_df = emp_df.drop(columns=['Protein ID', 'cf', 'mz'])
        lit_df = lit_df.drop(columns=['Protein ID', 'Sequence', 'cf', 'mz', 'abundances'])
        merged_df = pd.merge(emp_df, lit_df, how="outer", on=['id'], suffixes=('_emp', '_lit'))
        
        merged_df = merged_df.drop_duplicates()
        
        # remove any rows with modified peptides
        e_mask = merged_df['Sequence'].str.contains('|'.join(mods))
        merged_df = merged_df[~e_mask]
        
        # filter by abundance
        # sum abundances and take top 50%
        merged_df.loc[:, 'sum_abundances'] = merged_df['abundances'].apply(
            lambda v: sum([float(value) for value in v[1:-1].split(', ')]))

        emp_threshold = merged_df['sum_abundances'].quantile(0.50)

        # merged_df = merged_df[merged_df['sum_abundances'] >= emp_threshold]
        
        # filter out noise by setting a limit on n_value standard deviation
        merged_df = merged_df[merged_df['n_value_lit'] != "no valid time points"]
        merged_df = merged_df[merged_df['n_value_emp'] != "no valid time points"]
        merged_df['n_value_stddev'] = merged_df['n_value_stddev'].astype(float)
        merged_df['n_value_emp'] = merged_df['n_value_emp'].astype(float)
        merged_df = merged_df[merged_df['n_value_emp'] <= 100]
        filtered_df = merged_df[merged_df.loc[:, 'n_value_stddev'] <= 0.05]
        
        cols = ['n_value_emp', 'n_value_lit']
        for c in cols:
            merged_df[c] = merged_df[c].astype(float).round(2)
            filtered_df[c] = filtered_df[c].astype(float).round(2)
        
        x = merged_df['n_value_lit']
        y = merged_df['n_value_emp']
        f_x = filtered_df['n_value_lit']
        f_y = filtered_df['n_value_emp']
        
    #     plt.subplot(2, 2, count)
    #     plt.scatter(x, y, color=color, alpha=0.50, label='Unfiltered')
    #     coeffs1 = np.polyfit(x, y, 1)
    #     fit1 = np.poly1d(coeffs1)
    #     plt.plot(x, fit1(x), color="darkgrey", alpha=0.5, linewidth=2, label="Unfiltered Linear Fit")
    #
    #     distances = np.abs(y - x) / np.sqrt(2)
    #     std_distance = np.std(distances)
    #     std_distances[group] = std_distance
    #
    #     d_mask = distances <= std_distance
    #     filt_x = x[d_mask]
    #     filt_y = y[d_mask]
    #
    #     # plt.scatter(filt_x, filt_y, color=color2, alpha=0.2, label='Filtered')
    #     # coeffs2 = np.polyfit(filt_x, filt_y, 1)
    #     # fit2 = np.poly1d(coeffs2)
    #     # plt.plot(filt_x, fit2(filt_x), color=color2, linewidth=2, label="Filtered Linear Fit")
    #
    #     plt.scatter(f_x, f_y, color=color2, alpha=0.5, label='Filtered')
    #     coeffs2 = np.polyfit(f_x, f_y, 1)
    #     fit2 = np.poly1d(coeffs2)
    #     plt.plot(f_x, fit2(f_x), color=color2, linewidth=2, label="Filtered Linear Fit")
    #
    #     upper = [xs + (std_distance * np.sqrt(2)) for xs in list(range(0, 100))]
    #     lower = [xs - (std_distance * np.sqrt(2)) for xs in list(range(0, 100))]
    #
    #     # plt.plot([0, 100], [upper[0], upper[99]], color=color2, linestyle='--', alpha=0.5)
    #     # plt.plot([0, 100], [lower[0], lower[99]], color=color2, linestyle='--', alpha=0.5)
    #     plt.plot([0, 100], [0, 100], color='black', alpha=0.75, linestyle='--')
    #
    #     # plt.xlim(0, 100)
    #     # plt.ylim(0, 100)
    #     plt.title(f"{group} - Literature vs Empirical Peptide N-Values")
    #     plt.xlabel("Peptide Literature N-Values")
    #     plt.ylabel("Peptide Empirical N-Values")
    #     count += 1
    #
    # plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
    # plt.tight_layout()
    # plt.show()
    
    # with open("std_distances.txt", 'w') as outf:
    #     for g, value in std_distances.items():
    #         outf.write(f"{g}, {value}\n")
    
        no_mods = True
        if no_mods:
            amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W',
                           'Y']
        else:
            amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W',
                           'Y',
                           'm', 'q', '1', '2', '3', '4', 's', 't', 'y', 'k', 'c', 'o']
        
        # set plot colors
        if group == "Diet_A":
            color = "blue"
        elif group == "Diet_C":
            color = "red"
        elif group == "Diet_F":
            color = "green"
        else:
            color = "orange"
            
        # create numpy arrays with sequences as rows, and amino acid counts as columns
        emp_aa_matrix = np.zeros((len(f_x), len(amino_acids)), dtype=int)
        lit_aa_matrix = np.zeros((len(f_y), len(amino_acids)), dtype=int)
        
        for i, peptide in enumerate(filtered_df['Sequence'].values):
            for aa in peptide:
                if aa in amino_acids:
                    emp_aa_matrix[i, amino_acids.index(aa)] += 1
        
        for i, peptide in enumerate(filtered_df['Sequence'].values):
            for aa in peptide:
                if aa in amino_acids:
                    lit_aa_matrix[i, amino_acids.index(aa)] += 1
        
        emp_n_values = np.array(f_y, dtype=float)
        lit_n_values = np.array(f_x, dtype=float)
        
        # https://numpy.org/doc/stable/reference/generated/numpy.linalg.lstsq.html#numpy-linalg-lstsq
        # emp_amino_acid_values, emp_residuals, emp_rank, emp_s = np.linalg.lstsq(emp_aa_matrix, emp_n_values, rcond=None)
        # print("Empirical Amino Acid Values:", emp_amino_acid_values)

        lit_amino_acid_values, lit_residuals, lit_rank, lit_s = np.linalg.lstsq(lit_aa_matrix, lit_n_values, rcond=None)
        print("Literature Amino Acid Values:", lit_amino_acid_values)
        
        # this didn't really do anything, so I'm not using it.
        # emp_amino_acid_values = apply_constraint(emp_aa_matrix, emp_n_values, emp_amino_acid_values)
        # lit_amino_acid_values = apply_constraint(lit_aa_matrix, lit_n_values, lit_amino_acid_values)
        
        aa_df = pd.read_csv("aa_labeling_sites.tsv", sep='\t')
        aa_nv = list(aa_df.iloc[0, :].values)
        
        x0 = np.array(aa_nv[2:-12])
        
        def residuals(xs, A, b):
            return np.dot(A, xs) - b
        
        result = least_squares(residuals, x0, args=(emp_aa_matrix, emp_n_values))
        emp_amino_acid_values = result.x
        print('Empirical Amino Acid Values:', emp_amino_acid_values)
        
        if no_mods:
            lit_graph = graph_scatterplot(f"{group} Literature N-Values", "Literature N-Values", aa_nv[2:-12],
                                          "Estimated N-Values", lit_amino_acid_values, color, count)
            count += 1
            emp_graph = graph_scatterplot(f"{group} Empirical N-Values", "Literature N-Values", aa_nv[2:-12],
                                          "Empirical N-Values", emp_amino_acid_values, color, count)
        else:
            lit_graph = graph_scatterplot(f"{group} Literature N-Values", "Literature N-Values", aa_nv[2:],
                                          "Estimated N-Values", lit_amino_acid_values, color, count)
            count += 1
            emp_graph = graph_scatterplot(f"{group} Empirical N-Values", "Literature N-Values", aa_nv[2:],
                                          "Empirical N-Values", emp_amino_acid_values, color, count)
        count += 1
        results[f"{group}_empirical_n_value"] = [round(v, 2) for v in emp_amino_acid_values]
        results[f"{group}_literature_n_value"] = [round(v, 2) for v in lit_amino_acid_values]
        
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
    plt.tight_layout()
    # plt.show()
    
    amino_acid_names = [
        'Alanine',
        'Cysteine',
        'Aspartic Acid',
        'Glutamic Acid',
        'Phenylalanine',
        'Glycine',
        'Histidine',
        'Isoleucine',
        'Lysine',
        'Leucine',
        'Methionine',
        'Asparagine',
        'Proline',
        'Glutamine',
        'Arginine',
        'Serine',
        'Threonine',
        'Valine',
        'Tryptophan',
        'Tyrosine'
    ]
    
    df = pd.DataFrame(data=results, index=amino_acid_names)
    df.to_csv(path_or_buf="table_data/amino_acid_n_values", sep='\t')
    
    sys.exit()


if __name__ == "__main__":
    main()

# filter out high standard deviations
# make tablular form of the AA n-values
# recalculate peptide n-values from empirical AA n-values
# filter top 50% abundant peptides
