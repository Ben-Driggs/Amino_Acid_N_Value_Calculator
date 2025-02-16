import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def least_squares():
    pass


def graph_scatterplot(title, x_label, x_axis, y_label, y_axis):
    plt.scatter(x_axis, y_axis)
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    return plt


def main():
    args = sys.argv[1:]
    amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    
    # read .tsv files and remove unnecessary columns
    emp_df = pd.read_csv(args[0], sep='\t', low_memory=False)
    lit_df = pd.read_csv(args[1], sep='\t', low_memory=False)
    
    columns = ['Sequence', 'n_value']
    emp_df = emp_df.loc[:, columns]
    lit_df = lit_df.loc[:, columns]
    
    # drop any rows with string values
    emp_df = emp_df[emp_df.loc[:, 'n_value'] != "no valid time points"]
    lit_df = lit_df[lit_df.loc[:, 'n_value'] != "no valid time points"]
    
    # create numpy arrays with sequences as rows, and amino acid counts as columns
    emp_aa_matrix = np.zeros((len(emp_df['Sequence'].values), len(amino_acids)), dtype=int)
    lit_aa_matrix = np.zeros((len(lit_df['Sequence'].values), len(amino_acids)), dtype=int)
    
    for i, peptide in enumerate(emp_df['Sequence'].values):
        for aa in peptide:
            if aa in amino_acids:
                emp_aa_matrix[i, amino_acids.index(aa)] += 1
                
    for i, peptide in enumerate(lit_df['Sequence'].values):
        for aa in peptide:
            if aa in amino_acids:
                lit_aa_matrix[i, amino_acids.index(aa)] += 1
                
    emp_n_values = np.array(emp_df['n_value'].values, dtype=float)
    lit_n_values = np.array(lit_df['n_value'].values, dtype=float)
    
    emp_amino_acid_values, emp_residuals, emp_rank, emp_s = np.linalg.lstsq(emp_aa_matrix, emp_n_values, rcond=None)
    print("Empirical Amino Acid Values:", emp_amino_acid_values)
    
    lit_amino_acid_values, lit_residuals, lit_rank, lit_s = np.linalg.lstsq(lit_aa_matrix, lit_n_values, rcond=None)
    print("Literature Amino Acid Values:", lit_amino_acid_values)
    
    aa_df = pd.read_csv("aa_labeling_sites.tsv", sep='\t')
    aa_nv = list(aa_df.iloc[0, :].values)
    
    lit_graph = graph_scatterplot("Diet A Literature N-Values", "Literature N-Values", aa_nv[2:-12],
                                  "Calculated Literature N-Values", lit_amino_acid_values)
    lit_graph.show()
    
    emp_graph = graph_scatterplot("Diet A Empirical N-Values", "Literature N-Values", aa_nv[2:-12],
                                  "Calculated Empirical N-Values", emp_amino_acid_values)
    emp_graph.show()
    
    sys.exit()


if __name__ == "__main__":
    main()
    
# import numpy as np
# import pandas as pd
# from sklearn.decomposition import PCA
# from sklearn.preprocessing import StandardScaler
# import matplotlib.pyplot as plt
#
# # Example DataFrame
# data = pd.DataFrame({
#     'feature1': [1, 2, 3, 4, 5],
#     'feature2': [2, 3, 4, 5, 6],
#     'feature3': [3, 4, 5, 6, 7]
# })
#
# # Standardize the data before applying PCA
# scaler = StandardScaler()
# scaled_data = scaler.fit_transform(data)
#
# # Apply PCA
# pca = PCA(n_components=2)
# principal_components = pca.fit_transform(scaled_data)
#
# # Create a DataFrame with principal components
# pca_df = pd.DataFrame(data=principal_components, columns=['PCA1', 'PCA2'])
#
# # Create the scatterplot
# plt.scatter(pca_df['PCA1'], pca_df['PCA2'])
# plt.title("PCA Scatterplot")
# plt.xlabel("PCA1")
# plt.ylabel("PCA2")
# plt.show()

    