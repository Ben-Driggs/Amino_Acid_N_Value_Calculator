import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize


def apply_constraint(matrix, n_values, aa_n_values):
    def objective(x, a, B):
        return np.linalg.norm(a @ x - B)
    
    # Constraints function to ensure values are non-zero
    def constraint(x):
        return x - 1e-6  # Example constraint to ensure values are not zero
    
    # Define constraints dictionary
    constraints = {'type': 'ineq', 'fun': constraint}
    
    # Perform optimization
    result = minimize(objective, aa_n_values, args=(matrix, n_values), constraints=constraints)
    
    # Print the result
    print("Solution:", result.x)
    return aa_n_values


def graph_scatterplot(title, x_label, x_axis, y_label, y_axis, color, count):
    plt.subplot(4, 2, count)
    plt.scatter(x_axis, y_axis, color=color)
    plt.title(title)
    plt.xlabel(x_label, labelpad=5)
    plt.ylabel(y_label, labelpad=5)
    return plt


def main():
    args = sys.argv[1:]
    no_mods = True
    fig = plt.figure(figsize=(10, 10))
    
    if no_mods:
        amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W',
                       'Y']
    else:
        amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W',
                       'Y',
                       'm', 'q', '1', '2', '3', '4', 's', 't', 'y', 'k', 'c', 'o']
    
    count = 1
    for d in args:
        d = d.split(',')
        emp = d[0]
        lit = d[1]
        
        # read .tsv files and remove unnecessary columns
        emp_df = pd.read_csv(emp, sep='\t', low_memory=False)
        lit_df = pd.read_csv(lit, sep='\t', low_memory=False)
        
        # determine diet group
        group = emp[0:6]
        
        # set plot colors
        if group == "Diet_A":
            color = "blue"
        elif group == "Diet_C":
            color = "red"
        elif group == "Diet_F":
            color = "green"
        else:
            color = "orange"
        
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
        
        emp_Q, emp_R = np.linalg.qr(emp_aa_matrix)
        lit_Q, lit_R = np.linalg.qr(lit_aa_matrix)
        
        emp_Q_T_b = np.dot(emp_Q.T, emp_n_values)
        lit_Q_T_b = np.dot(lit_Q.T, lit_n_values)
        
        emp_amino_acid_values = np.linalg.solve(emp_R, emp_Q_T_b)
        lit_amino_acid_values = np.linalg.solve(lit_R, lit_Q_T_b)
        
        # this didn't really do anything, so I'm not using it.
        # lit_amino_acid_values = apply_constraint(lit_aa_matrix, lit_n_values, lit_amino_acid_values)
        
        aa_df = pd.read_csv("aa_labeling_sites.tsv", sep='\t')
        aa_nv = list(aa_df.iloc[0, :].values)
        
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
    
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
    plt.tight_layout()
    plt.show()
    
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

# def objective(x, A, b):
#     return np.linalg.norm(A @ x - b)
#
# # Constraints function to ensure values are non-zero
# def constraint(x):
#     return x - 1e-6  # Example constraint to ensure values are not zero
#
# # Example data
# A = np.array([[1, 2], [3, 4], [5, 6]])
# b = np.array([7, 8, 9])
#
# # Initial guess
# x0 = np.ones(A.shape[1])
#
# # Define constraints dictionary
# constraints = {'type': 'ineq', 'fun': constraint}
#
# # Perform optimization
# result = minimize(objective, x0, args=(A, b), constraints=constraints)
#
# # Print the result
# print("Solution:", result.x)

# import matplotlib.pyplot as plt
# import numpy as np
#
# # Example data
# x = np.array([1, 2, 3, 4, 5])
# y = np.array([2, 3, 5, 7, 11])
# yerr = np.array([0.5, 0.4, 0.6, 0.7, 0.3])  # Example error values
#
# plt.figure()
#
# # Create scatter plot with error bars
# plt.errorbar(x, y, yerr=yerr, fmt='o', capsize=5, capthick=1, ecolor='red')
#
# # Customize the plot
# plt.title("Scatter Plot with Error Bars")
# plt.xlabel("X-axis Label")
# plt.ylabel("Y-axis Label")
#
# # Display the plot
# plt.show()
