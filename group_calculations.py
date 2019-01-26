'''
Benjamin Rozonoyer
brozonoyer@brandeis.edu

This module contains methods dealing with correlation matrices, eigenvalues, and
matrix multiplication required over the course of the algorithm.
'''



import numpy as np, scipy as sp, pandas as pd, math
from scipy import linalg as LA

#----------------------------------------------------------------------------------#

def compute_factor_for_group(group):
    '''
    Braverman formula (3) p.125
    Takes a group of parameters A_l; creates a correlation matrix R_l = {(x_i, x_j)},
    where x_i and x_j are in A_l; extracts the eigenvector corresponding to the maximum
    eigenvalue, and computes the factor for that group as a linear combination of x_i
    parameter with the corresponding a_i element of the eigenvector.
    :param group: list of parameters
    :return: factor
    '''

    # create correlation matrix out of parameters in group
    group_matrix = np.array([x for (x_id, x) in group]).transpose()
    correlation_matrix = create_correlation_matrix(group)
    eigenvalues, eigenvectors = extract_maximum_eigenvalues(correlation_matrix, 1)

    # factor as linear combination of eigenvector's elements with group's parameters (equivalent to matrix multiplication)
    factor = np.nan_to_num(np.matmul(group_matrix, eigenvectors[:,-1]))

    # normalize to length 1
    #factor = factor / np.linalg.norm(factor)
    print(factor)
    return factor, eigenvectors[:,-1]

#----------------------------------------------------------------------------------#

def compute_functional_value_for_group(group, number_of_eigenvalues):
    '''
    For given group of parameters, extract n maximum eigenvalues from its correlation
    matrix and return its sum -- this is the component of the total functional J's value
    :param group: list of parameters in format [(parameter_id, parameter)]
    :param number_of_eigenvalues: number of highest eigenvalues to sum, specified in algorithm's initialization
    :return: value: sum of given number of highest eigenvalues
    '''

    # if empty group, its functional value is simply 0
    if not group:
        return 0

    # if group is not empty
    correlation_matrix = create_correlation_matrix(group)
    eigenvalues, eigenvectors = extract_maximum_eigenvalues(correlation_matrix, number_of_eigenvalues)
    return sum(eigenvalues)

#----------------------------------------------------------------------------------#

def extract_maximum_eigenvalues(C, number_of_eigenvalues):
    '''
    Inputs a correlation matrix and outputs a given number of its maximum eigenvalues
    with their corresponding eigenvectors
    :param C: correlation matrix
    :param number_of_eigenvalues: number of max eigenvalues to be returned
    :return: list of maximum eigenvalues
    '''

    #print(C)
    # how many maximum eigenvalues to use
    eigenvalue_range = (max(0, len(C) - number_of_eigenvalues), len(C) - 1)

    # extract eigenvalues and corresponding eigenvectors
    eigenvalues, eigenvectors = LA.eigh(C, eigvals=eigenvalue_range)
    return eigenvalues, eigenvectors

#----------------------------------------------------------------------------------#

def create_correlation_matrix(group):
    '''
    Creates Pearson correlation matrix for parameters in group, using Pandas DataFrame
    :param group: grouped parameters in format [(parameter_id, parameter)]
    :return: correlation matrix for parameters
    '''

    # create matrix (DataFrame) out of parameters in group
    group = pd.DataFrame([x for (x_id, x) in group], index=[x_id for (x_id, x) in group]).T.fillna(0)
    correlation_matrix = group.corr(method='pearson')
    return correlation_matrix.values