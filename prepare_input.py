'''
Benjamin Rozonoyer
brozonoyer@brandeis.edu

This module contains methods to read in parameters from a file and split them randomly
into the desired number of groups
'''

import pandas as pd, numpy as np, scipy as sp, random, math, group_calculations as gc, visualization as vis
import matplotlib.pyplot as plt

#----------------------------------------------------------------------------------#

def read_parameters(file):
    '''
    Scans in parameter vectors from file in csv format.
    It does not fill in missing values in vectors, because that may skew the results.
    Therefore all values in the parameter vectors must be filled.
    :param file: filename
    :return: parameters: parameter vectors in the format [(parameter_id, parameter)] list
    '''

    data = pd.read_csv(file, header=None)

    # scanning of parameters
    # NOTE: depending on the parameters' labels, you will want to add constraints to which items to scan.
    # For example, we only want to scan the "item" column in the "m255_v1.csv" file, etc.
    parameters = data[[x for x in data.columns]]#.fillna(0)# if x.startswith("item")]].fillna(0)

    # convert retained parameters into list of (id, vector) to be used in algorithm
    parameter_names = parameters.columns
    parameter_vectors = parameters.T.values.tolist()
    parameters = [(parameter_names[i], parameter_vectors[i]) for i in range(len(parameter_names))]

    return parameters

#----------------------------------------------------------------------------------#

def create_random_groupings(parameters, num_groups=2, random_seed=31, logfile="logfile"):
    '''
    Groups parameters randomly into equal groups
    :param parameters: original parameters
    :param num_groups: number of groups to split the parameters into (set for entire algorithm)
    :param random_seed: random seed to reproduce results
    :param logfile: name of file to record initialization
    :return: parameters: original parameters in (parameter_id, parameter) format
             groups: list of parameter groups in the format [(group_id, group)],
                     where a group is a list of parameters in the format [(parameter_id, parameter)]
    '''

    # WRITE TO LOGFILE
    #f = open(logfile, 'a')
    #f.write("\n\n****************************************************************************\n")
    #f.write("****************************************************************************\n\n")
    #f.write("\nInitial Grouping of Parameters")
    #f.write("\n" + str(num_groups) + " Groups")
    #f.write("\n(Random Seed = " + str(random_seed) + ")")

    # shuffle parameters
    random.seed(random_seed)
    parameters = random.sample(parameters, len(parameters))

    # split parameters into groups by ids
    parameter_ids = [x_id for (x_id, x) in parameters]
    split = np.array_split(parameter_ids, num_groups)

    # find parameter corresponding to id and place (id, parameter) tuples into group
    groups = []
    for i in range(len(split)):
        parameters_group = [(x_id, x) for (x_id, x) in parameters if x_id in split[i]]
        groups.append((i, parameters_group))

    return parameters, groups