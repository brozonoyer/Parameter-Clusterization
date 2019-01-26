'''
Benjamin Rozonoyer
brozonoyer@brandeis.edu

This module contains Braverman's modified algorithm, and also is the main module to
run the algorithm.
'''

import numpy as np
from numpy import linalg as LA
import pandas as pd
import scipy as sp
import random, math
import argparse, pickle
import csv
from scipy import linalg
import prepare_input as pi
import group_calculations as gc
import os
import visualization as vis
import matplotlib
import time

#----------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------#

class Train():

# The main algorithm class, which allows us to create an object to run Braverman's modified algorithm
# on given data, setting the relevant parameters (see initialization) as required. This class
# contains the algorithm itself, and calls on methods from the other modules.

#----------------------------------------------------------------------------------#

    def __init__(self, train_filename, num_groups=2, eigenvalues=1, logfile="logfile", roundoff=10):
        '''
        initialization method
        :param filename: file (in csv format) of parameter vectors, in "Input" folder
        :param num_groups: number of groups into which to split the parameters (default is 2)
        :param eigenvalues: how many of the maximum eigenvalues in a given group's correlation
                            matrix to sum when computing the functional value for each group
        :param logfile: text file to write computation and results
        :param roundoff: how many decimals to round off eigenvalue computations to (to eliminate
                         discrepancies when running same calculations in other platforms)
        '''

        # WRITE TO LOGFILE
        #self.logfile = open(logfile, 'a')
        #self.logfile.write("\n\n\n+++++++++++++++++++++++++++++++++++++++++++++++\n")
        #self.logfile.write("PARAMETERS FROM FILE \"" + filename + "\"")

        # decimals to round off functional's value
        self.roundoff = roundoff
        # how many eigenvalues to use when calculating functional value for group
        self.eigenvalues = eigenvalues
        # count of iterations (until convergence)
        self.iterations = 0
        # read in parameters from specified file
        self.sorted_parameters = pi.read_parameters(train_filename)
        self.alphabetized_parameter_ids = [id for (id, _) in self.sorted_parameters]

        # randomly groups parameters
        self.parameters, self.groups = pi.create_random_groupings(self.sorted_parameters, num_groups, logfile=logfile)

#----------------------------------------------------------------------------------#

    def iterate(self):
        '''
        Braverman iterative algorithm, pp.125-126
        Iterates over parameters and improves groupings until convergence
        :return: factors_with_groups: the resulting parameter groupings, and a factor corresponding
                                      to each of them, computed by Braverman formula (3)
        '''

        changes_made_during_cycle = True

        # iterate over parameters, regrouping them as long as regrouping increases functional value
        while changes_made_during_cycle:

            J = round(self.compute_total_functional_value(), self.roundoff)
            changes_made_during_cycle = False

            # WRITE TO LOGFILE
            #self.logfile.write('\n\n\n*************************************\n')
            #self.logfile.write("\nIteration " + str(self.iterations))
            #self.logfile.write("\nJ = " + str(J) + '\n')

            for tuple in self.parameters:

                # WRITE TO LOGFILE
                #self.logfile.write('\n---------------------------------------')
                #self.logfile.write("\nVariable \"" + tuple[0] + '\"')

                # for each parameter, obtain the group from which it was removed,
                # the group to which it was added, and true if transferred
                # (the result is None, None, and False if no change is made for parameter)
                (decreased_group_id, decreased_group),\
                (increased_group_id, increased_group),\
                                            changes_made_for_variable = self.transfer(tuple)

                # if parameter is transferred, then change has been made during iteration and groupings haven't yet converged
                if changes_made_for_variable:

                    changes_made_during_cycle = True

                    # update the groups according to changes
                    for i in range(len(self.groups)):
                        (id, group) = self.groups[i]
                        if id == decreased_group_id:
                            self.groups[i] = (decreased_group_id, decreased_group)
                        elif id == increased_group_id:
                            self.groups[i] = (increased_group_id, increased_group)

            self.iterations += 1

        # WRITE TO LOGFILE
        #self.logfile.write("\n\nCONVERGED\n")
        #self.logfile.write('\nIterations to Convergence:\t' + str(self.iterations) + '\n')


        ### AFTER CONVERGENCE ###


        factors_with_groups = []
        for (id, group) in self.groups:

            # WRITE TO LOGFILE
            #self.logfile.write("\n\nGroup " + str(id) + ": Length " + str(len(group)))
            #self.logfile.write('\n' + str(sorted([parameter_id for (parameter_id, parameter) in group])) + '\n')
            #self.logfile.write('\n' + str(sorted([parameter_id for (parameter_id, parameter) in group])) + '\n')

            if group:   # for each non-empty group, compute factor according to Braverman formula (3)
                factor, eigenvector = gc.compute_factor_for_group(group)
                factors_with_groups.append((id, factor, group))

                # WRITE TO LOGIFILE
                #self.logfile.write('\nGroup ' + str(id) + ' Factor\n' + str(factor))

        return factors_with_groups

#----------------------------------------------------------------------------------#

    def transfer(self, tuple):
        '''
        For given parameter x, find its group, then find the value of the functional's
        component (Braverman formula (2) p.125) for the group minus x. Then loop through the other
        groups. For each other group, find the value of the functional's component for that group
        plus x. If the net change in the functional (Braverman formula (1) p.125) is positive for any of the groups,
        make the transfer and return true; else, return false
        :param tuple: parameter that possibly needs to be transferred in the tuple format (id, parameter)
        :return: boolean indicating whether x has been transferred to another group
        '''

        # extract parameter and its id from the tuple
        x_id, x = tuple[0], tuple[1]

        # find group that x (current parameter) belongs to
        current_id, current_group = self.find_group(tuple)

        # value of component before x is removed
        SIGMA_1 = round(gc.compute_functional_value_for_group(current_group, self.eigenvalues), self.roundoff)

        # group with x removed
        current_group_minus_x = [(v_id, v) for (v_id, v) in current_group if v_id != x_id]

        # value of component after x is removed
        DELTA_1 = round(gc.compute_functional_value_for_group(current_group_minus_x, self.eigenvalues), self.roundoff)

        # decrease in component's value resulting from removal of x (negative if decreased)
        decrease = DELTA_1 - SIGMA_1

        # WRITE TO LOGFILE
        #self.logfile.write("\nCurrent Group\t\t\t" + str(len(current_group)))
        #self.logfile.write("\nFunctional Value for Current Group\t" + str(SIGMA_1))
        #self.logfile.write("\nCurrent Group Minus x\t" + str(len(current_group_minus_x)))
        #self.logfile.write("\nFunctional Value for Current Group Minus x\t" + str(DELTA_1))
        #self.logfile.write("\ndecrease\t" + str(decrease))

        # loop through other groups
        for (other_id, other_group) in self.groups:
            if other_id != current_id:

                # value of component before x is added
                DELTA_2 = round(gc.compute_functional_value_for_group(other_group, self.eigenvalues), self.roundoff)

                # group after x is added
                other_group_plus_x = other_group.copy()
                other_group_plus_x.append((x_id, x))

                # value of component after x is added
                SIGMA_2 = round(gc.compute_functional_value_for_group(other_group_plus_x, self.eigenvalues), self.roundoff)

                # increase in component's value after x is added
                increase = SIGMA_2 - DELTA_2

                # the net change in the functional's value
                NET_CHANGE = decrease + increase

                # WRITE TO LOGFILE
                #self.logfile.write("\nOther Group\t\t\t\t" + str(len(other_group)))
                #self.logfile.write("\nFunctional Value for Other Group\t" + str(DELTA_2))
                #self.logfile.write("\nFunctional Value for Other Group Plus x\t" + str(SIGMA_2))
                #self.logfile.write("\nincrease\t" + str(increase))
                #self.logfile.write("\nNet Change in J:\t\t" + str(NET_CHANGE))

                # if the transfer increases the total functional value, keep transfer
                if NET_CHANGE > 0:
                    return (current_id, current_group_minus_x),\
                           (other_id, other_group_plus_x),\
                           True # transfer has been made

        # otherwise don't make transfer
        return (None, None), (None, None), False # no transfers have been made

#----------------------------------------------------------------------------------#

    def find_group(self, tuple):
        '''
        Finds the grouping to which the parameter belongs at the current iteration
        :param tuple: parameter whose group needs to be found in tuple format (id, parameter)
        :return: group to which the parameter belongs
        '''

        # search through groups
        for (group_id, group) in self.groups:
            # search through parameters within each group
            for (x_id, x) in group:
                # return group if parameter_id matches
                if x_id == tuple[0]:
                    return group_id, group

#----------------------------------------------------------------------------------#

    def compute_total_functional_value(self):
        '''
        Computes the value of the functional J as the sum of sums of the indicated
        number of greatest eigenvalues for the correlation matrix of each group of
        parameters (c.f. corresponding method in group_calculations.py module):
        i.e. for k groups and n maximum eigenvalues to sum, the functional J is:
        J = (group_1_lambda_1 + ... + group1_lambda_n) + ... + (group_k_lambda_1 + ... + group_k_lambda_n)
        :return: value of J
        '''
        J = 0
        for (group_id, group) in self.groups:
            # sums the functional's component for each group of parameters
            J += gc.compute_functional_value_for_group(group, self.eigenvalues)
        return J

#----------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------#

def output_training_info(input_filename, output_filename, model_pickle_file, num_groups, eigenvalues, logfile, roundoff):
    """
    :param input_filename: training parameters file (in csv format)
    :param output_filename: where to write info about extremal groupings
    :param model_pickle_file: where to pickle group factors resulting from extremal groupings
    :param num_groups:
    :param eigenvalues:
    :param logfile:
    :param roundoff:
    :return:
    """

    # creates an instance of the Train class, with parameters read in from command line
    T = Train(train_filename=input_filename, num_groups=num_groups, eigenvalues=eigenvalues, logfile=logfile, roundoff=roundoff)
    factors_with_groups = T.iterate() # returns data structure with all the grouping information

    #print(factors_with_groups)

    # data structures for writing information to output file
    group2parameters = {triple[0]: tuple([tup[0] for tup in triple[2]]) for triple in factors_with_groups}
    parameter2group = {}
    for parameter_id in T.alphabetized_parameter_ids:
        for group_id in group2parameters:
            if parameter_id in group2parameters[group_id]:
                parameter2group[parameter_id] = group_id
    grouping_list = [parameter2group[parameter_id] for parameter_id in T.alphabetized_parameter_ids]
    group2factor = {group_id:group_factor for (group_id, group_factor, _) in factors_with_groups}
    factors = list(group2factor.values())
    weight_list = [np.inner(group2factor[parameter2group[parameter_id]], parameter) for (parameter_id, parameter) in T.sorted_parameters]
    #print(factors)
    # write to Extremal Groups Info File
    with open(output_filename, "a+") as f:
        f.write("Extremal Groups Info\nBuilt on features\n")
        f.write(",".join(T.alphabetized_parameter_ids))
        f.write("\nnParams," + str(len(T.alphabetized_parameter_ids)))
        f.write("\nnGroups," + str(num_groups))
        f.write("\ninGroup, group ids start from 0\n")
        f.write(",".join([str(x) for x in grouping_list]))
        f.write("\nWeights\n")
        f.write(",".join([str(x) for x in weight_list]))

    # pickle model (consisting of factors for the groups)
    with open(model_pickle_file, "wb") as f:
        pickle.dump(factors, f)

#----------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------#

def output_inference_info(test_filename, model, output_file):

    parameters = pi.read_parameters(test_filename)
    with open(model, "rb") as f:
        factors = pickle.load(f)
    print(factors)
    groups = [list() for factor in factors]

    # sort each parameter into a group corresponding to the factor with which it most closely correlates
    for (parameter_id, parameter) in parameters:
        best_factor_index = -1
        best_correlation = -float('inf')
        for i in range(len(factors)):
            print("i", i)
            print("parameter", parameter)
            new_correlation = np.inner(factors[i], parameter)
            print("correlation at i", new_correlation)
            if new_correlation > best_correlation:
                best_factor_index = i
                best_correlation = new_correlation
        groups[best_factor_index].append((parameter_id, parameter))

    print(len(groups[1]))

    # determine factors for each group of new parameters
    group_factors = []
    for group in groups:
        factor, eigenvector = gc.compute_factor_for_group(group)
        group_factors.append(factor)

    transposed_group_factors = list(map(list, zip(*group_factors)))

    # write as csv file
    with open(output_file, "w", newline='') as f:
        writer = csv.writer(f, delimiter=',')
        for row in transposed_group_factors:
            writer.writerow(row)

#----------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------#


if __name__ == '__main__':

    #start_time = time.time()

    # required and optional arguments to program
    parser = argparse.ArgumentParser()
    parser.add_argument("input_filename", help="input file (train if step 1, test if step 2)")
    parser.add_argument("step", type=int, help="step 1 = training\nstep 2 = inference")
    parser.add_argument("output_filename", help="output file")
    parser.add_argument("model_pickle_file", help="filepath to pickled model (empty if step 1)")
    parser.add_argument("--num_groups", type=int, help="number of parameter groupings")
    parser.add_argument("--eigenvalues", type=int, help="number of group correlation matrix eigenvalues to use for calculations of functional's value")
    parser.add_argument("--logfile", help="logfile to store information of run")
    parser.add_argument("--roundoff", type=int, help="roundoff for numbers in calculations (default set to 10)")

    # consumes command line parameters, and sets default values if optional parameters aren't provided
    args = parser.parse_args()
    input_filename = args.input_filename
    output_filename = args.output_filename
    model_pickle_file = args.model_pickle_file
    if args.step == 1: # if training, then set required parameters
        if args.num_groups:
            num_groups = args.num_groups
        else:
            num_groups = 2
        if args.eigenvalues:
            eigenvalues = args.eigenvalues
        else:
            eigenvalues = 1
        if args.roundoff:
            roundoff = args.roundoff
        else:
            roundoff = 1
        if args.logfile:
            logfile = args.logfile
        else:
            logfile = 1

        output_training_info(input_filename, output_filename, model_pickle_file, num_groups, eigenvalues, logfile, roundoff)

    elif args.step == 2: # if testing, then get output file
        output_inference_info(input_filename, model_pickle_file, output_filename)

    #end_time = time.time()
    #print("Time:", end_time - start_time)

