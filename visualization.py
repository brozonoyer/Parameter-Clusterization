import matplotlib.pyplot as plt
import numpy as np

#----------------------------------------------------------------------------------#

def histogram(observations):
    '''
    Plots histogram of observations for a given parameter vector
    :param observations: parameter vector (list), where each element is an observation
    '''

    plt.hist(observations)
    plt.title("Parameter Clusterization Histogram")
    plt.xlabel("Projection onto Factor")
    plt.ylabel("Frequency")
    plt.show()

#----------------------------------------------------------------------------------#

def scatterplot(X, Y):
    '''
    Given two factors, plots them as x and y coordinates
    :param X:
    :param Y:
    :return:
    '''

    plt.scatter(X,Y)
    plt.title("Parameter Clusterization Scatterplot")
    plt.xlabel("Projection onto Factor 1")
    plt.ylabel("Projection onto Factor 2")
    plt.show()

#----------------------------------------------------------------------------------#

def visualize(observations, plot="histogram"):

    if plot == "histogram":
        histogram(observations)
    elif plot == "scatterplot":
        pass