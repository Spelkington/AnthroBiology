from random import expovariate
from pgen import poidev
import plotly.plotly as py
import plotly.graph_objs as go

populationIntervals=[]
frequencies=[]

observedMutations = 82
experimentRepetitions = 100
frequency = 0.0

for j in range(4500,20000,100):
    for i in range(experimentRepetitions):

        mutationChance = 0.001      # Percent chance of mutation
        geneCopies = 77             # sample size (number of gene copies)
        twoN = j                 # population size (number of gene copies)

        treeDepth = 0.0            # age of last common ancestor in generations

        treeLength = 0.0           # current tree length

        # Each pass through loop deals with one coalescent interval.
        while geneCopies > 1:

            # determine hazard of a coalescent event
            hazard = geneCopies *(geneCopies - 1)/(2.0 * twoN)

            # determine time until next coalescent event
            timeInterval = expovariate(hazard)

            # depth: how much time has passed within the tree
            treeDepth += timeInterval
            # Record gene-hours(?) within the tree
            treeLength += timeInterval * geneCopies
            # Decrease the number of gene copies by 1
            geneCopies -= 1

        # Generate poisson distrib. with mean u * L
        treeMutations = poidev(mutationChance * treeLength)
        # If our experiments' mutations are less than or equal to the observed
        # mutations,
        if treeMutations <= observedMutations:
            frequency += 1  # Keep track of how many times that has happened.

    frequency /= float(experimentRepetitions)   # Calculate true frequency of
                                                # underperformance

    print "F[%d] = %6.3f for hypothesis: 2N=%i" % (observedMutations, frequency, j)

    populationIntervals.append(j)
    frequencies.append(frequency)

data = [go.Bar(
            x=populationIntervals,
            y=frequencies
        )]

py.iplot(data, filename='basic-bar')
