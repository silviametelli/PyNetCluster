import pynetclustering as pync
import os

## batch-clustering

#set your working directory
os.chdir("/Users/silvia/Desktop/PyNetCluster/pynetclustering/examples/")
os.getcwd()


pync.run_clustering("edgelist.txt", edgelist=True)


# if you need to see the sorted unique edge-list, use the following:
# pync.get_data("edgelist.txt") # > outputs a datafile called "sorted_edgelist.txt": unique list of nodes with associated connections
pync.run_clustering("connected_nodes.txt", edgelist=False) # equivalent to pync.run_clustering("edgelist.txt", edgelist=True)


# output: [('U1', 1), ('U2', 2), ('U3', 3), ('U4', 3), ('U5', 4), ('U6', 5)]