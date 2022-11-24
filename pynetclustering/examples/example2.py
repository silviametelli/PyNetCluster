import pynetclustering as pync
import os

## sequential-clustering

#set your working directory
os.chdir("/Users/silvia/Desktop/PyNetCluster/pynetclustering/examples/")
os.getcwd()


pync.run_clustering("edgelist.txt", edgelist=True, sequential=True)


# if you need to see the sorted unique edge-list, use the following:
pync.get_data("edgelist.txt") # > outputs a datafile called "sorted_edgelist.txt": unique list of nodes with associated connections
pync.run_clustering("connected_nodes.txt", edgelist=False, sequential=True) # equivalent to pync.run_clustering("edgelist.txt", edgelist=True, sequential=True))
