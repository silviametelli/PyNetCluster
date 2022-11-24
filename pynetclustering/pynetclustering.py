"""Created on Thu Sept 22 10:10:01 2022
@author: Silvia Metelli
"""
import operator
from math import log, lgamma
from collections import Counter
import sys
import os



def cluster_likel(n, x, alpha=1, beta=1):
    if x > n:
        sys.exit("Error: x>n\n")
    return lgamma(alpha + beta) + lgamma(alpha + x) + lgamma(beta + n - x) - lgamma(alpha) - \
           lgamma(beta) - lgamma(alpha + beta + n)


def cluster_likelihood(num_subjects, feature_counter):
    if len(feature_sizes) == 0 and len(feature_prior_counts) == 0:
        feature_count_freqs = Counter(feature_counter.values())
        n0 = num_features - len(feature_counter)
        likelihood = n0 * cluster_likel(num_subjects, 0) if n0 > 0 else 0
        for x in feature_count_freqs:
            likelihood += feature_count_freqs[x] * cluster_likel(num_subjects, x)
        return likelihood
    likelihood = 0
    n0 = num_features
    for f in feature_counter:
        feature_cluster_size = 1 if f not in feature_sizes else feature_sizes[f]
        feature_cluster_prior_count = int(
            num_subjects * feature_cluster_size / 2.0) if f not in feature_prior_counts else feature_prior_counts[f]
        p = (feature_cluster_prior_count + 1) / float((total_num_subjects + 2) * feature_cluster_size)
        feature_cluster_alpha, feature_cluster_beta = beta_parameters(num_features, p)
        likelihood += cluster_likel(num_subjects * feature_cluster_size, feature_counter[f], feature_cluster_alpha,
                                    feature_cluster_beta)
        n0 -= feature_cluster_size
    for f in feature_sizes:
        if f not in feature_counter:
            feature_cluster_prior_count = int(num_subjects / 2.0) if f not in feature_prior_counts else \
            feature_prior_counts[f]
            p = (feature_cluster_prior_count + 1) / float((total_num_subjects + 2) * feature_sizes[f])
            feature_cluster_alpha, feature_cluster_beta = beta_parameters(num_features, p)
            likelihood += cluster_likel(num_subjects * feature_sizes[f], 0, feature_cluster_alpha,
                                        feature_cluster_beta)
            n0 -= feature_sizes[f]
        else:
            for f in feature_prior_counts:
                if f not in feature_counter and f not in feature_sizes:
                    p = (feature_prior_counts[f] + 1) / float((total_num_subjects + 2))
                    feature_cluster_alpha, feature_cluster_beta = beta_parameters(num_features, p)
                    likelihood += cluster_likel(num_subjects * feature_prior_counts[f], 0, feature_cluster_alpha,
                                                feature_cluster_beta)
                    n0 -  1
    if n0 > 0:
        feature_cluster_prior_count = int(num_subjects * feature_cluster_size / 2.0)
        p = (feature_cluster_prior_count + 1) / float((total_num_subjects + 2) * feature_cluster_size)
        feature_cluster_alpha, feature_cluster_beta = beta_parameters(num_features, p)
        likelihood += n0 * cluster_likel(num_subjects, 0, feature_cluster_alpha, feature_cluster_beta)
    return likelihood


def cluster(cluster1, cluster2):
    return cluster_likelihood(len(cluster1.subjects) + len(cluster2.subjects),
                                  cluster1.feature_counts + cluster2.feature_counts
                              )

def beta_parameters(n, mu):
    s = ((mu * (1 - mu)) / float(n)) ** 0.5
    a = ((1 - mu) / float(s) - 1 / float(mu)) * (mu ** 2)
    b = a * (1 / float(mu) - 1)
    return a, b

########################################################################################
###############################   CLASS CLUSTER     ####################################
########################################################################################

class Cluster:
    count = 0
    lik = 0
    def __init__(self, subject, features=None):
        features = features or []
        try:
            self.subjects = [int(subject)]
        except:
            self.subjects = [subject]
        if all([":" in f for f in features]):
            self.feature_counts = Counter()
            for f in features:
                key_val = f.split(":")
                self.feature_counts[key_val[0]] += int(key_val[1])
        else:
            self.feature_counts = Counter(features)
        self.cluster_similarity = {}
        Cluster.count += 1
    def close(self):
        del self.subjects
        del self.feature_counts
        del self.lik
    def set_lik(self):
        self.lik = cluster_likelihood(len(self.subjects), self.feature_counts)
    def merge(self, other_cluster, sim=""):
        self.subjects += other_cluster.subjects
        self.feature_counts += other_cluster.feature_counts
        self.lik = sim + self.lik + other_cluster.lik
        other_cluster.close()
        del other_cluster
        Cluster.count = Cluster.count - 1

####################################################################### end of class

def initialise(Data_file, Feature_file=None, Clusters_size=None, informative_prior = True):
    global feature_sizes, feature_prior_counts, num_features, total_num_subjects
    feature_separator = ","
    cluster_list = []
    feature_totals = Counter()
    feature_sizes = {}
    feature_prior_counts = {}
    if Feature_file:
        F = open(Feature_file, "r")
        for line in F:
            d = line.strip().split()
            feature_sizes[d[0]] = int(d[1])
    total_num_subjects = 0
    D = open(Data_file, 'r')
    for line in D:
        key_val = line.strip().split("\t") #.replace("U", "").replace("C", "")
        try:
            cluster_list += [Cluster(key_val[0], key_val[1].split(feature_separator))]
            try:
                total_num_subjects = max(total_num_subjects, key_val[0])
            except:
                None
        except:
            cluster_list += [Cluster(Cluster.count, key_val[0].split(feature_separator))]
        for feature in key_val[-1].split(feature_separator):
            feature_totals[feature] += 1
    total_num_subjects = max(total_num_subjects, len(cluster_list))
    if Clusters_size:
        F = open(Clusters_size, "r")
        for line in F:
            d = line.strip().split()
            feature_prior_counts[d[0]] = int(d[1])
    else:
        if informative_prior:
            feature_prior_counts = feature_totals
    num_features = len(feature_totals)
    Cluster.likel = 0
    for cl in cluster_list:
        cl.set_lik()
        Cluster.lik += cl.lik
    for i in range(Cluster.count - 1):
        for j in range(i + 1, Cluster.count):
            sim = cluster_similarity(cluster_list[i], cluster_list[j])
            cluster_list[i].cluster_similarity[cluster_list[j]] = sim
            cluster_list[j].cluster_similarity[cluster_list[i]] = sim
    return cluster_list

def cluster_similarity(cluster1, cluster2):
    sim = cluster(cluster1, cluster2) - cluster1.lik - cluster2.lik
    return sim

def find_most_similar():
    max_similarity = -sys.float_info.max
    for i in range(Cluster.count - 1):
        for j in range(i + 1, Cluster.count):
            similarity = cluster_list[i].cluster_similarity[cluster_list[j]]
            if similarity > max_similarity:
                best_pair = [i, j]
                max_similarity = similarity
    return best_pair, max_similarity

def merge_clusters(i, j, sim=""):
    if i > j:
        temp = i
        i = j
        j = temp
    for cl in cluster_list:
        if cl != cluster_list[i] and cl != cluster_list[j]:
            cl.cluster_similarity.pop(cluster_list[j])
    cluster_list[i].merge(cluster_list[j], sim)
    cluster_list.pop(j)
    if i > j:
        i = i - 1
    for cl in cluster_list:
        if cl != cluster_list[i]:
            sim = cluster_similarity(cluster_list[i], cl)
            cluster_list[i].cluster_similarity[cl] = sim
            cl.cluster_similarity[cluster_list[i]] = sim

def get_clustering_from_cluster_list():
    clustering = {}
    counter = 0
    for cl in cluster_list:
        counter = counter + 1
        for sub in cl.subjects:
            clustering[sub] = counter
    return clustering

def user_clustering():
    max_l = Cluster.lik
    while Cluster.count > 1:
        pair, delta = find_most_similar()
        if delta > 0:
            merge_clusters(pair[0], pair[1], delta)
            Cluster.lik += delta
            if Cluster.lik > max_l:
                max_l = Cluster.lik
        else:
            break
    return get_clustering_from_cluster_list()


#####################################################################################
###################################### MAIN #########################################

def get_sorted(df):
        # first, sort df by node1, node2:
        with open(df, 'r') as DF:
            rows = DF.readlines()
            sorted_rows = sorted(rows, key=lambda x: (x.split(",")[0],x.split(",")[1]), reverse=False)
            filename = "sorted_edgelist.txt"
            with open(filename, 'w') as sorted_df:
                for row in sorted_rows:
                    sorted_df.write(row)
            return filename


def get_data(df, sep=","):

    """
    :param df: a txt data set containing an edge list as a pair of connected nodes. Example row: node1, node2
    :param sep: specifies the delimiter used to separate the connected nodes in each edge. Default comma separator (',').
    :return: sorted list with keys unique nodes and values list of nodes connected to each key.
    """

    for d in ('connected_nodes.txt', 'sorted_edgelist.txt'):
        if os.path.exists(d): os.remove(d)
    filename = 'connected_nodes.txt'
    old_node1 = ""
    sorted_df = get_sorted(df)
    with open(sorted_df, 'r') as DF:
        for line in DF:
            [node1, node2] = line.strip().split(sep)
            if node1 != old_node1:
                if old_node1 != "":
                    with open(filename, 'a') as f:
                        print(old_node1 + "\t" + ",".join(str(c) for c in sorted(connctd_nodes)), file=f)
                connctd_nodes = [(node2)]
                old_node1 = node1
            else:
                connctd_nodes += [(node2)]
        if old_node1 != "":
            with open(filename, 'a') as f:
                print(old_node1 + "\t" + ",".join(str(c) for c in sorted(connctd_nodes)), file=f)
    f.close()
    return filename


def run_clustering(df_filename, edgelist=True, sequential=False):
    """
    :param Reading_file:  data file containing edge list
    :param edgelist: if TRUE (default) reads from edgelist, else reads directly from sorted unique list of nodes with appended connections
    :param sequential: if TRUE sequential version is performed, otherwise batch-clustering
    :return: tuples of each item (network node) with associated cluster label
    """

    if edgelist:
        cnnctd_nodes = get_data(df_filename)
    else:
        cnnctd_nodes = df_filename

    if not sequential:
        global cluster_list
        if not cnnctd_nodes:
            raise ValueError("sorted edge-list not found: please do not remove 'sorted_edgelist.txt' or 'connected_nodes.txt' ")
        cluster_list = initialise(cnnctd_nodes, Feature_file=None, Clusters_size=None)
        clustering = user_clustering()
        print(sorted(clustering.items(), key=operator.itemgetter(0)))


