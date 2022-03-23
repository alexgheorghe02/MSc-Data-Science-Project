#This node2vec implementation uses the reference code provided by Grover and Leskovec.
#See https://github.com/aditya-grover/node2vec for the reference implementation.
#The reference code was adapted to work with string-labelled graphs.
#This node2vec implementation also includes code for hub attention-directed random walks.
#For the reference hub attention implementation,
#see https://github.com/Kombustor/submission-ecir2020-randomwalks/tree/submission
#Note that the original hub attention implementation by Schliski et al was modified to be compatible with
#the alias sampling implementation used by Grover and Leskovec.
#See node2vec: Scalable Feature Learning for Networks. Aditya Grover and Jure Leskovec. Knowledge Discovery and Data Mining, 2016.
#See Influence of Random Walk Parametrization on Graph Embeddings. Fabian Schliski, Jörg Schlötterer and Michael Granitzer. ECIR 2020.


import networkx as nx
import numpy as np
import random
from tqdm import tqdm
from gensim.models import Word2Vec
import os
os.environ['PYTHONHASHSEED']='123'

def read_graph(input, weighted = False):
    if weighted:
        G = nx.read_edgelist(input, nodetype=str, data=(('weight',float),), create_using=nx.Graph())
    else:
        G = nx.read_edgelist(input, nodetype=str, create_using=nx.Graph())
        for edge in G.edges():
            G[edge[0]][edge[1]]['weight'] = 1
    return G

def alias_setup(probs):
	'''
	Compute utility lists for non-uniform sampling from discrete distributions.
	Refer to https://hips.seas.harvard.edu/blog/2013/03/03/the-alias-method-efficient-sampling-with-many-discrete-outcomes/
	for details
	'''
	K = len(probs)
	q = np.zeros(K)
	J = np.zeros(K, dtype=np.int)

	smaller = []
	larger = []
	for kk, prob in enumerate(probs):
	    q[kk] = K*prob
	    if q[kk] < 1.0:
	        smaller.append(kk)
	    else:
	        larger.append(kk)

	while len(smaller) > 0 and len(larger) > 0:
	    small = smaller.pop()
	    large = larger.pop()

	    J[small] = large
	    q[large] = q[large] + q[small] - 1.0
	    if q[large] < 1.0:
	        smaller.append(large)
	    else:
	        larger.append(large)

	return J, q

def alias_draw(J, q):
	'''
	Draw sample from a non-uniform discrete distribution using alias sampling.
	'''
	K = len(J)

	kk = int(np.floor(np.random.rand()*K))
	if np.random.rand() < q[kk]:
	    return kk
	else:
	    return J[kk]


class Node2Vec():
    def __init__(self, graph, p, q):
        self.G = graph
        self.p = p
        self.q = q

    def node2vec_walk(self, walk_length, start_node):
        '''
		Simulate a random walk starting from start node.
		'''
        G = self.G
        alias_nodes = self.alias_nodes
        alias_edges = self.alias_edges

        walk = [start_node]

        while len(walk) < walk_length:
            cur = walk[-1]
            cur_nbrs = sorted(G.neighbors(cur))
            if len(cur_nbrs) > 0:
                if len(walk) == 1:
                    walk.append(cur_nbrs[alias_draw(alias_nodes[cur][0], alias_nodes[cur][1])])
                else:
                    prev = walk[-2]
                    next = cur_nbrs[alias_draw(alias_edges[(prev, cur)][0], alias_edges[(prev, cur)][1])]
                    walk.append(next)
            else:
                break

        return walk

    def simulate_walks(self, num_walks, walk_length):
        '''
        Repeatedly simulate random walks from each node.
        '''
        G = self.G
        walks = []
        nodes = list(G.nodes())
        for walk_iter in range(num_walks):
            random.shuffle(nodes)
            for node in tqdm(nodes, desc=f"Walk {walk_iter + 1}/{num_walks}"):
                walks.append(self.node2vec_walk(walk_length=walk_length, start_node=node))

        return walks

    def get_alias_edge(self, src, dst):
        '''
        Get the alias edge setup lists for a given edge.
        '''
        G = self.G
        p = self.p
        q = self.q

        unnormalized_probs = []
        for dst_nbr in sorted(G.neighbors(dst)):
            if dst_nbr == src:
                unnormalized_probs.append(G[dst][dst_nbr]['weight'] / p)
            elif G.has_edge(dst_nbr, src):
                unnormalized_probs.append(G[dst][dst_nbr]['weight'])
            else:
                unnormalized_probs.append(G[dst][dst_nbr]['weight'] / q)
        norm_const = sum(unnormalized_probs)
        normalized_probs = [float(u_prob) / norm_const for u_prob in unnormalized_probs]

        return alias_setup(normalized_probs)

    def preprocess_transition_probs(self):
        '''
        Preprocessing of transition probabilities for guiding the random walks.
        '''
        G = self.G

        alias_nodes = {}
        for node in tqdm(G.nodes()):
            unnormalized_probs = [G[node][nbr]['weight'] for nbr in sorted(G.neighbors(node))]
            norm_const = sum(unnormalized_probs)
            normalized_probs = [float(u_prob) / norm_const for u_prob in unnormalized_probs]
            alias_nodes[node] = alias_setup(normalized_probs)

        alias_edges = {}

        for edge in G.edges():
            alias_edges[edge] = self.get_alias_edge(edge[0], edge[1])
            alias_edges[(edge[1], edge[0])] = self.get_alias_edge(edge[1], edge[0])

        self.alias_nodes = alias_nodes
        self.alias_edges = alias_edges

        return


class HubsWalker(Node2Vec):

    def __init__(self, G, h):
        self.G = G
        self.h = h
        self.detect_hubs()

    def detect_hubs(self):
        self.hubs = []

        # [(node1, node2...), (degree1, degree2, ...)]
        node_degrees = list(zip(*self.G.degree()))
        degrees = node_degrees[1]
        nodes = node_degrees[0]

        avg_degree = np.mean(degrees)
        std_degree = np.std(degrees)

        for idx, node in enumerate(nodes):
            deg = degrees[idx]
            if(deg > (avg_degree + std_degree)):
                self.hubs.append(node)


    def is_hub(self, node):
        return node in self.hubs

    def get_alias_edge(self, src, dst):
        G = self.G
        h = self.h

        unnormalized_probs = []
        # for every (sorted) neighbor of destination node
        for dst_nbr in sorted(G.neighbors(dst)):
            # if the destination neighbor is a hub, modify probability with h (=> likelihood of revisiting a hub increased)
            if self.is_hub(dst_nbr):
                unnormalized_probs.append(G[dst][dst_nbr]['weight']/h)
            else:
                unnormalized_probs.append(G[dst][dst_nbr]['weight'])

        norm_const = sum(unnormalized_probs)
        normalized_probs = [float(u_prob) / norm_const for u_prob in unnormalized_probs]

        return alias_setup(normalized_probs)


def learn_embeddings(walks, dimensions, window, epochs, output):
    '''
    Learn embeddings by optimizing the Skipgram objective using SGD.
    '''
    model = Word2Vec(walks, vector_size=dimensions, window=window, min_count=0, sg=1, epochs=epochs, workers = 1, seed = 1)
    model.wv.save_word2vec_format(output)
    return model
