#This program verifies that the
#conditional fractional strong matching preclusion number of B3 is 4,
#and that B3 is conditionally fractional strongly super matched.
#This program runs in around a day, so file output is not necessary.

import networkx as nx
import itertools as itr
from scipy.optimize import linprog
import warnings
warnings.simplefilter('ignore')

#Return whether or not graph G is basic
#(contains an isolated vertex)
def is_basic(G):

    degrees = []
    for v in G.nodes():
        degrees.append(G.degree(v))
    
    if min(degrees) == 0:
        return True
    else:
        return False

#Return whether or not graph G is conditional basic
#(contains a set of three vertices where one is the only neighbor of the other two)
def is_cond_basic(G):
    degrees = []
    for v in G.nodes():
        degrees.append(G.degree(v))

    if is_basic(G):
        return False

    if degrees.count(1) == 2:
        nodes2 = []
        for v in G.nodes():
            if G.degree(v) == 1:
                nodes2.append(v)

        total = 0
        for i in nx.common_neighbors(G, nodes2[0], nodes2[1]):
            total+=1
        if total >= 1:
            return True
        else:
            return False
    else:
        return False

#Return whether or not graph P has a fractional perfect matching
#Done by creating a linear program which imposes the following bounds:
#1. The sum of the weights on the edges surrounding every vertex must be 1
#2. The weight of each edge is between 0 and 1
#Then returning whether or not the linear program found a solution (res.success)
def doesGraphHaveFPM(G):
    EdgeToNum = dict()
    iter = 0
    for e in G.edges():
        EdgeToNum[e] = iter
        EdgeToNum[e[::-1]] = iter
        iter+=1

    NodeToNum = dict()
    iter = 0
    for n in G.nodes():
        NodeToNum[n] = iter
        iter+=1

    c = [1 for e in G.edges()]

    A = [[0 for e in G.edges()] for n in G.nodes()]

    B = [1 for n in G.nodes()]

    for n in G.nodes():
        for e in G.edges(n):
            A[NodeToNum[n]][EdgeToNum[e]] = 1
    Bound = (0,1)

    res = linprog(c = c, A_eq = A, b_eq = B, bounds = Bound)
    return res.success


#Generate and return a burnt pancake graph of any dimension
def burntpancakegraph(n):
    A = nx.Graph()

    #All permutations of the first n natural numbers for which
    #each is either positive or negative
    nums = []
    for i in range(1, n+1):
        nums.append(i)
    naturalLabels = itr.permutations(nums, n) 
    binaryLabels = itr.product([0,1], repeat = n)

    naturalArray = list(naturalLabels)
    binaryArray = list(binaryLabels)

    for i in naturalArray:
        for j in binaryArray:
            newLabels = []
            for k in range(n):
                newLabels.append(list(i)[k])
                if (j[k] == 1):
                    newLabels[k]*=-1
            newLabelTuple = tuple(newLabels)
            A.add_node(newLabelTuple)

    
    #For the first k numbers in the compared vertices, the k string must be
    #reversed and negated, and the rest of the n-k numbers must be the same.
    C = A.copy()
    for u in A.nodes():
        B = A.copy()
        B.remove_node(u)
        for v in B.nodes():
            for k in range(1, n+1):
                
                flipTuple = u[::-1]
                compareTuple = u
                arrayRep = []
                for i in range(0, k):
                    arrayRep.append(flipTuple[n-k+i]*-1)
                for i in range(k, n):
                    arrayRep.append(compareTuple[i])
                newcomparetuple = tuple(arrayRep)
                if (newcomparetuple == v):
                    C.add_edge(u,v)
    
    return C


#Generate a B3 graph. This is the only time we call burntpancakegraph() because
#it is faster to copy this graph than generate a new one.
R = burntpancakegraph(3)

#To optimize using the vertex transitive property of B3,
#we remove the same vertex (1,2,3) every time the fault set
#has at least one vertex.
#To do this, we must iterate through all the vertices which
#aren't (1,2,3).
non123 = list(R.nodes())
non123.remove((1,2,3))

#Keep track of progress
#There are 1302609 cases in total
numTested = 0

#We remove i edges and 4-i vertices
#for a total of 4 faults
#(one vertex is removed later)
for i in range(0, 5):
    for k in itr.combinations(non123, max(0,3-i)):
        for m in itr.combinations(list(R.edges()), i):

            #nonTrivPreclusion changes to false as soon as we know that the current
            #fault set being tested can not be a non-trivial preclusion.
            nonTrivPreclusion = True

            #The 0 can be changed to however many cases have been tested so far
            #to pause and restart the program. A similar strategy can be used
            #to distribute the algorithm.
            if (numTested < 0):
                nonTrivPreclusion = False

            #Create a fault-ridden B3 and check if the fault set is trivial
            if (nonTrivPreclusion == True):

                P = R.copy()
                P.remove_edges_from(m)
                P.remove_nodes_from(k)

                if (i != 4):
                    P.remove_node((1,2,3))

                if (is_basic(P) or is_cond_basic(P)):
                    nonTrivPreclusion = False

            #If the fault set isn't trivial, check if it's a preclusion set
            if (nonTrivPreclusion == True):
                if (doesGraphHaveFPM(P)):
                    nonTrivPreclusion = False

            #If the fault set is a preclusion set, B3 fails the test.
            if (nonTrivPreclusion == True):
                #This code never runs because, as we find out,
                #B3 is conditionally fractional strongly super matched.
                print("Found a CFSMP set which does not trivially\
                      almost isolate two vertices.")
                print("vertices: " + str(k) + "  edges: " + str(m))

            #Print progress
            numTested+=1
            if (numTested%1000 == 0):
                print("Tested " + str(numTested) + " combos so far... (" +\
                str(round(100*numTested/1302609, 2)) + "%)")
                
print("Completed")
