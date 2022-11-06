"""
deBruijnGraph.py

Description: An implementation of a simplified version
of sequence assembly using de Bruijn graphs. Two modes 
can be used by the user, mode 'g' produces a filtered 
good_reads file that contains all reads that are consisted 
of k-mers which appear more than once; and mode 'c' produces
assembled non-ambiguious contigs 
    
Created by Etha Hua, March 14 2022
"""

import sys


# Class definition of a de-Bruijn-graph node
class graphNode:
    def __init__(self, value):
        self.val = value
        self.outgo = []
        self.ingo = []
        self.visited = False
        self.branch = False
        self.count = 1
    
    # Method for adding an outgoing node of a node
    def out_edge(self, value):
        self.outgo.append(value)

    # Method for adding an ingoing node of a node
    def in_edge(self, value):
        self.ingo.append(value)

# Class definition for a de Bruijn graph
class deBruijnGraph:
    def __init__(self):
        self.nodes = dict()

    # Method for building the graph
    # Parameters: fstream: the input file stream that contains
    #                      reads to be incorporated into the graph
    #             kval: k value, the length of a k-mer in the graph
    # Returns: N/A
    # Note: it updates the graph that has been called upon
    def buildGraph(self, fstream, kval):
        for line in fstream:
            numLoop = len(line) - kval
            for i in range(0, numLoop):
                kmer = line[i:i+kval]
                self.addNode(kmer)
                if i > 0 and i < numLoop - 1:
                    self.addInEdge(kmer, line[i - 1:i + kval - 1])
                    self.addOutEdge(kmer, line[i + 1:i + kval + 1])
                elif i == 0:
                    self.addOutEdge(kmer, line[i + 1:i + kval + 1])
                else: # i == numLoop - 1
                    self.addInEdge(kmer, line[i - 1:i + kval - 1])

    # Method for adding a node to the graph
    # Parameter: valNode: the value of the node to be added
    # Returns: N/A
    # Note: if that node value already exists in the graph,
    #       the node will be incremented its count by 1
    def addNode(self, valNode):
        if self.nodes.get(valNode) is not None:
            self.nodes[valNode].count += 1
        else:
            self.nodes[valNode] = graphNode(valNode)
    
    # Method for adding an ingoing edge to an existing node
    # Parameter: valNode: the value of the node for which an in 
    #            going edge is to be added;
    #            inValNode: the value of the node which the 
    #            ingoing edge is leading to
    # Note: it also keeps track of the branchiness of the node
    def addInEdge(self, valNode, inValNode):
        if inValNode not in self.nodes[valNode].ingo:
            self.nodes[valNode].ingo.append(inValNode)
        if len(self.nodes[valNode].ingo) > 1:
            self.nodes[valNode].branch = True
    
    # Method for adding an outgoing edge to an existing node
    # Parameter: valNode: the value of the node for which an out 
    #            going edge is to be added;
    #            outValNode: the value of the node which the 
    #            outgoing edge is leading to
    # Note: it also keeps track of the branchiness of the node
    def addOutEdge(self, valNode, outValNode):
        if outValNode not in self.nodes[valNode].outgo:
            self.nodes[valNode].outgo.append(outValNode)
        if len(self.nodes[valNode].outgo) > 1:
            self.nodes[valNode].branch = True
    
    # Method for removing all branching nodes in the graph,
    #        including their adjacent edges.
    def removeBranchingNodes(self):
        toremove = []
        for node in self.nodes:
            currNode = self.nodes[node]
            if currNode.branch:
                for inNode in currNode.ingo:
                    self.nodes[inNode].outgo.remove(currNode.val)
                for outNode in currNode.outgo:
                    self.nodes[outNode].ingo.remove(currNode.val)
                toremove.append(currNode.val)
        for bNode in toremove:
            self.nodes.pop(bNode)

    # Method for reading one contig off the de Bruijn graph, and mark
    # all the visited nodes as visited
    # Parameter: the node value of a node to start the path (contig)
    # Returns: a string representing the current contig read off the graph
    def readOneContig(self, start):
        currContig = start
        currNode = self.nodes[start]
        currNode.visited = True
        
        while(len(currNode.outgo) == 1 and not self.nodes[currNode.outgo[0]].visited):
            currNode = self.nodes[currNode.outgo[0]]
            currNode.visited = True
            currContig += currNode.val[-1]

        return currContig


    # Method for finding all the contigs in the de Bruijn Graph
    # Parameter: minLen: the minimum length of a contig that makes it 
    #                    an acceptable contig to be reported
    # Returns: an array of all contigs (strings)
    def findContigs(self, minLen):
        contigs = []
        for node in self.nodes:
            currNode = self.nodes[node]
            if len(currNode.ingo) == 0:
                currContig = self.readOneContig(currNode.val)
                if len(currContig) >= minLen:
                    contigs.append(currContig)
        for node in self.nodes:
            currNode = self.nodes[node]
            if not currNode.visited:
                currContig = self.readOneContig(currNode.val)
                if len(currContig) >= minLen:
                    contigs.append(currContig)
        return contigs

    # Method for printing all good reads, which is defined as 
    # containing all the k-mers that is not unique among all 
    # other reads
    # Parameter: reads: a fstream that contains all the original
    #                   reads sequences
    #            kval: the value of k, the length of a k-mer
    # Note: this function prints all the good reads to a file
    #       called "good_reads"
    def printGoodReads(self, reads, kval):
        original_stdout = sys.stdout 

        with open('good_reads', 'w') as f:
            sys.stdout = f 
            for read in reads:
                numLoop = len(read) - kval
                toPrint = True
                for i in range(0, numLoop):
                    kmer = read[i:i+kval]
                    if self.nodes[kmer].count <= 1:
                        toPrint = False
                        break
                if toPrint:
                    print(read[:-1])
        sys.stdout = original_stdout # restore the original stdout



# outputContigs
# Purpose: printing out contigs and lengths of contigs to 
#          "output_contigs" and "contig_lengths" respectively
# Parameter: contigs: an array of the contigs to be printed
def outputContigs(contigs):
    
    original_stdout = sys.stdout 

    with open('output_contigs', 'w') as f:
        sys.stdout = f 
        for contig in contigs:
            print(contig)

    with open('contig_lengths', 'w') as f:
        sys.stdout = f 
        for contig in contigs:
            print(len(contig))       
        
    sys.stdout = original_stdout # restore the original stdout
    

# Main function starts here
if __name__ == "__main__":

    if len(sys.argv) < 4:
        print("Too few arguments.")
        exit()
    elif len(sys.argv) > 5:
        print("Too many arguments.")
        exit()
    elif sys.argv[2] != "g" and sys.argv[2] != "c":
        print("Mode not supported.")
        exit()

    # set the length of a k-mer
    kval = int(sys.argv[3])

    fname = sys.argv[1]
    fstream = open(fname, "r")
    myGraph = deBruijnGraph()
    myGraph.buildGraph(fstream, kval)


    if sys.argv[2] == "g":
        fstream = open(fname, "r")
        myGraph.printGoodReads(fstream, kval)
    else: # sys.argv[2] == "c"
        # set the minimum contig length to filter out short contigs
        mincLength = 100
        if len(sys.argv) == 5:
            mincLength = int(sys.argv[4])
        myGraph.removeBranchingNodes()
        contigs = myGraph.findContigs(mincLength)
        outputContigs(contigs)