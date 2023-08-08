import math
import re
import sys 
import networkx as nx
import random
import numpy as np

class GraphTolls:
    
    def __init__( self, Path) -> None:
        self.Path = Path
        self.graph = self.Read_Graph()
        self.m = self.graph.number_of_edges()
        self.n = self.graph.number_of_nodes()
        self.adjency = self.graph.adj
        self.Node_list = {i : i for i in self.graph.nodes()}
        self.Degree = dict(self.graph.degree())
        self.DegCom = {}
        self.membership = { i : None for i in self.graph.nodes() }
        
    def Read_Graph( self):
        
        if self.Path[len(self.Path)-3: ] == 'txt' or self.Path[len(self.Path)-3: ] == 'dat':
            Graph = nx.read_edgelist(self.Path, nodetype = int, data = True)
            graph = nx.Graph(Graph)
        elif self.Path[len(self.Path)-3: ] == 'gml':
            Graph = nx.read_gml(self.Path,label = 'id')
            graph = nx.Graph(Graph)
        else :
            raise TypeError (" the type of graph is not suportable or not no such file or directory")

        return graph 
    
    
    def Read_GroundTruth( self, path):
        with open(path, "r") as file:
            lines = file.readlines()
            result = []
            for x in lines:
                x = x.rstrip()
                result.append(x.split()[1])

        true_partion = [int(x)for x in result]
        return true_partion
    
     
    def Is_Intersiction( self, communities):
        dupes = []
        flat = [item for sublist in communities for item in sublist]
        for f in flat:
            if flat.count(f) > 1:
                if f not in dupes:
                    dupes.append(f)

        if dupes:
            return True
        else:
            return False   
        
    def sum( self, arg):
        if len(arg) < 1:
            return None
        else:
            return sum(arg)

    def count( self, arg):
        return len(arg)
  
    def min( self, arg):
        if len(arg) < 1:
            return None
        else:
            return min(arg)
  
    def max( self, arg):
        if len(arg) < 1:
            return None
        else:
            return max(arg)
  
    def avg( self,arg):
        if len(arg) < 1:
            return None
        else:
            return sum(arg) / len(arg)   
  
    def median( self,arg):
        if len(arg) < 1:
            return None
        else:
            arg.sort()
            return  arg[len(arg) // 2]
  
    def stdev( self,arg):
        if len(arg) < 1 or len(arg) == 1:
            return None
        else:
            avg = self.avg(arg)
            sdsq = sum([(i - avg) ** 2 for i in arg])
            stdev = (sdsq / (len(arg) - 1)) ** .5
            return stdev
  
    def percentile( self, arg):
        if len(arg) < 1:
            value = None
        elif (arg >= 100):
            sys.stderr.write('ERROR: percentile must be < 100.  you supplied: %s\n'% arg)
            value = None
        else:
            element_idx = int(len(arg) * (arg / 100.0))
            self.arg.sort()
            value = self.arg[element_idx]
        return value  
    
    
    def keys_of_maximum_value( self, d):
        max_value = max( d.values())
        for key, value in d.items():
            if value == max_value:
                return key
        
        
    def ngh_node ( self, node, com):
        ngh_com = 0
        for ngh in self.graph[node]:  
            if  self.membership[ngh] == com:
                if ngh == node:
                    ngh_com += 1
                else:
                    ngh_com += 1 / 2.
                               
        return ngh_com
    
    def neigh_comm ( self, node):
        ngh_com = {}
        link = 1
        for ngh in self.adjency[node]:  
            if  self.membership[ngh] != None:
                com_id = self.membership[ngh]
                ngh_com[com_id] = ngh_com.get( com_id, 0) + link
                                             
        return ngh_com
    
    def com_ngh_com ( self, com_id):
        com_ngh = {}
        for node, com in self.membership.items():
            if  com == com_id:  
                link = 1
                for ngh in self.adjency[node]:  
                   comid = self.membership[ngh]
                   if com != comid :
                       com_ngh[comid] = com_ngh.get( comid, 0) + link                       
    
        return com_ngh
    
    def merge_com ( self, com_id1, com_id2 ):
        for node, com in self.membership.items():
            if com == com_id2:
                self.membership[node] = com_id1
                
    def sel_edge_btwc( self, com_id1, com_id2):
        Edg_betw = 0
        for node, com in self.membership.items():
            if com == com_id1:
                Edg_betw = Edg_betw + self.ngh_node( node, com_id2)
            
        return Edg_betw 
   
    def weighted_choice( self, objects, weights):
        weights = np.array( weights, dtype = np.float64)
        sum_of_weights = weights.sum()
        np.multiply( weights, 1 / sum_of_weights, weights)
        weights = weights.cumsum()
        x = random.random()
        for i in range(len(weights)):
            if x < weights[i]:
                return objects[i]
               
    def generate_random_not_in_list( self, my_list,min_val,max_val):
        while True:
            random_number = random.randint( min_val, max_val)
            if random_number not in my_list:
                return random_number            
    
    def interlink ( self, internals ):
        for node in self.graph.nodes():
            com = self.membership[node]
            inc = 0.        
            for neighbor in self.graph[node]:
                if self.membership[neighbor] == com:
                    if neighbor == node:
                        inc += 1
                    else:
                        inc += 1 / 2.
            internals[com] = internals.get( com, 0) + inc
                     
    def modularity( self):
        result = 0
        internals = {}
        self.interlink(internals)
        for com in set( self.membership.values()):
            in_degree = internals.get( com, 0.)
            degree = self.DegCom.get( com, 0.)
            result += in_degree / self.m - (( degree / ( 2. * self.m )) ** 2)
            
        return result
        