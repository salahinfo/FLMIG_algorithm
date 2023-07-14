import math
import re
import sys 
import networkx as nx
import random
import numpy as np

class GraphTolls:
    
    def __init__(self, Path) -> None:
        self.Path = Path
        self.graph = self.Read_Graph()
        self.m = self.graph.number_of_edges()
        self.n = self.graph.number_of_nodes()
        self.adjency = self.graph.adj
        self.Node_list = [i for i in self.graph.nodes()]
        self.Degree = self.graph.degree()
        
    def Read_Graph(self):
        
        if self.Path[len(self.Path)-3: ] == 'txt' or self.Path[len(self.Path)-3: ] == 'dat':
            Graph = nx.read_edgelist(self.Path, nodetype = int, data = True)
            graph = nx.Graph(Graph)
        elif self.Path[len(self.Path)-3: ] == 'gml':
            Graph = nx.read_gml(self.Path,label = 'id')
            graph = nx.Graph(Graph)
        else :
            raise TypeError (" the type of graph is not suportable or not no such file or directory")

        return graph 
    
    def Reve(self,x):
        sa = x.split()[::-1]
        l = []
        for i in sa:
            l.append(i)

        l=('  '.join(l))
        return l
    
       
    def Remove_Revers(self,path):
        with open(path, "r") as file:
            lines = file.readlines()
            result=[]
            for xa in lines:
                xa=re.sub(r'\s','  ',xa)
                if self.reve(xa) not in result:
                   result.append(xa.strip())   

        return(result)

    def Remove_Dublicate (self,path):
        pathw = path[ : -3]+'txt' 
        lines = self.remove_revers(path)
        with open(pathw, 'w') as f:
            for line in lines:
               f.write(line)
               f.write('\n')


    def Read_GroundTruth(self,path):
        with open(path, "r") as file:
            lines = file.readlines()
            result = []
            for x in lines:
                x = x.rstrip()
                result.append(x.split()[1])

        true_partion = [int(x)for x in result]
        return true_partion
    
     
    def Is_Intersiction(self,communities):
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
        
    def sum(self,arg):
        if len(arg) < 1:
            return None
        else:
            return sum(arg)

    def count(self,arg):
        return len(arg)
  
    def min(self,arg):
        if len(arg) < 1:
            return None
        else:
            return min(arg)
  
    def max(self,arg):
        if len(arg) < 1:
            return None
        else:
            return max(arg)
  
    def avg(self,arg):
        if len(arg) < 1:
            return None
        else:
            return sum(arg) / len(arg)   
  
    def median(self,arg):
        if len(arg) < 1:
            return None
        else:
            arg.sort()
            return  arg[len(arg) // 2]
  
    def stdev(self,arg):
        if len(arg) < 1 or len(arg) == 1:
            return None
        else:
            avg = self.avg(arg)
            sdsq = sum([(i - avg) ** 2 for i in arg])
            stdev = (sdsq / (len(arg) - 1)) ** .5
            return stdev
  
    def percentile(self, arg):
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
    
    def select_edge_betw(self,* arg):
        Edg_betw = 0
        for node in self.adjency[arg[0]]:
            if node in arg[1]:
                Edg_betw = Edg_betw + 1
        
        return Edg_betw 
    
    def is_edge_betw(self,vert,commu):
        for node in self.adjency[vert]:
            if node in commu:
                return True
            
        return False
    
    def select_edge_c(self,* arg):
        Edg_betw = 0
        for v in arg[0]:
            for node in self.adjency[v]:
                if node in arg[1]:
                    Edg_betw = Edg_betw + 1
        
        return Edg_betw 
   
    def weighted_choice(self,objects, weights):
        weights = np.array(weights, dtype = np.float64)
        sum_of_weights = weights.sum()
        np.multiply(weights, 1 / sum_of_weights, weights)
        weights = weights.cumsum()
        x = random.random()
        for i in range(len(weights)):
            if x < weights[i]:
                return objects[i]
