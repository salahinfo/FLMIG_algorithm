import math
import re
import sys 
import networkx as nx
import random
import numpy as np
from networkx.algorithms.community import is_partition

class GraphTolls:
    
    def __init__( self, Path , graph) -> None:
        self.Path = Path
        self.m = graph.size( weight= 'weight')
        self.n =  graph.number_of_nodes()
        #self.adjency = graph.adj
        self.Node_list = {i : i for i in graph.nodes()}
        self.Degree = dict(graph.degree())
        self.DegCom = {}
        #self.membership = { i : None for i in graph.nodes() }
        self.loops = {}
        self.internal = {}
        
    
    
    def Read_GroundTruth( self, path):
        with open(path, "r") as file:
            lines = file.readlines()
            result = []
            for x in lines:
                x = x.rstrip()
                result.append(x.split()[1])

        true_partion = [ int(x) for x in result]
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
  
    def avg( self, arg):
        if len(arg) < 1:
            return None
        else:
            return sum(arg) / len(arg)   
    
    def median( self, arg):
        if len(arg) < 1:
            return None
        else:
            arg.sort()
            return  arg[len(arg) // 2]
  
    def stdev( self, arg):
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
        
        
    def ngh_node ( self, membership,  node, com, graph ):
        ngh_com = 0.
        for ngh, datas in graph[node].items():
            link = datas.get('weight', 1.)  
            if  membership[ngh] == com:
                ngh_com += link
                             
        return ngh_com
    
    def neigh_comm ( self, membership, node, graph):
        ngh_com = {}
        for ngh, datas in graph[node].items():
            #print("weigh", link)            
            if  node != ngh:
                link = float(datas.get( 'weight', 1.))
                com_id = membership[ngh]
                ngh_com[com_id] = ngh_com.get( com_id, 0.) + link
                                                      
        return ngh_com
    
    def ngh_strcom ( self, membership, comid,  node, graph):
        ngh_com = {}
        for ngh, datas in graph[node].items():
            #print("weigh", link)            
            if  membership[ngh] != None  and  node != ngh and membership[ngh] != comid:
                link = float(datas.get( 'weight', 1.))
                com_id = membership[ngh]
                ngh_com[com_id] = ngh_com.get( com_id, 0.) + link
                                                      
        return ngh_com

    def boundry_node( self, membership, graph):
        b_node = {}
        for node in graph.nodes():
            com = membership[node]
            for ngh,data in graph[node].items():
                if membership[ngh]!= com:
                    b_node[node] = self.membership[ngh]
        
        return b_node
    
    def com_ngh_com ( self, membership, com_id, graph):
        com_ngh = {}
        for node in  graph.nodes(): 
            com = membership[node]
            if  com == com_id:  
                for ngh, datas in  graph[node].items():
                    link = float(datas.get('weight', 1.))
                    #print(link)
                    comid = membership[ngh]
                    if comid != com_id and node != ngh:
                        com_ngh[comid] = com_ngh.get( comid, 0.) + link                      
                    
        return com_ngh
    
    def merge_com ( self, membership, com_id1, com_id2 ):
        for node, com in membership.items():
            if com == com_id2:
                membership[node] = com_id1

        return membership        
                
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
               
    def generate_random_not_in_list( self, my_list):
        while True:
            random_number = random.randint( 0, self.m)
            if random_number not in my_list:
                return random_number            
    
                         
    def modularity( self, membership):
        result = 0
        for com in set( membership.values()):
            in_degree = self.internal.get( com, 0.)
            degree = self.DegCom.get( com, 0.)
            #print(self.m)
            result += in_degree / self.m - (( degree / ( 2. * self.m )) ** 2.)
            
        return result
    
    def induced_graph( self, p, graph, weight):
        ret = nx.Graph()
        ret.add_nodes_from(p.values())
        #print(p, graph.nodes())
        for node1, node2, datas in graph.edges( data = True ):
            edge_weight = datas.get(weight, 1)
            com1 = p[node1]
            com2 = p[node2]
            w_prec = ret.get_edge_data(com1, com2, {weight: 0}).get( weight, 1)
            ret.add_edge( com1, com2, **{weight: w_prec + edge_weight})

        return ret
    
    
    def check_connectivite( self, p, graph, weight= 'weight'):
        # FUNCTION TO CHECK THE CONNECTIVITE #
        com = set(p.values())
        check = []
        Ncom = {}
        for com_id in com:
            ret = nx.Graph()
            nodecom = [node for node in p if p[node] == com_id]

            ret.add_nodes_from(nodecom)
            ret = graph.subgraph(nodecom)
            if nx.is_connected(ret):
                check.append(True)
            else :
                Ncom[com_id] = list(nx.connected_components(ret))
                # RETURN DICONNECTED COMPONENTS #
            
        if len(check) == len(com):
            return True
        
        else :

            return Ncom
        
        
    def check_partition( self, p, graph, weight= 'weight'):
        # CHECK INTERSICTION OF PARTITION 
       
        partition = self.tarans( p)    
       
        if is_partition(graph, partition):
            
            return True        
        else : 
            return False    

    
            
       

    
    def modifie_status( self, graph, weight, solution = None ):
        """Initialize the status of a graph with every node in one community"""
        self.total_weight = 0
        self.Degree = {}
        self.DegCom = {}
        #self.adjency = graph.adj
        self.internal = {}
        membership = {}
        self.loops = {}
        self.m = graph.size(weight=weight)
        self.n = graph.number_of_nodes()
        if solution is None:
            for node in graph.nodes():
                membership[node] = node
                deg = float(graph.degree(node, weight='weight'))
                if deg < 0:
                    error = "Bad graph type ({})".format(type(graph))
                    raise ValueError(error)
                self.Degree[node] = deg
                self.DegCom[node] = deg
                edge_data = graph.get_edge_data( node, node, {'weight': 0})
                #print(edge_data)
                self.loops[node] = float(edge_data.get('weight', 1))
                self.internal[node] = self.loops[node]
       
        return membership
    
    def renumber( self , membership):
        """Renumber the values of the dictionary from 0 to n"""
        count = 0
        ret = membership.copy()
        new_values = dict([])
    
        for key in  membership.keys():
            value =  membership[key]
            new_value = new_values.get(value, -1)
            if new_value == -1:
                new_values[value] = count
                new_value = count
                count += 1
            membership[key] = new_value
        return  membership
        
    def delet_node( self, membership, node, com, weicom):
        
        self.DegCom[com] = float(self.DegCom.get(com, 0.) - self.Degree.get(node, 0.))
        self.internal[com] = float(self.internal.get(com, 0.) - weicom - self.loops.get(node, 0.))
        membership[node] = None

        return membership
        
    def insert_node( self, membership,  node, com, weicom):
            
        self.DegCom[com] = float(self.DegCom.get( com, 0.) + self.Degree.get( node, 0.))
        self.internal[com] = float(self.internal.get( com, 0.) + weicom + self.loops.get( node, 0.))
        membership[node] = com  

        return membership                    
    
    def init( self, graph, solution, weight = 'weight'):
        self.DegCom = {}
        self.internal = {}
        self.Degree= {}
        membership = {}
        self.m = graph.size( weight= 'weight')
        self.n =  graph.number_of_nodes()
        self.loops = {}    
        for node in  graph.nodes():
            com = solution[node]
            membership[node] = com
            deg = float( graph.degree( node, 'weight'))
            self.DegCom[com] = self.DegCom.get(com, 0.) + deg
            self.Degree[node] = deg
            edge_data = graph.get_edge_data(node, node, {weight: 0})
            self.loops[node] = float(edge_data.get(weight, 1))
            inc = 0.
            for neighbor, datas in graph[node].items():
                edge_weight = datas.get( 'weight', 1.)
                if edge_weight <= 0:
                    error = "Bad graph type ({})".format(type( graph))
                    raise ValueError(error)                 
                if solution[neighbor] == com:
                    if neighbor == node:
                        inc += float( edge_weight)
                    else:
                        inc += float( edge_weight) / 2.
                
            self.internal[com] = self.internal.get( com, 0.) + inc 

        return membership

    def best_sol( self, dend, level):
         
        #print(dend) 
        partition = dend[0].copy()
        #print("papap", partition)
        for index in range( 1, level + 1):
            for node, community in partition.items():
                partition[node] = dend[index][community]
        
        return partition
    


    def tarans ( self, clust_label):
        clust = set(clust_label.values())
        clusters = []
        for i in clust:
            cluste = set()
            for j in clust_label:
                if clust_label[j] == i:   
                    cluste.add(j)
            
            clusters.append(cluste)

        return clusters


def Read_Graph( Path):
        
    if Path[len(Path)-3: ] == 'txt' or Path[len(Path)-3: ] == 'dat':
        Graph = nx.read_edgelist( Path, nodetype = int , data = [('weight', float)])
        graph = nx.Graph(Graph)
    elif Path[len(Path)-3: ] == 'gml':
        Graph = nx.read_gml( Path, label = 'id')
        graph = nx.Graph(Graph)
    else :
        raise TypeError (" the type of graph is not suportable or not no such file or directory")

    return graph 

                