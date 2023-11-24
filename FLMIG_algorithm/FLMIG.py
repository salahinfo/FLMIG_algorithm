from cProfile import label
import random
from collections import deque,Counter
import copy
import time 
import sys
import math
from scipy.stats import expon
from sklearn.metrics.cluster import normalized_mutual_info_score
from GraphTools import GraphTolls, Read_Graph
import matplotlib.pyplot as  plt 
import networkx as nx



class Fast_local_Move_IG(GraphTolls) :
    def __init__( self, Nb, Beta,path, graph):    
        self.Nb = Nb
        self.Beta = Beta
        self.Mod_val = 0
        super(Fast_local_Move_IG, self).__init__(path, graph)

    def expon ( self , value , teta):

        value = math.log(value)
        x = value / teta
        p = math.exp(x)

        return p 

                          

    def Destruction( self, graph):
        
        s = self.renumber()
        super().init( graph, s, weight='weight') 
        cut_len = int(len(self.membership)* float(self.Beta)) 
        drop_node = random.sample( list(self.membership.keys()), cut_len )
        for al in drop_node:
            com_id = self.membership[al]
            wgh = super().neigh_comm(al, graph)    
            super().delet_node( al, com_id, wgh.get( com_id, 0.))
            if al not in self.membership.values():
                super().insert_node( al, al, wgh.get( al, 0.))
            else:
                com_id = super().generate_random_not_in_list(set(self.membership.values()))
                super().insert_node( al, com_id, wgh.get( com_id, 0.))     
   
        return  self.membership, drop_node 


    def Reconstruction( self,  graph, soltion, drop_node):
          
        self.__randomcom( graph, drop_node)
        #self.draw_communities(graph,soltion)
        soltion = self.con_dense( graph)
                
        return  soltion


    def flocalmove ( self, graph, resolution=1):
    
        Nodelist = deque(graph.nodes())
        random.shuffle(Nodelist)
        while Nodelist:
            node = Nodelist.popleft()
            com_node = self.membership[node]
            degc_totw = self.Degree.get(node, 0.) / (self.m * 2.)  # NOQA
            neigh_communities = super().neigh_comm( node, graph )
            super().delet_node(node, com_node, neigh_communities.get(com_node, 0.))
            best_com = com_node
            best_increase = 0
            for com, dnc in neigh_communities.items():
                Delat_Q = resolution * dnc - self.DegCom.get(com, 0.) * degc_totw
                    #print(incr)
                if Delat_Q > best_increase:
                    best_increase = Delat_Q
                    best_com = com
            super().insert_node(node, best_com, neigh_communities.get(best_com, 0.))
            if best_com != com_node:
                for veg in graph[node]:
                    if self.membership[veg] != best_com:    
                        Nodelist.append(veg)
           
        
        return self.membership


    def con_dense(self, graph):
        p_list = []
        n_mod = super().modularity()
        mod_graph = graph.copy()
        p = super().renumber()
        p_list.append(p)
        mod_graph = super().induced_graph( p, mod_graph, weight='weight')
        super().modifie_status( mod_graph, weight = 'weight')
        Q_val = n_mod
        while True :
            solution = self.flocalmove(mod_graph)
            #solution = self.FL_move(mod_graph)
            n_mod = super().modularity()
            if n_mod - Q_val < 0.0000001:
                break
            
            Q_val = n_mod
            p = super().renumber()
            p_list.append(p)
            mod_graph = super().induced_graph( p, mod_graph, weight='weight')
            super().modifie_status( mod_graph, weight='weight')
        
        #self.draw_communities( mod_graph, p)        
        P = super().best_sol( p_list, len(p_list)- 1)    
        super().init( graph, P, weight='weight')                                                                          

        return self.membership    
                
    def __randomcom ( self , graph, drop_node):
        random.shuffle(drop_node)
        for vsele in  drop_node:    
            #degree = self.Degree[vsele]
            qum={}
            check = False
            com_befor = self.membership[vsele]
            ngh_com = super().neigh_comm( vsele, graph)
            degc_totw = self.Degree.get(vsele, 0.) / (2.*self.m )
            #print(degc_totw)
            check = False
            com_n = com_befor
            Delat_Q = 0
            for com, dvcp  in ngh_com.items() :
                if com != com_befor:
                    Delat_Q = dvcp - self.DegCom.get(com, 0.) * degc_totw
                    if Delat_Q > 0:
                        #print(Delat_Q) 
                        qum[com] = Delat_Q
                        check = True

            if check :
                prb = [ self.expon( i, 0.01) for k,i in qum.items() ]
                com_n = super().weighted_choice( list(qum.keys()), prb)
                super().delet_node( vsele, com_befor, ngh_com.get( com_befor, 0.))
                super().insert_node( vsele, com_n, ngh_com.get( com_n, 0.))                         
          
    
    
    def draw_communities( self, G, community_map, alpha = 1):
         
        cmap = plt.get_cmap("magma")
        pos = nx.spring_layout(G)
        indexed = [community_map.get(node) for node in G]
        plt.axis("off")
        node_size = G.number_of_nodes()      
        nx.draw( G, with_labels=True, pos = pos, node_color = indexed, node_size = node_size*3 ,font_size = 10, font_color = 'white', font_weight = 'bold')
        labels = nx.get_edge_attributes(G, 'weight')
        nx.draw_networkx_edge_labels(G, pos, edge_labels=labels)
        plt.title('(a)')
        plt.show()
        
    
                
    def Run_FMLIG ( self, graph):
        start = time.time()
        super().modifie_status( graph, weight='weight') 
        soltion = self.flocalmove( graph)
        print( time.time() - start)
        #self.draw_communities( graph, soltion)
        best_solution = copy.copy( soltion)
        #print(" internal ", self.internal# )
        #super().init( graph, soltion)
        Q_best = super().modularity()
        print( "q1", Q_best)
        T_init = 0.025 * Q_best
        T = T_init
        nb_iter = 0
        status_list = []
        while nb_iter < self.Nb:
            Q1 = super().modularity()
            
            incumbent_solution = copy.copy(soltion)
        
            soltion, drop_nod = self.Destruction(graph)
            print(soltion)    
            soltion = self.Reconstruction( graph,soltion, drop_nod)
            #self.draw_communities(graph, soltion) 
            super().init(graph, soltion, weight = 'weight')  
        
            soltion = self.flocalmove(graph)
                                  
            Q2 = super().modularity()
            print("Q2", Q2, " number of iteration", nb_iter)
            #print("q2",Q2, "time", time.time() - start)        
            if Q2 > Q_best:
                best_solution = copy.copy(soltion)
                Q_best = Q2
                           
            P = random.random()
            if Q2 < Q1 and P > math.exp((Q2 - Q1)//T):
                soltion = copy.copy(incumbent_solution)
                T = T_init

            elif Q2  >=  Q1:
                T = T_init
              
            else:
                T = T*0.9
                           
            nb_iter = nb_iter + 1
        
        #print(status_list)
        end = time.time()
        t = end - start
 
        return Q_best, best_solution, t
        

def de_main():
    path = sys.argv[1]
    Number_iter = int(sys.argv[2])
    Beta = sys.argv[3]
    
    graph = Read_Graph(path)
    data = GraphTolls(path, graph)
    NMI_list = []
    Time_list = [] 
    Q_list = []
    nb_run = 0
    while nb_run < int(sys.argv[5]) :
        #print("rb",nb_run)
        communities = Fast_local_Move_IG( Number_iter, Beta, path, graph)
        mod,community,tim = communities.Run_FMLIG(graph) 
        #print(community)
        communities.draw_communities(graph, community)
        Q_list.append(mod)
        Time_list.append(tim)
        #label = communities.lebel_node(community)  
        if sys.argv[4]!= 'None':
            True_partition = data.Read_GroundTruth(sys.argv[4])
            #print(True_partition)
            community = dict(sorted(community.items()))
            NMI = normalized_mutual_info_score( True_partition, list(community.values()))
            NMI_list.append(NMI)
            
            
        elif sys.argv[4] == 'None':
            pass            
      
        nb_run = nb_run +1
            
    if sys.argv[4] != 'None':
        Q_avg = communities.avg(Q_list)
        Q_max = communities.max(Q_list)
        Q_std = communities.stdev(Q_list)
        NMI_max = communities.max(NMI_list)
        time_run = communities.avg(Time_list)
        return NMI_max, Q_max, Q_avg, Q_std, time_run     
    elif sys.argv[4] == 'None':
        Q_avg = communities.avg(Q_list)
        Q_max = communities.max(Q_list)
        Q_std = communities.stdev(Q_list)
        time_run = communities.avg(Time_list)
        return Q_max, Q_avg, Q_std, time_run

if __name__ == '__main__':

   
    if  sys.argv[4] != 'None':
        NMI_max, Q_Max,Q_avg , Q_std, time_run = de_main()
        print("NMI_max",NMI_max)
        print("the value of Q_max",Q_Max)
        print("the value of Q_avg",Q_avg)
        print("the value of Q_std",Q_std)
        print("time ",time_run)
    elif sys.argv[4] == 'None' :
        Q_max, Q_avg, Q_std, time_run = de_main()
        print("the value of Q_max",Q_max)
        print("the value of Q_avg",Q_avg)
        print("the value of Q_std",Q_std)
        print("the value of time ",time_run) 
     
