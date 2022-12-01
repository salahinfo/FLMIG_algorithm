from GraphTools import GraphTolls
import random
from collections import deque
import copy
import time 
from sklearn.metrics.cluster import normalized_mutual_info_score
import sys
import math


class Fast_local_Move_IG(GraphTolls) :
    def __init__(self,G,Nb,Beta):    
        self.G = G
        self.Nb = Nb
        self.Beta = Beta
        self.m = G.number_of_edges()
        self.n = G.number_of_nodes()
        self.Mod_val = 0
        self.Qv = deque([i for i in self.G.nodes()])
        
    def GCH(self):
        community = []
        vertex_list = [i for i in self.G.nodes()]
        node = random.choice(vertex_list)
        vertex_list.remove(node)
        community.append([node])
        while vertex_list != []:    
            node = random.choice(vertex_list)
            vertex_list.remove(node)
            MAX_Q = 0
            pos = -1
            for index,clusters in enumerate(community):    
                Kbv = super().select_edge_betw(self.G,node,clusters)
                db = sum([j for k,j in self.G.degree(clusters)])
                delta_Q = 1/self.m * Kbv -self.G.degree(node)/(2*self.m**2)*db
                if delta_Q > MAX_Q:
                    MAX_Q = delta_Q
                    pos = index
    
            if MAX_Q > 0:
                community[pos].append(node)
            else:
                community.append([node])
            
        return  community 

    def Destruction(self,community):
        vertex_list = [i for i in self.G.nodes()]
        drop_node = []
        cut_len = int(float(self.n)* float(self.Beta))
        random.shuffle(vertex_list)
        drop = vertex_list[self.n-cut_len:]
        for i in range(cut_len) :
            v = drop.pop()
            drop_node.append(v)
            for cluster in community:
                if v in cluster:
                    cluster.remove(v)
                    if len(cluster) == 0: 
                        del cluster 


        return  community , drop_node

    def Reconstruction(self,community,drop_node):

        random.shuffle(drop_node)
        for node in drop_node:
            MAX_Q = 0
            pos = -1
            for index,clusters in enumerate(community):    
                Kbv = super().select_edge_betw(self.G,node,clusters)
                db = sum([j for k,j in self.G.degree(clusters)])
                delta_Q = 1/self.m * Kbv - self.G.degree(node)/(2*self.m**2)*db
                if delta_Q > MAX_Q:
                    MAX_Q = delta_Q
                    pos = index
    
            if MAX_Q > 0:
                community[pos].append(node)
            else:
                community.append([node])

        while drop_node !=[]:
                vsele= random.choice(drop_node)
                degree = self.G.degree(vsele)
                qum=[]
                drop_node.remove(vsele)
                for clusters in community:
                    if vsele in clusters:
                        original_clusters = clusters
                        dvc = super().select_edge_betw(self.G,vsele,original_clusters)
                        devc = sum([j for k,j in self.G.degree(clusters)])
                        break 
                for index,clusters in enumerate(community) :
                    if vsele in clusters:
                        deq = 0
                        beforr = index
                        qum.append(deq)
                    else:
                        dvcp = super().select_edge_betw(self.G,vsele,clusters)
                        devcp = sum([j for k,j in self.G.degree(clusters)])
                        deq = (1/self.m)*(dvcp-dvc)+(degree/(2*self.m**2))*(devc-devcp-degree)
                        qum.append(deq)

                index_com = random.choice(qum)
                if index_com > 0:
                    poss = qum.index(index_com)
                    community[poss].append(vsele)
                    community[beforr].remove(vsele)
                    if community[beforr] == []:
                        del  community[beforr]
    
                
        return  community


    def FL_move(self,community):

        visted = [0 for i in range (self.G.number_of_nodes())]
        random.shuffle(self.Qv)
        while self.Qv:
            vsele = self.Qv.popleft()
            degree = self.G.degree(vsele)
            qum=[]
            for clusters in  community:
                if vsele in clusters:
                    original_clusters = clusters
                    dvc = super().select_edge_betw(self.G,vsele,original_clusters)
                    devc = sum([j for k,j in self.G.degree(clusters)])
                    break 
            for index,clusters in enumerate(community) :
                if vsele in clusters:
                    deq = 0
                    befor= index
                    qum.append(deq)
                else:
                    dvcp= super().select_edge_betw(self.G,vsele,clusters)
                    devcp = sum([j for k,j in self.G.degree(clusters)])
                    deq = (1/self.m)*(dvcp-dvc)+(degree/(2*self.m**2))*(devc-devcp-degree)
                    qum.append(deq)
         
            if max(qum) > 0:
                pos = qum.index(max(qum))
                community[pos].append(vsele)
                community[befor].remove(vsele)
                Neigh= list(self.G.neighbors(vsele))
                for index,veg in enumerate(Neigh):
                    if (veg not in community[pos]):
                        if (veg not in self.Qv) and (visted[index]==0):
                           self.Qv.append(veg)
                           visted[index]=1
                          

                if community[befor] == []:
                    del community[befor]
    

        return community

    def lebel_node (self,community):
        label = sorted([i for i in self.G.nodes()])
        for index,no in enumerate (label):
            for i in range(len(community)):
                if no in community[i]:
                    label[index] = i
        
        return label

    def Run_FMLIG (self):
        start = time.time()
        soltion = self.GCH()
        soltion = self.FL_move(soltion)
        best_solution = copy.deepcopy(soltion)
        T_init = 0.025*super().Modularity(self.G,soltion)
        T = T_init
        nb_iter = 0
        while nb_iter < self.Nb:
            Q1 = super().Modularity(self.G,soltion)
            incumbent_solution = copy.deepcopy(soltion)
            soltion,drop_nodes = self.Destruction(soltion)
            soltion = self.Reconstruction(soltion,drop_nodes)
            soltion = self.FL_move(soltion)            
            Q2 = super().Modularity(self.G,soltion)
            if Q2 > super().Modularity(self.G,best_solution):
                best_solution = copy.deepcopy(soltion)
                
            P = random.random()
            if Q2 < Q1 and P > math.exp((Q2 - Q1)//T):
                soltion = copy.deepcopy(incumbent_solution)
                T = T_init

            elif Q2  >=  Q1:
                T = T_init
              
            else:
                T = T*0.9
        
            nb_iter = nb_iter + 1
        
        self.Mod_val = super().Modularity(self.G,best_solution)
        end = time.time()
        t = end-start
        
        return self.Mod_val, best_solution,t
        

def de_main():
    path = sys.argv[1]
    Number_iter = int(sys.argv[2])
    Beta = sys.argv[3]
    data = GraphTolls()
    graph = data.Read_Graph(path)
    NMI_list = []
    Time_list = [] 
    Q_list = []
    nb_run = 0
    while nb_run <= int(sys.argv[5]) :
        communities = Fast_local_Move_IG(graph, Number_iter, Beta)
        mod,community,tim = communities.Run_FMLIG() 
        Q_list.append(mod)
        Time_list.append(tim)
        label = communities.lebel_node(community)  
        if sys.argv[4]!= 'None':
            True_partition = data.Read_GroundTruth(sys.argv[4])
            NMI = normalized_mutual_info_score(True_partition,label)
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
        return NMI_max, Q_max, Q_avg, Q_std,time_run     
    elif sys.argv[4] == 'None':
        Q_avg = communities.avg(Q_list)
        Q_max = communities.max(Q_list)
        Q_std = communities.stdev(Q_list)
        time_run = communities.avg(Time_list)
        return Q_max, Q_avg, Q_std,time_run

if __name__ == '__main__':

   
    if  sys.argv[4]!= 'None':
        NMI_max,Q_Max,Q_avg ,Q_std,time_run = de_main()
        print("NMI_max",NMI_max)
        print("Q_max",Q_Max)
        print("Q_avg",Q_avg)
        print("Q_std",Q_std)
        print("time ",time_run)
    elif sys.argv[4] == 'None' :
        Q_max, Q_avg, Q_std,time_run = de_main()
        print("Q_max",Q_max)
        print("Q_avg",Q_avg)
        print("Q_std",Q_std)
        print("time ",time_run) 
     