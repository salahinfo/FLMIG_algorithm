import random
from collections import deque
import copy
import time 
import sys
import math
import networkx.algorithms.community as nx_comm
from scipy.stats import expon
from sklearn.metrics.cluster import normalized_mutual_info_score
from GraphTools import GraphTolls



class Fast_local_Move_IG(GraphTolls) :
    def __init__(self,G,Nb,Beta):    
        self.G = G
        self.Nb = Nb
        self.Beta = Beta
        self.m = G.number_of_edges()
        self.n = G.number_of_nodes()
        self.Mod_val = 0

    def expon (self , value , teta):
        x = (1/teta)*(value)
        p =float(math.exp(x))
        return p      
    
    def GCH(self):
        community = []
        vertex_list = [i for i in self.G.nodes()]
        node = random.choice(vertex_list)
        vertex_list.remove(node)
        community.append({node})
        while vertex_list != []:    
            node = random.choice(vertex_list)
            vertex_list.remove(node)
            MAX_Q = 0
            pos = -1
            for index,clusters in enumerate(community):
                if super().is_edge_betw(self.G,node,clusters):        
                    Kbv = super().select_edge_betw(self.G,node,clusters)
                    db = sum([self.G.degree[j] for j in clusters])
                    #print(clusters , db)
                    delta_Q = 1/self.m * Kbv -self.G.degree[node]/(2*self.m**2)*db
                    if delta_Q > MAX_Q:
                        MAX_Q = delta_Q
                        pos = index
                else :
                    delta_Q = 0

            if MAX_Q > 0:
                community[pos].add(node)
            else:
                community.append({node})

            
        return  community 

    def Destruction(self,community):       
        drop_node = []
        merg_node = []
        graph_node = [i for i in self.G.nodes()]
        cut_len = int(len(community)* float(self.Beta))
        mov_len = int(self.n * float(self.Beta))
        moved_node = graph_node[:mov_len]
        index_community = random.sample(community,cut_len)
        for cluster in index_community:
            inedex_comu = community.index(cluster)
            for v in cluster :
                drop_node.append(v)
            del community[inedex_comu]
               
        for clust in community:
            for node in clust:
                merg_node.append(node)

              
        return  community, drop_node, merg_node
    
    def Reconstruction(self,community,drop_node, merg_node):
        
        random.shuffle(drop_node)
        
        for node in drop_node:
            MAX_Q = 0
            pos = -1
            for index,clusters in enumerate(community):  
                if super().is_edge_betw(self.G,node,clusters): 
                    Kbv = super().select_edge_betw(self.G,node,clusters)
                    db = sum([self.G.degree[j] for j in clusters])
                    delta_Q = 1/self.m * Kbv - self.G.degree[node]/(2*self.m**2)*db
                    if delta_Q > MAX_Q:
                        MAX_Q = delta_Q
                        pos = index
    
            if MAX_Q > 0:
                community[pos].add(node)
            else:
                community.append({node})
        
        
        cut_len = int(float(self.n)* float(self.Beta))
        while merg_node !=[]:
            vsele= random.choice(merg_node)
            degree = self.G.degree[vsele]
            qum=[]
            proba_ngh = []
            l = []
            chek = False
            merg_node.remove(vsele)
            for clusters in community:
                if vsele in clusters:
                    original_clusters = clusters
                    dvc = super().select_edge_betw(self.G,vsele,original_clusters)
                    devc = sum([self.G.degree[j] for j in clusters])
                    break 
            for index,clusters in enumerate(community) :
                if vsele in clusters:
                    deq = 0
                    beforr = index
                elif super().is_edge_betw(self.G,vsele,clusters): 
                    dvcp = super().select_edge_betw(self.G,vsele,clusters)
                    devcp = sum([self.G.degree[j] for j in clusters])
                    deq = (1/self.m)*(dvcp-dvc)-(degree/(2*self.m**2))*(devcp-devc+degree)
                    if deq > 0: 
                        qum.append(deq)
                        l.append(index)
                        
            if qum !=[]:    
                prb = [self.expon(i, 0.1) for i in qum]
                poss = super().weighted_choice(l , prb)
                community[poss].add(vsele)
                community[beforr].remove(vsele)
                if community[beforr] == []:
                    del  community[beforr]
        
        for community_condidate in community: 
            maxq=0
            pos=-1
            db = sum([self.G.degree[j] for j in community_condidate])
            for index,cluster in enumerate(community):

                if cluster == community_condidate:
                    maxq = 0 
                else:    
                    dbc = sum([self.G.degree[j] for j in cluster])
                    Kbv =  super().select_edge_c(self.G,cluster,community_condidate)
                    delta_Q =  Kbv - (dbc*db)/(2*self.m)
                    if delta_Q > maxq:
                            maxq=delta_Q
                            pos = index
                            
            if maxq > 0:
                community[pos].update(community_condidate)
                community.remove(community_condidate)
                
        return  community


    def FL_move(self,community): 
        Qv = deque([i for i in self.G.nodes()])
        random.shuffle(Qv)
        while Qv:
            vsele = Qv.popleft()
            degree = self.G.degree[vsele]
            qum=[]
            for clusters in community:
                if vsele in clusters:
                    original_clusters = clusters
                    dvc = super().select_edge_betw(self.G,vsele,original_clusters)
                    devc = sum([self.G.degree[j] for j in clusters])
                    break 
                
            for index,clusters in enumerate(community) :
                if vsele in clusters:
                    deq = 0
                    befor = index
                    qum.append(deq)
                elif super().is_edge_betw(self.G,vsele,clusters): 
                    dvcp = super().select_edge_betw(self.G,vsele,clusters)
                    devcp = sum([self.G.degree[j] for j in clusters])
                    deq = (1/self.m)*(dvcp-dvc)-(degree/(2*self.m**2))*(devcp-devc+degree)
                    qum.append(deq)
                else :
                    qum.append(0)
            
            if max(qum) > 0:
                pos = qum.index(max(qum))
                community[pos].add(vsele)
                community[befor].remove(vsele)
                Neigh = list(self.G.neighbors(vsele))
                for veg in Neigh:
                    if veg not in community[pos] and veg not in Qv:    
                           Qv.append(veg)       
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
        print(soltion)
        soltion = self.FL_move(soltion)
        best_solution = copy.deepcopy(soltion)
        Q_best = nx_comm.modularity(self.G, soltion)
        T_init = 0.025*Q_best
        T = T_init
        nb_iter = 0
        while nb_iter < self.Nb:
            Q1 = nx_comm.modularity(self.G, soltion)
            incumbent_solution = copy.deepcopy(soltion)
            soltion,drop_nodes,merg_node = self.Destruction(soltion)
            soltion = self.Reconstruction(soltion,drop_nodes,merg_node)
            soltion = self.FL_move(soltion) 
                       
            Q2 = nx_comm.modularity(self.G, soltion)
            if Q2 > Q_best:
                best_solution = copy.deepcopy(soltion)
                Q_best = Q2
                
            P = random.random()
            
            if Q2 < Q1 and P > math.exp((Q2 - Q1)//T):
                soltion = copy.deepcopy(incumbent_solution)
                T = T_init

            elif Q2  >=  Q1:
                T = T_init
              
            else:
                T = T*0.9
        
            nb_iter = nb_iter + 1
        
        end = time.time()
        t = end-start
        
        return Q_best, best_solution,t
        

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
    while nb_run < int(sys.argv[5]) :
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
        print("the value of Q_max",Q_Max)
        print("the value of Q_avg",Q_avg)
        print("the value of Q_std",Q_std)
        print("time ",time_run)
    elif sys.argv[4] == 'None' :
        Q_max, Q_avg, Q_std,time_run = de_main()
        print("the value of Q_max",Q_max)
        print("the value of Q_avg",Q_avg)
        print("the value of Q_std",Q_std)
        print("the value of time ",time_run) 
     
