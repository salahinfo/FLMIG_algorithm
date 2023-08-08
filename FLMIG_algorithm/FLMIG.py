import random
from collections import deque,Counter
import copy
import time 
import sys
import math
from scipy.stats import expon
from sklearn.metrics.cluster import normalized_mutual_info_score
from GraphTools import GraphTolls




class Fast_local_Move_IG(GraphTolls) :
    def __init__(self,Nb,Beta,path):    
        self.Nb = Nb
        self.Beta = Beta
        self.Mod_val = 0
        super(Fast_local_Move_IG, self).__init__(path)

    def expon (self , value , teta):
        x = (1/teta)*(value)
        p =float(math.exp(x))
        return p      
    
    def GCH(self):
        
        vertex_list = list(self.graph.nodes())
        node = random.choice(vertex_list)
        #print("nnnn",node)
        com_id = 0
        self.membership[node] = com_id
        self.DegCom[com_id] = self.Degree[node]
        #print(self.Degree)
        vertex_list.remove(node)
        for node in  vertex_list:    
            comm_ngh = super().neigh_comm( node )
            MAX_Q = 0
            pos = -1
            for com, Kbv in comm_ngh.items():    
                db = self.DegCom[com]
                delta_Q = 1/self.m * Kbv - self.Degree[node]/(2*self.m**2)*db
                if delta_Q > MAX_Q:
                    MAX_Q = delta_Q
                    pos = com
                else :
                    delta_Q = 0
            
            if MAX_Q > 0:
                self.membership[node] = pos
                self.DegCom[pos] = self.DegCom[pos] + self.Degree[node]
            else:
                com_id = com_id + 1
                self.membership[node] = com_id
                self.DegCom[com_id] = self.Degree[node]                     
        
        return  self.membership 

    def Destruction( self):       
        drop_node = []
        merg_node = []
        cut_len = int(len(self.DegCom)* float(self.Beta)) 
        index_community = random.sample( list(self.DegCom.keys()), cut_len )
        for com in index_community:
            for node, com_id in self.membership.items():
                if com_id == com:
                    self.membership[node] = None                    
                    drop_node.append(node)   
            del self.DegCom[com] 
        merg_node = [nod for nod in self.Node_list if nod not in drop_node] 
          
        return  self.membership, drop_node, merg_node
    
    def Reconstruction( self, drop_node, merg_node):
        random.shuffle(drop_node)
        for node in  drop_node:    
            comm_ngh = super().neigh_comm( node)
            MAX_Q = 0
            pos = -1
            for com, Kbv in comm_ngh.items():    
                db = self.DegCom[com]
                delta_Q = 1/self.m * Kbv - self.Degree[node]/(2*self.m**2)*db
                if delta_Q > MAX_Q:
                    MAX_Q = delta_Q
                    pos = com
                else :
                    delta_Q = 0
            
            if MAX_Q > 0:
                self.membership[node] = pos
                self.DegCom[pos] = self.DegCom[pos] + self.Degree[node]
            else:
                com_id = super().generate_random_not_in_list(list(self.DegCom.keys()),0,self.n)
                self.membership[node] = com_id
                self.DegCom[com_id] = self.Degree[node] 
                                    
        for vsele in merg_node:
            degree = self.Degree[vsele]
            qum={}
            com_befor = self.membership[vsele]
            ngh_com = super().neigh_comm( vsele)
            #print(ngh_com)
            dvc = super().ngh_node( vsele, com_befor)
            devc = self.DegCom[com_befor]              
            for com  in ngh_com :
                if com != com_befor:
                    dvcp = ngh_com[com]
                    devcp = self.DegCom[com]
                    deq = (1/self.m)*(dvcp-dvc)-(degree/(2*self.m**2))*(devcp-devc+degree)
                    if deq > 0: 
                        qum[com] = deq
            
            if len(qum) != 0:  
                prb = [ self.expon(i, 0.1) for k,i in qum.items() ]
                com_n = super().weighted_choice( list(qum.keys()), prb)
                self.membership[vsele] = com_n
                self.DegCom[com_n] = self.DegCom[com_n] + self.Degree[vsele]
                self.DegCom[com_befor] = self.DegCom[com_befor] - self.Degree[vsele]
                if self.DegCom[com_befor] <= 0:
                    del self.DegCom[com_befor]
        
        #print("degree community" ,self.DegCom)
        modified = True
        while modified :
            modified = False        
            community = list(self.DegCom.keys())
            #print("befor reconstruction", self.membership, community)
            for com1 in community:
                maxq = 0
                pos = -1
                ngh_comm = super().com_ngh_com(com1)
                #print(com1 , ngh_comm)
                for comm ,Kbv in ngh_comm.items():
                    delta_Q =  Kbv / self.m - ( self.DegCom[com1] * self.DegCom[comm] ) / ( 2 * self.m**2 )
                    #print("delta",delta_Q)
                    if delta_Q > maxq :
                        maxq= delta_Q
                        pos = comm
                            
                if maxq > 0:
                    modified = True
                    self.merge_com( com1, pos)              
                    self.DegCom[com1] = self.DegCom[com1] + self.DegCom[pos]
                    del self.DegCom[pos]
                                                  
        return  self.membership


    def FL_move( self):
        Qv = deque(self.Node_list)
        random.shuffle(Qv)
        while Qv:
            vsele = Qv.popleft()
            degree = self.Degree[vsele]
            #qum = 0
            com_befor = self.membership[vsele]
            #print(com_befor)
            ngh_com = super().neigh_comm(vsele)
            #print(ngh_com)
            dvc = super().ngh_node( vsele, com_befor)
            devc = self.DegCom[com_befor]
            maxq = 0              
            for com, dvcp  in ngh_com.items() :
                if com !=  com_befor: 
                    devcp = self.DegCom[com]
                    deq = ( 1/self.m )*( dvcp-dvc )-( degree/(2*self.m**2) )*( devcp-devc+degree )
                    if deq > maxq :
                        maxq = deq
                        m_com = com
                        
            if maxq > 0:
                self.membership[vsele] = m_com
                self.DegCom[m_com] = self.DegCom[m_com] + self.Degree[vsele]
                self.DegCom[com_befor] = self.DegCom[com_befor] - self.Degree[vsele]
                if self.DegCom[com_befor] <= 0:
                    del self.DegCom[com_befor]
                                         
                Neigh = list(self.adjency[vsele])
                for veg in Neigh:
                    if self.membership[veg] != self.membership[vsele] and veg not in Qv:    
                        Qv.append(veg)       
                                      
        return self.membership


    def Run_FMLIG (self):
        start = time.time()
        soltion = self.GCH()
        soltion = self.FL_move()
        #print("after fast local" , soltion)    
        best_solution = copy.copy(soltion)
        #print("best solution", best_solution)
        Q_best = super().modularity()
        #print("qqqqqqqqqqqq",Q_best)
        T_init = 0.025 * Q_best
        T = T_init
        nb_iter = 0
        while nb_iter < self.Nb:
            Q1 = super().modularity()
            print("q1", Q1)
            incumbent_solution = copy.copy(soltion)
            soltion,drop_nodes,merg_node = self.Destruction()
            soltion = self.Reconstruction( drop_nodes, merg_node)
            soltion = self.FL_move() 
                                   
            Q2 = super().modularity()
            print( "q2" , Q2)
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
        
        end = time.time()
        t = end-start
        
        return Q_best, best_solution, t
        

def de_main():
    path = sys.argv[1]
    Number_iter = int(sys.argv[2])
    Beta = sys.argv[3]
    data = GraphTolls(path)
    graph = data.Read_Graph()
    communities = Fast_local_Move_IG(Number_iter, Beta, path)
    NMI_list = []
    Time_list = [] 
    Q_list = []
    nb_run = 0
    while nb_run < int(sys.argv[5]) :
        print("rb",nb_run)
        mod,community,tim = communities.Run_FMLIG() 
        Q_list.append(mod)
        Time_list.append(tim)
        #label = communities.lebel_node(community)  
        if sys.argv[4]!= 'None':
            True_partition = data.Read_GroundTruth(sys.argv[4])
            NMI = normalized_mutual_info_score(True_partition,community)
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
     
