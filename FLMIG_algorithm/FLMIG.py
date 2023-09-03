import random
from collections import deque,Counter
import copy
import time 
import sys
import math
from scipy.stats import expon
from sklearn.metrics.cluster import normalized_mutual_info_score
from GraphTools import GraphTolls, Read_Graph




class Fast_local_Move_IG(GraphTolls) :
    def __init__( self, Nb, Beta,path, graph):    
        self.Nb = Nb
        self.Beta = Beta
        self.Mod_val = 0
        super(Fast_local_Move_IG, self).__init__(path, graph)

    def expon (self , value , teta):
        x = (1/teta)*(value)
        p =float(math.exp(x))
        return p      
    

                         
    def GCH( self, graph):
         
        vertex_list = list(graph.nodes())
        node = random.choice(vertex_list)
        #print("nnnn",node)
        com_id = 0
        self.membership[node] = com_id
        self.DegCom[com_id] = self.Degree[node]
        #print(self.Degree)
        vertex_list.remove(node)
        for node in  vertex_list:    
            comm_ngh = super().neigh_comm( node , graph)
            MAX_Q = 0
            pos = -1
            for com, Kbv in comm_ngh.items():    
                db = self.DegCom[com]
                delta_Q =  Kbv - self.Degree[node]/(2.*self.m)*db
                #print(delta_Q)
                if delta_Q > MAX_Q:
                    MAX_Q = delta_Q
                    pos = com
                else :
                    delta_Q = 0
            
            if MAX_Q > 0:
                super().insert_node( node, pos, comm_ngh.get( pos , 0))
            else:
                com_id = com_id + 1
                super().insert_node(node , com_id, comm_ngh.get(com_id, 0))
                                     
        
        return  self.membership 

    def Destruction( self, graph):       
        cut_len = int(len(self.membership)* float(self.Beta)) 
        drop_node = random.sample( list(self.membership.keys()), cut_len )
        #print("list", cut_len, index_community)
        #print("degcomunity", self.DegCom)
        #print("mm",self.membership)
        for al in drop_node:
            com_id = self.membership[al]
            wgh = super().neigh_comm(al, graph)    
            super().delet_node( al, com_id, wgh.get( com_id, 0.))
                                    
                      
        #merg_node = [ nod for nod in self.Node_list if nod not in index_community] 
        return  self.membership, drop_node
    
    def Reconstruction( self, drop_node, graph):
        p_list = []
        self.__affect_node( graph, drop_node)
        #print(solution)
        #print( self.internal)
        #print(self.DegCom)
        self.__randomcom( graph)
        s = self.__merge_community(graph)
        #print("befor internal", self.internal)
        n_mod = super().modularity()
        #print("befor", n_mod)
        mod_graph = graph.copy()
        p = super().renumber()
        p_list.append(p)
        mod_graph = super().induced_graph( p, mod_graph, weight='weight')
        #print(mod_graph)
        super().modifie_status( mod_graph, weight = 'weight')
        Q_val = n_mod
        while True :
            #print(len(s))
            solution = self.FL_move(mod_graph)
            #print(len(solution))
            #solution = self.__merge_community(mod_graph)
            #print(len(solution))
            n_mod = super().modularity()
            #print(n_mod)
            
            if n_mod - Q_val < 0.0000001:
                break
        
            Q_val = n_mod
            p = super().renumber()
            p_list.append(p)
            mod_graph = super().induced_graph( p, mod_graph, weight='weight')
            super().modifie_status( mod_graph, weight='weight')
        
        P = super().best_sol( p_list, len(p_list)- 1)    
        super().init( graph, P, weight='weight')                                                                          
        #print("after", super().modularity())
        #print("after internal",self.internal)    
        return  self.membership


    def FL_move( self, graph):
        Qv = deque([ i for i in graph.nodes()])
        random.shuffle(Qv)
        while Qv:
            vsele = Qv.popleft()
            #print(vsele)
            degree = self.Degree[vsele]
            #qum = 0
            #print(degree)
            #print(self.loops)
            com_befor = self.membership[vsele]
            ngh_com = super().neigh_comm(vsele, graph)
            #print("combefor", com_befor,"nghhhh", ngh_com)
            dvc = super().ngh_node(vsele, com_befor, graph)
            devc = self.DegCom[com_befor]
            maxq = 0
            #print("vsele", vsele)
            #print(self.m)
            #dd = (degree) / ( self.m * 2.)
            m_com = com_befor             
            for com, dvcp in ngh_com.items() : 
                devcp = self.DegCom[com]
                #deq = ( 1/self.m ) * ( dvcp - dvc ) - ( (degree)/(2.*self.m**2.)) *( devcp-devc+degree )
                #deq = dvcp - (devcp )*(degree / (self.m * 2.))
                deq = (1/self.m)*(dvcp-dvc)+(degree/(2.*self.m**2.))*(devc-devcp -degree)
                #print(deq)
                if deq > maxq :
                    maxq = deq
                    m_com = com
                        
            super().delet_node( vsele, com_befor, ngh_com.get( com_befor, 0.))            
            super().insert_node( vsele, m_com, ngh_com.get( m_com, 0.))            
            if m_com != com_befor:                                       
                for veg in graph[vsele]:
                    if self.membership[veg] != m_com and veg not in Qv:    
                        Qv.append(veg)       
                                      
        return self.membership
    
    def __affect_node ( self, graph, drop_node = None):
        random.shuffle(drop_node)
        for node in  drop_node:    
            comm_ngh = super().neigh_comm( node, graph)
            MAX_Q = 0
            pos = -1
            for com, Kbv in comm_ngh.items():    
                db = self.DegCom[com]
                delta_Q =  Kbv - self.Degree[node]/(2.*self.m)*db
                if delta_Q > MAX_Q:
                    MAX_Q = delta_Q
                    pos = com
                else :
                    delta_Q = 0
            
            if MAX_Q > 0:
                super().insert_node( node, pos, comm_ngh.get( pos, 0.))
            else:
                com_id = super().generate_random_not_in_list(set(self.membership.values()))
                super().insert_node( node, com_id, comm_ngh.get( com_id, 0.))
                
    def __randomcom ( self ,graph):
        for vsele in  graph.nodes():
            degree = self.Degree[vsele]
            qum={}
            com_befor = self.membership[vsele]
            ngh_com = super().neigh_comm( vsele, graph)
            #print( vsele, com_befor, ngh_com)
            dvc = super().ngh_node(vsele, com_befor, graph)
            devc = self.DegCom[com_befor]              
            for com  in ngh_com :
                if com != com_befor:
                    dvcp = ngh_com[com]
                    devcp = self.DegCom[com]
                    deq = (1/self.m)*(dvcp-dvc)+(degree/(2.*self.m**2.))*(devc-devcp -degree)
                    #deq = dvcp - (devcp )*(degree / (self.m * 2.))
                    #print(deq)
                    if deq > 0: 
                        qum[com] = deq
            
            #print("q in reconstruction", qum)
            if len(qum) > 0:  
                prb = [ self.expon( i, 0.1) for k,i in qum.items() ]
                com_n = super().weighted_choice( list(qum.keys()), prb)
                #print(com_n , ngh_com)
                super().delet_node( vsele, com_befor, ngh_com.get( com_befor, 0.))
                super().insert_node( vsele, com_n, ngh_com.get( com_n, 0.))        
        
    def __merge_community( self, graph):
        modified = True
        while modified :
            modified = False         
            community = set(self.membership.values())
            #print(community, self.membership)
            visited = { i : False for i in community }
            for com1 in community:
                if visited[com1] == False:
                    #print(com1)
                    maxqq = 0
                    poss = -1
                    ngh_comm = super().com_ngh_com(com1, graph)
                    #print(com1 , ngh_comm)
                    for comm ,Kbv in ngh_comm.items():
                        delta_Q =  Kbv/self.m  - ( self.DegCom[com1] * self.DegCom[comm] ) / ( 2. * self.m**2 )
                        #print("delta",delta_Q)
                        if delta_Q > maxqq :
                            maxqq= delta_Q
                            poss = comm
                            
                    if maxqq > 0:
                        modified = True
                        self.merge_com( com1, poss)              
                        self.DegCom[com1] = self.DegCom[com1] + self.DegCom[poss]
                        self.internal[com1] = self.internal[com1] + self.internal[poss]+ ngh_comm.get(poss, 0.)
                        del self.DegCom[poss]
                        del self.internal[poss]
                        visited[poss] = True
        
        return self.membership        
           
    def Run_FMLIG ( self, graph):
        start = time.time()
        soltion = self.GCH( graph)
        #print("solution", self.membership)
        #print(self.DegCom)
        #super().init( graph, soltion)
        soltion = self.FL_move( graph)
        best_solution = copy.copy(soltion)
        #print(" internal ", self.internal# )
        #super().init( graph, soltion)
        Q_best = super().modularity()
        #print("qqqqqqqqqqqq",Q_best)
        T_init = 0.025 * Q_best
        T = T_init
        nb_iter = 0
        status_list = []
        while nb_iter < self.Nb:
            Q1 = super().modularity()
            #print("q1q1qqqqqqq1", Q1)
            incumbent_solution = copy.copy(soltion)
            #print(self.membership)
            #soltion = self.FL_move( weight = 'weight')
            #print(super().modularity( weight ='weight'))
            #print(" befor internal", self.internal)
            soltion, drop_nodes = self.Destruction(graph)
            #print("after destruction", self.membership, drop_nodes , self.internal)
            soltion = self.Reconstruction( drop_nodes, graph)
            #print("after reconstruction",super().modularity())
            #super().init(graph, soltion, weight = 'weight')  
            soltion = self.FL_move(graph)   
            #print("after fast", self.internal)
            #print("after the fast loacl")
            #super().init(graph, soltion)                      
            Q2 = super().modularity()
            #print("Q2", Q2)
            #for data in c_graph.edges( data= True):
            #    print(data)
                
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
            
            soltion = super().renumber()
            super().init( graph, soltion, weight='weight')                 
            nb_iter = nb_iter + 1
        
        #print(status_list)
        end = time.time()
        t = end - start
        #partition = status_list[0].copy()
        #for index in range( 1, len(status_list)-1):
            #for node, community in partition.items():
                #partition[node] = status_list[index][community]
 
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
        print("rb",nb_run)
        communities = Fast_local_Move_IG( Number_iter, Beta, path, graph)
        mod,community,tim = communities.Run_FMLIG(graph) 
        #print(community)
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
        return NMI_max, Q_max, Q_avg, Q_std, time_run     
    elif sys.argv[4] == 'None':
        Q_avg = communities.avg(Q_list)
        Q_max = communities.max(Q_list)
        Q_std = communities.stdev(Q_list)
        time_run = communities.avg(Time_list)
        return Q_max, Q_avg, Q_std, time_run

if __name__ == '__main__':

   
    if  sys.argv[4] != 'None':
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
     
