from GraphTools import GraphTolls, Read_Graph
from FLMIG import Fast_local_Move_IG



path = '/home/yacine/Desktop/real_network/karate.gml'
g = Read_Graph(path)
#print(g.nodes())


d = Fast_local_Move_IG(10, 0.4, path, g)
sol = d.GCH(g)
d.modifie_status( g, weight='weight')
#print(d.membership)
sol = d.localmove(g)
#print(d.membership)
n_mod = d.modularity()
print("befor", n_mod)
p_list = []
mod_graph = g.copy()
p = d.renumber()
p_list.append(p)
mod_graph = d.induced_graph( p, mod_graph, weight = 'weight')
d.modifie_status( mod_graph, weight = 'weight')
Q_val = n_mod
while True :
    #print( mod_graph.nodes())
    solution = d.localmove(mod_graph)
    n_mod = d.modularity()
    print(n_mod)
    if n_mod - Q_val < 0.000000001:
        break
    
    #print(d.membership)

    Q_val = n_mod
    p = d.renumber()
    p_list.append(p)
    mod_graph = d.induced_graph( p, mod_graph, weight = 'weight')
    d.modifie_status( mod_graph, weight='weight')
    

P = d.best_sol(p_list, len(p_list) - 1)
d.init(g, P, weight='weight')








