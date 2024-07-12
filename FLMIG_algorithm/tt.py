import networkx as nx
import matplotlib.pyplot as plt
#
# Create a graph
G = nx.read_edgelist('/home/yacine/Desktop/test.txt',nodetype= int)
G = nx.Graph()

# Add nodes and specify colors
nodes = ['A', 'B', 'C', 'D']
node_colors = ['red', 'blue', 'green', 'purple']

for node, color in zip(nodes, node_colors):
    G.add_node(node)
    G.nodes[node]['color'] = color

# Add edges (you can customize this part as needed)

# Compute a layout
pos = nx.spring_layout(G)

# Extract node colors from node attributes
node_colors = [G.nodes[node]['color'] for node in G.nodes]

# Draw the graph with different node colors
nx.draw(G, pos, with_labels=True, node_size=700, node_color=node_colors)

plt.show()
plt.savefig("/home/yacine/Desktop/image eps/a.eps",format='eps')

def disconnected_compen( self, p, graph, com_id, weight= 'weight'):
        
        ret = nx.Graph()
        nodecom = [ node for node in p if p[node] == com_id ]
        
        if len(nodecom) <= 1:
            
            return True 
        
        ret.add_nodes_from(nodecom)
        ret = graph.subgraph(nodecom)
        if nx.is_connected(ret):
            
            return True
        
        else :

            return False  
        
#def Singelton_community( self, membership, graph , drop):
        singelton = []
        #membership = self.renumber(membership)
        print("BEFOR",membership, drop)
        membership = super().init( graph, membership, weight='weight') 
        for bcom, lst in  drop.items():
            drop_node = lst
            #print("drop node", drop_node)
            for al in drop_node:
                singelton.append(al)
                com_id = membership[al]
                wgh = super().neigh_comm(membership, al, graph)    
                membership = super().delet_node( membership, al, com_id, wgh.get( com_id, 0.))
                if al not in membership.values():
                    membership = super().insert_node( membership, al, al, wgh.get( al, 0.))
                    
                else:
                    com_id = super().generate_random_not_in_list(set(membership.values()))
                    membership = super().insert_node( membership, al, com_id, wgh.get( com_id, 0.))
                        
        print("after",membership)
        return  membership , singelton        