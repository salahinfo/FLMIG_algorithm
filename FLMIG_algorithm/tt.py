import networkx as nx
import matplotlib.pyplot as plt
#
# Create a graph
G = nx.read_edgelist('/home/yacine/Desktop/real_network/test.txt',nodetype= int)
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
