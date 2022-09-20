#### Constructing tree
import itertools as it
import os
from collections import defaultdict

import networkx as nx
import numpy as np
import pandas as pd
from pathogen_detection.televir_deploy_parameters import (
    Params_Illumina,
    Params_Nanopore,
)

tree = lambda: defaultdict(tree)


def make_tree(lst):
    d = tree()
    for x in lst:
        curr = d
        for item in x:
            curr = curr[item]
    return d


class pipeline_tree:
    def __init__(self):
        pass

    def param_input(self, technology):
        if technology == "nanopore":
            self.mod = Params_Nanopore
        elif technology == "illumina":
            self.mod = Params_Illumina

        self.SOFTWARE_LIST = list(self.mod.SOFTWARE.keys())
        self.params_lookup = {
            "PREPROCESS": self.mod.ARGS_QC,
            "ENRICHMENT": self.mod.ARGS_ENRICH,
            "ASSEMBLY": self.mod.ARGS_ASS,
            "CONTIG_CLASSIFICATION": self.mod.ARGS_ASS,
            "READ_CLASSIFICATION": self.mod.ARGS_CLASS,
            "REMAPPING": self.mod.ARGS_REMAP,
        }

    def fill_dict(self, ix, branch_left):
        """
        traverse dependcy tree
        """
        software_list = self.SOFTWARE_LIST
        if ix == (len(software_list)):
            return branch_left

        if branch_left:
            return {
                node: self.fill_dict(ix, branch) for node, branch in branch_left.items()
            }

        else:
            current = software_list[ix]

        soft_dict = {}
        suffix = ""

        for soft in self.mod.SOFTWARE[current]:
            param_names = []
            param_combs = []
            if soft in self.params_lookup[current].keys():
                params_dict = self.params_lookup[current][soft]

                for i, g in params_dict.items():
                    param_names.append(i + suffix)
                    param_combs.append([(i + suffix, x) for x in g])

                param_combs = list(it.product(*param_combs))
                param_tree = make_tree(param_combs)
                soft_dict[soft] = param_tree
            else:
                soft_dict[soft] = {(f"{soft}_ARGS", "None"): {}}

        return {
            (current, soft): self.fill_dict(ix + 1, g) for soft, g in soft_dict.items()
        }

    def tree_index(self, tree, root, node_index=[], edge_dict=[], leaves=[]):
        """ """
        if len(tree) == 0:
            leaves.append(root)
            return

        subix = {}
        for i, g in tree.items():

            ix = len(node_index)
            node_index.append([ix, i])
            edge_dict.append([root, ix])
            subix[i] = ix

        #
        td = [self.tree_index(g, subix[i]) for i, g in tree.items()]

        return node_index, edge_dict, leaves

    def create_pipe_tree(self):
        """ """
        self.pipeline_tree = self.fill_dict(0, {})
        node_index, edge_dict, leaves = self.tree_index(self.pipeline_tree, "root")
        self.node_index = node_index
        self.edge_dict = edge_dict
        self.leaves = leaves
        self.nodes = [x[0] for x in node_index]

    def node_scores(
        self,
        df,
        params,
        stat=[
            "coverage",
            "depth",
            "depthR",
            "complete",
            "runtime",
            "finished",
            "success",
        ],
    ):
        """ """
        if os.path.isfile("tree_weights.tsv"):
            with open("tree_weights.tsv", "r") as fl:
                lines = fl.readlines()
                lines = [x.strip().split() for x in lines]
                node_weights = {
                    x[0]: [float(y) for y in x[1].split(",")] for x in lines
                }
                for l in list(node_weights.keys()):
                    g = node_weights[l]
                    if l not in ["root", "sink"]:
                        node_weights.pop(l)
                        node_weights[int(l)] = g
                self.node_weights = node_weights
                self.dag_dict = {
                    z: [
                        self.edge_dict[x][1]
                        for x in range(len(self.edge_dict))
                        if self.edge_dict[x][0] == z
                    ]
                    for z in list(set(self.node_weights.keys()))
                }

        else:

            nid = pd.DataFrame(self.node_index, columns=["index", "node"])
            nid.set_index("index")

            def calculate_score(path, nid, stats):

                sub = []
                cols_addressed = []
                for ix in path:
                    if ix == "root":
                        continue
                    tpl = nid.loc[ix]["node"]
                    argname, param = tpl
                    # print(argname, argname in params.param)

                    if argname in params.param.unique():

                        if not len(sub):
                            sub = params.loc[
                                (params.param == argname) & (params.value == param)
                            ]

                        else:
                            sub = pd.merge(
                                sub,
                                params.loc[
                                    (params.param == argname) & (params.value == param)
                                ],
                                on="source",
                            )

                if len(sub) == 0:
                    score = [0] * len(stats)

                else:
                    sources = list(sub.source.unique())
                    summaries_source = df.loc[df.source.isin(sources)]
                    score = [np.median(summaries_source[x]) for x in stats]

                return score

            nodes = [x[0] for x in self.node_index]
            G = nx.DiGraph()  #  G is an empty Graph
            G.add_nodes_from(nodes)
            G.add_edges_from(self.edge_dict)
            G.remove_node("root")

            node_weights = {}
            for node in nodes[1:]:
                path = nx.all_simple_paths(G, source=0, target=node)
                path = list(path)[0]
                score = calculate_score(path, nid, stat)
                node_weights[node] = score

            node_weights[0] = calculate_score([], nid, stat)
            node_weights["root"] = calculate_score([], nid, stat)
            node_weights["sink"] = [0] * len(stat)

            lines = [
                "{}\t{}".format(x, ",".join(np.array(g, dtype=str)))
                for x, g in node_weights.items()
            ]
            with open("tree_weights.tsv", "w") as fl:
                fl.write("\n".join(lines))

            self.node_weights = node_weights
            self.dag_dict = {
                z: [
                    self.edge_dict[x][1]
                    for x in range(len(self.edge_dict))
                    if self.edge_dict[x][0] == z
                ]
                for z in list(set(node_weights.keys()))
            }

    def tree_scores(self, stats=[0, 1, 2]):
        """ """

        visited = {}
        dag_dict = self.dag_dict.copy()

        for nd, lt in dag_dict.items():
            visited[nd] = []

        path_to_leaf_weights = {x: [] for x in self.node_weights.keys()}

        for l in self.leaves:
            dag_dict[l] = ["sink"]

        path_to_leaf_weights = self.node_weights.copy()
        for nd in list(path_to_leaf_weights.keys()):

            g = [path_to_leaf_weights[nd][x] for x in stats]

            g = np.prod(np.array(g, dtype=float), axis=0)
            path_to_leaf_weights[nd] = g

        def tree_weights(node):

            for child in dag_dict[node]:
                if child not in self.leaves:
                    tree_weights(child)

            td = [path_to_leaf_weights[c] for c in dag_dict[node]]
            td = np.array(td)
            td[np.isnan(td)] = 0
            path_to_leaf_weights[node] = np.max(td)

        tree_weights("root")

        self.weights = path_to_leaf_weights

    def simplify_tree(self, links, root, party, nodes_compress=[], edge_keep=[]):
        """ """
        party.append(root)

        if root in self.leaves:
            nodes_compress.append([party[0], tuple(party)])
            return

        subix = {}

        if len(links) > 1:
            nnode = tuple(party)
            ori = nnode[0]
            nodes_compress.append([ori, nnode])

        for i, g in enumerate(links):

            if len(links) != 1:
                edge_keep.append([ori, g])
                party = []

            self.simplify_tree(
                self.dag_dict[g],
                g,
                party,
                nodes_compress=nodes_compress,
                edge_keep=edge_keep,
            )

        return nodes_compress, edge_keep

    def compress_tree(self):
        """ """
        nodes_compress, edges_compress = self.simplify_tree(
            self.dag_dict[0], 0, [], nodes_compress=[], edge_keep=[]
        )
        #
        self.nodes_compress = nodes_compress
        self.edge_compress = edges_compress


import random

import plotly.graph_objects as go


def plotly_network(G, pos):
    edge_x = []
    edge_y = []
    for edge in G.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.append(x0)
        edge_x.append(x1)
        edge_x.append(None)
        edge_y.append(y0)
        edge_y.append(y1)
        edge_y.append(None)

    edge_trace = go.Scatter(
        x=edge_x,
        y=edge_y,
        line=dict(width=0.5, color="#888"),
        hoverinfo="none",
        mode="lines",
    )

    node_x = []
    node_y = []
    for node in G.nodes():
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)

    node_trace = go.Scatter(
        x=node_x,
        y=node_y,
        mode="markers",
        hoverinfo="text",
        marker=dict(
            showscale=True,
            # colorscale options
            #'Greys' | 'YlGnBu' | 'Greens' | 'YlOrRd' | 'Bluered' | 'RdBu' |
            #'Reds' | 'Blues' | 'Picnic' | 'Rainbow' | 'Portland' | 'Jet' |
            #'Hot' | 'Blackbody' | 'Earth' | 'Electric' | 'Viridis' |
            colorscale="YlGnBu",
            reversescale=True,
            color=[],
            size=10,
            colorbar=dict(
                thickness=15,
                title="combine performace statistic",
                xanchor="left",
                titleside="right",
            ),
            line_width=2,
        ),
    )
    return node_trace, edge_trace


def plotly_node_text(G, node_trace, annotations, path_to_leaf_weights):
    node_adjacencies = []
    node_text = []
    for node_ix, adjacencies in enumerate(G.adjacency()):
        # print(node_ix, adjacencies)
        node = list(G.nodes)[node_ix]
        node_adjacencies.append(path_to_leaf_weights[node])
        if node in annotations.keys():
            node_text.append(annotations[node])

    node_trace.marker.color = node_adjacencies
    node_trace.text = node_text
    return node_trace


class tree_plot:
    def __init__(self):
        pass

    def graph(self, nodes, node_dict, edge_list, node_weights):
        """ """
        annotations = {
            x[0]: [node_dict[y] for y in x[1]]
            for x in nodes
            if x[0] not in ["root", "sink"]
        }
        annotations = {
            x: "<br>".join(["{}: {}".format(*y[1]) for y in g])
            for x, g in annotations.items()
        }

        nodes = [x[0] for x in nodes]
        self.annotations = annotations
        self.nodes = nodes
        self.edge_list = edge_list
        self.weights = node_weights

    def hierarchy_pos(
        self, G, root=None, width=1.0, vert_gap=0.2, vert_loc=0, xcenter=0.5
    ):

        """
        From Joel's answer at https://stackoverflow.com/a/29597209/2966723.
        Licensed under Creative Commons Attribution-Share Alike

        If the graph is a tree this will return the positions to plot this in a
        hierarchical layout.

        G: the graph (must be a tree)

        root: the root node of current branch
        - if the tree is directed and this is not given,
          the root will be found and used
        - if the tree is directed and this is given, then
          the positions will be just for the descendants of this node.
        - if the tree is undirected and not given,
          then a random choice will be used.

        width: horizontal space allocated for this branch - avoids overlap with other branches

        vert_gap: gap between levels of hierarchy

        vert_loc: vertical location of root

        xcenter: horizontal location of root
        """
        if not nx.is_tree(G):
            raise TypeError("cannot use hierarchy_pos on a graph that is not a tree")

        if root is None:
            if isinstance(G, nx.DiGraph):
                root = next(
                    iter(nx.topological_sort(G))
                )  # allows back compatibility with nx version 1.11
            else:
                root = random.choice(list(G.nodes))

        def _hierarchy_pos(
            G,
            root,
            width=1.0,
            vert_gap=0.2,
            vert_loc=0,
            xcenter=0.5,
            pos=None,
            parent=None,
        ):
            """
            see hierarchy_pos docstring for most arguments

            pos: a dict saying where all nodes go if they have been assigned
            parent: parent of this branch. - only affects it if non-directed

            """

            if pos is None:
                pos = {root: (xcenter, vert_loc)}
            else:
                pos[root] = (xcenter, vert_loc)
            children = list(G.neighbors(root))
            if not isinstance(G, nx.DiGraph) and parent is not None:
                children.remove(parent)
            if len(children) != 0:

                dx = width / len(children)
                nextx = xcenter - width / 2 - dx / 2
                for child in children:
                    nextx += dx
                    pos = _hierarchy_pos(
                        G,
                        child,
                        width=dx,
                        vert_gap=vert_gap,
                        vert_loc=vert_loc - vert_gap,
                        xcenter=nextx,
                        pos=pos,
                        parent=root,
                    )
            return pos

        return _hierarchy_pos(
            G, root, width=width, vert_gap=vert_gap, vert_loc=vert_loc, xcenter=xcenter
        )

    def generate_graph(self):

        Gsmall = nx.DiGraph()  #  G is an empty Graph
        Gsmall.add_nodes_from(self.nodes)
        Gsmall.add_edges_from(self.edge_list)

        values_short = [self.weights.get(x) for x in Gsmall.nodes]

        pos_short = self.hierarchy_pos(
            Gsmall,
            0,
        )

        self.graph = Gsmall
        self.pos = pos_short
        self.values = values_short

    def graph_plot(self):
        node_trace, edge_trace = plotly_network(self.graph, self.pos)
        node_trace = plotly_node_text(
            self.graph, node_trace, self.annotations, self.weights
        )
        fig = go.Figure(
            data=[edge_trace, node_trace],
            layout=go.Layout(
                title="<br>Pipeline graph",
                titlefont_size=16,
                showlegend=False,
                hovermode="closest",
                margin=dict(b=20, l=5, r=5, t=40),
                annotations=[
                    dict(
                        text="Python code: <a href='https://plotly.com/ipython-notebooks/network-graphs/'> https://plotly.com/ipython-notebooks/network-graphs/</a>",
                        showarrow=False,
                        xref="paper",
                        yref="paper",
                        x=0.005,
                        y=-0.002,
                    )
                ],
                xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            ),
        )
        return fig
