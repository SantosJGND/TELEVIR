#### Constructing tree
import itertools as it
import os
import sys
from collections import defaultdict
from wsgiref.simple_server import software_version

import networkx as nx
import numpy as np
import pandas as pd
from pathogen_identification.televir_deploy_parameters import (
    Params_Illumina,
    Params_Nanopore,
)
from settings.constants_settings import ConstantsSettings as CS

tree = lambda: defaultdict(tree)


def make_tree(lst):
    d = tree()
    for x in lst:
        curr = d
        for item in x:
            curr = curr[item]
    return d


class pipeline_tree:

    edge_dict: list
    leaves: list
    nodes: list
    node_index: pd.DataFrame

    def __init__(self):
        pass

    def param_input(self, technology):
        if technology == CS.TECHNOLOGY_minion:
            self.mod = Params_Nanopore
        elif technology == CS.TECHNOLOGY_illumina:
            self.mod = Params_Illumina

        self.SOFTWARE_LIST = list(self.mod.SOFTWARE.keys())
        self.params_lookup = {
            CS.PIPELINE_NAME_read_quality_analysis: self.mod.ARGS_QC,
            CS.PIPELINE_NAME_viral_enrichment: self.mod.ARGS_ENRICH,
            CS.PIPELINE_NAME_assembly: self.mod.ARGS_ASS,
            CS.PIPELINE_NAME_contig_classification: self.mod.ARGS_CLASS,
            CS.PIPELINE_NAME_read_classification: self.mod.ARGS_CLASS,
            CS.PIPELINE_NAME_remapping: self.mod.ARGS_REMAP,
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
                    param_combs.append([(i + suffix, x, "param") for x in g])

                param_combs = list(it.product(*param_combs))
                param_tree = make_tree(param_combs)
                soft_dict[soft] = param_tree
            else:
                soft_dict[soft] = {(f"{soft}_ARGS", "None", "param"): {}}

        return {
            (current, soft, "module"): self.fill_dict(ix + 1, g)
            for soft, g in soft_dict.items()
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
        td = [
            self.tree_index(
                g, subix[i], node_index=node_index, edge_dict=edge_dict, leaves=leaves
            )
            for i, g in tree.items()
        ]

        return node_index, edge_dict, leaves

    def create_pipe_tree(self):
        """ """
        self.pipeline_tree = self.fill_dict(0, {})
        node_index, edge_dict, leaves = self.tree_index(
            self.pipeline_tree, "root", node_index=[], edge_dict=[], leaves=[]
        )
        self.node_index = pd.DataFrame(node_index, columns=["index", "node"])
        self.node_index.set_index("index", inplace=True)

        self.edge_dict = edge_dict
        self.leaves = leaves
        self.nodes = [x[0] for x in node_index]

        self.dag_dict = {
            z: [
                self.edge_dict[x][1]
                for x in range(len(self.edge_dict))
                if self.edge_dict[x][0] == z
            ]
            for z in list(set(self.node_index.index))
        }

    def generate_graph_no_root(self):
        """generate networkx graph"""
        G = nx.DiGraph()  #  G is an empty Graph
        G.add_nodes_from(self.nodes)
        G.add_edges_from(self.edge_dict)
        G.remove_node("root")

        return G

    @staticmethod
    def path_collect_runs(path: list, nid: pd.DataFrame, params: pd.DataFrame):
        """
        follow path through tree and collect runids
        params: requires runid, module, software, value
        nid: requires index, node
        """
        sub = []

        dead = False

        module = None
        software = None

        for i, ix in enumerate(path):
            if ix == "root" or dead:
                continue
            tpl = nid.loc[ix]["node"]
            argname, param, col = tpl
            if param == "None":
                continue

            if col == "module":
                module = argname
                software = param
                new_subs = params[params["module"] == module][
                    params["software"] == software
                ].runid.unique()
            else:
                new_subs = params[params["module"] == module][
                    params["software"] == software
                ][params["value"] == str(param)].runid.unique()

            if not len(sub):
                sub = new_subs
            else:

                subs_keep = [x for x in sub if x in new_subs]
                if len(subs_keep) == 0:
                    dead = True
                    print(tpl)
                    print(path)
                    print([nid.loc[x]["node"] for x in path])
                    print(sub)
                    print("#####")
                sub = subs_keep

        return sub

    def get_node_path(self, node):
        """return path for a given node"""
        G = self.generate_graph_no_root()
        path = nx.all_simple_paths(G, source=0, target=node)
        path = list(path)[0]
        return path

    @staticmethod
    def get_node_path_from_graph(node, graph):
        """return path for a given node"""
        path = nx.all_simple_paths(graph, source=0, target=node)
        path = list(path)[0]
        return path

    def get_leaf_paths(self) -> dict:
        """ "return path list for each node"""
        node_paths = {}
        G = self.generate_graph_no_root()
        for node in self.leaves:
            path = self.get_node_path_from_graph(node, G)
            node_paths[node] = path

        return node_paths

    def collect_leaf_runs(self, params: pd.DataFrame) -> pd.DataFrame:
        """collect runids for each leaf node"""

        nid = self.node_index

        print("collecting leaf runs")
        params["value"] = params.value.astype(str)

        node_runs = {x: [] for x in self.nodes}
        leaf_paths = self.get_leaf_paths()
        for node, path in leaf_paths.items():
            node_runs[node] = self.path_collect_runs(path, nid, params)

        node_runs_df = pd.DataFrame(node_runs.items(), columns=["node", "runids"])
        node_runs_df = node_runs_df.explode("runids")
        node_runs_df = node_runs_df.dropna()

        return node_runs_df

    @staticmethod
    def calculate_node_score(stats: list, df: pd.DataFrame, sub: list):
        """calculate node scores
        stats: list of stats to calculate
        df: dataframe of all runids and stats
        sub: list of runids to calculate stats for
        """

        representation = 0
        score = []

        if len(sub) == 0:
            score = [0] * len(stats)

        else:

            sources = sub  # list(sub.source.unique())
            summaries_source = df[(df.runid.isin(sources))]

            if len(summaries_source) == 0:
                score = [0] * len(stats)

            else:
                score = [np.mean(summaries_source[x]) for x in stats]

            samples_present = list(summaries_source.source.unique())
            representation = len(samples_present) / df.source.nunique()

        score.append(representation)
        return score

    def calculate_node_scores_return_df(
        self,
        node_runs_df: pd.DataFrame,
        df,
        stat=[
            "coverage",
            "depth",
            "depthR",
            "complete",
            "ahelp",
            "runtime",
            "finished",
            "success_prop",
            "ref_proportion",
        ],
    ) -> pd.DataFrame:

        print("calculating node scores")
        print(df.shape)
        node_scores = {x: [0] * (len(stat) + 1) for x in self.nodes}
        for node, runids in node_runs_df.groupby("node"):
            node_scores[node] = self.calculate_node_score(
                stat, df, sub=runids.runids.to_list()
            )

        node_scores[0] = [0] * (len(stat) + 1)
        node_scores["root"] = [0] * (len(stat) + 1)
        node_scores["sink"] = [0] * (len(stat) + 1)

        nodes_df = pd.DataFrame(
            [[x] + g for x, g in node_scores.items()],
            columns=["node"] + stat + ["branch_coverage"],
        )

        return nodes_df

    def calculate_node_scores_with_source_subset_draw(
        self, node_runs_df: pd.DataFrame, df=pd.DataFrame, cv=10, subset=0.6
    ):
        """calculate node scores with source subset draw
        node_runs_df: dataframe of node and runids
        cv: number of cross validation folds
        subset: proportion of sources to use for each fold
        """

        print("calculating node scores with source subset draw")
        iter_scores = []
        df_sources = df.source.unique()

        for i in range(cv):
            sample_sources = np.random.choice(
                df_sources, int(len(df_sources) * subset), replace=False
            )
            sample_df = df[df.source.isin(sample_sources)]
            print(df.shape, sample_df.shape)
            nodes_iter = self.calculate_node_scores_return_df(node_runs_df, sample_df)
            iter_scores.append(nodes_iter)

        iter_scores = pd.concat(iter_scores)
        iter_scores = iter_scores.groupby("node").mean().reset_index()

        self.node_weights = iter_scores.set_index("node")
        self.dag_dict = {
            z: [
                self.edge_dict[x][1]
                for x in range(len(self.edge_dict))
                if self.edge_dict[x][0] == z
            ]
            for z in list(set(self.node_weights.index))
        }

        return iter_scores

    def calculate_node_scores(
        self,
        node_runs_df: pd.DataFrame,
        df,
        stat=[
            "coverage",
            "depth",
            "depthR",
            "complete",
            "ahelp",
            "runtime",
            "finished",
            "success_prop",
            "ref_proportion",
        ],
    ) -> pd.DataFrame:

        print("calculating node scores")
        node_scores = {x: [0] * (len(stat) + 1) for x in self.nodes}
        for node, runids in node_runs_df.groupby("node"):
            print(node)
            len(runids)
            node_scores[node] = self.calculate_node_score(
                stat, df, sub=runids.runids.to_list()
            )

        node_scores[0] = [0] * (len(stat) + 1)
        node_scores["root"] = [0] * (len(stat) + 1)
        node_scores["sink"] = [0] * (len(stat) + 1)

        nodes_df = pd.DataFrame(
            [[x] + g for x, g in node_scores.items()],
            columns=["node"] + stat + ["branch_coverage"],
        )

        self.node_weights = nodes_df.set_index("node")
        self.dag_dict = {
            z: [
                self.edge_dict[x][1]
                for x in range(len(self.edge_dict))
                if self.edge_dict[x][0] == z
            ]
            for z in list(set(self.node_weights.index))
        }

    def node_scores(
        self,
        df,
        params,
        stat=[
            "coverage",
            "depth",
            "depthR",
            "complete",
            "ahelp",
            "runtime",
            "finished",
            "success_prop",
            "ref_proportion",
        ],
    ):
        """ """
        if os.path.isfile("tree_weights.tsv"):

            node_df = pd.read_csv("tree_weights.tsv", sep="\t")
            node_df = node_df[~node_df.node.isin(["root", "sink"])]
            node_df["node"] = node_df["node"].astype(int)

            self.node_weights = node_df.set_index("node")

            self.dag_dict = {
                int(z): [
                    self.edge_dict[x][1]
                    for x in range(len(self.edge_dict))
                    if self.edge_dict[x][0] == (int(z))
                ]
                for z in list(set(self.node_weights.index))
            }

        else:

            nid = self.node_index
            print("calculating node scores")
            params["value"] = params.value.astype(str)

            lost = {}

            def calculate_score(path, nid, stats, df, params):

                sub = []

                dead = False

                module = None
                software = None

                for i, ix in enumerate(path):
                    if ix == "root" or dead:
                        continue
                    tpl = nid.loc[ix]["node"]
                    argname, param, col = tpl
                    if param == "None":
                        continue

                    if col == "module":
                        module = argname
                        software = param
                        new_subs = params[params["module"] == module][
                            params["software"] == software
                        ].runid.unique()
                    else:
                        new_subs = params[params["module"] == module][
                            params["software"] == software
                        ][params["value"] == str(param)].runid.unique()

                    if not len(sub):
                        sub = new_subs
                    else:

                        subs_keep = [x for x in sub if x in new_subs]
                        if len(subs_keep) == 0:
                            dead = True
                            print(tpl)
                            print(path)
                            print([nid.loc[x]["node"] for x in path])
                            print(sub)
                            print("#####")

                        sub = subs_keep

                representation = 0

                if len(sub) == 0:
                    score = [0] * len(stats)

                else:

                    sources = sub  # list(sub.source.unique())
                    summaries_source = df[(df.runid.isin(sources))]

                    if len(summaries_source) == 0:
                        score = [0] * len(stats)

                    else:
                        score = [np.median(summaries_source[x]) for x in stats]

                    samples_present = list(summaries_source.source.unique())
                    representation = len(samples_present) / df.source.nunique()

                score.append(representation)
                return score

            G = nx.DiGraph()  #  G is an empty Graph
            G.add_nodes_from(self.nodes)
            G.add_edges_from(self.edge_dict)
            G.remove_node("root")

            node_weights = {x: [0] * (len(stat) + 1) for x in self.nodes}
            for node in self.leaves:
                path = nx.all_simple_paths(G, source=0, target=node)
                path = list(path)[0]
                score = calculate_score(path, nid, stat, df, params)

                node_weights[node] = score

            node_weights[0] = calculate_score([], nid, stat, df, params)
            node_weights["root"] = calculate_score([], nid, stat, df, params)
            node_weights["sink"] = [0] * (len(stat) + 1)

            lines = [
                "{}\t{}".format(x, ",".join(np.array(g, dtype=str)))
                for x, g in node_weights.items()
            ]

            nodes_df = pd.DataFrame(
                [[x] + g for x, g in node_weights.items()],
                columns=["node"] + stat + ["branch_coverage"],
            )

            nodes_df.to_csv("tree_weights.tsv", sep="\t", index=False)

            self.node_weights = nodes_df.set_index("node")
            self.dag_dict = {
                z: [
                    self.edge_dict[x][1]
                    for x in range(len(self.edge_dict))
                    if self.edge_dict[x][0] == z
                ]
                for z in list(set(self.node_weights.index))
            }

            ##

    def tree_scores(self, stats=["coverage", "depth"], node_color=None):
        """ """

        if isinstance(node_color, int):
            node_color = [node_color]

        visited = {}
        dag_dict = self.dag_dict.copy()

        for nd, lt in dag_dict.items():
            visited[nd] = []

        for l in self.leaves:
            dag_dict[l] = ["sink"]

        path_to_leaf_weights = self.node_weights.copy()
        for nd in path_to_leaf_weights.index:

            g = [path_to_leaf_weights.loc[nd][x] for x in stats]
            g = np.array(g, dtype=float)
            g = g[~np.isnan(g)]

            g = np.prod(g, axis=0)

            if node_color is not None:
                if nd in node_color:
                    g = g * 10
                else:
                    g = 0

            path_to_leaf_weights.loc[nd] = g

        def tree_weights(node):

            for child in dag_dict[node]:
                if child not in self.leaves:
                    tree_weights(child)

            td = [path_to_leaf_weights.loc[c] for c in dag_dict[node]]
            td = np.array(td)
            td[np.isnan(td)] = 0

            path_to_leaf_weights.loc[node] = np.max(td)

        tree_weights(0)

        self.weights = path_to_leaf_weights

    def simplify_tree(self, links, root, party: list, nodes_compress=[], edge_keep=[]):
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

    def split_modules(self):

        if self.nodes_compress is None:
            self.compress_tree()

        nodes_compress = self.nodes_compress
        edge_compress = self.edge_compress

        new_nodes = []
        new_edges = []

        for node in nodes_compress:
            internal_nodes = []
            internal_edges = []

            internal_splits = [0]

            for ix, internal_node in enumerate(node[1]):
                internal_name = self.node_index.loc[internal_node].node
                is_module = internal_name[2] == "module"

                if is_module and ix != 0:
                    internal_splits.append(ix)

            if len(internal_splits) > 1:
                internal_splits.append(len(node[1]))
                for ix, split in enumerate(internal_splits[:-1]):
                    internal_nodes.append(
                        [node[1][split], (node[1][split : internal_splits[ix + 1]])]
                    )

                for ix, split in enumerate(internal_nodes[:-1]):
                    internal_edges.append([split[0], internal_nodes[ix + 1][0]])

                new_nodes.extend(internal_nodes)
                edge_compress.extend(internal_edges)
            else:
                new_nodes.append(node)

        self.nodes_compress = new_nodes

        self.compress_dag_dict = {
            z: [
                self.edge_compress[x][1]
                for x in range(len(self.edge_compress))
                if self.edge_compress[x][0] == z
            ]
            for z in list(set([x[0] for x in self.nodes_compress]))
        }

    def same_module_children(self, node, party, branches=[]):
        """ """
        children = self.compress_dag_dict[node]

        if len(children) == 0:
            return branches

        for child in children:
            child_name = self.node_index.loc[child].node
            new_party = party.copy()

            if child_name[2] != "module":
                new_party.append(child)
                self.same_module_children(child, new_party, branches=branches)
            else:
                if len(new_party) > 0:
                    branches.append(new_party)
                new_party = []

        return branches

    def get_module_tree(self):

        if self.nodes_compress is None:
            self.compress_tree()

            nodes_compress = self.nodes_compress
            edge_compress = self.edge_compress

            self.split_modules()

        nodes_df = pd.DataFrame(self.nodes_compress, columns=["node", "branch"])
        edge_df = pd.DataFrame(self.edge_compress, columns=["parent", "child"])

        def edit_branches(branches, nodes_df, edge_df, parent_node):
            """ """
            new_edges = []
            new_nodes = []
            for branch in branches:
                if len(branch) == 1:
                    continue

                nodes_df = nodes_df[~nodes_df.node.isin(branch)]
                edge_df = edge_df[~edge_df.parent.isin(branch[:-1])]
                edge_df = edge_df[~edge_df.child.isin(branch)]

                new_edges.append([parent_node, branch[-1]])
                new_nodes.append([branch[-1], tuple(branch)])

            return new_nodes, new_edges, nodes_df, edge_df

        def merge_new_branches(new_nodes, new_edges, nodes_df, edge_df):

            new_nodes = pd.DataFrame(new_nodes, columns=["node", "branch"])
            new_edges = pd.DataFrame(new_edges, columns=["parent", "child"])

            nodes_df = pd.concat([nodes_df, new_nodes], ignore_index=True)
            edge_df = pd.concat([edge_df, new_edges], ignore_index=True)

            return nodes_df, edge_df

        for node in self.nodes_compress:
            node_name = self.node_index.loc[node[0]].node
            if not node_name[2] == "module":
                continue

            parent_node = edge_df[edge_df.child == node[0]].parent.values

            if len(parent_node) == 0:
                parent_node = [0]
            parent_node = parent_node[0]

            same_module_branches = self.same_module_children(
                node[0], [node[0]], branches=[]
            )
            new_nodes, new_edges, nodes_df, edge_df = edit_branches(
                same_module_branches, nodes_df, edge_df, parent_node
            )
            nodes_df, edge_df = merge_new_branches(
                new_nodes, new_edges, nodes_df, edge_df
            )

        self.nodes_compress = nodes_df.to_numpy().tolist()
        self.edge_compress = edge_df.to_numpy().tolist()


import random

# import plotly.graph_objects as go


def plotly_network(G, pos, node=""):
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

    if node:
        node_trace = go.Scatter(
            x=node_x,
            y=node_y,
            mode="markers",
            hoverinfo="text",
            marker=dict(
                showscale=True,
                colorscale="YlGnBu",
                reversescale=True,
                color=[],
                size=10,
                colorbar=dict(
                    thickness=15,
                    title="Node Connections",
                    xanchor="left",
                    titleside="right",
                ),
                line_width=2,
            ),
        )

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
        if isinstance(node_weights, pd.DataFrame):
            node_weights = {x: max(node_weights.loc[x]) for x in node_weights.index}

        annotations = {
            x[0]: [node_dict[y] for y in x[1]]
            for x in nodes
            if x[0] not in ["root", "sink"]
        }
        annotations = {
            x: "<br>".join(["{}: {}".format(*y[1]) for y in g])
            for x, g in annotations.items()
        }
        annotations = {
            x: "<br>".join([f"node: {x}", f"w: {round(node_weights[x], 5)}", g])
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

        values_short = [self.weights[x] for x in Gsmall.nodes]

        pos_short = self.hierarchy_pos(
            Gsmall,
            0,
        )

        self.graph = Gsmall
        self.pos = pos_short
        self.values = values_short

    def graph_plot(self, node=""):
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
