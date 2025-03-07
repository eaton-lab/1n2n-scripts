#!/usr/bin/env python

"""Write a table with triplets for every orthogroup.

This is used to analyze trees produced by orthofinder and then reconciled
by generax to produce rooted gene trees. In our analysis each tip in each
gene tree represents a gene. We want to extract all unique triplets from
each gene tree representing an ingroup-sister-outgroup trio. Some trees
contain many and other trees contain none.

Usage
-----
# call on a single newick file
$ python og_tree_to_table.py NWK -i C_frag -s C_prot -o G_dry > table.tsv

# call on many newick files at once to create one big table
$ python og_tree_to_table.py DIR/*.nwk -i C_frag -s C_prot -o G_dry > table.tsv

"""

from pathlib import Path
from argparse import ArgumentParser
import toytree
import pandas as pd


def get_combinatorial_triplets(ogid: str, tree: toytree.ToyTree, ingroup_prefix: str, sister_prefix: str, outg_prefix: str, relabel: bool) -> pd.DataFrame | None:
    """Return minimum spanning taxon triplets from an orthogroup tree

    rep    OG   i-clade   mean-group                 i            s               o
      0   001         1            0      C_fragilis_a  C_protusa_a  C_dryopteris_a
      1   001         1            1      C_fragilis_a  C_protusa_a  C_dryopteris_b
      2   001         2            0      C_fragilis_b  C_protusa_b  C_dryopteris_c
      3   001         3            0      C_fragilis_c  C_protusa_c  C_dryopteris_d
      4   001         3            1      C_fragilis_c  C_protusa_d  C_dryopteris_d

    Relabel tips to be like this:
      C_fragilis_CfragiEVm038212t1  ->  Cfragilis_CfragiEVm038212t1
    """
    data = []
    for node in tree[:tree.ntips]:

        # only focus on ingroup nodes
        if not node.name.startswith(ingroup_prefix):
            continue

        # store results for this node
        sisters = []
        outgroups = []

        # get a sister
        trace = node
        while trace.up:
            candidates = trace.get_sisters()
            for cnode in candidates:
                for tip in cnode.get_leaf_names():
                    if tip.startswith("C_prot"):
                        sisters.append(tip)
            trace = trace.up
            if sisters:
                break

        # get outgroup
        while trace.up:
            candidates = trace.get_sisters()
            for cnode in candidates:
                for tip in cnode.get_leaf_names():
                    if tip.startswith("G_dry"):
                        outgroups.append(tip)
            trace = trace.up
            if outgroups:
                break

        # if valid triplet
        if sisters and outgroups:
            for sister in sisters:
                for outg in outgroups:
                    # relabel
                    if relabel:
                        names = [i.replace("_", "", 1) for i in [node.name, sister, outg]]
                    else:
                        names = [node.name, sister, outg]
                    data.append([ogid] + names)

    # If no valid trios in the tree return None
    if not data:
        return None

    # convert to a dataframe
    data = pd.DataFrame(data, columns=["OG", "ingroup", "sister", "outgroup"])

    # add column grouping those that share the same sister+outgroup
    data['same_sister'] = data.groupby(['sister']).ngroup()
    data['same_sister_and_og'] = data.groupby(['sister', 'outgroup']).ngroup()
    return data.sort_values(by=["sister", "outgroup", "ingroup"]).reset_index(drop=True)


def get_parser():
    parser = ArgumentParser("og_tree_to_table")
    parser.add_argument("newicks", type=Path, help="path to one or more newick files (regex allowed)")
    parser.add_argument("-i", type=str, required=True, help="prefix for ingroup gene names")
    parser.add_argument("-s", type=str, required=True, help="prefix for sister gene names")
    parser.add_argument("-o", type=str, required=True, help="prefix for outgroup gene names")
    parser.add_argument("--relabel", action="store_true", help="remove first underscore in gene names")
    return parser


def parse_newicks_as_one_or_more_paths(path: Path) -> list[Path]:
    if path.is_file():
        return [path]
    if any(char in path for char in "*?[]"):
        matched_files = list(path.parent.glob(path.name))
        return matched_files if matched_files else []
    return []


def main():
    # get command line arguments
    parser = get_parser()
    args = parser.parse_args()

    # get a dataframe for every newick file
    dfs = []
    for newick in parse_newicks_as_one_or_more_paths(args.newicks):

        ogid = newick.stem
        tree = toytree.tree(newick)
        df = get_combinatorial_triplets(ogid, tree, args.i, args.s, args.o, args.relabel)
        dfs.append(df)

    return pd.concat(dfs).to_csv(sep="\t")


if __name__ == "__main__":
    main()
