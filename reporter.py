
import argparse
import tempfile
from os.path import join
import os
import asterid as ad
from random import sample
import pandas as pd
import dendropy
from dendropy.calculate.treecompare \
    import false_positives_and_negatives
import sys
 
def compare_trees(tr1, tr2):
    # Find leaf labels that are in both trees
    lb1 = set([l.taxon.label for l in tr1.leaf_nodes()])
    lb2 = set([l.taxon.label for l in tr2.leaf_nodes()])
    com = lb1.intersection(lb2)
 
    # Restrict trees to shared leaf set
    if com != lb1 or com != lb2:
        com = list(com)
        tns = dendropy.TaxonNamespace(com)
 
        tr1.retain_taxa_with_labels(com)
        tr1.migrate_taxon_namespace(tns)
 
        tr2.retain_taxa_with_labels(com)
        tr2.migrate_taxon_namespace(tns)
    com = list(com)
 
    # Update tree bipartitions
    tr1.update_bipartitions()
    tr2.update_bipartitions()
 
    # Compute number of leaves and number of internal edges
    nl = len(com)
    ei1 = len(tr1.internal_edges(exclude_seed_edge=True))
    ei2 = len(tr2.internal_edges(exclude_seed_edge=True))
 
    # Compute number of false positives and false negatives
    [fp, fn] = false_positives_and_negatives(tr1, tr2)
 
    # Compute symmetric difference rate
    sd = float(fp + fn) / (ei1 + ei2)
 
    # Compute Robinson-Foulds error rate
    rf = float(fp + fn) / (2 * nl - 6)
 
    return(nl, ei1, ei2, fp, fn, sd, rf)

def comparetreepath(tr1, tr2):
    tax = dendropy.TaxonNamespace()
    tr1 = dendropy.Tree.get(path=tr1,
                            schema='newick',
                            rooting='force-unrooted',
                            taxon_namespace=tax)
 
    tr2 = dendropy.Tree.get(path=tr2,
                            schema='newick',
                            rooting='force-unrooted',
                            taxon_namespace=tax)
 
    # Unroot trees
    tr1.collapse_basal_bifurcation(set_as_unrooted_tree=True)
    tr2.collapse_basal_bifurcation(set_as_unrooted_tree=True)
 
    # Compute RF distance
    [nl, ei1, ei2, fp, fn, sd, rf] = compare_trees(tr1, tr2)
    return [nl, ei1, ei2, fp, fn, sd, rf]

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str,
                        help="Input tree list file", required=True)
# parser.add_argument("-o", "--output", type=str, required=True)
args = parser.parse_args()
df = pd.read_csv(args.input)
nrf = []

for i, r in df.iterrows():
    nrf.append(comparetreepath(r.streepath, r.inputpath)[-1])
df["rf"] = nrf
dfa = df.set_index(['condition', 'k', 'method'])[['rf']]
with pd.option_context("display.max_rows", 1000):
    mdfa = dfa.groupby(level=list(range(3))).mean()
    print(mdfa)
    print(mdfa.to_latex(escape=False,float_format="{:0.3f}".format))
