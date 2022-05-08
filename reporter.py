import argparse
import tempfile
from os.path import join, isfile
import os
# import asterid as ad

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon
import numpy as np
from random import sample
import pandas as pd
import dendropy
from dendropy.calculate.treecompare \
    import false_positives_and_negatives
import sys
import json
from subprocess import Popen, PIPE
import re
import shlex
localization = {}

if isfile("localization.json"):
    with open("localization.json") as fh: localization = json.load(fh)


ASTRALP = "/home/baqiaol2/scratch/atsume/methods/ASTRAL/astral.5.7.7.jar"
def compute_mqsst(gtree, stree):
    process = Popen(shlex.split(f"java -jar {ASTRALP} -i {gtree} -q {stree}"), stdout=PIPE, stderr=PIPE)
    output, err = process.communicate()
    m = re.search("Final normalized quartet score is: (\d\.\d*)", err.decode('ascii'))
    return float(m.group(1))

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
    return(nl, ei1, ei2, fp / ei2, fn / ei1, sd, rf)

def comparetreepath(tr1, tr2):
    tax = dendropy.TaxonNamespace()
    tr1 = dendropy.Tree.get(path=tr1,
                            schema='newick',
                            rooting='force-unrooted',
                            taxon_namespace=tax)
    import treeswift as ts
    with open(tr2) as fh: tr2ts = ts.read_tree_newick(fh.read())
    for n in tr2ts.traverse_postorder(True, True):
        if n.label and "_" in n.label:
            n.label = n.label.split("_")[0]
    tr2 = dendropy.Tree.get(data=tr2ts.newick(),
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
args = parser.parse_args()
df = pd.read_csv(args.input)
nrf = []
nfp = []
nfn = []
time = []
mem = []
mqsst = []
for i, r in df.iterrows():
    print(f"comparing {r.streepath} {r.inputpath}")
    #nrf.append(
    ctp = comparetreepath(r.streepath, r.inputpath)
    nrf.append(ctp[-1])
    nfp.append(ctp[3])
    nfn.append(ctp[4])
    #)
    #mqsst.append(compute_mqsst(r.gtreepath, r.inputpath))
    mqsst.append(1)
    metafile = r.inputpath + ".meta"
    if os.path.isfile(metafile):
        with open(metafile) as fh:
            s = fh.read().split(";")
            time.append(float(s[0]))
            mem.append(float(s[1]) / 1e6)
    else:
        time.append(float('nan'))
        mem.append(float('nan'))
df["rf"] = nrf
df["fp"] = nfp
df["fn"] = nfn
df["time"] = time
df["mem"] = mem
df["mqsst"] = mqsst

df.to_csv("rcsv.csv")

dfa = df.set_index(['condition', 'k', 'method'])[['rf', 'time', 'mem']]

with pd.option_context("display.max_rows", 1000):
    ddf = dfa.groupby(level=list(range(3)))
    # import code; code.interact(local=dict(globals(), **locals()))
    mdfa = ddf.mean()
    print(mdfa)
    print(mdfa.to_latex(escape=False,float_format="{:0.3f}".format))

def f7(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]
order = f7(df.method)
first_method = df.method[0]

for c in set(df.condition):
    for k in set(df.k):
        for m in set(df.method):
            #import code; code.interact(local=dict(globals(), **locals()))
            myrf = dfa.loc[(c, k, m)].rf
            herrf = dfa.loc[(c, k, first_method)].rf
            if (len(myrf) != len(herrf)):
                continue
            if np.allclose(myrf, herrf):
                continue
            print(c, k, m)
            pvalue = wilcoxon(myrf, herrf).pvalue
            if pvalue < 0.05:
                print("Significant!")
            print("{:.4f}".format(pvalue))# rfs = df[(df.condition == c) ]
toremove = []
for c in set(df.condition):
    for k in set(df.k):
        for m in set(df.method):
            myrf = dfa.loc[(c, k, m)].rf
            herrf = dfa.loc[(c, k, first_method)].rf
            if len(myrf) < len(herrf):
                toremove.append((c,k,m))

for c,k,m in toremove:
    i = dfa.loc[(c,k, m)].index
    # print(i)
    df = df[df.apply(lambda x: not (x.condition == c and x.k == k and x.method == m), axis=1)]
    # df.drop(i, inplace=True)
dfa = df.set_index(['condition', 'k', 'method'])[['rf', 'time', 'mem']]


import seaborn as sns
sns.set_style("whitegrid")
c = set(df.condition)
for a in ["rf", "time", "mem"]:
    results = {}
    for x in c:
        # we should normalize across all conditions
        with sns.plotting_context(font_scale=1.2):
            g = sns.catplot(x="k", y=a, hue="method",
                        data=df[df.condition == x], kind="bar", ci="sd",hue_order=order)
            legend = g._legend
            legend.set_title("Method")
            axes = g.axes.flatten()
            axes[0].set_title(localization.get(x, x))
            # g.set_titles()
            # for t in legend.texts:
                # t.set_text({
                    # "fastral": "FASTRAL",
                    # "astral": "ASTRAL",
                    # "astrid": "ASTRID",
                    # "fastral2ftc": "FASTRAL-J",
                    # "astralj2": "ASTRAL-J",
                # }[t.get_text()])
            g.set_axis_labels("# Genes",
            {"time": "Time (s)", "rf": "nRF", "mem": "Peak Memory (GB)"}[a])
            results[f"{x}_{a}.pdf"] = g
    maxylim = max(results[k].axes[0,0].get_ylim()[1] for k in results)
    for k in results:
        results[k].set(ylim=(0, maxylim))
        results[k].savefig(k)
        # plt.close()
# sns.factorplot(x='k', y ='rf', hue='method', data = df, kind='bar')
