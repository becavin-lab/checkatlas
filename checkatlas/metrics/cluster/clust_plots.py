from typing import Collection

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from . import clust_compute as compute


def hist_external(
    adata,
    partition_keys: Collection[str],
    obsm_representation: Collection[str] = None,
):
    metrics_dict = dict(
        Davies_Bouldin=compute.davies_bouldin,
        silouhette=compute.silhouette,
        calinski=compute.calinski_harabasz,
    )
    # ,
    # DBCV=compute.dbcv)  # Import dbcv...
    # f = plt.figure()
    i = 1
    for name, func in metrics_dict.items():
        plt.subplot(1, len(metrics_dict), i)
        to_plot_array = []
        for key in partition_keys:
            score = func(adata, partition_key=key, obsm_representation=None)
            method_representation = f"{key}_original"
            metric = name
            to_plot_array.append([score, method_representation, metric])
            if obsm_representation:
                for rep in obsm_representation:
                    score = func(
                        adata, partition_key=key, obsm_representation=rep
                    )
                    method_representation = f"{key}_{rep}"
                    metric = name
                    to_plot_array.append(
                        [score, method_representation, metric]
                    )
        to_plot_df = pd.DataFrame(
            to_plot_array,
            columns=["score", "clustering_representation", "metric"],
        )
        sns.barplot(
            x="metric",
            y="score",
            hue="clustering_representation",
            ax=plt.gca(),
            data=to_plot_df,
        )
        i += 1


def hist_internal(adata, partition_keys: Collection[str], reference):
    metrics_dict = dict(
        rand=compute.rand,
        FowlkesMallows=compute.fowlkes_mallows,
        NMI=compute.nmi,
        vmeasure=compute.vmeasure,
    )
    # f = plt.figure()
    to_plot_array = []
    i = 0
    for key in partition_keys:
        for name, func in metrics_dict.items():
            score = func(adata, key, reference)
            clustering = key
            metric = name
            to_plot_array.append([score, clustering, metric])
            i += 1
    to_plot_df = pd.DataFrame(
        to_plot_array, columns=["score", "clustering", "metric"]
    )
    out = sns.barplot(x="metric", y="score", hue="clustering", data=to_plot_df)
    return out
