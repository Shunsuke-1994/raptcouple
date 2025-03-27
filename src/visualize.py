import matplotlib.pyplot as plt
import seaborn as sns
import sys
sys.path.append("./")
from src.plmc import read_params


def plot_from_param_file(
        paramfile,
        figsize = (10, 10),
        use_gap = False,
        cmap = "coolwarm",
        height_ratios = [1, 10],
        title = "", 
        cbar_labels = (
            "Positional field",
            "Frobenius Norm of APC coupling"
            ),
        square = False
        ):


    params = read_params(paramfile)

    pos_table = params["hi"].T
    score_table = params["FN_apc"]
    sequence = params["target_seq"]

    if not use_gap:
        col_match, sequence = zip(*[(i,s) for i,s in enumerate(params["target_seq"]) if s!="."])
        pos_table = pos_table[:, col_match]
        score_table = score_table[:, col_match][col_match, :]



    fig, axes = plt.subplots(2, 1, figsize = figsize, gridspec_kw = {"height_ratios":height_ratios})

    # plot positional field
    sns.heatmap(
        pos_table,
        cmap = cmap,
        linewidths=0.5,
        center=0,
        square=square,
        ax=axes[0],
        cbar_kws = {"label":cbar_labels[0]}
        )
    axes[0].set_xticks([i+0.5 for i in range(len(sequence))], list(sequence), rotation = 0);
    axes[0].set_yticks([i+0.5 for i in range(5)], list(params["alphabet"]), rotation = 270);
    axes[0].set_title(title)

    # plot coupling
    sns.heatmap(
        score_table,
        cmap = cmap,
        linewidths=0.5,
        center=0,
        square=square,
        ax=axes[1],
        cbar_kws = {"label":cbar_labels[1]}
        )
    axes[1].set_xticks(
        [i+0.5 for i in range(len(sequence))],
        list(sequence),
        rotation=0
    );
    axes[1].set_yticks(
        [i+0.5 for i in range(len(sequence))],
        list(sequence),
        rotation=270
    );

    return fig, axes