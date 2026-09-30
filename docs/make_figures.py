"""Draw the figures used in README.md.

    python docs/make_figures.py

Writes docs/pipeline.png (the two ways to run RaptCouple) and docs/example_output.png
(couplings and the structure they support, from the model in example/Ishida2020/).
Needs matplotlib; the second figure also reads the bundled model file.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from src.plmc import read_params
from src.structure import fold

HERE = os.path.dirname(__file__)
BLUE, GREEN, GREY, TEXT = "#0072B2", "#009E73", "#EEEEEE", "#222222"
plt.rcParams.update({"font.family": "sans-serif",
                     "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"]})


def box(ax, x, y, w, h, label, colour, text_colour=TEXT, size=9, weight="normal"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.06",
                                linewidth=0, facecolor=colour))
    ax.text(x + w / 2, y + h / 2, label, ha="center", va="center", fontsize=size,
            color=text_colour, weight=weight, linespacing=1.45)


def arrow(ax, start, end, colour=TEXT):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=11,
                                 linewidth=1.1, color=colour, shrinkA=2, shrinkB=2))


def pipeline():
    fig, ax = plt.subplots(figsize=(9.6, 3.5))
    ax.set_xlim(0, 10); ax.set_ylim(0, 3.6); ax.axis("off")

    box(ax, 0.1, 1.45, 1.35, 0.7, "SELEX\nreads", GREY)
    box(ax, 1.75, 1.45, 1.5, 0.7, "preprocess\n(trim, count, merge)", GREY, size=8)

    # the two ways to obtain a query
    box(ax, 3.55, 2.35, 2.5, 0.75,
        "you supply a query\ntarget-guided mode", BLUE, "white", 8.5, "bold")
    box(ax, 3.55, 0.5, 2.5, 0.75,
        "seeds from the pool itself\nquery-free mode (STREME / vsearch)", GREEN, "white", 8.5, "bold")

    box(ax, 6.35, 1.45, 1.5, 0.7, "jackhmmer\nalignment", GREY, size=8)
    box(ax, 8.15, 1.45, 1.7, 0.7, "Potts model\n(plmc)", GREY, size=8)

    arrow(ax, (1.45, 1.8), (1.75, 1.8))
    arrow(ax, (3.25, 1.8), (3.55, 2.72))
    arrow(ax, (3.25, 1.8), (3.55, 0.87))
    arrow(ax, (6.05, 2.72), (6.35, 1.8))
    arrow(ax, (6.05, 0.87), (6.35, 1.8))
    arrow(ax, (7.85, 1.8), (8.15, 1.8))

    # what the trained model gives you
    outputs = ["sequence motifs", "coupling-supported structure", "mutation effects"]
    ax.text(9.0, 1.2, "\n".join("• " + o for o in outputs), ha="center", va="top",
            fontsize=8, color=TEXT, linespacing=1.6)

    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "pipeline.png"), dpi=200, bbox_inches="tight",
                facecolor="white")
    plt.close(fig)


def example_output():
    model = os.path.join(HERE, "..", "example", "Ishida2020", "outputs",
                         "Ishida2020-6R-1-2626-55264.43.model_params")
    params = read_params(model)
    sequence, _, structure, _ = fold(params, min_loop_length=3, threshold=3)

    columns = [i for i, base in enumerate(params["target_seq"]) if base in "AUGC"]
    scores = np.asarray(params["FN_apc"], float)[np.ix_(columns, columns)]
    upper = np.triu(scores, k=4)
    values = upper[upper != 0]
    z = (scores - values.mean()) / values.std()
    n = len(sequence)
    separation = np.abs(np.subtract.outer(np.arange(n), np.arange(n)))
    z[separation < 4] = np.nan

    fig, (top, bottom) = plt.subplots(2, 1, figsize=(6.4, 6.8),
                                      gridspec_kw={"height_ratios": [1, 9], "hspace": 0.16})
    top.axis("off")
    top.set_xlim(-0.5, n - 0.5)
    for i, (base, bracket) in enumerate(zip(sequence, structure)):
        top.text(i, 0.62, base, ha="center", va="center", fontsize=6.5, color=TEXT)
        top.text(i, 0.12, bracket, ha="center", va="center", fontsize=7,
                 color=BLUE if bracket != "." else "#BBBBBB")
    top.set_ylim(-0.1, 0.95)

    image = bottom.imshow(z, cmap="magma", vmin=-1, vmax=8, interpolation="nearest")
    pairs, stack = [], []
    for i, bracket in enumerate(structure):
        if bracket == "(":
            stack.append(i)
        elif bracket == ")":
            pairs.append((stack.pop(), i))
    for i, j in pairs:
        bottom.plot(j, i, marker="s", markersize=5, markerfacecolor="none",
                    markeredgecolor="#56B4E9", markeredgewidth=1.1)
    ticks = [0] + list(range(9, n, 10))
    bottom.set_xticks(ticks, [t + 1 for t in ticks], fontsize=7)
    bottom.set_yticks(ticks, [t + 1 for t in ticks], fontsize=7)
    bottom.set_xlabel("position", fontsize=8)
    bar = fig.colorbar(image, ax=bottom, fraction=0.045, pad=0.03)
    bar.set_label("coupling z-score", fontsize=8)
    bar.ax.tick_params(labelsize=7)
    top.set_title("couplings, and the base pairs they support (squares)", fontsize=9, pad=4)

    fig.savefig(os.path.join(HERE, "example_output.png"), dpi=200, bbox_inches="tight",
                facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    pipeline()
    example_output()
    print("wrote docs/pipeline.png and docs/example_output.png")
