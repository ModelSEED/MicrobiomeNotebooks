"""pheatmap-style heatmaps in matplotlib (the original Figure 1 and Figure S4 were drawn with pheatmap).

draw_heatmap(blocks, ...) draws one or more row blocks side by side that share the column order and the column
dendrogram; each block has left annotation bars (row annotations), the matrix, and right row labels. Legends go on the
right. Sizes are in inches.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
from scipy.cluster.hierarchy import dendrogram

FONT = "Liberation Sans"
plt.rcParams["font.family"] = FONT
BWR = LinearSegmentedColormap.from_list("pheat_bwr", ["#1029ff", "#8fa0ff", "#ffffff", "#ff8f87", "#ff1c0d"])


class Canvas:
    def __init__(self, w, h):
        self.w, self.h = w, h
        self.fig = plt.figure(figsize=(w, h))

    def ax(self, x, y, w, h, **kw):
        """x, y from the top-left corner, in inches."""
        a = self.fig.add_axes([x / self.w, 1 - (y + h) / self.h, w / self.w, h / self.h], **kw)
        return a


def blank(a):
    a.set_xticks([]); a.set_yticks([])
    for s in a.spines.values():
        s.set_visible(False)


def draw_dendro(a, Z, n):
    dendrogram(Z, ax=a, no_labels=True, color_threshold=0, above_threshold_color="#000000",
               link_color_func=lambda k: "#000000")
    a.set_xlim(0, 10 * n)
    blank(a)
    for c in a.collections:
        c.set_linewidth(1.0)


def color_bar(a, colors, vertical):
    n = len(colors)
    for i, c in enumerate(colors):
        if vertical:
            a.add_patch(Rectangle((0, n - i - 1), 1, 1, facecolor=c, edgecolor="none"))
        else:
            a.add_patch(Rectangle((i, 0), 1, 1, facecolor=c, edgecolor="none"))
    if vertical:
        a.set_xlim(0, 1); a.set_ylim(0, n)
    else:
        a.set_xlim(0, n); a.set_ylim(0, 1)
    blank(a)


def legend_block(c, x, y, title, items, fs=8.5, box=0.13, step=0.19):
    """Draw a pheatmap-style legend at (x, y) inches; returns the y below it."""
    c.fig.text(x / c.w, 1 - y / c.h, title, fontsize=fs + 1, weight="bold", va="top")
    y += 0.24
    for lab, col in items:
        a = c.ax(x, y, box, box)
        a.add_patch(Rectangle((0, 0), 1, 1, facecolor=col, edgecolor="none")); a.set_xlim(0, 1); a.set_ylim(0, 1)
        blank(a)
        c.fig.text((x + box + 0.07) / c.w, 1 - (y + box / 2) / c.h, lab, fontsize=fs, va="center")
        y += step
    return y + 0.12


def colorbar(c, x, y, h, cmap, vmin, vmax, ticks, fs=8.5, title=None):
    if title:
        c.fig.text(x / c.w, 1 - y / c.h, title, fontsize=fs + 1, weight="bold", va="top")
        y += 0.24
    a = c.ax(x, y, 0.16, h)
    grad = np.linspace(vmax, vmin, 256)[:, None]
    a.imshow(grad, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax, extent=[0, 1, vmin, vmax])
    a.set_xticks([]); a.yaxis.tick_right(); a.set_yticks(ticks)
    a.set_yticklabels([("%g" % t).replace("-", "−") for t in ticks], fontsize=fs)
    a.tick_params(axis="y", length=0, pad=3)
    for s in a.spines.values():
        s.set_visible(False)
    return y + h + 0.3
