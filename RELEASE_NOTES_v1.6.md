# Release v1.6 (3 October 2026)

Figure layout and figure provenance. No analysis and no result table changes.

## Changes

- `make_manuscript_figures.py` (new): rebuilds manuscript Figure 1 (five panels) and Figure 4 (three panels) from the deposited tables. Earlier, no script in the repository produced these two composites; `make_figures.py` writes the older single-panel versions. Panels are placed in fixed positions with a label zone for each panel, and the script stops with an error if any axis label, tick label or in-panel text would enter another panel or leave the canvas.
- `make_graphical_abstract.py`: panels are now lettered a to c, positioned explicitly, and checked with the same layout test. The axis label of panel b is shortened and the threshold label of panel c is placed inside the plot. Plotted numbers are unchanged.
- Layout fixes in Figure 1 (long labels of panel b no longer run into panel a), Figure 4 (labels of panel c no longer run into panel b, the threshold label is no longer clipped, the legend of panel a is inside the plot) and the graphical abstract (panel labels no longer overlap).
- `figures/`: the five submitted images, numbered as in the paper (`Figure_1.png` to `Figure_4.png`, `Graphical_Abstract.tif`).
- `README.md`: describes the figure scripts and the numbering.

## Not changed

Figures 2 and 3 are stored as images and are pixel-identical to `figures/Figure3.png` and `figures/Figure4.png`; no script in this repository draws them.
