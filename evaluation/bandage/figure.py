#!/usr/bin/env python3
"""Appendix figures: the recorded Bandage panels of the forward-strand top-depth set of every dataset (rows) and tool
(columns HARICOT, Rawsamble, Flye, wtdbg2), as two separate PDFs for a US-letter manuscript: OUTPREFIX_bacteria.pdf
(D1-D6) and OUTPREFIX_eukaryotes.pdf (D7-D14), each one grid trimmed to its content (no page, no margins), 6.5 in wide
(the text width of a letter page with 1 in margins; include with width=\\textwidth), at most 7.9 in tall (a letter page
with its caption), plus a PNG preview of each.

Inputs (read only): ROOT/<key>/bandage/{ours,rawsamble,flye,wtdbg2}.png (panels.sh), TABLE (table.py: contigs, QUAST
contigs >= 1 kb, for the per-panel count) and ROOT/<key>/rawsamble/miniasm_protocol.txt (PLAIN / DEVIATION =
density-filtered retry, marked with a dagger / FAIL = no assembly). Colours (panels.sh, color_gfa.py): each segment's
hue is the reference midpoint of its longest primary minimap2 (asm10) block, 0 -> 0.85 along the main contig; light
grey = not aligned, dark grey = aligned to another contig. HARICOT's unitig graph is coloured by the same rule from its
evaluation sequence (anchor_seq.py), so the four columns share one colouring method. Rawsamble's panel is its miniasm
graph (asm.gfa). Each panel is cropped to its drawing (white margins removed) and centred in its cell; nothing else is
changed (native 2000 px, embedded unresampled). Text is vector (Liberation Sans, Arial metrics, TrueType-embedded).

Usage: figure.py ROOT TABLE OUTPREFIX   -> OUTPREFIX_{bacteria,eukaryotes}.pdf, OUTPREFIX_{bacteria,eukaryotes}.png
"""
import argparse
import colorsys
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib import font_manager  # noqa: E402
from PIL import Image, ImageOps  # noqa: E402

FONT = "Liberation Sans"
font_manager.findfont(FONT, fallback_to_default=False)  # fail rather than silently fall back
matplotlib.rcParams.update({"font.family": FONT, "pdf.fonttype": 42, "font.size": 7})

DATASETS = [  # (label, key, species, chromosome)
    ("D1", "cj_fwd", "C. jejuni", ""),
    ("D2", "lm_fwd", "L. monocytogenes", ""),
    ("D3", "vp1_fwd", "V. parahaemolyticus", "chr. 1"),
    ("D4", "vp2_fwd", "V. parahaemolyticus", "chr. 2"),
    ("D5", "se_fwd", "S. enterica", ""),
    ("D6", "sd_fwd", "S. dysgalactiae", ""),
    ("D7", "at2_fwd", "A. thaliana", "chr. 2"),
    ("D8", "at3_fwd", "A. thaliana", "chr. 3"),
    ("D9", "at5_fwd", "A. thaliana", "chr. 5"),
    ("D10", "dm2L_fwd", "D. melanogaster", "2L"),
    ("D11", "dm2R_fwd", "D. melanogaster", "2R"),
    ("D12", "dm3L_fwd", "D. melanogaster", "3L"),
    ("D13", "dmX_fwd", "D. melanogaster", "X"),
    ("D14", "hg18_fwd", "H. sapiens HG002", "chr. 18"),
]
FIGURES = [("bacteria", "bacteria", DATASETS[0:6]), ("eukaryotes", "eukaryotes", DATASETS[6:14])]
TOOLS = [("HARICOT", "ours", "ours_seq"), ("Rawsamble", "rawsamble", "rawsamble"), ("Flye", "flye", "flye"),
         ("wtdbg2", "wtdbg2", "wtdbg2")]  # (column title, panel name, table tool)

# geometry, inches; savefig trims to the drawing plus TRIM_PAD on each side, so the PDF is FIG_W + 2 TRIM_PAD = 6.5 wide
TRIM_PAD = 0.01
FIG_W, MAX_FIG_H = 6.5 - 2 * TRIM_PAD, 7.9 - 2 * TRIM_PAD
HEAD_H, NOTE_H, GAP = 0.17, 0.12, 0.035  # column-title band; contig-count line under each panel; space between rows
LEG_GAP, LEG_H = 0.06, 0.13  # space above the legend line; its height
LABEL_PAD = 0.07  # between the longest row label and the first panel
# font sizes, pt
FS_HEAD, FS_ROW, FS_SP, FS_NOTE, FS_LEG = 8, 8, 7, 6.5, 6.5


def text_w(s, size, **kw):
    """width of s in inches when rendered at size pt"""
    fig = plt.figure(figsize=(2, 2))
    t = fig.text(0, 0, s, fontsize=size, **kw)
    w = t.get_window_extent(fig.canvas.get_renderer()).width / fig.dpi
    plt.close(fig)
    return w


LABEL_W = max(max(text_w(lab, FS_ROW, fontweight="bold"), text_w(sp, FS_SP, style="italic"), text_w(ch, FS_SP))
              for lab, _, sp, ch in DATASETS) + LABEL_PAD
PANEL_W = (FIG_W - LABEL_W) / len(TOOLS)


def contigs(table):
    out = {}
    with open(table) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        ki, ti, qi = head.index("key"), head.index("tool"), head.index("contigs")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out[(f[ki], f[ti])] = f[qi]
    return out


def panel(path, aspect):
    """the panel cropped to its drawing (white margins removed), padded white to the cell's aspect (w / h)"""
    im = Image.open(path).convert("RGB")
    box = ImageOps.invert(im).getbbox()
    if box:
        pad = 20
        im = im.crop((max(box[0] - pad, 0), max(box[1] - pad, 0), min(box[2] + pad, im.width),
                      min(box[3] + pad, im.height)))
    w, h = im.size
    tw, th = (w, round(w / aspect)) if w / h >= aspect else (round(h * aspect), h)
    canvas = Image.new("RGB", (tw, th), "white")
    canvas.paste(im, ((tw - w) // 2, (th - h) // 2))
    return np.asarray(canvas)


def grid(rows, root, nc, rs_proto, rs_missing):
    n = len(rows)
    fixed = HEAD_H + LEG_GAP + LEG_H - GAP
    ph = min(PANEL_W * 3 / 4, (MAX_FIG_H - fixed) / n - NOTE_H - GAP)  # Bandage's 4:3 where the height allows
    row_h = ph + NOTE_H + GAP
    fig_h = fixed + n * row_h
    fig = plt.figure(figsize=(FIG_W, fig_h))

    def X(x):  # inches from the left edge -> figure fraction
        return x / FIG_W

    def Y(t):  # inches from the top edge -> figure fraction (from the bottom)
        return 1 - t / fig_h

    for j, (name, _, _) in enumerate(TOOLS):
        fig.text(X(LABEL_W + (j + 0.5) * PANEL_W), Y(0.5 * HEAD_H), name, ha="center", va="center",
                 fontsize=FS_HEAD, fontweight="bold")
    for i, (lab, key, sp, chrom) in enumerate(rows):
        top = HEAD_H + i * row_h
        lines = [(lab, dict(fontsize=FS_ROW, fontweight="bold")), (sp, dict(fontsize=FS_SP, style="italic"))]
        if chrom:
            lines.append((chrom, dict(fontsize=FS_SP)))
        lead = 0.125  # line spacing, inches
        y0 = top + 0.5 * ph - 0.5 * lead * (len(lines) - 1)
        for k, (s, kw) in enumerate(lines):
            fig.text(X(0), Y(y0 + k * lead), s, va="center", ha="left", **kw)
        for j, (name, suffix, tkey) in enumerate(TOOLS):
            ax = fig.add_axes([X(LABEL_W + j * PANEL_W), Y(top + ph), PANEL_W / FIG_W, ph / fig_h])
            ax.set_xticks([])
            ax.set_yticks([])
            for s in ax.spines.values():
                s.set_linewidth(0.4)
                s.set_color("#999999")
            p = os.path.join(root, key, "bandage", f"{suffix}.png")
            proto = rs_proto.get(key, "") if name == "Rawsamble" else "PLAIN"
            if proto == "FAIL" or not os.path.exists(p):
                ax.text(0.5, 0.5, rs_missing if proto == "FAIL" else "no panel", ha="center", va="center",
                        fontsize=FS_SP, style="italic", color="#666666", transform=ax.transAxes)
                note = ""
            else:
                ax.imshow(panel(p, PANEL_W / ph), interpolation="none", aspect="auto")
                c = nc.get((key, tkey), "")
                note = f"{c} contig{'s' if c != '1' else ''}" + (" †" if proto == "DEVIATION" else "")
            fig.text(X(LABEL_W + (j + 0.5) * PANEL_W), Y(top + ph + 0.5 * NOTE_H), note, ha="center",
                     va="center", fontsize=FS_NOTE)
    legend(fig, X, Y, fig_h - LEG_H)
    return fig, ph, fig_h


def legend(fig, X, Y, top):
    """one line under the grid, from the first panel's left edge: hue ramp (reference start -> end), the two greys"""
    yc = top + 0.5 * LEG_H
    sw_h = 0.075  # swatch / ramp height, inches
    x = LABEL_W
    for s in ("Position on the reference:", "start"):
        fig.text(X(x), Y(yc), s, fontsize=FS_LEG, va="center", ha="left")
        x += text_w(s, FS_LEG) + 0.06
    ramp_w = 1.3
    ax = fig.add_axes([X(x), Y(yc + 0.5 * sw_h), ramp_w / FIG_W, sw_h / fig.get_figheight()])
    ramp = np.array([[colorsys.hsv_to_rgb(0.85 * v, 0.85, 0.95) for v in np.linspace(0, 1, 256)]])
    ax.imshow(ramp, aspect="auto", extent=(0, 1, 0, 1), interpolation="bilinear")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_linewidth(0.4)
    x += ramp_w + 0.06
    fig.text(X(x), Y(yc), "end", fontsize=FS_LEG, va="center", ha="left")
    x += text_w("end", FS_LEG) + 0.3
    for col, s in (("#cccccc", "not aligned"), ("#555555", "aligned to another contig")):
        a = fig.add_axes([X(x), Y(yc + 0.5 * sw_h), 0.13 / FIG_W, sw_h / fig.get_figheight()])
        a.add_patch(plt.Rectangle((0, 0), 1, 1, color=col))
        a.axis("off")
        fig.text(X(x + 0.18), Y(yc), s, fontsize=FS_LEG, va="center", ha="left")
        x += 0.18 + text_w(s, FS_LEG) + 0.3


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root")
    ap.add_argument("table")
    ap.add_argument("out")
    a = ap.parse_args()
    nc = contigs(a.table)
    rs_missing = "failed"  # Rawsamble left no assembly
    rs_proto = {}
    for _, key, _, _ in DATASETS:
        p = os.path.join(a.root, key, "rawsamble", "miniasm_protocol.txt")
        if os.path.exists(p):
            rs_proto[key] = open(p).read().split(":")[0].strip()
    print(f"font {font_manager.findfont(FONT)}; label column {LABEL_W:.3f} in; panel width {PANEL_W:.3f} in")
    for suffix, title, rows in FIGURES:
        fig, ph, fig_h = grid(rows, a.root, nc, rs_proto, rs_missing)
        out = f"{a.out}_{suffix}"
        fig.savefig(out + ".pdf", bbox_inches="tight", pad_inches=TRIM_PAD, metadata={
            "Title": f"Bandage panels, forward-strand {title} datasets at the top depth",
            "Subject": f"Bandage panels of {a.root}",
            "Creator": "HARICOT evaluation/bandage/figure.py"})
        fig.savefig(out + ".png", bbox_inches="tight", pad_inches=TRIM_PAD, dpi=200)
        plt.close(fig)
        print(f"{out}.pdf  panel {PANEL_W:.3f} x {ph:.3f} in; layout {FIG_W} x {fig_h:.3f} in before trimming")


if __name__ == "__main__":
    main()
