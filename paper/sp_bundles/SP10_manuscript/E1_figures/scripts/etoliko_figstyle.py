"""etoliko_figstyle.py -- manuscript (MDPI Microorganisms) version of the thesis house style.

Drop-in replacement for workspace_2/scripts/etoliko_figstyle.py: same names, so every copied
fig*.py imports this one unchanged. Differences from the thesis module:
  - English; sans-serif (Arial, else Liberation Sans / DejaVu Sans)
  - white background
  - decimal point "."
  - DOMAIN colours Bacteria #1F5FA6, Archaea #D9622B
  - caption() draws nothing (captions live in captions.md)
  - widths: FULL_W_IN = 17 cm, HALF_W_IN = 8.5 cm
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

SURFACE   = "#ffffff"
INK       = "#0b0b0b"
INK_2     = "#3d3c3a"
INK_MUTED = "#6b6a66"
GRID      = "#e3e2df"
RULE      = "#b9b8b4"

PHYLUM = {
    "Others": "#acadaf", "Actinomycetota": "#71427c", "Atribacterota": "#205578",
    "Spirochaetota": "#336c8a", "Pseudomonadota": "#e1662c", "Chloroflexota": "#ed9f47",
    "Cyanobacteriota": "#f0b532", "Bacillota": "#bf66c5", "Acidobacteriota": "#932284",
    "Thermodesulfobacteriota": "#a689a0", "Bacteroidota": "#836389", "Unassigned": "#f5e030",
}
DOMAIN    = {"Bacteria": "#1F5FA6", "Archaea": "#D9622B"}
DOMAIN_EL = {"Bacteria": "Bacteria", "Archaea": "Archaea"}   # name kept so the scripts run unchanged

HILITE = "#932284"
STATUS = {"good": "#1baf7a", "warning": "#eda100", "critical": "#e34948"}
PHYLUM_ORDER = ["Pseudomonadota", "Atribacterota", "Cyanobacteriota", "Acidobacteriota",
                "Chloroflexota", "Actinomycetota", "Bacillota", "Spirochaetota",
                "Bacteroidota", "Thermodesulfobacteriota", "Unassigned", "Others"]

DECIMAL = "."

def num(v, dp=4):
    return ("%.*f" % (dp, v)).replace(".", DECIMAL)

# >= 7 pt everywhere at final size
PT_TITLE, PT_LABEL, PT_TICK, PT_VALUE, PT_CAPTION = 8, 8, 7, 7, 7
CM = 1 / 2.54
FULL_W_IN, HALF_W_IN = 17 * CM, 8.5 * CM

_SANS_PREF = ("Arial", "Helvetica", "Liberation Sans", "DejaVu Sans")

def _pick_font():
    have = {f.name for f in font_manager.fontManager.ttflist}
    for c in _SANS_PREF:
        if c in have:
            return c
    return "DejaVu Sans"

FONT = _pick_font()

def check_greek(verbose=True):
    if verbose:
        print("   font: %s (manuscript style, English)" % FONT)
    return True

def apply(verify_greek=True):
    check_greek()
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "font.family": "sans-serif", "font.sans-serif": [FONT, "DejaVu Sans"], "font.size": PT_LABEL,
        "text.color": INK, "axes.labelcolor": INK_2,
        "axes.edgecolor": RULE, "axes.linewidth": 0.6,
        "axes.titlesize": PT_TITLE, "axes.titleweight": "bold", "axes.labelsize": PT_LABEL,
        "xtick.color": INK_2, "ytick.color": INK_2,
        "xtick.labelsize": PT_TICK, "ytick.labelsize": PT_TICK,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "legend.frameon": False, "legend.fontsize": PT_TICK,
        "grid.color": GRID, "grid.linewidth": 0.6,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    })

def despine(ax, keep=("left", "bottom")):
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(s in keep)

def recessive_grid(ax, axis="x"):
    ax.grid(axis=axis, zorder=0)
    ax.set_axisbelow(True)

def caption(fig, text, y=0.012, width=112):
    """Manuscript figures carry no caption. The text is still laid out (some scripts
    measure it to size the figure) but is invisible and removed again by save()."""
    import textwrap
    CAPTIONS.append(text)
    t = fig.text(0.012, y, textwrap.fill(text, width), ha="left", va="bottom",
                 fontsize=PT_CAPTION, color=INK_2, linespacing=1.5)
    t.set_gid("caption"); t.set_visible(False)
    return t

CAPTIONS = []

def _strip_titles(fig):
    """No title or caption text inside the image: drop suptitles and axes titles
    (panel letters a, b, ... are kept), and figure-level free text."""
    import re
    for t in list(fig.texts):
        if t.get_gid() == "caption":
            t.remove()
    if fig._suptitle is not None:
        fig._suptitle.set_text("")
    for ax in fig.axes:
        for loc in ("center", "left", "right"):
            t = ax.get_title(loc=loc)
            if t and not re.fullmatch(r"\(?[a-z]\)?\.?", t.strip()):
                ax.set_title("", loc=loc)

def save(fig, stem, outdir=".", dpi=600):
    os.makedirs(outdir, exist_ok=True)
    _strip_titles(fig)
    w, h = fig.get_size_inches()
    if w > FULL_W_IN * 1.02:          # scale to 17 cm: fonts were authored for the thesis block
        print("   note: authored at %.1f cm wide; exported at that size, fonts checked separately" % (w * 2.54))
    made = []
    for ext, kw in (("pdf", {}), ("png", {"dpi": dpi})):
        p = os.path.join(outdir, "%s.%s" % (stem, ext))
        fig.savefig(p, bbox_inches="tight", pad_inches=0.03, **kw)
        made.append(p)
    plt.close(fig)
    for p in made:
        print("   wrote %s (%.0f KB)" % (p, os.path.getsize(p) / 1024))
    return made
