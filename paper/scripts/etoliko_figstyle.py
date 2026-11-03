"""etoliko_figstyle.py -- one house style for every thesis figure.

Greek thesis, Times. Imported by every fig*.py, so a font or colour changed here
changes every panel together.

FONT. Times New Roman if the machine has it; otherwise Liberation Serif, which is
metric-compatible with Times and carries the full Greek range. The list is tried in
order and the choice is REPORTED, because a silent fallback to a font without Greek
prints tofu boxes that look like a rendering glitch rather than a missing font.
`check_greek()` refuses to continue if the chosen face cannot draw Greek.

DECIMAL SEPARATOR. Greek convention is the comma: 0,8936. Set DECIMAL = "." if the
thesis uses the point instead -- it is one line, and it is the sort of thing an
examiner notices.

PALETTE. Taken from the thesis's existing donut figure so the chapter reads as one
document -- see the block below for exactly which values were kept, which was nudged
and why, and what the original palette can and cannot be asked to do.
Colour encodes DOMAIN only in these panels. Nothing is identified by colour alone.
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

# ---------------------------------------------------------------- palette
# TAKEN FROM THE THESIS'S EXISTING FIGURE so the chapter reads as one document.
# The twelve phylum colours below are sampled straight out of the published donut
# (inner ring) and are unchanged.
#
# HONEST NOTE ON WHAT THAT PALETTE CAN AND CANNOT DO. Run through the six-check
# validator on a #fcfcfb surface, the twelve-slot set FAILS four checks:
#     lightness band   #f0b532 (0.81) and #f5e030 (0.90) sit above the band
#     chroma floor     #acadaf, #205578, #336c8a, #a689a0, #836389 read as grey
#     CVD separation   #205578 vs #71427c  dE 2.4 protan  -- indistinguishable
#     normal vision    #f0b532 vs #ed9f47  dE 6.1         -- hard for anyone
#     contrast         four slots below 3:1 against the surface
# It works in the donut anyway because every wedge is DIRECTLY LABELLED with its
# percentage and the legend sits alongside: colour is not carrying identity alone.
# Any figure here that uses PHYLUM must do the same -- direct labels, or a table.
# Where colour must carry identity by itself, use DOMAIN below instead.
SURFACE   = "#fcfcfb"
INK       = "#0b0b0b"
INK_2     = "#52514e"
INK_MUTED = "#7d7c78"
GRID      = "#e3e2df"
RULE      = "#c9c8c4"

PHYLUM = {                       # exact values from the existing figure
    "Others":                 "#acadaf",
    "Actinomycetota":         "#71427c",
    "Atribacterota":          "#205578",
    "Spirochaetota":          "#336c8a",
    "Pseudomonadota":         "#e1662c",
    "Chloroflexota":          "#ed9f47",
    "Cyanobacteriota":        "#f0b532",
    "Bacillota":              "#bf66c5",
    "Acidobacteriota":        "#932284",
    "Thermodesulfobacteriota":"#a689a0",
    "Bacteroidota":           "#836389",
    "Unassigned":             "#f5e030",
}

# DOMAIN is drawn FROM that palette so the two figure families match, but the pair
# has to survive on its own -- in a bar chart colour IS the identity.
#   Archaea  #e1662c  = Pseudomonadota, unchanged.
#   Bacteria #145d90  = Atribacterota #205578 nudged up in chroma. The original
#            fails the chroma floor at 0.08 ("reads grey"); #145d90 is the smallest
#            step that clears 0.10 and is visually the same dark blue.
# Validated as a pair on #fcfcfb:
#   lightness PASS · chroma PASS · contrast >=3:1 PASS
#   CVD separation dE 21.0 (protan)   -- target >=8
#   normal vision  dE 33.1            -- floor  >=15
# The middle ring's own domain colours (#face4b Bacteria, #f5bd72 Archaea) are NOT
# used here: as a standalone pair they score dE 6.4 to normal vision, which is
# below the floor -- they only work in the donut because of the direct labels.
DOMAIN    = {"Bacteria": "#145d90", "Archaea": "#e1662c"}
DOMAIN_EL = {"Bacteria": "Βακτήρια", "Archaea": "Αρχαία"}

# annotation / callout ink, from the same palette (Acidobacteriota)
HILITE = "#932284"
STATUS = {"good": "#1baf7a", "warning": "#eda100", "critical": "#e34948"}

# ordered list for figures that need several categories, worst adjacent pairs first
# separated; use ONLY with direct labels (see the note above).
PHYLUM_ORDER = ["Pseudomonadota", "Atribacterota", "Cyanobacteriota", "Acidobacteriota",
                "Chloroflexota", "Actinomycetota", "Bacillota", "Spirochaetota",
                "Bacteroidota", "Thermodesulfobacteriota", "Unassigned", "Others"]

# ---------------------------------------------------------------- numbers
DECIMAL = ","          # Greek convention. Set to "." to use the point.

def num(v, dp=4):
    """Format a number in the thesis's decimal convention."""
    return ("%.*f" % (dp, v)).replace(".", DECIMAL)

# ---------------------------------------------------------------- typography
PT_TITLE, PT_LABEL, PT_TICK, PT_VALUE, PT_CAPTION = 10.5, 9, 8, 7.5, 7.5

_SERIF_PREF = ("Times New Roman", "Liberation Serif", "Nimbus Roman",
               "Nimbus Roman No9 L", "FreeSerif", "DejaVu Serif")
GREEK_PROBE = "ΑΒΓΔΕΖΗΘΙΚΛΜΝΞΟΠΡΣΤΥΦΧΨΩαβγδεζηθικλμνξοπρστυφχψωάέίόύήώϊϋΐΰς"

def _pick_font():
    have = {f.name for f in font_manager.fontManager.ttflist}
    for c in _SERIF_PREF:
        if c in have:
            return c
    return "DejaVu Serif"

FONT = _pick_font()

def check_greek(verbose=True):
    """Refuse to render Greek in a face that cannot draw it.

    A missing glyph in matplotlib is a hollow rectangle. In a 29-row figure that is
    easy to mistake for a layout bug, and easy to miss altogether in a caption.
    """
    path = font_manager.findfont(font_manager.FontProperties(family=FONT),
                                 fallback_to_default=False)
    missing = []
    try:
        from fontTools.ttLib import TTFont
        t = TTFont(path, fontNumber=0, lazy=True)
        cmap = set()
        for tb in t["cmap"].tables:
            cmap |= set(tb.cmap.keys())
        t.close()
        missing = [ch for ch in GREEK_PROBE if ord(ch) not in cmap]
    except ImportError:
        if verbose:
            print("   [--]  fontTools not installed; Greek coverage NOT verified")
        return True
    if missing:
        raise SystemExit(
            "FONT ERROR: %s (%s) is missing %d Greek glyph(s): %s\n"
            "            Install Liberation Serif or Times New Roman, or edit\n"
            "            _SERIF_PREF in etoliko_figstyle.py."
            % (FONT, os.path.basename(path), len(missing), "".join(missing[:12])))
    if verbose:
        print("   font: %s  (%s) — Greek coverage complete"
              % (FONT, os.path.basename(path)))
    return True

def apply(verify_greek=True):
    if verify_greek:
        check_greek()
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "font.family": "serif", "font.serif": [FONT], "font.size": PT_LABEL,
        "text.color": INK, "axes.labelcolor": INK_2,
        "axes.edgecolor": RULE, "axes.linewidth": 0.6,
        "axes.titlesize": PT_TITLE, "axes.titleweight": "bold",
        "axes.labelsize": PT_LABEL,
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
    import textwrap
    fig.text(0.012, y, textwrap.fill(text, width), ha="left", va="bottom",
             fontsize=PT_CAPTION, color=INK_2, linespacing=1.5)

def save(fig, stem, outdir=".", dpi=600):
    os.makedirs(outdir, exist_ok=True)
    made = []
    for ext, kw in (("pdf", {}), ("png", {"dpi": dpi})):
        p = os.path.join(outdir, "%s.%s" % (stem, ext))
        fig.savefig(p, bbox_inches="tight", pad_inches=0.02, **kw)
        made.append(p)
    plt.close(fig)
    for p in made:
        print("   wrote %s (%.0f KB)" % (p, os.path.getsize(p) / 1024))
    return made
