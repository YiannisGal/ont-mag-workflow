# -*- coding: utf-8 -*-
"""etoliko_tree.py -- Newick parsing and rectangular phylogram drawing.

Written rather than imported. ete3 and Bio.Phylo would each do this, but this
project's rule is that nothing new goes into the machine's conda environments,
several of which are already broken, and a tree renderer is thirty lines of
layout plus a parser. Every figure that draws a tree uses this module, so a
change to branch styling or support filtering changes all of them together --
the same contract etoliko_figstyle.py has for colour and type.

WHAT IT DELIBERATELY DOES NOT DO. No ultrametric or circular layouts, no
automatic clade collapsing by taxonomy, no re-estimation of anything. A figure
must draw the tree that was inferred, and a renderer that quietly reshapes it
is a way of publishing a different tree from the one the methods describe.
"""
import re


class Node(object):
    __slots__ = ("name", "length", "support", "children", "parent", "x", "y")

    def __init__(self, name="", length=0.0, support=None):
        self.name = name
        self.length = length
        self.support = support
        self.children = []
        self.parent = None
        self.x = 0.0
        self.y = 0.0

    @property
    def is_tip(self):
        return not self.children

    def walk(self):
        yield self
        for c in self.children:
            for n in c.walk():
                yield n

    def tips(self):
        return [n for n in self.walk() if n.is_tip]


_NUM = re.compile(r"^[\d.eE+-]+$")


def _label(tok):
    """split 'name:0.123' or '95:0.123' or '95/100:0.12' into (name, length)"""
    if ":" in tok:
        nm, _, ln = tok.rpartition(":")
    else:
        nm, ln = tok, ""
    try:
        length = float(ln)
    except ValueError:
        length = 0.0
    return nm.strip().strip("'\""), length


def parse(text):
    """Newick -> root Node. Internal labels are kept as `support` when they
    look numeric and as `name` otherwise; IQ-TREE writes ultrafast bootstrap
    there, FastTree writes SH-like values, and a tree with neither must be
    reported as unsupported rather than drawn as if it were supported.

    The label that follows a ')' belongs to the node that was just CLOSED, not
    to the cursor -- the cursor has already moved to its parent. A first
    version of this parser assigned it to the parent and silently lost one
    support value per internal node; the unit test at the bottom of this file
    is what caught it, which is why it is there.
    """
    s = text.strip()
    if ";" in s:
        s = s[:s.rindex(";")]
    root = Node()
    cur = root
    buf = ""
    closed = None

    def flush():
        if not buf.strip():
            return
        if closed is not None:
            _assign(closed, buf)
        else:
            nm, ln = _label(buf)
            t = Node(nm, ln)
            t.parent = cur
            cur.children.append(t)

    for ch in s:
        if ch == "(":
            n = Node()
            n.parent = cur
            cur.children.append(n)
            cur = n
            buf = ""
            closed = None
        elif ch == ",":
            flush()
            buf = ""
            closed = None
        elif ch == ")":
            flush()
            buf = ""
            closed = cur
            cur = cur.parent
        else:
            buf += ch
    flush()
    while len(root.children) == 1 and not root.name and root.support is None:
        root = root.children[0]
        root.parent = None
    return root


def _assign(node, tok):
    """apply a 'label:length' token to an existing node"""
    nm, ln = _label(tok)
    node.length = ln
    if nm:
        if _NUM.match(nm):
            try:
                node.support = float(nm.split("/")[0])
            except ValueError:
                node.name = nm
        else:
            node.name = nm


def ladderize(node, ascending=True):
    """order children by clade size -- the only reordering applied, and it
    changes no topology"""
    for c in node.children:
        ladderize(c, ascending)
    node.children.sort(key=lambda c: len(c.tips()), reverse=not ascending)
    return node


def reroot_on(root, node):
    """place the root on the branch leading to `node`, so an outgroup clade
    can be put at one end. Branch lengths on the reversed path are preserved."""
    if node is root or node.parent is None:
        return root
    path = []
    n = node
    while n.parent is not None:
        path.append(n)
        n = n.parent
    old_root = n
    new_root = Node()
    half = (node.length or 0.0) / 2.0
    parent = node.parent
    parent.children.remove(node)
    node.length = half
    # reverse the chain from `parent` up to the old root
    prev, prev_len = node, half
    cur = parent
    chain = []
    while cur is not None:
        chain.append(cur)
        cur = cur.parent
    for i, nd in enumerate(chain):
        nxt = chain[i + 1] if i + 1 < len(chain) else None
        if nxt is not None:
            nxt.children.remove(nd)
    acc = half
    cur = parent
    carried = acc
    pieces = []
    for nd in chain:
        pieces.append((nd, nd.length or 0.0))
    # rebuild: node and the reversed spine both hang off new_root
    new_root.children = [node]
    node.parent = new_root
    spine = parent
    spine.length = carried
    new_root.children.append(spine)
    spine.parent = new_root
    prev = spine
    for i in range(1, len(chain)):
        nd = chain[i]
        nd.length = pieces[i - 1][1]
        nd.parent = prev
        prev.children.append(nd)
        prev = nd
    # drop any node left with exactly one child (the old root usually)
    _suppress_unifurcations(new_root)
    return new_root


def _suppress_unifurcations(node):
    for c in list(node.children):
        _suppress_unifurcations(c)
    if node.parent is not None and len(node.children) == 1:
        only = node.children[0]
        only.length = (only.length or 0.0) + (node.length or 0.0)
        only.parent = node.parent
        idx = node.parent.children.index(node)
        node.parent.children[idx] = only


def layout(root, cladogram=False):
    """x = distance from root (or depth), y = tip index. Returns the tip list
    in drawing order, top to bottom."""
    tips = []

    def setx(n, acc):
        n.x = acc + (0.0 if cladogram else (n.length or 0.0))
        if cladogram:
            n.x = acc + 1.0
        for c in n.children:
            setx(c, n.x)
    root.length = 0.0
    setx(root, 0.0)

    def sety(n):
        if n.is_tip:
            n.y = float(len(tips))
            tips.append(n)
        else:
            for c in n.children:
                sety(c)
            n.y = (n.children[0].y + n.children[-1].y) / 2.0
    sety(root)
    return tips


def draw(ax, root, color=None, lw=0.7, support_min=None, support_fmt=None,
         support_size=6.0, support_color="#52514e", support_halo="#fcfcfb"):
    """Rectangular phylogram. `color(node)` may return a per-branch colour.

    support_min: nodes at or above this get a printed value. Passing None
    prints nothing -- which is the correct default, because a tree drawn with
    some support values shown and no stated threshold invites the reader to
    assume the unlabelled nodes were tested and failed.
    """
    for n in root.walk():
        if n.parent is not None:
            c = color(n) if color else "#0b0b0b"
            ax.plot([n.parent.x, n.x], [n.y, n.y], color=c, lw=lw,
                    solid_capstyle="butt", zorder=2)
        if n.children:
            c = color(n) if color else "#0b0b0b"
            ax.plot([n.x, n.x], [n.children[0].y, n.children[-1].y],
                    color=c, lw=lw, solid_capstyle="butt", zorder=2)
        if (support_min is not None and n.support is not None
                and n.children and n.support >= support_min):
            # a halo in the surface colour, because in a dense clade the value
            # otherwise prints on top of a branch and becomes unreadable --
            # which is worse than not showing it, since a half-legible number
            # still looks like a number
            eff = None
            if support_halo:
                from matplotlib import patheffects
                eff = [patheffects.withStroke(linewidth=1.9,
                                              foreground=support_halo)]
            ax.annotate(support_fmt(n.support) if support_fmt else "%.0f" % n.support,
                        (n.x, n.y), textcoords="offset points", xytext=(-2.5, 1.8),
                        ha="right", va="bottom", fontsize=support_size,
                        color=support_color, zorder=5, annotation_clip=False,
                        path_effects=eff)


def mrca(root, names):
    """smallest clade containing every named tip; None if any name is absent"""
    want = set(names)
    best = None
    for n in root.walk():
        got = {t.name for t in n.tips()}
        if want <= got:
            best = n
    return best


def scale_bar(ax, x0, y0, length, label, color="#52514e", lw=0.9, size=7.0):
    ax.plot([x0, x0 + length], [y0, y0], color=color, lw=lw, zorder=3,
            clip_on=False)
    ax.annotate(label, (x0 + length / 2.0, y0), textcoords="offset points",
                xytext=(0, -3), ha="center", va="top", fontsize=size,
                color=color, zorder=3, annotation_clip=False)
