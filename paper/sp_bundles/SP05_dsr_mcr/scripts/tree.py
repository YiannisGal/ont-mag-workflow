"""Minimal Newick reader for IQ-TREE output: leaves, patristic distances, node supports, bipartitions."""
import re
class T:
    def __init__(self, s):
        s = s.strip().rstrip(';'); self.par, self.bl, self.name, self.sup = [], [], [], []
        stack = []; cur = self._new(None); self.root = cur; i = 0
        while i < len(s):
            c = s[i]
            if c == '(':
                n = self._new(cur); stack.append(cur); cur = n; i += 1
            elif c == ',':
                cur = self._new(stack[-1]); i += 1
            elif c == ')':
                cur = stack.pop(); i += 1
            else:
                m = re.match(r"([^:,();]*)(?::([-\d.eE]+))?", s[i:])
                lab = m.group(1)
                if lab:
                    if self._children(cur): self.sup[cur] = lab
                    else: self.name[cur] = lab
                if m.group(2): self.bl[cur] = float(m.group(2))
                i += max(len(m.group(0)), 1)
        self.kids = {}
        for n, p in enumerate(self.par):
            if p is not None: self.kids.setdefault(p, []).append(n)
    def _new(self, p):
        self.par.append(p); self.bl.append(0.0); self.name.append(None); self.sup.append(None); return len(self.par) - 1
    def _children(self, n): return any(p == n for p in self.par)
    def leaves(self): return {self.name[i]: i for i in range(len(self.par)) if self.name[i]}
    def anc(self, x):
        d, t = {}, 0.0
        while x is not None: d[x] = t; t += self.bl[x]; x = self.par[x]
        return d
    def dist(self, a, b):
        A = self.anc(a); y, t = b, 0.0
        while y not in A: t += self.bl[y]; y = self.par[y]
        return t + A[y], y
    def below(self, n):
        out, st = set(), [n]
        while st:
            x = st.pop()
            if self.name[x]: out.add(self.name[x])
            st.extend(self.kids.get(x, []))
        return out
    def splits(self):
        """every internal edge as (leafset_below, support)"""
        return [(self.below(n), self.sup[n]) for n in range(len(self.par)) if n != self.root and self.kids.get(n)]
