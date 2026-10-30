"""Pour actions xi(t) = (x_p(t), Q(t), h_p(t)) that a milk pitcher can actually pour.

A pour is written as a sequence of moves of one pitcher (`Pitcher`).  Every channel -- spout
position x, y [mm], flow Q [ml/s], spout height h [m], wobble amplitude a [mm] and frequency
f [Hz], spout direction th [rad from +y] -- is a keyframe track interpolated with a monotone C1
cubic (PCHIP): the pitcher never teleports, the flow never jumps, the spout is lifted or lowered
at a finite rate, and a channel held between two keys stays exactly constant (no overshoot).
The wobble is a lateral sine a(t) sin(phi(t)), phi' = 2 pi f(t), perpendicular to the spout.

Coordinates: cup centre is the origin, "forward" (B -> C -> D in the videos) is +y, the pitcher
stands on the -y side.  The pattern pours follow the on-screen instructions and the timing of the
reference videos: the stream is never interrupted between the first touch-down and the final
lift-off.  Between strokes the barista lifts the pitcher and thins the stream while moving back,
then lowers it again; every pattern ends with a pull-through at low flow with the spout lifted,
the flow running out while the pitcher still moves.  Each stroke (and each lift between strokes)
carries its own passive batch tag.
"""
from dataclasses import dataclass
import json
import math
import numpy as np
from scipy.interpolate import PchipInterpolator


CH = ("x", "y", "Q", "h", "a", "f", "th")
DT_GRID = 1e-3          # [s] tabulation step of the tracks

H_DRAW = 0.010          # spout just above the surface while drawing [m]
H_LIFT = 0.050          # lifted spout between strokes and in the pull-through [m]
Q_THIN = 3.0            # thin stream while lifted [ml/s] (task-book minimum)
Q_PULL = 4.0            # pull-through flow [ml/s]


@dataclass
class Sample:
    x: np.ndarray       # spout position [m]
    v: np.ndarray       # spout velocity [m/s]
    Q: float            # [m^3/s]
    h: float            # [m]
    e: np.ndarray       # spout forward unit vector
    tag: int            # batch tag index or -1


class Pitcher:
    """Builder: `move(T, to=..., Q=..., h=..., a=..., f=..., th=...)` brings the named channels to
    their new values in T seconds; channels not named are held.  `stroke(label)` starts a new
    batch tag, `mark(label)` records a time for the gates."""

    def __init__(self, name, at, h=H_LIFT, Q=Q_THIN, th=0.0):
        self.name = name
        self.t = 0.0
        self.k = {c: [] for c in CH}
        self._key(0.0, dict(x=float(at[0]), y=float(at[1]), Q=float(Q), h=float(h), a=0.0, f=0.0, th=float(th)))
        self.strokes = []
        self.marks = {}
        self.meta = {}

    def last(self, c):
        return self.k[c][-1][1]

    def at(self):
        return np.array([self.last("x"), self.last("y")])

    def _key(self, t, vals):
        for c in CH:
            v = float(vals[c]) if c in vals else self.last(c)
            if self.k[c] and abs(self.k[c][-1][0] - t) < 1e-9:
                self.k[c][-1] = (t, v)
            else:
                self.k[c].append((t, v))

    def move(self, T, to=None, **ch):
        bad = set(ch) - set(CH)
        if bad:
            raise ValueError(f"unknown channels {bad}")
        if to is not None:
            ch["x"], ch["y"] = float(to[0]), float(to[1])
        self.t += float(T)
        self._key(self.t, ch)
        return self

    def path(self, T, pts, **ch):
        """Run through the points `pts` in T seconds (time proportional to arc length)."""
        pts = [self.at()] + [np.asarray(p, float) for p in pts]
        L = np.cumsum([0.0] + [np.linalg.norm(b - a) for a, b in zip(pts[:-1], pts[1:])])
        t0 = self.t
        for p, l in zip(pts[1:], L[1:]):
            self.t = t0 + T * l / L[-1]
            self._key(self.t, dict(x=p[0], y=p[1], **(ch if l == L[-1] else {})))
        return self

    def stroke(self, label):
        self.strokes.append([label, self.t])
        return self

    def mark(self, label):
        self.marks[label] = self.t
        return self

    def build(self, tail=1.6, T_end=None, continuous=True, notes=""):
        return Pour(self.name, self.k, self.strokes or [["A", 0.0]], self.marks,
                    T_end if T_end is not None else self.t + tail, self.meta, continuous, notes)


# --------------------------------------------------------------------------- moves shared by the patterns
def touch_down(p, Q, T=0.5, h=H_DRAW):
    """Lower the spout onto the surface while opening the flow."""
    return p.move(T, Q=Q, h=h)


def lift_back(p, back, Q, T=1.5, h=H_LIFT):
    """Between two strokes: lift the spout and thin the stream while moving back, lower again at
    `back` and open the flow to Q (the stream is not interrupted)."""
    x0, x1 = p.at(), np.asarray(back, float)
    p.move(0.30 * T, to=x0 + 0.25 * (x1 - x0), Q=Q_THIN, h=h)
    p.move(0.40 * T, to=x0 + 0.90 * (x1 - x0))
    return p.move(0.30 * T, to=x1, Q=Q, h=H_DRAW)


def pull_through(p, to, T, Q=Q_PULL, h=H_LIFT, T_lift=0.4, T_off=0.3):
    """Finishing stroke: thin the stream and lift while starting forward, cross the pattern, and
    let the flow run out while the pitcher still moves and rises (the lift-off)."""
    x0, x1 = p.at(), np.asarray(to, float)
    d = x1 - x0
    p.move(T_lift, to=x0 + d * 0.6 * T_lift / T, Q=Q, h=h)
    p.move(T - T_lift - T_off, to=x0 + d * (1 - T_off / T))
    p.move(T_off, to=x1, Q=0.0, h=h + 0.02)
    return p.move(0.4, to=x1 + 0.1 * d)          # the hand coasts to rest after the stream is off


# --------------------------------------------------------------------------- built pour
class Pour:
    def __init__(self, name, keys, strokes, marks, T_end, meta=None, continuous=True, notes=""):
        self.name = name
        self.keys = {c: [tuple(k) for k in keys[c]] for c in CH}
        self.strokes = [list(s) for s in strokes]
        self.marks = dict(marks)
        self.T_end = float(T_end)
        self.meta = dict(meta or {})
        self.continuous = continuous
        self.notes = notes
        self.n_tags = len(self.strokes)
        self._tabulate()

    def _tabulate(self):
        n = int(math.ceil(self.T_end / DT_GRID)) + 3
        self.tg = np.arange(n) * DT_GRID
        self.val, self.der = {}, {}
        for c in CH:
            kt = np.array([k[0] for k in self.keys[c]])
            kv = np.array([k[1] for k in self.keys[c]])
            if len(kt) == 1:
                self.val[c], self.der[c] = np.full(n, kv[0]), np.zeros(n)
                continue
            ip = PchipInterpolator(kt, kv)
            tt = np.clip(self.tg, kt[0], kt[-1])
            self.val[c] = ip(tt)
            d = ip.derivative()(tt)
            d[(self.tg < kt[0]) | (self.tg > kt[-1])] = 0.0
            self.der[c] = d
        self.val["Q"] = np.maximum(self.val["Q"], 0.0)
        f = self.val["f"]
        self.phi = 2 * math.pi * np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * DT_GRID)])
        self.t_str = np.array([s[1] for s in self.strokes])
        self.t_flow_end = float(self.tg[np.nonzero(self.val["Q"] > 0)[0][-1]]) if (self.val["Q"] > 0).any() else 0.0

    # ---- evaluation
    def _w(self, t):
        u = min(max(t, 0.0), self.tg[-2]) / DT_GRID
        i = int(u)
        return i, u - i

    def _g(self, arr, i, w):
        return arr[i] + w * (arr[i + 1] - arr[i])

    def tag_at(self, t):
        return max(int(np.searchsorted(self.t_str, t, side="right")) - 1, 0)

    def sample(self, t, dt_fd=None):
        i, w = self._w(t)
        g = lambda c: self._g(self.val[c], i, w)
        d = lambda c: self._g(self.der[c], i, w)
        th = g("th")
        e = np.array([math.sin(th), math.cos(th)])
        nrm = np.array([e[1], -e[0]])
        a, phi = g("a"), self._g(self.phi, i, w)
        wob = a * math.sin(phi)
        dwob = d("a") * math.sin(phi) + a * 2 * math.pi * g("f") * math.cos(phi)
        dnrm = np.array([-e[0], -e[1]]) * d("th")
        x = (np.array([g("x"), g("y")]) + wob * nrm) * 1e-3
        v = (np.array([d("x"), d("y")]) + dwob * nrm + wob * dnrm) * 1e-3
        Q = max(g("Q"), 0.0) * 1e-6
        return Sample(x, v, Q, g("h"), e, self.tag_at(t) if Q > 0 else -1)

    def pos_mm(self, t):
        return self.sample(t).x * 1e3

    def total_volume_ml(self):
        return float(np.trapezoid(self.val["Q"], self.tg))

    # ---- strokes / tags
    def tag(self, label):
        return [s[0] for s in self.strokes].index(label)

    def start(self, label):
        return self.strokes[self.tag(label)][1]

    def end(self, label):
        k = self.tag(label)
        return self.strokes[k + 1][1] if k + 1 < len(self.strokes) else self.t_flow_end

    def labels(self):
        return [s[0] for s in self.strokes]

    def to_json(self):
        return json.dumps({"name": self.name, "T_end": self.T_end, "notes": self.notes, "strokes": self.strokes,
                           "marks": self.marks, "meta": self.meta, "keys": self.keys}, indent=1, default=float)

    def mirrored(self):
        """x -> -x mirror image of the action (G0.5)."""
        return MirroredAction(self)

    # ---- variants used by the gates (same hand motion, one thing changed)
    def _with(self, name=None, keys=None, strokes=None, marks=None, T_end=None, meta=None):
        return Pour(name or self.name, keys or self.keys, strokes or self.strokes,
                    self.marks if marks is None else marks, T_end or self.T_end, meta or self.meta,
                    self.continuous, self.notes)

    def with_flow(self, scale, name=None):
        k = dict(self.keys)
        k["Q"] = [(t, v * scale) for t, v in self.keys["Q"]]
        return self._with(name or f"{self.name}_Qx{scale:g}", keys=k)

    def with_height_cap(self, hmax, name=None):
        k = dict(self.keys)
        k["h"] = [(t, min(v, hmax)) for t, v in self.keys["h"]]
        return self._with(name or f"{self.name}_hcap", keys=k)

    def time_warped(self, t0, t1, factor, name=None):
        """The interval [t0, t1] runs `factor` times faster; everything after it moves earlier."""
        sh = (t1 - t0) * (1 - 1 / factor)

        def m(t):
            return t if t <= t0 else (t0 + (t - t0) / factor if t <= t1 else t - sh)
        k = {c: [(m(t), v) for t, v in self.keys[c]] for c in CH}
        return self._with(name or f"{self.name}_warp", keys=k, strokes=[[l, m(t)] for l, t in self.strokes],
                          marks={l: m(t) for l, t in self.marks.items()}, T_end=self.T_end - sh)

    def perturbed(self, rng, eps=0.10, name=None):
        """+-eps random scaling of the timing, the flow and the positions."""
        ts, q, sc = 1 + rng.uniform(-eps, eps), 1 + rng.uniform(-eps, eps), 1 + rng.uniform(-eps, eps)
        k = {c: [(t * ts, v * (q if c == "Q" else sc if c in ("x", "y") else 1.0)) for t, v in self.keys[c]]
             for c in CH}
        return self._with(name or f"{self.name}_pert", keys=k, strokes=[[l, t * ts] for l, t in self.strokes],
                          marks={l: t * ts for l, t in self.marks.items()}, T_end=self.T_end * ts)

    def truncated(self, t_stop, T_end, ramp=0.45, name=None):
        """Same pour up to t_stop, then the flow runs out in `ramp` seconds with the pitcher held."""
        i, w = self._w(t_stop)
        k = {}
        for c in CH:
            k[c] = [(t, v) for t, v in self.keys[c] if t < t_stop - 1e-9] + [(t_stop, self._g(self.val[c], i, w))]
        k["Q"].append((t_stop + ramp, 0.0))
        return self._with(name or f"{self.name}_stop", keys=k, strokes=[s for s in self.strokes if s[1] < t_stop],
                          marks={l: t for l, t in self.marks.items() if t <= t_stop}, T_end=T_end)


class MirroredAction:
    def __init__(self, base):
        self.base = base
        self.name = base.name + "_mirror"
        self.T_end = base.T_end
        self.n_tags = base.n_tags
        self.strokes, self.marks, self.meta, self.notes = base.strokes, base.marks, base.meta, base.notes

    def sample(self, t, dt_fd=None):
        s = self.base.sample(t)
        f = np.array([-1.0, 1.0])
        return Sample(s.x * f, s.v * f, s.Q, s.h, s.e * f, s.tag)

    def pos_mm(self, t):
        return self.base.pos_mm(t) * np.array([-1.0, 1.0])

    def total_volume_ml(self):
        return self.base.total_volume_ml()

    def to_json(self):
        return self.base.to_json()


# --------------------------------------------------------------------------- simple pours for the mechanism gates
def dwell(Q, T_inj=3.0, T_end=None, x=(0, -10), h=0.012, ramp=None, name=None):
    """Pour at a fixed point: the flow opens and closes within the pitcher's rate limit."""
    ramp = ramp or max(0.25, 1.6 * Q / LIMITS["dQdt"])
    p = Pitcher(name or f"dwell_Q{Q:g}", at=x, h=h, Q=0.0)
    p.move(ramp, Q=Q).move(T_inj - 2 * ramp).move(ramp, Q=0.0)
    return p.build(T_end=T_end or T_inj, continuous=False)


# --------------------------------------------------------------------------- the six patterns
def heart(cut_to=30.0):
    """Video: inject at B and level the cup (0-3 s), the foam forms a half-disc (3-4.5 s), move the
    pitcher slowly to the centre (4.5-6.5 s), lower the flow and lift to finish (6.5-8.5 s)."""
    p = Pitcher("heart", at=(0, -20))
    p.stroke("A")
    touch_down(p, Q=10.0)
    p.move(3.5)                                   # pour at B: the white half-disc
    p.move(2.5, to=(0, -4))                       # slowly towards the centre
    p.stroke("cut")
    pull_through(p, to=(0, cut_to), T=1.6)
    return p.build(notes="continuous pour; lift-off pull-through at h 5 cm, Q 4 ml/s")


def push_heart():
    """Video: 1st part at B, gently pushed towards D (0-4 s); raise the spout and thin the flow
    (4.5-5.5 s); 2nd part back at B, large flow pushes to C (6-9.5 s); back to B again (10-11 s);
    large flow for a heart (11.5-13 s); small flow lifted to finish (13.5-14.5 s)."""
    p = Pitcher("push_heart", at=(0, -22))
    p.stroke("A")
    touch_down(p, Q=10.0)
    p.move(1.5)                                   # first layer at B
    p.move(2.0, to=(0, -16), Q=11.0)              # keep pouring: the pattern is pushed towards D
    p.stroke("lift1")
    lift_back(p, back=(0, -24), Q=11.0)
    p.stroke("B")
    p.move(1.0)
    p.move(2.4, to=(0, -19), Q=12.0)              # large flow pushes the layers to C
    p.stroke("lift2")
    lift_back(p, back=(0, -26), Q=11.0)
    p.stroke("C")
    p.move(1.6)                                   # the last heart at B
    p.stroke("cut")
    pull_through(p, to=(0, 30), T=1.5)
    return p.build(notes="three parts at B separated by lifts with a thin stream; lift-off pull-through")


def layered_heart():
    """Video: inject at B (0-1.5 s); wobble once white lines appear (1.5-4 s); when the wobbled
    lines form a half-disc, move forward (4.5-7.5 s); lift with a small flow to finish (7.5-10 s)."""
    p = Pitcher("layered_heart", at=(0, -20))
    p.stroke("A")
    touch_down(p, Q=11.0)
    p.move(1.0)                                   # white lines appear
    p.move(0.4, a=8.0, f=2.5)                     # start the wobble
    p.move(2.4)                                   # wobble at B: half-disc of ribbons
    p.move(3.0, to=(0, -6))                       # move forward while wobbling
    p.move(0.3, a=0.0)
    p.stroke("cut")
    pull_through(p, to=(0, 30), T=1.6)
    return p.build(notes="wobble at B, forward while wobbling, lift-off pull-through")


def tulip():
    """Video: wobble at B once white lines appear (0-3.5 s), move forward (4-4.5 s), let the lines
    wrap back (5-6 s); 2nd part back at B (6.5-8.5 s), wobble while pushing forward (9-10 s), press
    the base pattern open (10.5-11.5 s); 3rd part back at B (12-14 s), a heart (14.5-15.5 s), lift."""
    p = Pitcher("tulip", at=(0, -20))
    p.stroke("A")
    touch_down(p, Q=10.0)
    p.move(0.4, a=7.0, f=3.0)
    p.move(2.4)                                   # wobble at B
    p.move(1.2, to=(0, -8))                       # forward while wobbling
    p.move(0.3, a=0.0)
    p.move(1.0)                                   # keep pouring: the lines wrap back
    p.stroke("lift1")
    lift_back(p, back=(0, -22), Q=10.0)
    p.stroke("B")
    p.move(0.6)
    p.move(0.3, a=6.0, f=3.0)
    p.move(1.4, to=(0, -14))                      # wobble while pushing forward
    p.move(0.3, a=0.0)
    p.move(1.0, Q=12.0)                           # large flow presses the base pattern open
    p.stroke("lift2")
    lift_back(p, back=(0, -24), Q=11.0)
    p.stroke("C")
    p.move(1.4)                                   # the heart
    p.stroke("cut")
    pull_through(p, to=(0, 30), T=1.4)
    return p.build(notes="three parts, lifts with a thin stream between them; lift-off pull-through")


def leaf():
    """Video: inject at B (0-1 s), keep the flow and start the wobble (1-2.5 s), two wobbles forward
    once the half-disc forms (3-5 s), retreat quickly to the cup edge (5.5-8.5 s), lower the flow
    and lift to finish (9-10.5 s)."""
    p = Pitcher("leaf", at=(0, -18))
    p.stroke("A")
    touch_down(p, Q=11.0)
    p.move(0.4, a=6.0, f=3.5)
    p.move(1.6)                                   # wobble at B
    p.move(1.0, to=(0, -8))                       # two wobbles forward
    p.mark("retreat0")
    p.move(3.0, to=(0, -36), Q=10.0)              # retreat to the edge while wobbling
    p.mark("retreat1")
    p.move(0.3, a=0.0)
    p.stroke("cut")
    pull_through(p, to=(0, 32), T=1.6)
    return p.build(notes="wobble at B, forward, retreat to the edge while wobbling, lift-off pull-through")


def swan(ang_deg=22.0):
    """Video: the wing is a leaf with the cup turned slightly (0-9 s); small flow, lifted, finish
    along the edge of the pattern (9.5-12 s); lower the pitcher again and move quickly back towards
    B (12.5-16 s, the neck); large flow for a white heart, lift to finish (16.5-17 s, head + beak)."""
    ang = math.radians(ang_deg)
    ax = np.array([math.sin(ang), math.cos(ang)])
    P0 = np.array([5.0, -4.0])
    P1, P2 = P0 + 8.0 * ax, P0 - 26.0 * ax
    c, r0, r1, th0, th1 = np.array([-2.0, 0.0]), 24.0, 21.0, 0.5 * math.pi, 1.22 * math.pi
    p = Pitcher("swan", at=P0)
    p.stroke("wing")
    p.move(0.5, Q=11.0, h=H_DRAW, th=ang)         # touch down with the cup turned by ang
    p.move(0.4, a=5.5, f=3.5)
    p.move(1.4)                                   # wobble at the start point
    p.move(0.8, to=P1)                            # forward
    p.move(2.6, to=P2, Q=10.0)                    # retreat to the edge while wobbling
    p.move(0.3, a=0.0)
    p.stroke("edge")                              # lifted thin stream along the wing's edge
    n0 = c + r0 * np.array([math.cos(th0), math.sin(th0)])
    p.move(0.4, to=P2 + np.array([-4.0, 2.0]), Q=Q_THIN, h=H_LIFT, th=0.0)
    p.path(1.4, [(-13.0, -8.0), (-10.0, 14.0), n0 + np.array([2.0, -2.0])])
    p.stroke("neck")
    p.move(0.4, to=n0, Q=6.0, h=0.015)            # lower the pitcher again
    ths = np.linspace(th0, th1, 7)[1:]
    p.path(0.9, [c + (r0 + (r1 - r0) * (t - th0) / (th1 - th0)) * np.array([math.cos(t), math.sin(t)])
                 for t in ths])                   # quickly back towards B: the neck
    hd = p.at()
    p.stroke("head")
    p.move(0.3, Q=12.0, h=H_DRAW)
    p.move(0.8)                                   # large flow: the head
    p.stroke("beak")
    p.move(0.35, to=hd + np.array([0.0, 4.0]), Q=Q_PULL, h=0.04)
    p.move(0.25, to=hd + np.array([0.0, 10.0]), Q=0.0, h=0.06)
    p.move(0.4, to=hd + np.array([0.0, 12.0]))
    p.meta.update(wing_P0=P0.tolist(), wing_ang_deg=ang_deg, wing_s=(-28.0, 10.0),
                  neck_len_mm=float(abs(th1 - th0) * 0.5 * (r0 + r1)))
    return p.build(tail=2.0, notes="leaf wing with the cup turned, lifted edge run, neck, head, lift-off")


ALL = {"heart": heart, "push_heart": push_heart, "layered_heart": layered_heart,
       "tulip": tulip, "leaf": leaf, "swan": swan}

# physical feasibility box (section 1.4) and the hand limits of a real pitcher
LIMITS = {"Q_ml": (3.0, 25.0), "h_m": (0.005, 0.15), "v": (0.0, 0.3), "a_mm": (0.0, 12.0),
          "f": (0.0, 6.0), "V_ml": 150.0,
          "dQdt": 40.0,     # [ml/s^2]  thin stream (3) -> full pour (12) takes >= ~0.25 s of tilting
          "dhdt": 0.25,     # [m/s]     the spout is lifted / lowered by 5 cm in >= 0.2 s
          "acc": 3.0,       # [m/s^2]   acceleration of the guided path (the wobble is bounded by a, f)
          "ramp": 0.6}      # [s]       flow below Q_min only while opening / closing the stream


def check_feasible(act, dt=2e-3):
    """Return the list of violations (empty = a pour a pitcher can make)."""
    act = getattr(act, "base", act)
    bad = []
    tg = np.arange(0.0, act.T_end, dt)
    idx = np.minimum((tg / DT_GRID).round().astype(int), len(act.tg) - 1)
    Q, h = act.val["Q"][idx], act.val["h"][idx]
    on = Q > 1e-9
    if Q.max() > LIMITS["Q_ml"][1]:
        bad.append(("Q_max", float(Q.max())))
    # contiguous pouring intervals; the stream is thin (< Q_min) only right after opening / before closing
    edges = np.diff(np.concatenate([[0], on.astype(int), [0]]))
    runs = list(zip(np.nonzero(edges == 1)[0], np.nonzero(edges == -1)[0]))
    if act.continuous and len(runs) > 1:
        bad.append(("flow interrupted", [(round(tg[a], 2), round(tg[b - 1], 2)) for a, b in runs]))
    for a, b in runs:
        ok = np.nonzero(Q[a:b] >= LIMITS["Q_ml"][0] - 1e-9)[0]
        if len(ok) == 0:
            bad.append(("Q below minimum", float(tg[a])))
            continue
        lo, hi = a + ok[0], a + ok[-1]
        if (Q[lo:hi + 1] < LIMITS["Q_ml"][0] - 1e-6).any():
            bad.append(("Q dips below minimum", float(tg[lo + np.argmin(Q[lo:hi + 1])])))
        if (lo - a) * dt > LIMITS["ramp"] or (b - 1 - hi) * dt > LIMITS["ramp"]:
            bad.append(("slow open / close", float(tg[a])))
    hon = h[on]
    if hon.size and (hon.min() < LIMITS["h_m"][0] or hon.max() > LIMITS["h_m"][1]):
        bad.append(("h", float(hon.min()), float(hon.max())))
    for c, lim in (("Q", "dQdt"), ("h", "dhdt")):
        r = float(np.abs(act.der[c]).max())
        if r > LIMITS[lim] * (1 + 1e-9):
            bad.append((lim, r))
    acc = np.hypot(np.gradient(act.der["x"], DT_GRID), np.gradient(act.der["y"], DT_GRID))[idx][on] * 1e-3
    if acc.size and acc.max() > LIMITS["acc"]:
        bad.append(("acc", float(acc.max()), float(tg[on][np.argmax(acc)])))
    vmax = max((float(np.linalg.norm(act.sample(t).v)) for t in tg[on]), default=0.0)
    if vmax > LIMITS["v"][1]:
        bad.append(("v", vmax))
    if act.val["a"].max() > LIMITS["a_mm"][1] or act.val["f"].max() > LIMITS["f"][1]:
        bad.append(("wobble", float(act.val["a"].max()), float(act.val["f"].max())))
    if act.total_volume_ml() > LIMITS["V_ml"]:
        bad.append(("V", act.total_volume_ml()))
    return bad


def describe(act):
    """One-line summary of a pour (duration, volume, peak rates) for the reports."""
    act = getattr(act, "base", act)
    on = act.val["Q"] > 1e-9
    acc = np.hypot(np.gradient(act.der["x"], DT_GRID), np.gradient(act.der["y"], DT_GRID))[on] * 1e-3
    return dict(T_pour=float(act.tg[on][-1] - act.tg[on][0]) if on.any() else 0.0,
                V_ml=round(act.total_volume_ml(), 1), Q_max=round(float(act.val["Q"].max()), 2),
                Q_min_on=round(float(act.val["Q"][on].min()), 2) if on.any() else 0.0,
                dQdt_max=round(float(np.abs(act.der["Q"]).max()), 1),
                dhdt_max=round(float(np.abs(act.der["h"]).max()), 3),
                acc_max=round(float(acc.max()), 2), strokes=act.labels())
