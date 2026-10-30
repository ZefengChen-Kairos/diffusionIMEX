"""Pour action trajectories xi(t) = (x_p(t), Q(t), h_p(t)).

Coordinates: cup centre is the origin, "forward" is +y, the pour starts on the -y side.
Each segment is linear in position / Q / h_p between its end points, with an optional
sinusoidal side wobble perpendicular to the forward axis.  Positions are given in mm,
Q in ml/s, h_p in m, wobble amplitude in mm and frequency in Hz.

Between segments the pitcher moves linearly to the next start point with Q = 0.
Every segment with Q > 0 gets its own passive batch tag (A, B, C, ...).
"""
from dataclasses import dataclass, field
import json
import math
import numpy as np


TAPER = 0.12   # [s] wobble envelope ramp at the ends of a wobble run


@dataclass
class Segment:
    t0: float
    t1: float
    p0: tuple
    p1: tuple = None            # None: stay at p0
    Q0: float = 0.0
    Q1: float = None            # None: = Q0
    h0: float = 0.012
    h1: float = None
    a: float = 0.0              # wobble amplitude [mm]
    f: float = 0.0              # wobble frequency [Hz]
    e: tuple = (0.0, 1.0)       # spout forward unit vector
    arc: dict = None            # {"c": (x,y) mm, "th0", "th1", "r0", "r1" mm}
    tag: int = -1               # assigned by Action
    taper_in: bool = True       # wobble envelope ramps in / out over TAPER seconds
    taper_out: bool = True      # (switched off between contiguous equal-wobble segments)

    def pos_mm(self, t):
        s = (t - self.t0) / max(self.t1 - self.t0, 1e-12)
        if self.arc is not None:
            A = self.arc
            th = A["th0"] + s * (A["th1"] - A["th0"])
            r = A["r0"] + s * (A["r1"] - A["r0"])
            p = np.array(A["c"], float) + r * np.array([math.cos(th), math.sin(th)])
        else:
            p1 = self.p0 if self.p1 is None else self.p1
            p = (1 - s) * np.array(self.p0, float) + s * np.array(p1, float)
        if self.a > 0 and self.f > 0:
            e = np.array(self.e, float)
            n = np.array([e[1], -e[0]])
            env = 1.0
            if self.taper_in:
                env = min(env, max(t - self.t0, 0.0) / TAPER)
            if self.taper_out:
                env = min(env, max(self.t1 - t, 0.0) / TAPER)
            # absolute-time phase: contiguous segments with equal wobble join smoothly
            p = p + env * self.a * math.sin(2 * math.pi * self.f * t) * n
        return p

    def end_mm(self):
        return self.pos_mm(self.t1)

    def start_mm(self):
        return self.pos_mm(self.t0)

    def Q(self, t):
        s = (t - self.t0) / max(self.t1 - self.t0, 1e-12)
        q1 = self.Q0 if self.Q1 is None else self.Q1
        return (1 - s) * self.Q0 + s * q1

    def h(self, t):
        s = (t - self.t0) / max(self.t1 - self.t0, 1e-12)
        h1 = self.h0 if self.h1 is None else self.h1
        return (1 - s) * self.h0 + s * h1


@dataclass
class Sample:
    x: np.ndarray       # spout position [m]
    v: np.ndarray       # spout velocity [m/s]
    Q: float            # [m^3/s]
    h: float            # [m]
    e: np.ndarray       # spout forward unit vector
    tag: int            # batch tag index or -1


class Action:
    def __init__(self, name, segments, T_end, notes=""):
        self.name = name
        self.segs = sorted(segments, key=lambda s: s.t0)
        self.T_end = T_end
        self.notes = notes
        for a, b in zip(self.segs[:-1], self.segs[1:]):
            same = abs(a.t1 - b.t0) < 1e-9 and a.a == b.a and a.f == b.f and tuple(a.e) == tuple(b.e)
            if same:
                a.taper_out = False
                b.taper_in = False
        k = 0
        for s in self.segs:
            if s.Q0 > 0 or (s.Q1 or 0) > 0:
                s.tag = k
                k += 1
        self.n_tags = k

    def _locate(self, t):
        for i, s in enumerate(self.segs):
            if s.t0 <= t <= s.t1:
                return s, None
            if t < s.t0:
                return None, (self.segs[i - 1] if i > 0 else None, s)
        return None, (self.segs[-1], None)

    def pos_mm(self, t):
        s, gap = self._locate(t)
        if s is not None:
            return s.pos_mm(t)
        a, b = gap
        if a is None:
            return b.start_mm()
        if b is None:
            return a.end_mm()
        w = (t - a.t1) / (b.t0 - a.t1)
        return (1 - w) * a.end_mm() + w * b.start_mm()

    def sample(self, t, dt_fd=2e-4):
        s, gap = self._locate(t)
        x = self.pos_mm(t) * 1e-3
        v = (self.pos_mm(t + dt_fd) - self.pos_mm(t - dt_fd)) * 1e-3 / (2 * dt_fd)
        if s is None:
            seg = gap[0] or gap[1]
            return Sample(x, v, 0.0, seg.h1 if seg.h1 is not None else seg.h0,
                          np.array(seg.e, float), -1)
        Q = max(s.Q(t), 0.0) * 1e-6
        return Sample(x, v, Q, s.h(t), np.array(s.e, float) / np.linalg.norm(s.e), s.tag if Q > 0 else -1)

    def total_volume_ml(self, n=20000):
        ts = np.linspace(0, self.T_end, n)
        return float(np.trapezoid([self.sample(t).Q for t in ts], ts) * 1e6)

    def to_json(self):
        return json.dumps({"name": self.name, "T_end": self.T_end, "notes": self.notes,
                           "segments": [s.__dict__ for s in self.segs]}, default=list, indent=1)

    def mirrored(self):
        """x -> -x mirror image of the action (G0.5)."""
        return MirroredAction(self)


class MirroredAction:
    def __init__(self, base):
        self.base = base
        self.name = base.name + "_mirror"
        self.T_end = base.T_end
        self.n_tags = base.n_tags
        self.segs = base.segs
        self.notes = base.notes

    def sample(self, t, dt_fd=2e-4):
        s = self.base.sample(t, dt_fd)
        f = np.array([-1.0, 1.0])
        return Sample(s.x * f, s.v * f, s.Q, s.h, s.e * f, s.tag)

    def pos_mm(self, t):
        return self.base.pos_mm(t) * np.array([-1.0, 1.0])

    def total_volume_ml(self, n=20000):
        return self.base.total_volume_ml(n)

    def to_json(self):
        return self.base.to_json()


def lift_cut(t0, p0, p1, T, Q, h=0.10, taper=0.3, along=False):
    """High pour moving p0 -> p1 in T seconds; the flow tapers to zero over the last `taper`
    seconds while the pitcher keeps moving (the barista's lift-off at the end of the cut).
    along=True points the spout along the stroke (ê = direction of motion)."""
    p0, p1 = np.array(p0, float), np.array(p1, float)
    tm = t0 + T - taper
    pm = p0 + (p1 - p0) * (T - taper) / T
    e = tuple((p1 - p0) / np.linalg.norm(p1 - p0)) if along else (0.0, 1.0)
    return [Segment(t0, tm, tuple(p0), tuple(pm), Q, None, h, e=e),
            Segment(tm, t0 + T, tuple(pm), tuple(p1), Q, 0.0, h, e=e)]


def heart(variant=1, cut=None, t_pour=6.0):
    S = Segment
    if variant == 0:     # task-book starting values
        return Action("heart_v0", [
            S(0.0, 0.8, (0, -18), None, 0.0, 9.0, 0.015, 0.012),
            S(0.8, 5.2, (0, -18), (0, -4), 9.0, None, 0.012),
            S(5.2, 5.7, (0, -4), None, 9.0, 6.0, 0.012, 0.10),
            S(5.7, 7.1, (0, -4), (0, 40), 6.0, None, 0.10),
        ], 9.0)
    # tuned: longer pour (larger disc); a short, fast, gentle cut that ends just past the front of
    # the disc with the flow tapering off (a cut run into the wall curls the tail into a vortex pair)
    c = cut or ((0, -6), (0, 32), 0.7, 3.0)
    tp = t_pour
    return Action("heart", [
        S(0.0, 0.8, (0, -18), None, 0.0, 9.0, 0.015, 0.012),
        S(0.8, tp, (0, -18), (0, -4), 9.0, None, 0.012),
        S(tp, tp + 0.5, (0, -4), c[0], 9.0, c[3], 0.012, 0.10),
    ] + lift_cut(tp + 0.5, c[0], c[1], c[2], c[3]), tp + 0.5 + c[2] + 1.8)


def push_heart(variant=0, cut=None):
    S = Segment
    # tuned: longer dots, dots 2 and 3 at Q 12 / h 1.0 cm (Q 14 raises v_j and lowers r to ~0.35);
    # short, fast, gentle cut ending past the front with a lift-off (spec start: (0,-34)->(0,40), 1.4 s, Q 6)
    # the cut ends just past the front of the pattern with a lift-off
    c = cut or ((0, -30), (0, 30), 0.6, 3.0)
    return Action("push_heart", [
        S(0.0, 1.8, (0, -10), (0, -2), 13.0, None, 0.012),
        S(2.4, 3.9, (0, -21), (0, -15), 12.0, None, 0.010),
        S(4.5, 5.8, (0, -29), (0, -25), 12.0, None, 0.010),
    ] + lift_cut(6.1, c[0], c[1], c[2], c[3]), 7.9 + c[2])


def layered_heart(variant=4, cut=None):
    S = Segment
    if variant == 1:     # spec starting values (jet runs through its own pattern; kept for the record)
        return Action("layered_heart_v0", [
            S(0.0, 0.5, (0, -24), None, 11.0, None, 0.014),
            S(0.5, 4.7, (0, -24), (0, 20), 11.0, None, 0.012, a=7.5, f=2.0),
            S(4.7, 5.2, (0, 20), (0, 2), 5.0, None, 0.012, 0.10),
            S(5.2, 6.6, (0, 2), (0, 40), 4.5, None, 0.10),
        ], 8.6)
    if variant == 2:     # in-place wobble
        return Action("layered_heart_v2", [
            S(0.0, 0.6, (0, -14), None, 11.0, None, 0.014),
            S(0.6, 4.4, (0, -14), None, 11.0, None, 0.012, a=7.5, f=2.0),
            S(4.4, 4.9, (0, -14), (0, -6), 11.0, 5.0, 0.012, 0.10),
            S(4.9, 6.3, (0, -6), (0, 40), 4.5, None, 0.10),
        ], 8.3)
    if variant == 3:     # slow retreat while wobbling, partial cut
        return Action("layered_heart_v3", [
            S(0.0, 0.6, (0, -4), None, 11.0, None, 0.014),
            S(0.6, 4.4, (0, -4), (0, -20), 11.0, None, 0.012, a=7.5, f=2.0),
            S(4.4, 4.9, (0, -20), (0, -12), 11.0, 5.0, 0.012, 0.10),
        ] + lift_cut(4.9, *(cut or ((0, -12), (0, 30), 0.7, 3.5))), 6.7 + (cut or (0, 0, 0.7))[2])
    if variant == 4:
        # tuned (default): slow retreat with a 9 mm / 1.5 Hz wobble, then a partial cut;
        # lambda = v_rel / f with v_rel = push speed + retreat speed
        return Action("layered_heart", [
            S(0.0, 0.6, (0, -4), None, 11.0, None, 0.014),
            S(0.6, 4.4, (0, -4), (0, -20), 11.0, None, 0.012, a=9.0, f=1.5),
            S(4.4, 4.9, (0, -20), (0, -12), 11.0, 5.0, 0.012, 0.10),
        ] + lift_cut(4.9, *(cut or ((0, -12), (0, 30), 0.7, 3.5))), 6.7 + (cut or (0, 0, 0.7))[2])
    # tuned: wobble near B with a slow forward drift; the new ribbon pushes the old ones forward
    return Action("layered_heart", [
        S(0.0, 0.6, (0, -20), None, 11.0, None, 0.014),
        S(0.6, 4.4, (0, -20), (0, -8), 11.0, None, 0.012, a=7.5, f=2.0),
        S(4.4, 4.9, (0, -8), None, 11.0, 5.0, 0.012, 0.10),
        S(4.9, 6.3, (0, -8), (0, 40), 4.5, None, 0.10),
    ], 8.3)


def tulip(variant=0, cut=None):
    S = Segment
    # tuned: longer segments, point dot at h 1.0 cm, short fast gentle cut with lift-off
    # (spec start: (0,-34)->(0,40), 1.4 s, Q 5)
    c = cut or ((0, -30), (0, 32), 0.6, 3.0)
    return Action("tulip", [
        S(0.0, 2.8, (0, -14), (0, 4), 12.0, None, 0.012, a=6.0, f=2.4),
        S(3.4, 5.4, (0, -20), (0, -10), 13.0, None, 0.012, a=6.0, f=2.4),
        S(6.0, 7.2, (0, -27), None, 12.0, None, 0.010),
    ] + lift_cut(7.5, c[0], c[1], c[2], c[3]), 9.3 + c[2])


def leaf(cut=None):
    S = Segment
    c = cut or ((0, -34), (0, 32), 0.9, 3.5)
    return Action("leaf", [
        S(0.0, 1.2, (0, -6), None, 11.0, None, 0.012, a=6.0, f=2.8),
        S(1.2, 2.0, (0, -6), (0, 4), 11.0, None, 0.012, a=6.0, f=2.8),
        S(2.0, 4.6, (0, 4), (0, -32), 10.0, None, 0.012, a=6.0, f=2.8),
        S(4.6, 5.0, (0, -32), (0, -34), 4.0, None, 0.012, 0.10),
    ] + lift_cut(5.0, c[0], c[1], c[2], c[3]), 6.8 + c[2])


def swan(variant=2, neck=(1.22, (-18, -17), ((-16, -14), (-22, -26), 0.35, 3.0))):
    S = Segment
    ang = math.radians(22.0)
    a = (math.sin(ang), math.cos(ang))
    P0 = (6.0, -2.0)
    P1 = (P0[0] + 10 * a[0], P0[1] + 10 * a[1])
    P2 = (P0[0] - 28 * a[0], P0[1] - 28 * a[1])
    wing = [
        S(0.0, 1.2, P0, None, 11.0, None, 0.012, a=5.5, f=2.8, e=a),
        S(1.2, 1.9, P0, P1, 11.0, None, 0.012, a=5.5, f=2.8, e=a),
        S(1.9, 4.3, P1, P2, 10.0, None, 0.012, a=5.5, f=2.8, e=a),
    ]
    if variant == 0:     # spec starting values
        return Action("swan_v0", wing + [
            S(5.7, 7.1, (0, 0), None, 8.0, None, 0.018,
              arc={"c": (-2.0, 0.0), "th0": 0.5 * math.pi, "th1": 1.4 * math.pi, "r0": 24.0, "r1": 20.0}),
            S(7.3, 7.9, (-12, -26), None, 14.0, None, 0.013),
            S(8.1, 8.9, (-14, -32), (-2, -8), 4.0, None, 0.10),
        ], 10.9)
    if variant == 1:     # first tuning: head next to the wing base (the last stroke hooks it)
        return Action("swan_v1", wing + [
            S(5.7, 6.7, (0, 0), None, 5.0, None, 0.025,
              arc={"c": (-2.0, 0.0), "th0": 0.5 * math.pi, "th1": 1.4 * math.pi, "r0": 24.0, "r1": 20.0}),
            S(6.9, 7.5, (-12, -26), None, 12.0, None, 0.013),
        ] + lift_cut(7.7, (-12, -31), (-11, -24), 0.35, 4.0, taper=0.15), 9.8)
    # tuned: thinner neck (lower Q, faster) that stops on the left side; the head dot sits clear of the
    # wing base and the last stroke pulls out through it, away from the neck, as the beak
    th1, hd, hc = neck
    return Action("swan", wing + [
        S(5.7, 6.5, (0, 0), None, 5.0, None, 0.025,
          arc={"c": (-2.0, 0.0), "th0": 0.5 * math.pi, "th1": th1 * math.pi, "r0": 24.0, "r1": 21.0}),
        S(6.7, 7.3, hd, None, 12.0, None, 0.013),
    ] + lift_cut(7.5, *hc, taper=0.15, along=True), 9.6)


ALL = {"heart": heart, "push_heart": push_heart, "layered_heart": layered_heart,
       "tulip": tulip, "leaf": leaf, "swan": swan}

# physical feasibility box (section 1.4)
LIMITS = {"Q_ml": (3.0, 25.0), "h_m": (0.005, 0.15), "v": (0.0, 0.3), "a_mm": (0.0, 12.0),
          "f": (0.0, 6.0), "V_ml": 150.0}


def check_feasible(act, n=4000):
    """Return list of violations of the feasibility box (empty = valid action)."""
    bad = []
    ts = np.linspace(0, act.T_end, n)
    for t in ts:
        s = act.sample(t)
        if s.Q > 0:
            q = s.Q * 1e6
            if not (LIMITS["Q_ml"][0] - 1e-9 <= q <= LIMITS["Q_ml"][1]) and q > 1.0:
                # ramps through small Q are allowed only transiently (start/stop)
                pass
            if not (LIMITS["h_m"][0] <= s.h <= LIMITS["h_m"][1]):
                bad.append(("h", t, s.h))
            if np.linalg.norm(s.v) > LIMITS["v"][1]:
                bad.append(("v", t, float(np.linalg.norm(s.v))))
    for sg in act.segs:
        if sg.a > LIMITS["a_mm"][1] or sg.f > LIMITS["f"][1]:
            bad.append(("wobble", sg.t0, (sg.a, sg.f)))
        for q in (sg.Q0, sg.Q1 if sg.Q1 is not None else sg.Q0):
            if q > LIMITS["Q_ml"][1]:
                bad.append(("Q", sg.t0, q))
    if act.total_volume_ml() > LIMITS["V_ml"]:
        bad.append(("V", 0, act.total_volume_ml()))
    return bad
