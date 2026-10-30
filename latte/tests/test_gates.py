"""pytest entry points.  Fast: G0 (no refinement run).  Full gates: LATTE_FULL=1 pytest -q tests."""
import os
import sys
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FULL = os.environ.get("LATTE_FULL") == "1"


def _assert(R):
    bad = [c for c in R.checks if not c["PASS"]]
    assert not bad, "failed checks: " + "; ".join(f"{c['id']} {c['desc']} = {c['value']}" for c in bad)


def test_g0_fast():
    from gates.g0 import run_g0
    _assert(run_g0(full=False))


@pytest.mark.skipif(not FULL, reason="set LATTE_FULL=1 for the full gate suite")
@pytest.mark.parametrize("gate", ["G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "CF"])
def test_gate(gate):
    from gates.run_gates import run
    _assert(run(gate))
