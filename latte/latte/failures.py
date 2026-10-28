"""Structured failure codes (V3 Sec. 15)."""

CODES = (
    "INPUT_INVALID", "GEOMETRY_INVALID", "RESOLUTION_INSUFFICIENT",
    "PRESSURE_FAILED", "VISCOSITY_FAILED", "CONSERVATION_FAILED",
    "BOUNDS_FAILED", "TIME_MAP_UNKNOWN", "OUT_OF_SCOPE",
)

# Codes that may be retried with a smaller time step; everything else is a hard stop
# (missing input, invalid geometry and unsealed physics are never retried).
RETRYABLE = ("PRESSURE_FAILED", "VISCOSITY_FAILED", "CONSERVATION_FAILED", "BOUNDS_FAILED")


class SimFailure(Exception):
    def __init__(self, code, reason, **context):
        if code not in CODES:
            raise ValueError(f"unknown failure code {code}")
        super().__init__(f"{code}: {reason}")
        self.code = code
        self.reason = reason
        self.context = context

    @property
    def retryable(self):
        return self.code in RETRYABLE
