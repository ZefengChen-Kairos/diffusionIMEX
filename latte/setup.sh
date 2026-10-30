#!/usr/bin/env bash
# Environment setup for AutoDL (single GPU, Linux, no sudo) or any Linux box.
#   bash setup.sh            # default PyPI
#   PIP_MIRROR=tuna bash setup.sh   # Tsinghua mirror (mainland China)
set -e
cd "$(dirname "$0")"
PIP="python3 -m pip"
IDX=""
case "${PIP_MIRROR:-}" in
  tuna)   IDX="-i https://pypi.tuna.tsinghua.edu.cn/simple" ;;
  aliyun) IDX="-i https://mirrors.aliyun.com/pypi/simple" ;;
esac
$PIP install $IDX -r requirements.txt || $PIP install -i https://pypi.tuna.tsinghua.edu.cn/simple -r requirements.txt
# phase-two GPU backend (optional now): pinned Warp; the first run JIT-compiles kernels (~1 min, normal)
if [ "${WITH_WARP:-0}" = "1" ]; then
  $PIP install $IDX "warp-lang==1.5.1" || echo "WARNING: warp-lang install failed; CPU backends still work"
fi
OUT=/root/autodl-tmp/latte_out
[ -d /root/autodl-tmp ] || OUT="$(pwd)/out"
mkdir -p "$OUT"
echo "outputs go to $OUT"
python3 - <<'PY'
import subprocess
try:
    print(subprocess.run(["nvidia-smi", "--query-gpu=name,memory.total,driver_version", "--format=csv,noheader"],
                         capture_output=True, text=True).stdout.strip() or "nvidia-smi: no GPU listed")
except FileNotFoundError:
    print("WARNING: nvidia-smi not found: no GPU. Phase one runs on the CPU (numba backend).")
try:
    import warp as wp
    wp.init()
    print("warp devices:", wp.get_devices())
except Exception as e:
    print("warp not available (%s): phase-one CPU backends only" % e.__class__.__name__)
import numba
print("numba", numba.__version__, "ok")
PY
