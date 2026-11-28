#!/bin/bash
# 编译并在后台跑,自动定线程数,脱离终端(关掉 ssh 也继续)。
#   ./run.sh              # 默认:论文 §4.6 的 1000x1000 blast
#   ./run.sh tend=1e-4    # 先试跑约 50 步
# 传给脚本的参数原样转给程序。
set -eu
SRC=Final_emb_omp_clean.cpp
BIN=mhd
LOG=run_$(date +%Y%m%d_%H%M%S).log

# ---- 可用核数:优先读 cgroup 配额,nproc 在容器里常常报宿主机的 ----
cores() {
  local q p
  if [ -r /sys/fs/cgroup/cpu.max ]; then                       # cgroup v2
    read -r q p < /sys/fs/cgroup/cpu.max || true
    [ "${q:-max}" != "max" ] && [ "${p:-0}" -gt 0 ] && { echo $(( (q + p - 1) / p )); return; }
  fi
  if [ -r /sys/fs/cgroup/cpu/cpu.cfs_quota_us ]; then          # cgroup v1
    q=$(cat /sys/fs/cgroup/cpu/cpu.cfs_quota_us); p=$(cat /sys/fs/cgroup/cpu/cpu.cfs_period_us)
    [ "$q" -gt 0 ] && [ "$p" -gt 0 ] && { echo $(( (q + p - 1) / p )); return; }
  fi
  nproc 2>/dev/null || echo 1
}
N=$(cores); [ "$N" -lt 1 ] && N=1

# ---- 编译(源码没变就跳过)----
if [ ! -x "$BIN" ] || [ "$SRC" -nt "$BIN" ]; then
  EIG=$(for d in /usr/include/eigen3 /usr/local/include/eigen3 ./eigen3 ./eigen; do [ -d "$d" ] && echo "$d" && break; done)
  [ -z "${EIG:-}" ] && { echo "找不到 Eigen 3。装:apt-get install -y libeigen3-dev"; exit 1; }
  echo "编译中 (Eigen: $EIG) ..."
  g++ -O3 -march=native -fopenmp -std=c++17 -I"$EIG" "$SRC" -o "$BIN"
fi

echo "线程数 $N   日志 $LOG"
OMP_NUM_THREADS=$N OMP_WAIT_POLICY=passive OMP_PROC_BIND=close OMP_PLACES=cores \
  setsid nohup ./"$BIN" "$@" > "$LOG" 2>&1 < /dev/null &
PID=$!
echo "$PID" > run.pid
echo "PID $PID  (已脱离终端,ssh 断开不影响)"
echo
echo "看进度:   tail -f $LOG"
echo "看步数:   grep -c '^step ' $LOG"
echo "停掉:     kill \$(cat run.pid)"
