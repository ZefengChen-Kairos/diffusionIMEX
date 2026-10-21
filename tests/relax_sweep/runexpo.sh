#!/bin/bash
R=/home/user/diffusionIMEX/run/expo; BIN=$R/pc
one() { IFS='|' read -r n a <<< "$1"; cd $R || return 1
  [ -s "$n.log" ] && grep -q "^case" "$n.log" && return 0
  OMP_NUM_THREADS=1 "$BIN" $a out="$R/$n.dat" > "$R/$n.log" 2>&1; }
export -f one; export R BIN
xargs -P ${NP:-3} -I{} bash -c 'one "$@"' _ {} < "$1"
echo ALLDONE
