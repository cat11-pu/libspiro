#!/bin/sh
# Build repro against the library sources. Two flavours:
#   default  -> mirrors ./configure with no --enable-test_inputs
#   checked  -> mirrors --enable-test_inputs (CHECK_INPUT_FINITENESS on)
set -e
cd "$(dirname "$0")/.."

cc -g -O0 -I. -Irepro  \
   spiro.c bezctx.c spiroentrypoints.c repro/repro.c -lm -o repro/repro_default

cc -g -O0 -I. -Irepro  \
   -DCHECK_INPUT_FINITENESS=1 \
   spiro.c bezctx.c spiroentrypoints.c repro/repro.c -lm -o repro/repro_checked
echo "built repro/repro_default and repro/repro_checked"
