#!/usr/bin/env bash
# Vendor the msprime C library, the tskit C library and kastore into src/,
# applying the small CRAN-compliance patches (no abort()/stderr/stdout).
#
# Usage: tools/sync_c_sources.sh [path/to/msprime] [path/to/tskit]
#        (defaults: ../msprime and ../tskit)
set -euo pipefail

PKG_DIR="$(cd "$(dirname "$0")/.." && pwd)"
MSP_DIR="$(cd "${1:-$PKG_DIR/../msprime}" && pwd)"
TSKIT_DIR="$(cd "${2:-$PKG_DIR/../tskit}" && pwd)"

SRC="$PKG_DIR/src"
mkdir -p "$SRC/tskit" "$SRC/kastore" "$SRC/lib"

for f in core convert genotypes haplotype_matching stats tables trees; do
    cp "$TSKIT_DIR/c/tskit/$f.c" "$TSKIT_DIR/c/tskit/$f.h" "$SRC/tskit/"
done
cp "$TSKIT_DIR/c/tskit.h" "$SRC/tskit.h"
cp "$TSKIT_DIR/c/subprojects/kastore/kastore.c" "$TSKIT_DIR/c/subprojects/kastore/kastore.h" \
    "$SRC/kastore/"

for f in msprime.c msprime.h mutgen.c likelihood.c likelihood.h rate_map.c rate_map.h \
         fenwick.c fenwick.h avl.c avl.h object_heap.c object_heap.h util.c util.h; do
    cp "$MSP_DIR/lib/$f" "$SRC/lib/"
done

# ISO C forbids returning a void expression (flagged by -Wpedantic).
perl -pi -e 's/^(\s*)return avl_clear_tree\(avltree\);/$1avl_clear_tree(avltree);/' \
    "$SRC/lib/avl.c"

# Route bug asserts to R's error handler instead of fprintf(stderr)+abort().
perl -0pi -e 's/fprintf\(stderr, "Bug detected in %s at line %d\. %s\\n", __FILE__, __LINE__,\s*\\\n\s*TSK_BUG_ASSERT_MESSAGE\);\s*\\\n\s*abort\(\);/msprimer_bug_assert_failed(__FILE__, __LINE__);/' \
    "$SRC/tskit/core.h"
perl -0pi -e 's/fprintf\(stderr, "Bug detected in %s at line %d\. %s\\n", __FILE__, __LINE__,\s*\\\n\s*KAS_BUG_ASSERT_MESSAGE\);\s*\\\n\s*abort\(\);/msprimer_bug_assert_failed(__FILE__, __LINE__);/' \
    "$SRC/kastore/kastore.h"
perl -0pi -e 's/(#define tsk_bug_assert\(condition\))/void msprimer_bug_assert_failed(const char *file, int line);\n$1/' \
    "$SRC/tskit/core.h"
perl -0pi -e 's/(#define kas_bug_assert\(condition\))/void msprimer_bug_assert_failed(const char *file, int line);\n$1/' \
    "$SRC/kastore/kastore.h"
# The debug stream is never used (TSK_DEBUG is never passed).
perl -pi -e 's/^#define TSK_DEFAULT_DEBUG_STREAM stdout/#define TSK_DEFAULT_DEBUG_STREAM NULL/' \
    "$SRC/tskit/core.h"
# The msprime library includes headers as <tskit.h>; it lives in src/.
for pat in "msprimer_bug_assert_failed(__FILE__" "TSK_DEFAULT_DEBUG_STREAM NULL"; do
    grep -q "$pat" "$SRC/tskit/core.h" || { echo "patch failed: $pat" >&2; exit 1; }
done
grep -q "msprimer_bug_assert_failed(__FILE__" "$SRC/kastore/kastore.h" \
    || { echo "kastore patch failed" >&2; exit 1; }
if grep -n "stdout\|stderr\|abort()" "$SRC"/lib/*.c "$SRC"/lib/*.h | grep -v "^\S*:\s*\*\|//\|printf(\"\|check the stderr"; then
    echo "WARNING: msprime C sources reference stdout/stderr/abort" >&2
fi

{
    echo "msprime C library: $(git -C "$MSP_DIR" describe --tags --always 2>/dev/null || echo unknown)"
    echo "msprime git commit: $(git -C "$MSP_DIR" rev-parse HEAD 2>/dev/null || echo unknown)"
    echo "tskit C version: $(cat "$TSKIT_DIR/c/VERSION.txt")"
    echo "kastore version: $(cat "$TSKIT_DIR/c/subprojects/kastore/VERSION.txt")"
    echo "tskit git commit: $(git -C "$TSKIT_DIR" rev-parse HEAD 2>/dev/null || echo unknown)"
} > "$PKG_DIR/inst/C_SOURCES_VERSION"

echo "Vendored C sources into $SRC"
cat "$PKG_DIR/inst/C_SOURCES_VERSION"
