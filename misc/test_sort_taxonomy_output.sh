#!/usr/bin/env bash
set -euo pipefail

test_dir="$(mktemp -d "${TMPDIR:-/tmp}/sort-taxonomy.XXXXXX")"
trap 'rm -rf "$test_dir"' EXIT
program="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)/sort_taxonomy_output"

cat > "$test_dir/nodes.dmp" <<'EOF'
1 | 1 | no rank |
2 | 1 | kingdom |
3 | 2 | species |
4 | 2 | species |
5 | 1 | kingdom |
6 | 5 | species |
EOF

cat > "$test_dir/input.lca" <<'EOF'
# command line
queryid	seq	len	naln	nspec	dustscore	gc	lca	taxa_path
read-6	AC	2	1	1	0	0	6	root;B;B1
read-4	AC	2	1	1	0	0	4	root;A;A2
read-3	AC	2	1	1	0	0	3	root;A;A1
read-5	AC	2	1	1	0	0	5	root;B
EOF
gzip -c "$test_dir/input.lca" > "$test_dir/input.lca.gz"

"$program" --input "$test_dir/input.lca.gz" \
  --nodes "$test_dir/nodes.dmp" --output "$test_dir/output.lca.gz"

actual="$(gunzip -c "$test_dir/output.lca.gz" | tail -n 4 | cut -f8)"
expected=$'3\n4\n5\n6'
test "$actual" = "$expected"

printf '%s\n' \
  $'taxid\tA\tq\tc\tphi\tllh' \
  $'6\t0\t0\t0\t0\t0' \
  $'3\t0\t0\t0\t0\t0' \
  $'4\t0\t0\t0\t0\t0' > "$test_dir/input.dfit"
gzip -c "$test_dir/input.dfit" > "$test_dir/input.dfit.gz"
"$program" --input "$test_dir/input.dfit.gz" --nodes "$test_dir/nodes.dmp" --output "$test_dir/output.dfit.gz"
test "$(gunzip -c "$test_dir/output.dfit.gz" | tail -n 3 | cut -f1)" = $'3\n4\n6'

printf '%s\n' \
  $'taxid\tname\trank\tnalign\tnreads' \
  $'6\tB1\tspecies\t1\t1' \
  $'3\tA1\tspecies\t1\t1' \
  $'4\tA2\tspecies\t1\t1' > "$test_dir/input.stat"
"$program" --input "$test_dir/input.stat" --nodes "$test_dir/nodes.dmp" --output "$test_dir/output.stat"
test "$(tail -n 3 "$test_dir/output.stat" | cut -f1)" = $'3\n4\n6'

awk 'BEGIN { print "queryid\tseq\tlen\tnaln\tnspec\tdustscore\tgc\tlca\ttaxa_path" }
     { taxid = (NR % 2 ? 6 : 3); print "read" NR "\tAC\t2\t1\t1\t0\t0\t" taxid "\troot" }' \
    < <(seq 1 100001) | gzip -c > "$test_dir/large.lca.gz"
"$program" --input "$test_dir/large.lca.gz" --nodes "$test_dir/nodes.dmp" --output "$test_dir/large.sorted.lca.gz"
test "$(gunzip -c "$test_dir/large.sorted.lca.gz" | sed -n '2p' | cut -f8)" = "3"
test "$(gunzip -c "$test_dir/large.sorted.lca.gz" | tail -n 1 | cut -f8)" = "6"

if "$program" --input "$test_dir/not-a.bdamage.gz" --nodes "$test_dir/nodes.dmp" >/dev/null 2>&1; then
    echo "Expected .bdamage.gz rejection" >&2
    exit 1
fi

echo "sort_taxonomy_output test passed"
