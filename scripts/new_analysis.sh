#!/usr/bin/env bash
# Stamp out a dated analysis directory.
#   usage:  bash scripts/new_analysis.sh park-compatibility
set -euo pipefail

SLUG="${1:?usage: new_analysis.sh <short-slug>}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DIR="$ROOT/analyses/$(date +%Y-%m)_${SLUG}"

[[ -e "$DIR" ]] && { echo "error: $DIR exists" >&2; exit 1; }
mkdir -p "$DIR"/{src,figures,results}
touch "$DIR"/figures/.gitkeep "$DIR"/results/.gitkeep

cat > "$DIR/CLAUDE.md" <<EOF
# Analysis: ${SLUG}

**Question.** <one sentence — what this analysis decides>

**Protocol.** See \`notes/sciphy_notes.md\` §<section>.

**Inputs.** <paths on scratch; do not copy data into this directory>

**State.** not started

## Conventions here
- Scripts in \`src/\`, figures in \`figures/\`, small tables in \`results/\`.
- Record findings in \`README.md\` as you go, not at the end.
EOF

cat > "$DIR/README.md" <<EOF
# ${SLUG}

**Question.**

**Answer.** <fill in when you have one — this is the line future-you reads>

## What was run

## Findings

## Caveats
EOF

echo "created $DIR"
echo "-> add a row to the analysis index in CLAUDE.md"
