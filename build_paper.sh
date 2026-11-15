#!/usr/bin/env bash
# Build paper.pdf from paper.md in an ACM-like single-column style.
set -euo pipefail
cd "$(dirname "$0")"
pandoc paper.md --citeproc --pdf-engine=xelatex \
    --shift-heading-level-by=-1 -H paper-style.tex -o paper.pdf
