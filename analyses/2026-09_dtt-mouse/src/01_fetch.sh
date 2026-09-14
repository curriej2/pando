#!/usr/bin/env bash
# Discover and fetch the DTT mouse (Yu, Kim, Seidel et al. 2026) public data.
# ⚠ Must run on a COMPUTE node: ncbi.nlm.nih.gov does not resolve from the login
# nodes, though github.com and shendure-web.gs.washington.edu do.
set -uo pipefail
DEST=/data1/choij10/justin/pando/data/dtt_mouse
mkdir -p "$DEST"/{geo,trees,meta}
GEO="GSE341627"; TOK="gbwpuyeqhrypxyt"
echo "=== connectivity ==="
for h in www.ncbi.nlm.nih.gov ftp.ncbi.nlm.nih.gov nextcell.pages.dev shendure-web.gs.washington.edu; do
  printf "  %-38s " "$h"; curl -sI --max-time 20 "https://$h/" -o /dev/null -w "%{http_code}\n" || echo "FAIL"
done
echo; echo "=== GEO series page ==="
curl -sL --max-time 60 "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${GEO}&token=${TOK}" \
     -o "$DEST/geo/series.html" && echo "  series.html $(wc -c < "$DEST/geo/series.html") bytes"
echo; echo "=== GEO supplementary listing (FTP over HTTPS) ==="
for sub in suppl matrix; do
  u="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE341nnn/${GEO}/${sub}/"
  printf "  %-64s " "$sub/"; curl -sL --max-time 40 "$u" -o "$DEST/geo/list_${sub}.html" -w "%{http_code}\n"
done
echo; echo "=== trees from nextcell ==="
for f in nextcell_backbone_tree_dated.nwk.gz nextcell_full_tree_dated.nwk.gz; do
  printf "  %-44s " "$f"
  curl -sL --max-time 900 "https://nextcell.pages.dev/data/$f" -o "$DEST/trees/$f" \
       -w "%{http_code} %{size_download} bytes\n"
done
echo; echo "=== cell metadata (150 MB) ==="
curl -sL --max-time 1800 \
  "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/NextCell/sc_transcriptome/cell_metadata.annotation.txt" \
  -o "$DEST/meta/cell_metadata.annotation.txt" -w "  %{http_code} %{size_download} bytes\n"
echo; echo "=== what landed ==="; find "$DEST" -type f -printf "  %-70p %10s\n" | sort
