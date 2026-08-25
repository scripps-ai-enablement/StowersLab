#!/bin/bash
# Install the VNO receptor RNA-seq pipeline as a shared Environment Modules
# module, so end users run it with no clone and no pip install:
#
#     module load vno-receptor
#     vr_analyze.py --results <nfcore results> --gtf <genes.gtf>
#
# Two directories are written:
#   $PREFIX/vno-receptor/$VERSION            the payload
#   $PREFIX/modules/vno-receptor/$VERSION    the modulefile
#
# At many sites those have DIFFERENT permissions -- the payload area is often
# group-writable while the modulefile area is owner-only. This script installs
# the payload, then tells you exactly what to hand the modules owner if it
# cannot write the second one.
#
# Usage:
#   PREFIX=/path/to/software \
#   bash install_shared_module.sh vno-receptor-rnaseq-repo.tar.gz
#
# Optional overrides: VERSION, PYTHON_MODULE, SAMTOOLS_MODULE
set -euo pipefail

TARBALL="${1:?usage: PREFIX=/path/to/software bash install_shared_module.sh <repo tarball>}"
PREFIX="${PREFIX:?set PREFIX to the software root for your site, e.g. /opt/shared/software}"
VERSION="${VERSION:-1.0.0}"
PYTHON_MODULE="${PYTHON_MODULE:-python/3.11.4}"
SAMTOOLS_MODULE="${SAMTOOLS_MODULE:-samtools/1.19}"

DEST="$PREFIX/vno-receptor/$VERSION"
MODDIR="$PREFIX/modules/vno-receptor"
HERE="$(cd "$(dirname "$0")" && pwd)"

echo "== 1. payload -> $DEST"
mkdir -p "$DEST"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
tar -xzf "$TARBALL" -C "$TMP"
SRC="$TMP/vno-receptor-rnaseq"
[[ -f "$SRC/bin/vr_analyze.py" ]] || { echo "unexpected tarball layout" >&2; exit 1; }

for d in bin config ref docs; do
  rm -rf "$DEST/$d"; cp -r "$SRC/$d" "$DEST/$d"
done
cp "$SRC"/*.md "$SRC/kernel.py" "$SRC/SKILL.md" "$SRC/install.py" "$DEST/" 2>/dev/null || true
chmod -R g+rX,o+rX "$DEST"
chmod +x "$DEST/bin"/*.py "$DEST/bin"/*.sh
echo "   $(find "$DEST" -type f | wc -l) files installed"

echo "== 2. modulefile -> $MODDIR/$VERSION"
MODFILE="$TMP/modulefile"
sed -e "s|^set vno_path .*|set vno_path      $DEST|" \
    -e "s|^set python_module .*|set python_module $PYTHON_MODULE|" \
    -e "s|^set samtools_module .*|set samtools_module $SAMTOOLS_MODULE|" \
    "$HERE/vno-receptor.modulefile" > "$MODFILE"

if mkdir -p "$MODDIR" 2>/dev/null && cp "$MODFILE" "$MODDIR/$VERSION" 2>/dev/null; then
  # default-version file, so `module load vno-receptor` needs no version suffix
  cat > "$MODDIR/.version" <<VEOF
#%Module
set ModulesVersion "$VERSION"
VEOF
  chmod g+rX,o+rX "$MODDIR/$VERSION" "$MODDIR/.version"
  echo "   modulefile installed"
  echo
  echo "== verify (in a fresh shell) =="
  echo "   module use $PREFIX/modules   # if not already on MODULEPATH"
  echo "   module load vno-receptor && vr_analyze.py --selftest"
else
  OUT="$HOME/vno-receptor-$VERSION.modulefile"
  cp "$MODFILE" "$OUT"
  echo "   !! cannot write $MODDIR -- it is likely owner-only."
  echo "      The payload IS installed. Send this to whoever owns the"
  echo "      modulefile directory, asking them to copy it to"
  echo "        $MODDIR/$VERSION"
  echo "      Ready, with paths already filled in:"
  echo "        $OUT"
fi
