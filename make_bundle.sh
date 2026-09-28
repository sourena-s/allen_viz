#!/usr/bin/env bash
# Packs the meshes, atlas and sample table into assets/bundle.js so index.html
# also works when opened directly from disk (file://), where fetch() is blocked.
# Re-run after changing SampleAnnot_all.csv or anything in assets/.
set -euo pipefail
cd "$(dirname "$0")"
b64() { base64 -w0 "$1"; }
{
  echo "window.BRAIN_BUNDLE = {"
  echo "  \"assets/BrainMesh_ICBM152.lh.mz3\": \"$(b64 assets/BrainMesh_ICBM152.lh.mz3)\","
  echo "  \"assets/AllenAtlas.nii.gz\": \"$(b64 assets/AllenAtlas.nii.gz)\","
  echo "  \"assets/mni152.nii.gz\": \"$(b64 assets/mni152.nii.gz)\","
  echo "  \"SampleAnnot_all.csv\": \"$(b64 SampleAnnot_all.csv)\""
  echo "};"
} > assets/bundle.js
echo "wrote assets/bundle.js ($(du -h assets/bundle.js | cut -f1))"
