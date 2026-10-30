#!/usr/bin/env bash
# Package the website and the WordPress plugin into dist/:
#   dist/degenprime-netlify-site-<version>.zip      drag onto https://app.netlify.com/drop
#   dist/degenprime-wordpress-plugin-<version>.zip  Plugins > Add New > Upload Plugin
# Uses the committed web/assets/degenprime.wasm; pass --wasm to rebuild it first.
# Both include the logo and icons from web/assets (made by scripts/make_icon.py).
set -euo pipefail
cd "$(dirname "$0")/.."
[ "${1:-}" = "--wasm" ] && bash scripts/build_wasm.sh
VERSION="$(grep -m1 '^version' Cargo.toml | cut -d'"' -f2)"
[ -s web/assets/degenprime.wasm ] || { echo "web/assets/degenprime.wasm is missing: run scripts/build_wasm.sh" >&2; exit 1; }
cp data/example_anfD.clust web/assets/example_anfD.clust
rm -rf dist && mkdir -p dist/netlify dist/wordpress/degenprime/assets
cp -r web/. dist/netlify/
cp wordpress/degenprime/degenprime.php wordpress/degenprime/readme.txt wordpress/degenprime/degenprime-embed.css dist/wordpress/degenprime/
cp LICENSE dist/wordpress/degenprime/LICENSE.txt
cp LICENSE dist/netlify/LICENSE.txt
# The built-in aligner's licence notices travel with every copy of the engine.
cp src/mafft/LICENSE-MAFFT.txt dist/wordpress/degenprime/LICENSE-MAFFT.txt
cp src/mafft/LICENSE-MAFFT.txt dist/netlify/LICENSE-MAFFT.txt
cp -r web/assets/. dist/wordpress/degenprime/assets/
# Browser-tab and home-screen icons belong to the website, not to a WordPress theme.
rm -f dist/wordpress/degenprime/assets/{favicon.ico,icon-192.png,icon-512.png,apple-touch-icon.png}
cp wordpress/degenprime/assets.htaccess dist/wordpress/degenprime/assets/.htaccess
zipdir() { # zipdir <source dir> <archive> <entry>
    python3 - "$1" "$2" "$3" <<'PY'
import os, sys, zipfile
src, out, entry = sys.argv[1:]
with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
    for root, _, files in os.walk(os.path.join(src, entry)):
        for f in sorted(files):
            p = os.path.join(root, f)
            z.write(p, os.path.relpath(p, src))
PY
}
zipdir dist/netlify "dist/degenprime-netlify-site-$VERSION.zip" .
zipdir dist/wordpress "dist/degenprime-wordpress-plugin-$VERSION.zip" degenprime
ls -l dist/*.zip
