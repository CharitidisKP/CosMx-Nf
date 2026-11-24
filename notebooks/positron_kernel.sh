#!/usr/bin/env bash
# Runs Positron's R kernel in the pipeline's image, with this workspace setting
#   "positron.r.kernel.path": "<Pipeline>/notebooks/positron_kernel.sh"
set -euo pipefail
PIPE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Positron's own ark sits in resources/ark, beside the resources/scripts/startup.R it passes
ARK=""
prev=""
for a in "$@"; do
  if [ "$prev" = "--startup-file" ]; then ARK="${a%/*/*}/ark/ark"; fi
  prev=$a
done
[ -x "$ARK" ] || { echo "positron_kernel.sh: no ark beside --startup-file" >&2; exit 1; }

# Passes Positron's own variables through --cleanenv as APPTAINERENV_ ones
for v in POSITRON POSITRON_VERSION POSITRON_LONG_VERSION POSITRON_MODE RUST_LOG RUST_BACKTRACE; do
  if [ -n "${!v:-}" ]; then export "APPTAINERENV_$v=${!v}"; fi
done
# The kernel's R profile, which sets Giotto's Python
export APPTAINERENV_R_PROFILE_USER="$PIPE/notebooks/kernel.Rprofile"

# The image and settings the stages use
# ark loads R's shared library itself, so it takes the library path from R's own ldpaths
exec "$PIPE/Environment/in_image.sh" sh -c '. "$R_HOME/etc/ldpaths" && exec "$@"' ark "$ARK" "$@"
