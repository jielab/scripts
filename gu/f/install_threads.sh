#!/usr/bin/env bash
# Compatibility entry point; installation is implemented in ../install.sh.
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
exec bash "$ROOT/install.sh" --threads "$@"
