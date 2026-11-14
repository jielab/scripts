#!/usr/bin/env bash
# Use the same OS table and existing scientific runtimes as le8/f/runtime_env.sh.
# Source after 0f/0os.sh and workflow_args; never activate or install environments.
if [[ $submit == hpc ]]; then
    export DIET_RSCRIPT=${DIET_RSCRIPT:-${R_BIN:-$LE8_RSCRIPT}}
    export DIET_PYTHON=${DIET_PYTHON:-${PYTHON_BIN:-$LE8_PYTHON}}
else
    source "$HELPERS/local_runtime.sh"
    export DIET_RSCRIPT=${DIET_RSCRIPT:-${R_BIN:-${LE8_RSCRIPT:-Rscript}}}
    export DIET_PYTHON=${DIET_PYTHON:-${PYTHON_BIN:-${LE8_PYTHON:-$(workflow_local_runtime le8 python3 python3)}}}
fi
# Select runtime before launch, including explicit CLI overrides. Keep these
# arguments for the Python parser and record the selected paths in job plans.
diet_runtime_args=("$@")
for ((diet_i=0; diet_i<${#diet_runtime_args[@]}; diet_i++)); do
    case ${diet_runtime_args[diet_i]} in
        --rscript|--r-bin|--python|--python-bin)
            ((diet_i+1<${#diet_runtime_args[@]})) || { echo 'Missing interpreter path' >&2; return 2; }
            diet_key=${diet_runtime_args[diet_i]}; ((diet_i+=1))
            case $diet_key in --rscript|--r-bin) export DIET_RSCRIPT=${diet_runtime_args[diet_i]} ;; *) export DIET_PYTHON=${diet_runtime_args[diet_i]} ;; esac ;;
    esac
done
unset diet_runtime_args diet_i diet_key
if [[ $submit == hpc ]]; then
    # Match LE8 isolation from login-shell SAIGE/TaiYi R settings.
    unset R_HOME R_LIBS R_LIBS_USER R_LIBS_SITE
    export R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null
    export PATH="$(dirname -- "$DIET_RSCRIPT"):$(dirname -- "$DIET_PYTHON"):$PATH"
    diet_runtime_lib="$(dirname -- "$DIET_PYTHON")/../lib"
    if [[ -r $diet_runtime_lib/libstdc++.so.6 ]]; then
        export LD_LIBRARY_PATH="$diet_runtime_lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
    fi
    unset diet_runtime_lib
fi
export PYTHONDONTWRITEBYTECODE=1

# Optional project-only additions; never replace the selected LE8 environment.
export DIET_R_LIBRARY=${DIET_R_LIBRARY:-$DIR_SOFTWARE/diet/R-library}
export DIET_PYTHON_LIBRARY=${DIET_PYTHON_LIBRARY:-$DIR_SOFTWARE/diet/python}
if [[ -d $DIET_R_LIBRARY ]]; then
    export R_LIBS_USER="$DIET_R_LIBRARY${R_LIBS_USER:+:$R_LIBS_USER}"
fi
if [[ -d $DIET_PYTHON_LIBRARY ]]; then
    export PYTHONPATH="$DIET_PYTHON_LIBRARY${PYTHONPATH:+:$PYTHONPATH}"
fi

# Optional private system libraries prepared by the explicit installer.
export DIET_SYSTEM_PREFIX=${DIET_SYSTEM_PREFIX:-$DIR_SOFTWARE/diet/system-libs}
if [[ -d $DIET_SYSTEM_PREFIX/lib ]]; then
    export LD_LIBRARY_PATH="$DIET_SYSTEM_PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
    export PKG_CONFIG_PATH="$DIET_SYSTEM_PREFIX/lib/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
    export PROJ_DATA=${PROJ_DATA:-$DIET_SYSTEM_PREFIX/share/proj}
    export GDAL_DATA=${GDAL_DATA:-$DIET_SYSTEM_PREFIX/share/gdal}
    export PATH="$DIET_SYSTEM_PREFIX/bin:$PATH"
fi
if [[ -d $DIET_SYSTEM_PREFIX/lib/R/library ]]; then
    export R_LIBS_USER="${R_LIBS_USER:+$R_LIBS_USER:}$DIET_SYSTEM_PREFIX/lib/R/library"
fi

# CPU-only cluster nodes need the project-scoped OSMesa backend.
export DIET_GRAPHICS_PREFIX=${DIET_GRAPHICS_PREFIX:-$DIR_SOFTWARE/diet/graphics}
if [[ -z ${DISPLAY:-} && -f $DIET_GRAPHICS_PREFIX/lib/libOSMesa.so.8 &&
      ( -z ${VTK_DEFAULT_OPENGL_WINDOW:-} || ${VTK_DEFAULT_OPENGL_WINDOW:-} == vtkOSOpenGLRenderWindow ) ]]; then
    export LD_LIBRARY_PATH="$DIET_GRAPHICS_PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
    export VTK_DEFAULT_OPENGL_WINDOW=vtkOSOpenGLRenderWindow
fi
