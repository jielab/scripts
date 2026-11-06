"""Shared, stdlib-only selection of the main analysis and optional IBDmix."""

STAGES = ("primary", "archaic", "summary", "spatial", "windows", "alleles", "time",
          "selection", "candidates", "genotypes", "ld", "report")
DEFAULT_STAGES = tuple(name for name in STAGES if name != "archaic")


def boolean_option(value):
    if isinstance(value, bool):
        return value
    value = str(value).upper()
    if value not in {"TRUE", "FALSE"}:
        raise ValueError("Expected TRUE or FALSE")
    return value == "TRUE"


def select_stages(command="pipeline", stages=None, ibdmix=False):
    enabled = boolean_option(ibdmix)
    requested = [name.strip() for name in stages.split(",") if name.strip()] if stages else None
    chosen = requested if requested is not None else list(STAGES if enabled else DEFAULT_STAGES) if command == "pipeline" else [command]
    if not chosen or len(chosen) != len(set(chosen)) or set(chosen) - set(STAGES):
        raise ValueError("Select unique stages from " + ",".join(STAGES))
    if "archaic" in chosen and not enabled:
        raise ValueError("IBDmix is disabled by default; exploratory archaic analysis requires --ibdmix TRUE")
    return [name for name in STAGES if name in chosen]
