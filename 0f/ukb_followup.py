"""Resolve the UKB administrative limit from the canonical R definition."""
from datetime import date
from functools import lru_cache
import hashlib
import os
from pathlib import Path
import re
import subprocess

_SOURCE = Path(__file__).resolve().with_name('phenotype.R')
_R = r'''
args <- commandArgs(trailingOnly=TRUE)
e <- new.env(parent=globalenv())
invisible(capture.output(source(args[1], local=e, echo=FALSE)))
if (!exists("date_follow_end", envir=e, inherits=FALSE))
  stop("Missing shared date_follow_end")
x <- get("date_follow_end", envir=e, inherits=FALSE)
if (!inherits(x, "Date") || length(x) != 1L || is.na(x))
  stop("Shared date_follow_end must be one non-missing Date")
cat("JIELAB_FOLLOWUP\t", format(x, "%Y-%m-%d"), "\n", sep="")
'''


@lru_cache(maxsize=4)
def _read(source, source_sha256, rscript):
    try:
        result = subprocess.run([rscript, '--vanilla', '-e', _R, source],
                                check=True, capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.SubprocessError) as exc:
        detail = getattr(exc, 'stderr', '') or str(exc)
        raise RuntimeError(f'Cannot resolve shared date_follow_end using {rscript}: {detail}') from exc
    prefix = 'JIELAB_FOLLOWUP\t'
    values = [line[len(prefix):] for line in result.stdout.splitlines() if line.startswith(prefix)]
    if len(values) != 1:
        raise ValueError('Cannot resolve one shared date_follow_end from phenotype.R')
    value = values[0]
    if date.fromisoformat(value).isoformat() != value:
        raise ValueError('Shared follow-up date is not an ISO date')
    if hashlib.sha256(Path(source).read_bytes()).hexdigest() != source_sha256:
        raise RuntimeError('phenotype.R changed while resolving date; retry')
    return value


def resolve_followup(*, rscript=None, environ=None):
    env = os.environ if environ is None else environ
    if str(env.get('DATE_FOLLOW_END', '')).strip():
        raise ValueError('DATE_FOLLOW_END is no longer supported; unset it. '
                         'The date comes from 0f/phenotype.R::date_follow_end')
    source = str(_SOURCE)
    sha = hashlib.sha256(_SOURCE.read_bytes()).hexdigest()
    value = _read(source, sha, str(rscript or env.get('R_BIN') or 'Rscript'))
    return dict(date=value, source=source, variable='date_follow_end', source_sha256=sha)


def validate_record(config, section, key, followup):
    """Accept a saved resolved value only when it agrees with the current source.

    The source hash is provenance, not a date-only cache invalidation key.
    This function never treats a saved value as a date override.
    """
    data = config.get(section, {})
    if not isinstance(data, dict):
        raise ValueError(f'{section} must be a configuration mapping')
    if section == 'phenotypes':
        outcomes = data.get('outcomes', {})
        if not isinstance(outcomes, dict):
            raise ValueError('phenotypes.outcomes must be a mapping')
        for outcome, fields in outcomes.items():
            if isinstance(fields, dict) and 'end_date' in fields:
                raise ValueError(f'phenotypes.outcomes.{outcome}.end_date is not supported; '
                                 'edit 0f/phenotype.R::date_follow_end')
    record = config.get('followup_provenance')
    if key not in data and record is None:
        return
    if (not isinstance(record, dict) or
            not re.fullmatch(r'[0-9a-f]{64}', str(record.get('source_sha256', '')))):
        raise ValueError(f'{section}.{key} cannot override 0f/phenotype.R::date_follow_end; '
                         'remove the date from user configuration')
    if (record.get('variable') != followup['variable'] or
            record.get('source') != followup['source'] or
            record.get('date') != followup['date'] or data.get(key) != followup['date']):
        raise ValueError('Stale resolved follow-up plan: regenerate the affected plan from '
                         '0f/phenotype.R::date_follow_end')
