"""Public score names and the existing engines that can produce them.

Names are the XXX in diet.XXX.sum, not manuscript/study names. Engines may
calculate a batch of related scores; publication selects the requested totals.
"""
UKB_COMPARATORS = (
    "medito", "medi24", "dash", "mind", "hpdi", "ahei", "phdi", "modern",
    "digm", "maha", "hlcd", "hlfd", "redii", "fds",
)
MAHA_VARIANTS = ("maha_bal", "maha_strict", "maha_nodairy", "maha_noprotein")
DOSE_SCORES = ("dose", "dose11", "dose9", "dose9half")
SCORE_ENGINES = {
    "ukb": {"maha": (*UKB_COMPARATORS, *MAHA_VARIANTS, "upf"),
            "dose": (*UKB_COMPARATORS, *DOSE_SCORES)},
    "nhanes": {"maha": ("maha", *MAHA_VARIANTS, "maha_processed", "maha_nomeat", "dash", "medi", "mind"),
               "dose": ("dose", "dose9half", "dash", "mind", "amed")},
    "chns": {"maha": ("chdi", "maha", *MAHA_VARIANTS, "maha_density", "maha_residual",
                       "dash", "dash_fung", "dash_density", "mind", "medi")},
    "charls": {"charls": ("dose9half", "maha", "dash", "medi")},
}
SCORES = tuple(dict.fromkeys(s for engines in SCORE_ENGINES.values() for scores in engines.values() for s in scores))


def score_names(value):
    if value is None or value == "all":
        return list(SCORES)
    names = [s.strip() for s in value.split(",")]
    unknown = set(names) - set(SCORES)
    if unknown:
        raise ValueError("Unknown --score name(s): " + ", ".join(sorted(unknown)) +
                         "; choose all or " + ",".join(SCORES))
    return list(dict.fromkeys(names))
