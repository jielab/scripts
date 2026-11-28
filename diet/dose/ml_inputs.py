"""Validate the FoodWAS feature contract before launching ML."""
import csv
from pathlib import Path

META = {'eid', 'region', 'region_code', 'cvd_stroke_i.t2e', 'cvd_stroke_i.Yt2e'}


def selected_features(dose_output):
    root = Path(dose_output)
    path = root / 'ygb/03ml/ML_input.csv'
    selected = root / 'ygb/02foodwas/group/Need_Food_lst.txt'
    with path.open(newline='') as handle:
        reader = csv.reader(handle)
        header = next(reader, [])
        if not META <= set(header) or len(header) != len(set(header)):
            raise ValueError(f'Invalid ML input header: {path}')
        if next(reader, None) is None:
            raise ValueError(f'Empty ML cohort: {path}')
    features = [name for name in header if name not in META]
    expected = [line.strip() for line in selected.read_text().splitlines() if line.strip()]
    if len(expected) != len(set(expected)) or set(features) != set(expected):
        raise ValueError('ML input features disagree with the current FoodWAS selection; rerun foodwas groups')
    return features
