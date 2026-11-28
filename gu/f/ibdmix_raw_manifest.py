"""Save interpretation settings, never transformed calls, during IBDmix calling."""
import argparse
import json
from pathlib import Path
from ibdmix_layout import cohort_for, unit_metadata_name


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--run",type=Path,required=True)
    for name in ("unit","chrom","build","locus"):p.add_argument("--"+name,required=True)
    for name in ("core-start","core-end"):p.add_argument("--"+name,type=int,default=0)
    a=p.parse_args()
    cohort=cohort_for(a.run)
    name=unit_metadata_name(a.unit,a.run.name)
    target=(cohort/"metadata"/name) if cohort else a.run/"metadata"/(a.unit+".json")
    target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(dict(schema="ibdmix-raw-v1",genome_build=a.build,locus_id=a.locus,chromosome=a.chrom,
        native_coordinate_system="1-based-exclusive-end",analysis_coordinate_system="0-based-half-open",
        raw_interpretation=dict(core_start=a.core_start,core_end=a.core_end)),indent=2)+"\n")

if __name__=="__main__":main()
