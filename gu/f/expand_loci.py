#!/usr/bin/env python3
"""Expand normalized BED loci while retaining an auditable core/analysis map."""
from __future__ import annotations
import argparse
import re
from pathlib import Path
from phyml import CHROM_LENGTHS

def parse_bp(value:str)->int:
    m=re.fullmatch(r'\s*([0-9]+)\s*(bp|k|kb|m|mb)?\s*',value,re.I)
    if not m:raise argparse.ArgumentTypeError('use a non-negative size such as 0, 100kb, or 1mb')
    scale={None:1,'bp':1,'k':1000,'kb':1000,'m':1000000,'mb':1000000}[m.group(2).lower() if m.group(2) else None]
    return int(m.group(1))*scale

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--input',required=True,type=Path); ap.add_argument('--output',required=True,type=Path)
    ap.add_argument('--map',required=True,type=Path); ap.add_argument('--flank',default='100kb',type=parse_bp)
    ap.add_argument('--build',required=True,choices=['37','38']); args=ap.parse_args()
    rows=[]
    for line_no,line in enumerate(args.input.read_text().splitlines(),1):
        if not line.strip() or line.lstrip().startswith('#'):continue
        fields=line.split('\t');
        if len(fields)<4:raise SystemExit(f'ERROR: normalized BED line {line_no} has fewer than four columns')
        chrom,start,end,locus=fields[:4]; start=int(start); end=int(end); chrom_len=CHROM_LENGTHS[args.build][chrom]
        analysis_start=max(0,start-args.flank); analysis_end=min(chrom_len,end+args.flank)
        rows.append((chrom,start,end,analysis_start,analysis_end,locus,args.flank))
    if not rows:raise SystemExit('ERROR: no loci to expand')
    args.output.parent.mkdir(parents=True,exist_ok=True); args.map.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(''.join(f'{c}\t{a0}\t{a1}\t{locus}\n' for c,_,_,a0,a1,locus,_ in rows))
    args.map.write_text('chr\tcore_start\tcore_end\tanalysis_start\tanalysis_end\tlocus_id\tflank_bp\n'+
                        ''.join(f'{c}\t{s}\t{e}\t{a0}\t{a1}\t{locus}\t{flank}\n' for c,s,e,a0,a1,locus,flank in rows))
    print(args.flank)

if __name__=='__main__':main()
