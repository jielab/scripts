#!/usr/bin/env python3


# 🚩 Imports, paths, and shared inputs
import argparse
import re
import sys
import urllib.error
import urllib.request
from urllib.parse import urljoin

def read_url(url):
    req = urllib.request.Request(url, headers={'User-Agent': 'pfi-pipeline/1.0'})
    with urllib.request.urlopen(req, timeout=120) as r:
        return r.read().decode('utf-8', errors='ignore')

def deprecated_pmc_url(url):
    marker = '/pub/pmc/'
    deprecated_marker = '/pub/pmc/deprecated/'
    if marker in url and deprecated_marker not in url:
        return url.replace(marker, deprecated_marker, 1)
    return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--url', required=True)
    ap.add_argument('--pattern', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--limit', type=int, default=0)
    args = ap.parse_args()

    url = args.url
    try:
        html = read_url(url)
    except urllib.error.HTTPError as e:
        fallback_url = deprecated_pmc_url(url) if e.code == 404 else None
        if not fallback_url:
            print(f'ERROR: cannot open {url}: {e}', file=sys.stderr)
            raise
        print(f'WARNING: {url} returned 404; trying {fallback_url}', file=sys.stderr)
        url = fallback_url
        try:
            html = read_url(url)
        except Exception as fallback_error:
            print(f'ERROR: cannot open {url}: {fallback_error}', file=sys.stderr)
            raise
    except Exception as e:
        print(f'ERROR: cannot open {url}: {e}', file=sys.stderr)
        raise

    hrefs = re.findall(r'href=["\']([^"\']+)["\']', html, flags=re.I)
    rx = re.compile(args.pattern)
    urls = []
    for h in hrefs:
        name = h.split('/')[-1]
        if rx.search(name):
            urls.append(urljoin(url, h))
    urls = sorted(set(urls))
    if args.limit and args.limit > 0:
        urls = urls[:args.limit]

    with open(args.out, 'w') as f:
        for u in urls:
            f.write(u + '\n')
    print(f'Found {len(urls)} URLs -> {args.out}')

if __name__ == '__main__':
    main()
