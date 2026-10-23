"""Run the FINAL workflow's synthetic scientific and integration tests."""
import argparse
import json
from pathlib import Path
import sys
import time
import unittest
from threadpoolctl import threadpool_limits


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output')
    args = p.parse_args(argv)
    root = Path(__file__).resolve().parents[1]
    started = time.monotonic()
    suite = unittest.TestSuite([
        unittest.TestLoader().discover(str(root / 'tests'), pattern="test_*.py"),
        unittest.TestLoader().discover(str(root / 'tests/abm'), pattern="test_*.py"),
    ])
    with threadpool_limits(4):
        result = unittest.TextTestRunner(verbosity=1, buffer=True).run(suite)
    report = dict(status='passed' if result.wasSuccessful() else 'failed',
                  synthetic_only=True, tests=result.testsRun, failures=len(result.failures),
                  errors=len(result.errors), skipped=len(result.skipped),
                  passed=result.testsRun-len(result.failures)-len(result.errors)-len(result.skipped),
                  skip_reasons=[reason for _, reason in result.skipped],
                  seconds=time.monotonic()-started)
    if args.output:
        Path(args.output).write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))
    if not result.wasSuccessful():
        sys.exit(1)


if __name__ == '__main__':
    main()
