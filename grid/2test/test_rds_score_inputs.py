"""Published RDS provenance must not prevent downstream score consumption."""
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "f"))
from shared_cohort import table
from twin import read_table

RSCRIPT = shutil.which("Rscript") or str(Path(sys.prefix) / "bin/Rscript")


@unittest.skipUnless(Path(RSCRIPT).is_file(), "Native Rscript is required for RDS provenance")
class RdsScoreInputsTests(unittest.TestCase):
    def test_raw_provenance_preserves_identifiers_scores_and_missing_values(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = root / "scores.rds"
            script = root / "write.R"
            script.write_text('''args <- commandArgs(TRUE)
x <- data.frame(eid=c("001", "002"), csx.EUR=c(0.25, NA_real_),
                group=c("01", "02"), check.names=FALSE)
attr(x, "table_exchange") <- list(bytes=charToRaw("source table bytes"))
saveRDS(x, args[1])
''')
            subprocess.run([RSCRIPT, "--vanilla", str(script), str(path)], check=True,
                           capture_output=True, text=True)
            native = read_table(path)
            self.assertEqual(native.eid.tolist(), ["001", "002"])
            self.assertEqual(native.group.tolist(), ["01", "02"])
            self.assertEqual(native["csx.EUR"].iloc[0], .25)
            self.assertTrue(np.isnan(native["csx.EUR"].iloc[1]))
            self.assertEqual(table(path).id.tolist(), ["001", "002"])


if __name__ == "__main__":
    unittest.main()
