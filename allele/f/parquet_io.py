"""Bounded Parquet batches and temporary outputs committed only on success."""

from pathlib import Path
import os
import tempfile
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


def batches(path, columns=None, size=60000):
    for batch in pq.ParquetFile(path).iter_batches(batch_size=size, columns=columns):
        yield batch.to_pandas()


def selected_rows(path, ids, column="site_id", columns=None):
    ids = set(ids)
    parts = []
    for part in batches(path, columns):
        keep = part[part[column].isin(ids)]
        if len(keep):
            parts.append(keep)
    if not parts:
        return next(batches(path, columns)).iloc[:0]
    return pd.concat(parts, ignore_index=True)


class ParquetSink:
    def __init__(self, path):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        fd, name = tempfile.mkstemp(
            prefix=self.path.name + ".", suffix=".tmp", dir=self.path.parent
        )
        os.close(fd)
        self.temporary = Path(name)
        self.writer = None

    def write(self, data):
        if len(data) == 0:
            return
        table = pa.Table.from_pandas(data, preserve_index=False)
        if self.writer is None:
            self.writer = pq.ParquetWriter(
                self.temporary, table.schema, compression="zstd"
            )
        self.writer.write_table(table)

    def commit(self):
        if self.writer is None:
            raise ValueError("No rows to write")
        self.writer.close()
        self.writer = None
        self.temporary.replace(self.path)

    def close(self):
        if self.writer is not None:
            self.writer.close()
            self.writer = None
        self.temporary.unlink(missing_ok=True)
