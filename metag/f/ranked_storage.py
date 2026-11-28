"""Publish completed directories on the same filesystem without symlinks."""
from contextlib import contextmanager
import ctypes
import errno
import os
from pathlib import Path
import shutil
import tempfile


def remove_entry(path):
    if path.is_symlink() or path.is_file():
        path.unlink()
    elif path.exists():
        shutil.rmtree(path)


def recover_directory(dest):
    """Restore an interrupted rotation on filesystems without directory exchange."""
    journal = dest.with_name('.' + dest.name + '.metag-rotation')
    if not journal.exists():
        return
    marker = journal / 'destination'
    if not marker.is_file() or marker.read_text() != dest.name:
        raise RuntimeError('Unrecognized rotation directory: ' + str(journal))
    backup = journal / 'previous'
    if backup.exists() or backup.is_symlink():
        if not (dest.exists() or dest.is_symlink()):
            backup.rename(dest)
        else:
            remove_entry(backup)
    marker.unlink()
    journal.rmdir()


def publish_directory(stage, dest):
    # Caller holds the destination lock. Both paths live on the destination FS.
    recover_directory(dest)
    if not (dest.exists() or dest.is_symlink()):
        stage.rename(dest)
        return
    exchange = getattr(ctypes.CDLL(None, use_errno=True), 'renameat2', None)
    if exchange is not None:
        if exchange(-100, os.fsencode(stage), -100, os.fsencode(dest), 2) == 0:
            return  # The caller's staging cleanup removes the old entry.
        code = ctypes.get_errno()
        if code not in (errno.ENOSYS, errno.EINVAL, errno.EOPNOTSUPP, errno.EXDEV):
            raise OSError(code, os.strerror(code), str(dest))
    journal = dest.with_name('.' + dest.name + '.metag-rotation')
    journal.mkdir()
    (journal / 'destination').write_text(dest.name)
    try:
        dest.rename(journal / 'previous')
        stage.rename(dest)
    finally:
        recover_directory(dest)



@contextmanager
def staging(parent, prefix):
    parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(dir=parent, prefix=prefix))
    try:
        yield stage
    finally:
        remove_entry(stage)
