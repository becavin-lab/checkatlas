from types import ModuleType
from typing import Optional


def _patch_resource_tracker():
    # Python 3.13 resource_tracker doesn't recognise joblib/loky resource
    # types (file, folder, semlock) and floods stderr with ValueError.
    try:
        import multiprocessing.resource_tracker as _rt
        import os
        import signal
        import sys

        _KNOWN = frozenset(("semlock", "file", "folder", "shared_memory"))

        def _safe_main():
            keep_running = True
            locked_files: set = set()

            def _cleanup(signum=-1, frame=None):
                nonlocal keep_running
                keep_running = False

            signal.signal(signal.SIGINT, _cleanup)
            signal.signal(signal.SIGTERM, _cleanup)

            fd = int(sys.argv[1]) if len(sys.argv) > 1 else None
            if fd is None:
                return

            while keep_running:
                try:
                    line = os.read(fd, 300)
                    if not line:
                        break
                    line = line.decode("ascii", errors="replace").strip()
                    if not line:
                        break
                    parts = line.split()
                    if len(parts) < 3:
                        continue
                    cmd, name, rtype = parts[0], parts[1], parts[2]
                    if rtype not in _KNOWN:
                        continue
                    if cmd == "REGISTER":
                        locked_files.add(name)
                    elif cmd == "UNREGISTER":
                        locked_files.discard(name)
                except (OSError, EOFError):
                    break

            for name in list(locked_files):
                try:
                    if os.path.exists(name):
                        if os.path.isdir(name):
                            os.rmdir(name)
                        else:
                            os.unlink(name)
                except OSError:
                    pass

        _rt.main = _safe_main
    except Exception:
        pass


_patch_resource_tracker()

from . import atlas, cellranger, check

seurat: Optional[ModuleType] = None
try:
    from . import seurat
except (ImportError, OSError) as e:
    import logging

    logging.getLogger("checkatlas").warning(
        f"Could not import seurat module (rpy2/R not available): {e}"
    )

__all__ = ["atlas", "seurat", "check", "cellranger"]
