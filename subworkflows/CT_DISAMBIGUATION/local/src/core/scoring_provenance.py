"""Content identities for observed/null scoring compatibility (no inference)."""
import hashlib
import json
from pathlib import Path


def digest(path):
    if not path:
        return None
    h=hashlib.sha256()
    with Path(path).open('rb') as fh:
        for chunk in iter(lambda:fh.read(1048576),b''):
            h.update(chunk)
    return h.hexdigest()


def scientific(tree,model,threshold,taxid,weights):
    return dict(tree_sha256=digest(tree),model=model,posterior_threshold=threshold,
                taxid_sha256=digest(taxid),weights_sha256=digest(weights))


def save(path,record):
    Path(path).write_text(json.dumps(record,sort_keys=True,indent=2))
