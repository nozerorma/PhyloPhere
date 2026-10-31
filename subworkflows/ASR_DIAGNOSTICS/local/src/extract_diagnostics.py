"""Select one observed contract from a gzip tar archive without extracting paths or links."""
import argparse
from pathlib import Path, PurePosixPath
import shutil
import tarfile


def extract_contract(archive, destination):
    names = {'caas_convergence_master.csv', 'observed_scoring.json'}
    with tarfile.open(archive, 'r:gz') as source:
        selected = {}
        seen = set()
        for member in source:
            path = PurePosixPath(member.name)
            if path.is_absolute() or '..' in path.parts or member.issym() or member.islnk():
                raise ValueError(f'Unsafe archive member: {member.name}')
            if not member.isfile() and not member.isdir():
                raise ValueError(f'Unsupported archive member: {member.name}')
            if member.isfile():
                if path in seen:
                    raise ValueError(f'Duplicate archive member: {member.name}')
                seen.add(path)
                if path.name in names:
                    selected.setdefault(path.parent, {})[path.name] = member
        candidates = [items for items in selected.values() if set(items) == names]
        if len(candidates) != 1 or len(selected) != 1:
            raise ValueError('Archive must contain exactly one master CSV and its completion record together')
        destination = Path(destination)
        destination.mkdir(parents=True, exist_ok=False)
        for name, member in candidates[0].items():
            with source.extractfile(member) as inp, (destination / name).open('wb') as out:
                shutil.copyfileobj(inp, out)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive')
    parser.add_argument('destination')
    args = parser.parse_args()
    extract_contract(args.archive, args.destination)
