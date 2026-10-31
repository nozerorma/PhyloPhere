#!/usr/bin/env python3
"""Task entry point: parameters are JSON data, never shell or Python source."""
import json
from pathlib import Path
import subprocess
import sys


def main():
    config = json.loads(Path(sys.argv[1]).read_text())
    subprocess.run([sys.executable, 'observed_b0_main.py', *config['scoring']], check=True)
    subprocess.run([sys.executable, 'contract_main.py', *config['contract']], check=True)


if __name__ == '__main__':
    main()
