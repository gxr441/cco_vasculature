#! /usr/bin/env python

import os
from pathlib import Path

p = Path(__file__).parent

for f in sorted(p.iterdir()):
    if f.suffix == ".py" and f.name[0].isdigit() \
                 and f.name[1].isdigit() and f.name[2] == '-':
        print(f)
        os.system(f'python3 {f.as_posix()}')
