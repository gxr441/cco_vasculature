#this is global config file for exp_solid
from pathlib import Path
import os
import platform
import sys
import re

class Config:
    def __init__(self):
        self.use_implicit_builtins = "--implicit" in sys.argv

        self.enable_pickle_cache = True
        self.pickle_cache_dir = self.get_pickle_cache_dir()

        self.openscad_library_paths = self.get_openscad_library_paths()

    def get_openscad_library_paths(self):
        """
        Return system-dependent OpenSCAD library paths or paths defined in
        os.environ['OPENSCADPATH'] """
        paths = [Path('.')]

        user_path = os.environ.get('OPENSCADPATH')
        if user_path:
            for s in re.split(r'\s*[;:]\s*', user_path):
                paths.append(Path(s))

        #user wide path
        default_paths = {
            'Linux':   Path.home() / '.local/share/OpenSCAD/libraries',
            'Darwin':  Path.home() / 'Documents/OpenSCAD/libraries',
            'Windows': Path('My Documents/OpenSCAD/libraries')
        }

        paths.append(default_paths[platform.system()])

        #system wide paths
        if platform.system() == 'Linux':
            #sorry, but I've no clue what the paths are on other operating systems
            paths.append(Path("/usr/share/openscad/libraries"))

        return paths

    def get_pickle_cache_dir(self):
        default_paths = {
            'Linux':   Path.home() / '.local/share/expSolidPython/pickle_cache',
            'Darwin':  Path.home() / 'Documents/expSolidPython/pickle_cache',
            'Windows': Path('My Documents/expSolidPython/pickle_cache')
        }
        return default_paths[platform.system()]


config = Config()

