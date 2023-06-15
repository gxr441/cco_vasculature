from solid2.core.object_base import OpenSCADConstant as _OpenSCADConstant
from solid2.core.scad_import import extra_scad_include as _extra_scad_include
from pathlib import Path as _Path

from .bosl2_base import Bosl2Base as _Bosl2Base

_extra_scad_include(f"{_Path(__file__).parent.parent / '../libs/BOSL2/constants.scad'}", use_not_include=False)

_UNDEF = _OpenSCADConstant('_UNDEF')
INCH = _OpenSCADConstant('INCH')
LEFT = _OpenSCADConstant('LEFT')
RIGHT = _OpenSCADConstant('RIGHT')
FRONT = _OpenSCADConstant('FRONT')
FWD = _OpenSCADConstant('FWD')
FORWARD = _OpenSCADConstant('FORWARD')
BACK = _OpenSCADConstant('BACK')
BOTTOM = _OpenSCADConstant('BOTTOM')
BOT = _OpenSCADConstant('BOT')
DOWN = _OpenSCADConstant('DOWN')
TOP = _OpenSCADConstant('TOP')
UP = _OpenSCADConstant('UP')
CENTER = _OpenSCADConstant('CENTER')
CTR = _OpenSCADConstant('CTR')
CENTRE = _OpenSCADConstant('CENTRE')
SEGMENT = _OpenSCADConstant('SEGMENT')
RAY = _OpenSCADConstant('RAY')
LINE = _OpenSCADConstant('LINE')
IDENT = _OpenSCADConstant('IDENT')
class get_slop(_Bosl2Base):
    def __init__(self, **kwargs):
       super().__init__("get_slop", {**kwargs})

