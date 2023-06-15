from solid2.core.object_base import OpenSCADConstant as _OpenSCADConstant
from solid2.core.scad_import import extra_scad_include as _extra_scad_include
from pathlib import Path as _Path

from .bosl2_base import Bosl2Base as _Bosl2Base

_extra_scad_include(f"{_Path(__file__).parent.parent / '../libs/BOSL2/walls.scad'}", use_not_include=False)

class sparse_wall(_Bosl2Base):
    def __init__(self, h=None, l=None, thick=None, maxang=None, strut=None, max_bridge=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("sparse_wall", {"h" : h, "l" : l, "thick" : thick, "maxang" : maxang, "strut" : strut, "max_bridge" : max_bridge, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class corrugated_wall(_Bosl2Base):
    def __init__(self, h=None, l=None, thick=None, strut=None, wall=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("corrugated_wall", {"h" : h, "l" : l, "thick" : thick, "strut" : strut, "wall" : wall, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class thinning_wall(_Bosl2Base):
    def __init__(self, h=None, l=None, thick=None, ang=None, braces=None, strut=None, wall=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("thinning_wall", {"h" : h, "l" : l, "thick" : thick, "ang" : ang, "braces" : braces, "strut" : strut, "wall" : wall, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class thinning_triangle(_Bosl2Base):
    def __init__(self, h=None, l=None, thick=None, ang=None, strut=None, wall=None, diagonly=None, center=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("thinning_triangle", {"h" : h, "l" : l, "thick" : thick, "ang" : ang, "strut" : strut, "wall" : wall, "diagonly" : diagonly, "center" : center, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class narrowing_strut(_Bosl2Base):
    def __init__(self, w=None, l=None, wall=None, ang=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("narrowing_strut", {"w" : w, "l" : l, "wall" : wall, "ang" : ang, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

