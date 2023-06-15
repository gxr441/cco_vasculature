from solid2.core.object_base import OpenSCADConstant as _OpenSCADConstant
from solid2.core.scad_import import extra_scad_include as _extra_scad_include
from pathlib import Path as _Path

from .bosl2_base import Bosl2Base as _Bosl2Base

_extra_scad_include(f"{_Path(__file__).parent.parent / '../libs/BOSL2/screw_drive.scad'}", use_not_include=False)

class _phillips_shaft(_Bosl2Base):
    def __init__(self, x=None, **kwargs):
       super().__init__("_phillips_shaft", {"x" : x, **kwargs})

class _ph_bot_angle(_Bosl2Base):
    def __init__(self, **kwargs):
       super().__init__("_ph_bot_angle", {**kwargs})

class _ph_side_angle(_Bosl2Base):
    def __init__(self, **kwargs):
       super().__init__("_ph_side_angle", {**kwargs})

class phillips_depth(_Bosl2Base):
    def __init__(self, size=None, d=None, **kwargs):
       super().__init__("phillips_depth", {"size" : size, "d" : d, **kwargs})

class phillips_diam(_Bosl2Base):
    def __init__(self, size=None, depth=None, **kwargs):
       super().__init__("phillips_diam", {"size" : size, "depth" : depth, **kwargs})

class torx_diam(_Bosl2Base):
    def __init__(self, size=None, **kwargs):
       super().__init__("torx_diam", {"size" : size, **kwargs})

class _torx_inner_diam(_Bosl2Base):
    def __init__(self, size=None, **kwargs):
       super().__init__("_torx_inner_diam", {"size" : size, **kwargs})

class torx_depth(_Bosl2Base):
    def __init__(self, size=None, **kwargs):
       super().__init__("torx_depth", {"size" : size, **kwargs})

class _torx_tip_radius(_Bosl2Base):
    def __init__(self, size=None, **kwargs):
       super().__init__("_torx_tip_radius", {"size" : size, **kwargs})

class _torx_rounding_radius(_Bosl2Base):
    def __init__(self, size=None, **kwargs):
       super().__init__("_torx_rounding_radius", {"size" : size, **kwargs})

class phillips_mask(_Bosl2Base):
    def __init__(self, size=None, _fn=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("phillips_mask", {"size" : size, "_fn" : _fn, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class torx_mask(_Bosl2Base):
    def __init__(self, size=None, l=None, center=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("torx_mask", {"size" : size, "l" : l, "center" : center, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class torx_mask2d(_Bosl2Base):
    def __init__(self, size=None, **kwargs):
       super().__init__("torx_mask2d", {"size" : size, **kwargs})

class robertson_mask(_Bosl2Base):
    def __init__(self, size=None, extra=None, ang=None, **kwargs):
       super().__init__("robertson_mask", {"size" : size, "extra" : extra, "ang" : ang, **kwargs})

