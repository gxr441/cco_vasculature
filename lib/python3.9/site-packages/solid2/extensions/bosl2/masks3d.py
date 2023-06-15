from solid2.core.object_base import OpenSCADConstant as _OpenSCADConstant
from solid2.core.scad_import import extra_scad_include as _extra_scad_include
from pathlib import Path as _Path

from .bosl2_base import Bosl2Base as _Bosl2Base

_extra_scad_include(f"{_Path(__file__).parent.parent / '../libs/BOSL2/masks3d.scad'}", use_not_include=False)

class chamfer_edge_mask(_Bosl2Base):
    def __init__(self, l=None, chamfer=None, excess=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("chamfer_edge_mask", {"l" : l, "chamfer" : chamfer, "excess" : excess, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class chamfer_corner_mask(_Bosl2Base):
    def __init__(self, chamfer=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("chamfer_corner_mask", {"chamfer" : chamfer, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class chamfer_cylinder_mask(_Bosl2Base):
    def __init__(self, r=None, chamfer=None, d=None, ang=None, from_end=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("chamfer_cylinder_mask", {"r" : r, "chamfer" : chamfer, "d" : d, "ang" : ang, "from_end" : from_end, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class rounding_edge_mask(_Bosl2Base):
    def __init__(self, l=None, r=None, r1=None, r2=None, d=None, d1=None, d2=None, excess=None, anchor=None, spin=None, orient=None, h=None, **kwargs):
       super().__init__("rounding_edge_mask", {"l" : l, "r" : r, "r1" : r1, "r2" : r2, "d" : d, "d1" : d1, "d2" : d2, "excess" : excess, "anchor" : anchor, "spin" : spin, "orient" : orient, "h" : h, **kwargs})

class rounding_corner_mask(_Bosl2Base):
    def __init__(self, r=None, d=None, style=None, excess=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("rounding_corner_mask", {"r" : r, "d" : d, "style" : style, "excess" : excess, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class rounding_angled_edge_mask(_Bosl2Base):
    def __init__(self, h=None, r=None, r1=None, r2=None, d=None, d1=None, d2=None, ang=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("rounding_angled_edge_mask", {"h" : h, "r" : r, "r1" : r1, "r2" : r2, "d" : d, "d1" : d1, "d2" : d2, "ang" : ang, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class rounding_angled_corner_mask(_Bosl2Base):
    def __init__(self, r=None, ang=None, d=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("rounding_angled_corner_mask", {"r" : r, "ang" : ang, "d" : d, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class rounding_cylinder_mask(_Bosl2Base):
    def __init__(self, r=None, rounding=None, d=None, **kwargs):
       super().__init__("rounding_cylinder_mask", {"r" : r, "rounding" : rounding, "d" : d, **kwargs})

class rounding_hole_mask(_Bosl2Base):
    def __init__(self, r=None, rounding=None, excess=None, d=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("rounding_hole_mask", {"r" : r, "rounding" : rounding, "excess" : excess, "d" : d, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class teardrop_edge_mask(_Bosl2Base):
    def __init__(self, l=None, r=None, angle=None, excess=None, d=None, **kwargs):
       super().__init__("teardrop_edge_mask", {"l" : l, "r" : r, "angle" : angle, "excess" : excess, "d" : d, **kwargs})

class teardrop_corner_mask(_Bosl2Base):
    def __init__(self, r=None, angle=None, excess=None, d=None, **kwargs):
       super().__init__("teardrop_corner_mask", {"r" : r, "angle" : angle, "excess" : excess, "d" : d, **kwargs})

class chamfer_edge_mask(_Bosl2Base):
    def __init__(self, l=None, chamfer=None, excess=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("chamfer_edge_mask", {"l" : l, "chamfer" : chamfer, "excess" : excess, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class chamfer_corner_mask(_Bosl2Base):
    def __init__(self, chamfer=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("chamfer_corner_mask", {"chamfer" : chamfer, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class chamfer_cylinder_mask(_Bosl2Base):
    def __init__(self, r=None, chamfer=None, d=None, ang=None, from_end=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("chamfer_cylinder_mask", {"r" : r, "chamfer" : chamfer, "d" : d, "ang" : ang, "from_end" : from_end, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class rounding_edge_mask(_Bosl2Base):
    def __init__(self, l=None, r=None, r1=None, r2=None, d=None, d1=None, d2=None, excess=None, anchor=None, spin=None, orient=None, h=None, **kwargs):
       super().__init__("rounding_edge_mask", {"l" : l, "r" : r, "r1" : r1, "r2" : r2, "d" : d, "d1" : d1, "d2" : d2, "excess" : excess, "anchor" : anchor, "spin" : spin, "orient" : orient, "h" : h, **kwargs})

class rounding_corner_mask(_Bosl2Base):
    def __init__(self, r=None, d=None, style=None, excess=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("rounding_corner_mask", {"r" : r, "d" : d, "style" : style, "excess" : excess, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class rounding_angled_edge_mask(_Bosl2Base):
    def __init__(self, h=None, r=None, r1=None, r2=None, d=None, d1=None, d2=None, ang=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("rounding_angled_edge_mask", {"h" : h, "r" : r, "r1" : r1, "r2" : r2, "d" : d, "d1" : d1, "d2" : d2, "ang" : ang, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class rounding_angled_corner_mask(_Bosl2Base):
    def __init__(self, r=None, ang=None, d=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("rounding_angled_corner_mask", {"r" : r, "ang" : ang, "d" : d, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class rounding_cylinder_mask(_Bosl2Base):
    def __init__(self, r=None, rounding=None, d=None, **kwargs):
       super().__init__("rounding_cylinder_mask", {"r" : r, "rounding" : rounding, "d" : d, **kwargs})

class rounding_hole_mask(_Bosl2Base):
    def __init__(self, r=None, rounding=None, excess=None, d=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("rounding_hole_mask", {"r" : r, "rounding" : rounding, "excess" : excess, "d" : d, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class teardrop_edge_mask(_Bosl2Base):
    def __init__(self, l=None, r=None, angle=None, excess=None, d=None, **kwargs):
       super().__init__("teardrop_edge_mask", {"l" : l, "r" : r, "angle" : angle, "excess" : excess, "d" : d, **kwargs})

class teardrop_corner_mask(_Bosl2Base):
    def __init__(self, r=None, angle=None, excess=None, d=None, **kwargs):
       super().__init__("teardrop_corner_mask", {"r" : r, "angle" : angle, "excess" : excess, "d" : d, **kwargs})

