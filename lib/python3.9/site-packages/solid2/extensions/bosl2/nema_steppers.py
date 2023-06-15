from solid2.core.object_base import OpenSCADConstant as _OpenSCADConstant
from solid2.core.scad_import import extra_scad_include as _extra_scad_include
from pathlib import Path as _Path

from .bosl2_base import Bosl2Base as _Bosl2Base

_extra_scad_include(f"{_Path(__file__).parent.parent / '../libs/BOSL2/nema_steppers.scad'}", use_not_include=False)

class nema_motor_width(_Bosl2Base):
    def __init__(self, size=None, **kwargs):
       super().__init__("nema_motor_width", {"size" : size, **kwargs})

class nema_motor_plinth_height(_Bosl2Base):
    def __init__(self, size=None, **kwargs):
       super().__init__("nema_motor_plinth_height", {"size" : size, **kwargs})

class nema_motor_plinth_diam(_Bosl2Base):
    def __init__(self, size=None, **kwargs):
       super().__init__("nema_motor_plinth_diam", {"size" : size, **kwargs})

class nema_motor_screw_spacing(_Bosl2Base):
    def __init__(self, size=None, **kwargs):
       super().__init__("nema_motor_screw_spacing", {"size" : size, **kwargs})

class nema_motor_screw_size(_Bosl2Base):
    def __init__(self, size=None, **kwargs):
       super().__init__("nema_motor_screw_size", {"size" : size, **kwargs})

class nema_motor_screw_depth(_Bosl2Base):
    def __init__(self, size=None, **kwargs):
       super().__init__("nema_motor_screw_depth", {"size" : size, **kwargs})

class nema11_stepper(_Bosl2Base):
    def __init__(self, h=None, shaft=None, shaft_len=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("nema11_stepper", {"h" : h, "shaft" : shaft, "shaft_len" : shaft_len, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class nema14_stepper(_Bosl2Base):
    def __init__(self, h=None, shaft=None, shaft_len=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("nema14_stepper", {"h" : h, "shaft" : shaft, "shaft_len" : shaft_len, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class nema17_stepper(_Bosl2Base):
    def __init__(self, h=None, shaft=None, shaft_len=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("nema17_stepper", {"h" : h, "shaft" : shaft, "shaft_len" : shaft_len, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class nema23_stepper(_Bosl2Base):
    def __init__(self, h=None, shaft=None, shaft_len=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("nema23_stepper", {"h" : h, "shaft" : shaft, "shaft_len" : shaft_len, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class nema34_stepper(_Bosl2Base):
    def __init__(self, h=None, shaft=None, shaft_len=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("nema34_stepper", {"h" : h, "shaft" : shaft, "shaft_len" : shaft_len, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class nema_mount_holes(_Bosl2Base):
    def __init__(self, size=None, depth=None, l=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("nema_mount_holes", {"size" : size, "depth" : depth, "l" : l, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class nema11_mount_holes(_Bosl2Base):
    def __init__(self, depth=None, l=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("nema11_mount_holes", {"depth" : depth, "l" : l, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class nema14_mount_holes(_Bosl2Base):
    def __init__(self, depth=None, l=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("nema14_mount_holes", {"depth" : depth, "l" : l, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class nema17_mount_holes(_Bosl2Base):
    def __init__(self, depth=None, l=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("nema17_mount_holes", {"depth" : depth, "l" : l, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class nema23_mount_holes(_Bosl2Base):
    def __init__(self, depth=None, l=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("nema23_mount_holes", {"depth" : depth, "l" : l, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class nema34_mount_holes(_Bosl2Base):
    def __init__(self, depth=None, l=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("nema34_mount_holes", {"depth" : depth, "l" : l, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

