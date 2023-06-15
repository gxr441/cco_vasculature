from solid2.core.object_base import OpenSCADConstant as _OpenSCADConstant
from solid2.core.scad_import import extra_scad_include as _extra_scad_include
from pathlib import Path as _Path

from .bosl2_base import Bosl2Base as _Bosl2Base

_extra_scad_include(f"{_Path(__file__).parent.parent / '../libs/BOSL2/gears.scad'}", use_not_include=False)

class spur_gear(_Bosl2Base):
    def __init__(self, pitch=None, teeth=None, thickness=None, shaft_diam=None, hide=None, pressure_angle=None, clearance=None, backlash=None, helical=None, slices=None, interior=None, mod=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("spur_gear", {"pitch" : pitch, "teeth" : teeth, "thickness" : thickness, "shaft_diam" : shaft_diam, "hide" : hide, "pressure_angle" : pressure_angle, "clearance" : clearance, "backlash" : backlash, "helical" : helical, "slices" : slices, "interior" : interior, "mod" : mod, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class spur_gear2d(_Bosl2Base):
    def __init__(self, pitch=None, teeth=None, hide=None, pressure_angle=None, clearance=None, backlash=None, interior=None, mod=None, anchor=None, spin=None, **kwargs):
       super().__init__("spur_gear2d", {"pitch" : pitch, "teeth" : teeth, "hide" : hide, "pressure_angle" : pressure_angle, "clearance" : clearance, "backlash" : backlash, "interior" : interior, "mod" : mod, "anchor" : anchor, "spin" : spin, **kwargs})

class rack(_Bosl2Base):
    def __init__(self, pitch=None, teeth=None, thickness=None, height=None, pressure_angle=None, backlash=None, clearance=None, helical=None, mod=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("rack", {"pitch" : pitch, "teeth" : teeth, "thickness" : thickness, "height" : height, "pressure_angle" : pressure_angle, "backlash" : backlash, "clearance" : clearance, "helical" : helical, "mod" : mod, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class rack2d(_Bosl2Base):
    def __init__(self, pitch=None, teeth=None, height=None, pressure_angle=None, backlash=None, clearance=None, mod=None, anchor=None, spin=None, **kwargs):
       super().__init__("rack2d", {"pitch" : pitch, "teeth" : teeth, "height" : height, "pressure_angle" : pressure_angle, "backlash" : backlash, "clearance" : clearance, "mod" : mod, "anchor" : anchor, "spin" : spin, **kwargs})

class bevel_gear(_Bosl2Base):
    def __init__(self, pitch=None, teeth=None, face_width=None, pitch_angle=None, mate_teeth=None, hide=None, pressure_angle=None, clearance=None, backlash=None, cutter_radius=None, spiral_angle=None, left_handed=None, slices=None, interior=None, mod=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("bevel_gear", {"pitch" : pitch, "teeth" : teeth, "face_width" : face_width, "pitch_angle" : pitch_angle, "mate_teeth" : mate_teeth, "hide" : hide, "pressure_angle" : pressure_angle, "clearance" : clearance, "backlash" : backlash, "cutter_radius" : cutter_radius, "spiral_angle" : spiral_angle, "left_handed" : left_handed, "slices" : slices, "interior" : interior, "mod" : mod, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class worm(_Bosl2Base):
    def __init__(self, pitch=None, d=None, l=None, starts=None, left_handed=None, pressure_angle=None, backlash=None, clearance=None, mod=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("worm", {"pitch" : pitch, "d" : d, "l" : l, "starts" : starts, "left_handed" : left_handed, "pressure_angle" : pressure_angle, "backlash" : backlash, "clearance" : clearance, "mod" : mod, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class worm_gear(_Bosl2Base):
    def __init__(self, pitch=None, teeth=None, worm_diam=None, worm_starts=None, worm_arc=None, crowning=None, left_handed=None, pressure_angle=None, backlash=None, clearance=None, mod=None, slices=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("worm_gear", {"pitch" : pitch, "teeth" : teeth, "worm_diam" : worm_diam, "worm_starts" : worm_starts, "worm_arc" : worm_arc, "crowning" : crowning, "left_handed" : left_handed, "pressure_angle" : pressure_angle, "backlash" : backlash, "clearance" : clearance, "mod" : mod, "slices" : slices, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class _gear_tooth_profile(_Bosl2Base):
    def __init__(self, pitch=None, teeth=None, pressure_angle=None, clearance=None, backlash=None, interior=None, valleys=None, center=None, mod=None, **kwargs):
       super().__init__("_gear_tooth_profile", {"pitch" : pitch, "teeth" : teeth, "pressure_angle" : pressure_angle, "clearance" : clearance, "backlash" : backlash, "interior" : interior, "valleys" : valleys, "center" : center, "mod" : mod, **kwargs})

class circular_pitch(_Bosl2Base):
    def __init__(self, pitch=None, mod=None, **kwargs):
       super().__init__("circular_pitch", {"pitch" : pitch, "mod" : mod, **kwargs})

class diametral_pitch(_Bosl2Base):
    def __init__(self, pitch=None, mod=None, **kwargs):
       super().__init__("diametral_pitch", {"pitch" : pitch, "mod" : mod, **kwargs})

class pitch_value(_Bosl2Base):
    def __init__(self, mod=None, **kwargs):
       super().__init__("pitch_value", {"mod" : mod, **kwargs})

class module_value(_Bosl2Base):
    def __init__(self, pitch=None, **kwargs):
       super().__init__("module_value", {"pitch" : pitch, **kwargs})

class _adendum(_Bosl2Base):
    def __init__(self, pitch=None, mod=None, **kwargs):
       super().__init__("_adendum", {"pitch" : pitch, "mod" : mod, **kwargs})

class _dedendum(_Bosl2Base):
    def __init__(self, pitch=None, clearance=None, mod=None, **kwargs):
       super().__init__("_dedendum", {"pitch" : pitch, "clearance" : clearance, "mod" : mod, **kwargs})

class pitch_radius(_Bosl2Base):
    def __init__(self, pitch=None, teeth=None, mod=None, **kwargs):
       super().__init__("pitch_radius", {"pitch" : pitch, "teeth" : teeth, "mod" : mod, **kwargs})

class outer_radius(_Bosl2Base):
    def __init__(self, pitch=None, teeth=None, clearance=None, interior=None, mod=None, **kwargs):
       super().__init__("outer_radius", {"pitch" : pitch, "teeth" : teeth, "clearance" : clearance, "interior" : interior, "mod" : mod, **kwargs})

class _root_radius(_Bosl2Base):
    def __init__(self, pitch=None, teeth=None, clearance=None, interior=None, mod=None, **kwargs):
       super().__init__("_root_radius", {"pitch" : pitch, "teeth" : teeth, "clearance" : clearance, "interior" : interior, "mod" : mod, **kwargs})

class _base_radius(_Bosl2Base):
    def __init__(self, pitch=None, teeth=None, pressure_angle=None, mod=None, **kwargs):
       super().__init__("_base_radius", {"pitch" : pitch, "teeth" : teeth, "pressure_angle" : pressure_angle, "mod" : mod, **kwargs})

class bevel_pitch_angle(_Bosl2Base):
    def __init__(self, teeth=None, mate_teeth=None, drive_angle=None, **kwargs):
       super().__init__("bevel_pitch_angle", {"teeth" : teeth, "mate_teeth" : mate_teeth, "drive_angle" : drive_angle, **kwargs})

class worm_gear_thickness(_Bosl2Base):
    def __init__(self, pitch=None, teeth=None, worm_diam=None, worm_arc=None, crowning=None, clearance=None, mod=None, **kwargs):
       super().__init__("worm_gear_thickness", {"pitch" : pitch, "teeth" : teeth, "worm_diam" : worm_diam, "worm_arc" : worm_arc, "crowning" : crowning, "clearance" : clearance, "mod" : mod, **kwargs})

class _gear_polar(_Bosl2Base):
    def __init__(self, r=None, t=None, **kwargs):
       super().__init__("_gear_polar", {"r" : r, "t" : t, **kwargs})

class _gear_iang(_Bosl2Base):
    def __init__(self, r1=None, r2=None, **kwargs):
       super().__init__("_gear_iang", {"r1" : r1, "r2" : r2, **kwargs})

class _gear_q6(_Bosl2Base):
    def __init__(self, b=None, s=None, t=None, d=None, **kwargs):
       super().__init__("_gear_q6", {"b" : b, "s" : s, "t" : t, "d" : d, **kwargs})

class _gear_q7(_Bosl2Base):
    def __init__(self, f=None, r=None, b=None, r2=None, t=None, s=None, **kwargs):
       super().__init__("_gear_q7", {"f" : f, "r" : r, "b" : b, "r2" : r2, "t" : t, "s" : s, **kwargs})

class spur_gear(_Bosl2Base):
    def __init__(self, pitch=None, teeth=None, thickness=None, shaft_diam=None, hide=None, pressure_angle=None, clearance=None, backlash=None, helical=None, slices=None, interior=None, mod=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("spur_gear", {"pitch" : pitch, "teeth" : teeth, "thickness" : thickness, "shaft_diam" : shaft_diam, "hide" : hide, "pressure_angle" : pressure_angle, "clearance" : clearance, "backlash" : backlash, "helical" : helical, "slices" : slices, "interior" : interior, "mod" : mod, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class spur_gear2d(_Bosl2Base):
    def __init__(self, pitch=None, teeth=None, hide=None, pressure_angle=None, clearance=None, backlash=None, interior=None, mod=None, anchor=None, spin=None, **kwargs):
       super().__init__("spur_gear2d", {"pitch" : pitch, "teeth" : teeth, "hide" : hide, "pressure_angle" : pressure_angle, "clearance" : clearance, "backlash" : backlash, "interior" : interior, "mod" : mod, "anchor" : anchor, "spin" : spin, **kwargs})

class rack(_Bosl2Base):
    def __init__(self, pitch=None, teeth=None, thickness=None, height=None, pressure_angle=None, backlash=None, clearance=None, helical=None, mod=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("rack", {"pitch" : pitch, "teeth" : teeth, "thickness" : thickness, "height" : height, "pressure_angle" : pressure_angle, "backlash" : backlash, "clearance" : clearance, "helical" : helical, "mod" : mod, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class rack2d(_Bosl2Base):
    def __init__(self, pitch=None, teeth=None, height=None, pressure_angle=None, backlash=None, clearance=None, mod=None, anchor=None, spin=None, **kwargs):
       super().__init__("rack2d", {"pitch" : pitch, "teeth" : teeth, "height" : height, "pressure_angle" : pressure_angle, "backlash" : backlash, "clearance" : clearance, "mod" : mod, "anchor" : anchor, "spin" : spin, **kwargs})

class bevel_gear(_Bosl2Base):
    def __init__(self, pitch=None, teeth=None, face_width=None, pitch_angle=None, mate_teeth=None, shaft_diam=None, hide=None, pressure_angle=None, clearance=None, backlash=None, cutter_radius=None, spiral_angle=None, left_handed=None, slices=None, interior=None, mod=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("bevel_gear", {"pitch" : pitch, "teeth" : teeth, "face_width" : face_width, "pitch_angle" : pitch_angle, "mate_teeth" : mate_teeth, "shaft_diam" : shaft_diam, "hide" : hide, "pressure_angle" : pressure_angle, "clearance" : clearance, "backlash" : backlash, "cutter_radius" : cutter_radius, "spiral_angle" : spiral_angle, "left_handed" : left_handed, "slices" : slices, "interior" : interior, "mod" : mod, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class worm(_Bosl2Base):
    def __init__(self, pitch=None, d=None, l=None, starts=None, left_handed=None, pressure_angle=None, backlash=None, clearance=None, mod=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("worm", {"pitch" : pitch, "d" : d, "l" : l, "starts" : starts, "left_handed" : left_handed, "pressure_angle" : pressure_angle, "backlash" : backlash, "clearance" : clearance, "mod" : mod, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class worm_gear(_Bosl2Base):
    def __init__(self, pitch=None, teeth=None, worm_diam=None, worm_starts=None, worm_arc=None, crowning=None, left_handed=None, pressure_angle=None, backlash=None, slices=None, clearance=None, mod=None, shaft_diam=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("worm_gear", {"pitch" : pitch, "teeth" : teeth, "worm_diam" : worm_diam, "worm_starts" : worm_starts, "worm_arc" : worm_arc, "crowning" : crowning, "left_handed" : left_handed, "pressure_angle" : pressure_angle, "backlash" : backlash, "slices" : slices, "clearance" : clearance, "mod" : mod, "shaft_diam" : shaft_diam, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class _gear_tooth_profile(_Bosl2Base):
    def __init__(self, pitch=None, teeth=None, pressure_angle=None, backlash=None, clearance=None, interior=None, valleys=None, center=None, mod=None, **kwargs):
       super().__init__("_gear_tooth_profile", {"pitch" : pitch, "teeth" : teeth, "pressure_angle" : pressure_angle, "backlash" : backlash, "clearance" : clearance, "interior" : interior, "valleys" : valleys, "center" : center, "mod" : mod, **kwargs})

