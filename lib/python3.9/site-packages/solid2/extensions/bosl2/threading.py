from solid2.core.object_base import OpenSCADConstant as _OpenSCADConstant
from solid2.core.scad_import import extra_scad_include as _extra_scad_include
from pathlib import Path as _Path

from .bosl2_base import Bosl2Base as _Bosl2Base

_extra_scad_include(f"{_Path(__file__).parent.parent / '../libs/BOSL2/threading.scad'}", use_not_include=False)

class threaded_rod(_Bosl2Base):
    def __init__(self, d=None, l=None, pitch=None, left_handed=None, bevel=None, bevel1=None, bevel2=None, starts=None, internal=None, d1=None, d2=None, higbee=None, higbee1=None, higbee2=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("threaded_rod", {"d" : d, "l" : l, "pitch" : pitch, "left_handed" : left_handed, "bevel" : bevel, "bevel1" : bevel1, "bevel2" : bevel2, "starts" : starts, "internal" : internal, "d1" : d1, "d2" : d2, "higbee" : higbee, "higbee1" : higbee1, "higbee2" : higbee2, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class threaded_nut(_Bosl2Base):
    def __init__(self, od=None, id=None, h=None, pitch=None, starts=None, left_handed=None, bevel=None, bevel1=None, bevel2=None, id1=None, id2=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("threaded_nut", {"od" : od, "id" : id, "h" : h, "pitch" : pitch, "starts" : starts, "left_handed" : left_handed, "bevel" : bevel, "bevel1" : bevel1, "bevel2" : bevel2, "id1" : id1, "id2" : id2, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class trapezoidal_threaded_rod(_Bosl2Base):
    def __init__(self, d=None, l=None, pitch=None, thread_angle=None, thread_depth=None, left_handed=None, bevel=None, bevel1=None, bevel2=None, starts=None, internal=None, higbee=None, higbee1=None, higbee2=None, d1=None, d2=None, center=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("trapezoidal_threaded_rod", {"d" : d, "l" : l, "pitch" : pitch, "thread_angle" : thread_angle, "thread_depth" : thread_depth, "left_handed" : left_handed, "bevel" : bevel, "bevel1" : bevel1, "bevel2" : bevel2, "starts" : starts, "internal" : internal, "higbee" : higbee, "higbee1" : higbee1, "higbee2" : higbee2, "d1" : d1, "d2" : d2, "center" : center, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class trapezoidal_threaded_rod(_Bosl2Base):
    def __init__(self, d=None, l=None, pitch=None, thread_angle=None, thread_depth=None, left_handed=None, bevel=None, bevel1=None, bevel2=None, starts=None, internal=None, higbee=None, higbee1=None, higbee2=None, d1=None, d2=None, center=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("trapezoidal_threaded_rod", {"d" : d, "l" : l, "pitch" : pitch, "thread_angle" : thread_angle, "thread_depth" : thread_depth, "left_handed" : left_handed, "bevel" : bevel, "bevel1" : bevel1, "bevel2" : bevel2, "starts" : starts, "internal" : internal, "higbee" : higbee, "higbee1" : higbee1, "higbee2" : higbee2, "d1" : d1, "d2" : d2, "center" : center, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class acme_threaded_rod(_Bosl2Base):
    def __init__(self, d=None, l=None, tpi=None, pitch=None, starts=None, left_handed=None, bevel=None, bevel1=None, bevel2=None, internal=None, higbee=None, higbee1=None, higbee2=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("acme_threaded_rod", {"d" : d, "l" : l, "tpi" : tpi, "pitch" : pitch, "starts" : starts, "left_handed" : left_handed, "bevel" : bevel, "bevel1" : bevel1, "bevel2" : bevel2, "internal" : internal, "higbee" : higbee, "higbee1" : higbee1, "higbee2" : higbee2, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class acme_threaded_nut(_Bosl2Base):
    def __init__(self, od=None, id=None, h=None, tpi=None, pitch=None, starts=None, left_handed=None, bevel=None, bevel1=None, bevel2=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("acme_threaded_nut", {"od" : od, "id" : id, "h" : h, "tpi" : tpi, "pitch" : pitch, "starts" : starts, "left_handed" : left_handed, "bevel" : bevel, "bevel1" : bevel1, "bevel2" : bevel2, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class npt_threaded_rod(_Bosl2Base):
    def __init__(self, size=None, left_handed=None, bevel=None, bevel1=None, bevel2=None, hollow=None, internal=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("npt_threaded_rod", {"size" : size, "left_handed" : left_handed, "bevel" : bevel, "bevel1" : bevel1, "bevel2" : bevel2, "hollow" : hollow, "internal" : internal, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class buttress_threaded_rod(_Bosl2Base):
    def __init__(self, d=None, l=None, pitch=None, left_handed=None, bevel=None, bevel1=None, bevel2=None, internal=None, higbee=None, higbee1=None, higbee2=None, d1=None, d2=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("buttress_threaded_rod", {"d" : d, "l" : l, "pitch" : pitch, "left_handed" : left_handed, "bevel" : bevel, "bevel1" : bevel1, "bevel2" : bevel2, "internal" : internal, "higbee" : higbee, "higbee1" : higbee1, "higbee2" : higbee2, "d1" : d1, "d2" : d2, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class buttress_threaded_nut(_Bosl2Base):
    def __init__(self, od=None, id=None, h=None, pitch=None, left_handed=None, bevel=None, bevel1=None, bevel2=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("buttress_threaded_nut", {"od" : od, "id" : id, "h" : h, "pitch" : pitch, "left_handed" : left_handed, "bevel" : bevel, "bevel1" : bevel1, "bevel2" : bevel2, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class square_threaded_rod(_Bosl2Base):
    def __init__(self, d=None, l=None, pitch=None, left_handed=None, bevel=None, bevel1=None, bevel2=None, starts=None, internal=None, higbee=None, higbee1=None, higbee2=None, d1=None, d2=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("square_threaded_rod", {"d" : d, "l" : l, "pitch" : pitch, "left_handed" : left_handed, "bevel" : bevel, "bevel1" : bevel1, "bevel2" : bevel2, "starts" : starts, "internal" : internal, "higbee" : higbee, "higbee1" : higbee1, "higbee2" : higbee2, "d1" : d1, "d2" : d2, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class square_threaded_nut(_Bosl2Base):
    def __init__(self, od=None, id=None, h=None, pitch=None, left_handed=None, bevel=None, bevel1=None, bevel2=None, starts=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("square_threaded_nut", {"od" : od, "id" : id, "h" : h, "pitch" : pitch, "left_handed" : left_handed, "bevel" : bevel, "bevel1" : bevel1, "bevel2" : bevel2, "starts" : starts, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class ball_screw_rod(_Bosl2Base):
    def __init__(self, d=None, l=None, pitch=None, ball_diam=None, ball_arc=None, starts=None, left_handed=None, internal=None, bevel=None, bevel1=None, bevel2=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("ball_screw_rod", {"d" : d, "l" : l, "pitch" : pitch, "ball_diam" : ball_diam, "ball_arc" : ball_arc, "starts" : starts, "left_handed" : left_handed, "internal" : internal, "bevel" : bevel, "bevel1" : bevel1, "bevel2" : bevel2, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class generic_threaded_rod(_Bosl2Base):
    def __init__(self, d=None, l=None, pitch=None, profile=None, left_handed=None, bevel=None, bevel1=None, bevel2=None, starts=None, internal=None, d1=None, d2=None, higbee=None, higbee1=None, higbee2=None, center=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("generic_threaded_rod", {"d" : d, "l" : l, "pitch" : pitch, "profile" : profile, "left_handed" : left_handed, "bevel" : bevel, "bevel1" : bevel1, "bevel2" : bevel2, "starts" : starts, "internal" : internal, "d1" : d1, "d2" : d2, "higbee" : higbee, "higbee1" : higbee1, "higbee2" : higbee2, "center" : center, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class generic_threaded_nut(_Bosl2Base):
    def __init__(self, od=None, id=None, h=None, pitch=None, profile=None, left_handed=None, starts=None, bevel=None, bevel1=None, bevel2=None, id1=None, id2=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("generic_threaded_nut", {"od" : od, "id" : id, "h" : h, "pitch" : pitch, "profile" : profile, "left_handed" : left_handed, "starts" : starts, "bevel" : bevel, "bevel1" : bevel1, "bevel2" : bevel2, "id1" : id1, "id2" : id2, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class thread_helix(_Bosl2Base):
    def __init__(self, d=None, pitch=None, thread_depth=None, flank_angle=None, turns=None, profile=None, starts=None, left_handed=None, internal=None, d1=None, d2=None, higbee=None, higbee1=None, higbee2=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("thread_helix", {"d" : d, "pitch" : pitch, "thread_depth" : thread_depth, "flank_angle" : flank_angle, "turns" : turns, "profile" : profile, "starts" : starts, "left_handed" : left_handed, "internal" : internal, "d1" : d1, "d2" : d2, "higbee" : higbee, "higbee1" : higbee1, "higbee2" : higbee2, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class threaded_rod(_Bosl2Base):
    def __init__(self, d=None, l=None, pitch=None, left_handed=None, bevel=None, bevel1=None, bevel2=None, starts=None, internal=None, d1=None, d2=None, higbee=None, higbee1=None, higbee2=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("threaded_rod", {"d" : d, "l" : l, "pitch" : pitch, "left_handed" : left_handed, "bevel" : bevel, "bevel1" : bevel1, "bevel2" : bevel2, "starts" : starts, "internal" : internal, "d1" : d1, "d2" : d2, "higbee" : higbee, "higbee1" : higbee1, "higbee2" : higbee2, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class threaded_nut(_Bosl2Base):
    def __init__(self, od=None, id=None, h=None, pitch=None, starts=None, left_handed=None, bevel=None, bevel1=None, bevel2=None, id1=None, id2=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("threaded_nut", {"od" : od, "id" : id, "h" : h, "pitch" : pitch, "starts" : starts, "left_handed" : left_handed, "bevel" : bevel, "bevel1" : bevel1, "bevel2" : bevel2, "id1" : id1, "id2" : id2, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class trapezoidal_threaded_rod(_Bosl2Base):
    def __init__(self, d=None, l=None, pitch=None, thread_angle=None, thread_depth=None, left_handed=None, bevel=None, bevel1=None, bevel2=None, starts=None, internal=None, higbee=None, higbee1=None, higbee2=None, d1=None, d2=None, center=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("trapezoidal_threaded_rod", {"d" : d, "l" : l, "pitch" : pitch, "thread_angle" : thread_angle, "thread_depth" : thread_depth, "left_handed" : left_handed, "bevel" : bevel, "bevel1" : bevel1, "bevel2" : bevel2, "starts" : starts, "internal" : internal, "higbee" : higbee, "higbee1" : higbee1, "higbee2" : higbee2, "d1" : d1, "d2" : d2, "center" : center, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class trapezoidal_threaded_nut(_Bosl2Base):
    def __init__(self, od=None, id=None, h=None, pitch=None, thread_angle=None, thread_depth=None, left_handed=None, starts=None, bevel=None, bevel1=None, bevel2=None, id1=None, id2=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("trapezoidal_threaded_nut", {"od" : od, "id" : id, "h" : h, "pitch" : pitch, "thread_angle" : thread_angle, "thread_depth" : thread_depth, "left_handed" : left_handed, "starts" : starts, "bevel" : bevel, "bevel1" : bevel1, "bevel2" : bevel2, "id1" : id1, "id2" : id2, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class acme_threaded_rod(_Bosl2Base):
    def __init__(self, d=None, l=None, tpi=None, pitch=None, starts=None, left_handed=None, bevel=None, bevel1=None, bevel2=None, internal=None, higbee=None, higbee1=None, higbee2=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("acme_threaded_rod", {"d" : d, "l" : l, "tpi" : tpi, "pitch" : pitch, "starts" : starts, "left_handed" : left_handed, "bevel" : bevel, "bevel1" : bevel1, "bevel2" : bevel2, "internal" : internal, "higbee" : higbee, "higbee1" : higbee1, "higbee2" : higbee2, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class acme_threaded_nut(_Bosl2Base):
    def __init__(self, od=None, id=None, h=None, tpi=None, pitch=None, starts=None, left_handed=None, bevel=None, bevel1=None, bevel2=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("acme_threaded_nut", {"od" : od, "id" : id, "h" : h, "tpi" : tpi, "pitch" : pitch, "starts" : starts, "left_handed" : left_handed, "bevel" : bevel, "bevel1" : bevel1, "bevel2" : bevel2, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class npt_threaded_rod(_Bosl2Base):
    def __init__(self, size=None, left_handed=None, bevel=None, bevel1=None, bevel2=None, hollow=None, internal=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("npt_threaded_rod", {"size" : size, "left_handed" : left_handed, "bevel" : bevel, "bevel1" : bevel1, "bevel2" : bevel2, "hollow" : hollow, "internal" : internal, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class buttress_threaded_rod(_Bosl2Base):
    def __init__(self, d=None, l=None, pitch=None, left_handed=None, bevel=None, bevel1=None, bevel2=None, internal=None, higbee=None, higbee1=None, higbee2=None, d1=None, d2=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("buttress_threaded_rod", {"d" : d, "l" : l, "pitch" : pitch, "left_handed" : left_handed, "bevel" : bevel, "bevel1" : bevel1, "bevel2" : bevel2, "internal" : internal, "higbee" : higbee, "higbee1" : higbee1, "higbee2" : higbee2, "d1" : d1, "d2" : d2, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class buttress_threaded_nut(_Bosl2Base):
    def __init__(self, od=None, id=None, h=None, pitch=None, left_handed=None, bevel=None, bevel1=None, bevel2=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("buttress_threaded_nut", {"od" : od, "id" : id, "h" : h, "pitch" : pitch, "left_handed" : left_handed, "bevel" : bevel, "bevel1" : bevel1, "bevel2" : bevel2, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class square_threaded_rod(_Bosl2Base):
    def __init__(self, d=None, l=None, pitch=None, left_handed=None, bevel=None, bevel1=None, bevel2=None, starts=None, internal=None, higbee=None, higbee1=None, higbee2=None, d1=None, d2=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("square_threaded_rod", {"d" : d, "l" : l, "pitch" : pitch, "left_handed" : left_handed, "bevel" : bevel, "bevel1" : bevel1, "bevel2" : bevel2, "starts" : starts, "internal" : internal, "higbee" : higbee, "higbee1" : higbee1, "higbee2" : higbee2, "d1" : d1, "d2" : d2, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class square_threaded_nut(_Bosl2Base):
    def __init__(self, od=None, id=None, h=None, pitch=None, left_handed=None, bevel=None, bevel1=None, bevel2=None, starts=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("square_threaded_nut", {"od" : od, "id" : id, "h" : h, "pitch" : pitch, "left_handed" : left_handed, "bevel" : bevel, "bevel1" : bevel1, "bevel2" : bevel2, "starts" : starts, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class ball_screw_rod(_Bosl2Base):
    def __init__(self, d=None, l=None, pitch=None, ball_diam=None, ball_arc=None, starts=None, left_handed=None, internal=None, bevel=None, bevel1=None, bevel2=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("ball_screw_rod", {"d" : d, "l" : l, "pitch" : pitch, "ball_diam" : ball_diam, "ball_arc" : ball_arc, "starts" : starts, "left_handed" : left_handed, "internal" : internal, "bevel" : bevel, "bevel1" : bevel1, "bevel2" : bevel2, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class generic_threaded_rod(_Bosl2Base):
    def __init__(self, d=None, l=None, pitch=None, profile=None, left_handed=None, bevel=None, bevel1=None, bevel2=None, starts=None, internal=None, d1=None, d2=None, higbee=None, higbee1=None, higbee2=None, center=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("generic_threaded_rod", {"d" : d, "l" : l, "pitch" : pitch, "profile" : profile, "left_handed" : left_handed, "bevel" : bevel, "bevel1" : bevel1, "bevel2" : bevel2, "starts" : starts, "internal" : internal, "d1" : d1, "d2" : d2, "higbee" : higbee, "higbee1" : higbee1, "higbee2" : higbee2, "center" : center, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class generic_threaded_nut(_Bosl2Base):
    def __init__(self, od=None, id=None, h=None, pitch=None, profile=None, left_handed=None, starts=None, bevel=None, bevel1=None, bevel2=None, id1=None, id2=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("generic_threaded_nut", {"od" : od, "id" : id, "h" : h, "pitch" : pitch, "profile" : profile, "left_handed" : left_handed, "starts" : starts, "bevel" : bevel, "bevel1" : bevel1, "bevel2" : bevel2, "id1" : id1, "id2" : id2, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class thread_helix(_Bosl2Base):
    def __init__(self, d=None, pitch=None, thread_depth=None, flank_angle=None, turns=None, profile=None, starts=None, left_handed=None, internal=None, d1=None, d2=None, higbee=None, higbee1=None, higbee2=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("thread_helix", {"d" : d, "pitch" : pitch, "thread_depth" : thread_depth, "flank_angle" : flank_angle, "turns" : turns, "profile" : profile, "starts" : starts, "left_handed" : left_handed, "internal" : internal, "d1" : d1, "d2" : d2, "higbee" : higbee, "higbee1" : higbee1, "higbee2" : higbee2, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

