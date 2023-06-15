from solid2.core.object_base import OpenSCADConstant as _OpenSCADConstant
from solid2.core.scad_import import extra_scad_include as _extra_scad_include
from pathlib import Path as _Path

from .bosl2_base import Bosl2Base as _Bosl2Base

_extra_scad_include(f"{_Path(__file__).parent.parent / '../libs/BOSL2/hingesnaps.scad'}", use_not_include=False)

class apply_folding_hinges_and_snaps(_Bosl2Base):
    def __init__(self, thick=None, foldangle=None, hinges=None, snaps=None, sockets=None, snaplen=None, snapdiam=None, hingegap=None, layerheight=None, **kwargs):
       super().__init__("apply_folding_hinges_and_snaps", {"thick" : thick, "foldangle" : foldangle, "hinges" : hinges, "snaps" : snaps, "sockets" : sockets, "snaplen" : snaplen, "snapdiam" : snapdiam, "hingegap" : hingegap, "layerheight" : layerheight, **kwargs})

class folding_hinge_mask(_Bosl2Base):
    def __init__(self, l=None, thick=None, layerheight=None, foldangle=None, hingegap=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("folding_hinge_mask", {"l" : l, "thick" : thick, "layerheight" : layerheight, "foldangle" : foldangle, "hingegap" : hingegap, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class snap_lock(_Bosl2Base):
    def __init__(self, thick=None, snaplen=None, snapdiam=None, layerheight=None, foldangle=None, hingegap=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("snap_lock", {"thick" : thick, "snaplen" : snaplen, "snapdiam" : snapdiam, "layerheight" : layerheight, "foldangle" : foldangle, "hingegap" : hingegap, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class snap_socket(_Bosl2Base):
    def __init__(self, thick=None, snaplen=None, snapdiam=None, layerheight=None, foldangle=None, hingegap=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("snap_socket", {"thick" : thick, "snaplen" : snaplen, "snapdiam" : snapdiam, "layerheight" : layerheight, "foldangle" : foldangle, "hingegap" : hingegap, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

