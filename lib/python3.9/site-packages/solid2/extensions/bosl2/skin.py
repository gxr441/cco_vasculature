from solid2.core.object_base import OpenSCADConstant as _OpenSCADConstant
from solid2.core.scad_import import extra_scad_include as _extra_scad_include
from pathlib import Path as _Path

from .bosl2_base import Bosl2Base as _Bosl2Base

_extra_scad_include(f"{_Path(__file__).parent.parent / '../libs/BOSL2/skin.scad'}", use_not_include=False)

_MAP_DIAG = _OpenSCADConstant('_MAP_DIAG')
_MAP_LEFT = _OpenSCADConstant('_MAP_LEFT')
_MAP_UP = _OpenSCADConstant('_MAP_UP')
class skin(_Bosl2Base):
    def __init__(self, profiles=None, slices=None, refine=None, method=None, sampling=None, caps=None, closed=None, z=None, style=None, anchor=None, cp=None, spin=None, orient=None, atype=None, **kwargs):
       super().__init__("skin", {"profiles" : profiles, "slices" : slices, "refine" : refine, "method" : method, "sampling" : sampling, "caps" : caps, "closed" : closed, "z" : z, "style" : style, "anchor" : anchor, "cp" : cp, "spin" : spin, "orient" : orient, "atype" : atype, **kwargs})

class linear_sweep(_Bosl2Base):
    def __init__(self, region=None, height=None, center=None, twist=None, scale=None, shift=None, slices=None, maxseg=None, style=None, cp=None, atype=None, h=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("linear_sweep", {"region" : region, "height" : height, "center" : center, "twist" : twist, "scale" : scale, "shift" : shift, "slices" : slices, "maxseg" : maxseg, "style" : style, "cp" : cp, "atype" : atype, "h" : h, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class _taperfunc(_Bosl2Base):
    def __init__(self, x=None, **kwargs):
       super().__init__("_taperfunc", {"x" : x, **kwargs})

class _ss_polygon_r(_Bosl2Base):
    def __init__(self, N=None, theta=None, **kwargs):
       super().__init__("_ss_polygon_r", {"N" : N, "theta" : theta, **kwargs})

class spiral_sweep(_Bosl2Base):
    def __init__(self, poly=None, h=None, r=None, turns=None, higbee=None, center=None, r1=None, r2=None, d=None, d1=None, d2=None, higbee1=None, higbee2=None, internal=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("spiral_sweep", {"poly" : poly, "h" : h, "r" : r, "turns" : turns, "higbee" : higbee, "center" : center, "r1" : r1, "r2" : r2, "d" : d, "d1" : d1, "d2" : d2, "higbee1" : higbee1, "higbee2" : higbee2, "internal" : internal, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class path_sweep(_Bosl2Base):
    def __init__(self, shape=None, path=None, method=None, normal=None, closed=None, twist=None, twist_by_length=None, symmetry=None, last_normal=None, tangent=None, uniform=None, relaxed=None, caps=None, style=None, transforms=None, anchor=None, cp=None, spin=None, orient=None, atype=None, **kwargs):
       super().__init__("path_sweep", {"shape" : shape, "path" : path, "method" : method, "normal" : normal, "closed" : closed, "twist" : twist, "twist_by_length" : twist_by_length, "symmetry" : symmetry, "last_normal" : last_normal, "tangent" : tangent, "uniform" : uniform, "relaxed" : relaxed, "caps" : caps, "style" : style, "transforms" : transforms, "anchor" : anchor, "cp" : cp, "spin" : spin, "orient" : orient, "atype" : atype, **kwargs})

class path_sweep2d(_Bosl2Base):
    def __init__(self, shape=None, path=None, closed=None, caps=None, quality=None, style=None, anchor=None, cp=None, spin=None, orient=None, atype=None, **kwargs):
       super().__init__("path_sweep2d", {"shape" : shape, "path" : path, "closed" : closed, "caps" : caps, "quality" : quality, "style" : style, "anchor" : anchor, "cp" : cp, "spin" : spin, "orient" : orient, "atype" : atype, **kwargs})

class _ofs_vmap(_Bosl2Base):
    def __init__(self, ofs=None, closed=None, **kwargs):
       super().__init__("_ofs_vmap", {"ofs" : ofs, "closed" : closed, **kwargs})

class _ofs_face_edge(_Bosl2Base):
    def __init__(self, face=None, firstlen=None, second=None, **kwargs):
       super().__init__("_ofs_face_edge", {"face" : face, "firstlen" : firstlen, "second" : second, **kwargs})

class sweep(_Bosl2Base):
    def __init__(self, shape=None, transforms=None, closed=None, caps=None, style=None, anchor=None, cp=None, spin=None, orient=None, atype=None, **kwargs):
       super().__init__("sweep", {"shape" : shape, "transforms" : transforms, "closed" : closed, "caps" : caps, "style" : style, "anchor" : anchor, "cp" : cp, "spin" : spin, "orient" : orient, "atype" : atype, **kwargs})

class subdivide_and_slice(_Bosl2Base):
    def __init__(self, profiles=None, slices=None, numpoints=None, method=None, closed=None, **kwargs):
       super().__init__("subdivide_and_slice", {"profiles" : profiles, "slices" : slices, "numpoints" : numpoints, "method" : method, "closed" : closed, **kwargs})

class slice_profiles(_Bosl2Base):
    def __init__(self, profiles=None, slices=None, closed=None, **kwargs):
       super().__init__("slice_profiles", {"profiles" : profiles, "slices" : slices, "closed" : closed, **kwargs})

class _closest_angle(_Bosl2Base):
    def __init__(self, alpha=None, beta=None, **kwargs):
       super().__init__("_closest_angle", {"alpha" : alpha, "beta" : beta, **kwargs})

class _smooth(_Bosl2Base):
    def __init__(self, data=None, len=None, closed=None, angle=None, **kwargs):
       super().__init__("_smooth", {"data" : data, "len" : len, "closed" : closed, "angle" : angle, **kwargs})

class rot_resample(_Bosl2Base):
    def __init__(self, rotlist=None, n=None, twist=None, scale=None, smoothlen=None, long=None, turns=None, closed=None, method=None, **kwargs):
       super().__init__("rot_resample", {"rotlist" : rotlist, "n" : n, "twist" : twist, "scale" : scale, "smoothlen" : smoothlen, "long" : long, "turns" : turns, "closed" : closed, "method" : method, **kwargs})

class _dp_distance_array(_Bosl2Base):
    def __init__(self, small=None, big=None, abort_thresh=None, **kwargs):
       super().__init__("_dp_distance_array", {"small" : small, "big" : big, "abort_thresh" : abort_thresh, **kwargs})

class _dp_distance_row(_Bosl2Base):
    def __init__(self, small=None, big=None, small_ind=None, tdist=None, **kwargs):
       super().__init__("_dp_distance_row", {"small" : small, "big" : big, "small_ind" : small_ind, "tdist" : tdist, **kwargs})

class _dp_extract_map(_Bosl2Base):
    def __init__(self, map=None, **kwargs):
       super().__init__("_dp_extract_map", {"map" : map, **kwargs})

class _skin_distance_match(_Bosl2Base):
    def __init__(self, poly1=None, poly2=None, **kwargs):
       super().__init__("_skin_distance_match", {"poly1" : poly1, "poly2" : poly2, **kwargs})

class _skin_aligned_distance_match(_Bosl2Base):
    def __init__(self, poly1=None, poly2=None, **kwargs):
       super().__init__("_skin_aligned_distance_match", {"poly1" : poly1, "poly2" : poly2, **kwargs})

class _skin_tangent_match(_Bosl2Base):
    def __init__(self, poly1=None, poly2=None, **kwargs):
       super().__init__("_skin_tangent_match", {"poly1" : poly1, "poly2" : poly2, **kwargs})

class _find_one_tangent(_Bosl2Base):
    def __init__(self, curve=None, edge=None, curve_offset=None, closed=None, **kwargs):
       super().__init__("_find_one_tangent", {"curve" : curve, "edge" : edge, "curve_offset" : curve_offset, "closed" : closed, **kwargs})

class associate_vertices(_Bosl2Base):
    def __init__(self, polygons=None, split=None, curpoly=None, **kwargs):
       super().__init__("associate_vertices", {"polygons" : polygons, "split" : split, "curpoly" : curpoly, **kwargs})

class _get_texture(_Bosl2Base):
    def __init__(self, tex=None, n=None, m=None, **kwargs):
       super().__init__("_get_texture", {"tex" : tex, "n" : n, "m" : m, **kwargs})

class textured_linear_sweep(_Bosl2Base):
    def __init__(self, path=None, texture=None, tex_size=None, h=None, counts=None, inset=None, rot=None, tscale=None, caps=None, col_wrap=None, twist=None, scale=None, shift=None, style=None, reverse=None, l=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("textured_linear_sweep", {"path" : path, "texture" : texture, "tex_size" : tex_size, "h" : h, "counts" : counts, "inset" : inset, "rot" : rot, "tscale" : tscale, "caps" : caps, "col_wrap" : col_wrap, "twist" : twist, "scale" : scale, "shift" : shift, "style" : style, "reverse" : reverse, "l" : l, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class textured_revolution(_Bosl2Base):
    def __init__(self, path=None, texture=None, tex_size=None, tscale=None, inset=None, rot=None, caps=None, wrap=None, shift=None, style=None, reverse=None, counts=None, **kwargs):
       super().__init__("textured_revolution", {"path" : path, "texture" : texture, "tex_size" : tex_size, "tscale" : tscale, "inset" : inset, "rot" : rot, "caps" : caps, "wrap" : wrap, "shift" : shift, "style" : style, "reverse" : reverse, "counts" : counts, **kwargs})

class textured_cylinder(_Bosl2Base):
    def __init__(self, h=None, r=None, texture=None, tex_size=None, counts=None, tscale=None, inset=None, rot=None, caps=None, style=None, reverse=None, shift=None, l=None, r1=None, r2=None, d=None, d1=None, d2=None, chamfer=None, chamfer1=None, chamfer2=None, rounding=None, rounding1=None, rounding2=None, **kwargs):
       super().__init__("textured_cylinder", {"h" : h, "r" : r, "texture" : texture, "tex_size" : tex_size, "counts" : counts, "tscale" : tscale, "inset" : inset, "rot" : rot, "caps" : caps, "style" : style, "reverse" : reverse, "shift" : shift, "l" : l, "r1" : r1, "r2" : r2, "d" : d, "d1" : d1, "d2" : d2, "chamfer" : chamfer, "chamfer1" : chamfer1, "chamfer2" : chamfer2, "rounding" : rounding, "rounding1" : rounding1, "rounding2" : rounding2, **kwargs})

class skin(_Bosl2Base):
    def __init__(self, profiles=None, slices=None, refine=None, method=None, sampling=None, caps=None, closed=None, z=None, style=None, convexity=None, anchor=None, cp=None, spin=None, orient=None, atype=None, **kwargs):
       super().__init__("skin", {"profiles" : profiles, "slices" : slices, "refine" : refine, "method" : method, "sampling" : sampling, "caps" : caps, "closed" : closed, "z" : z, "style" : style, "convexity" : convexity, "anchor" : anchor, "cp" : cp, "spin" : spin, "orient" : orient, "atype" : atype, **kwargs})

class linear_sweep(_Bosl2Base):
    def __init__(self, region=None, height=None, center=None, twist=None, scale=None, shift=None, slices=None, maxseg=None, style=None, convexity=None, cp=None, atype=None, h=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("linear_sweep", {"region" : region, "height" : height, "center" : center, "twist" : twist, "scale" : scale, "shift" : shift, "slices" : slices, "maxseg" : maxseg, "style" : style, "convexity" : convexity, "cp" : cp, "atype" : atype, "h" : h, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class spiral_sweep(_Bosl2Base):
    def __init__(self, poly=None, h=None, r=None, turns=None, higbee=None, center=None, r1=None, r2=None, d=None, d1=None, d2=None, higbee1=None, higbee2=None, internal=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("spiral_sweep", {"poly" : poly, "h" : h, "r" : r, "turns" : turns, "higbee" : higbee, "center" : center, "r1" : r1, "r2" : r2, "d" : d, "d1" : d1, "d2" : d2, "higbee1" : higbee1, "higbee2" : higbee2, "internal" : internal, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class path_sweep(_Bosl2Base):
    def __init__(self, shape=None, path=None, method=None, normal=None, closed=None, twist=None, twist_by_length=None, symmetry=None, last_normal=None, tangent=None, uniform=None, relaxed=None, caps=None, style=None, convexity=None, anchor=None, cp=None, spin=None, orient=None, atype=None, profiles=None, width=None, **kwargs):
       super().__init__("path_sweep", {"shape" : shape, "path" : path, "method" : method, "normal" : normal, "closed" : closed, "twist" : twist, "twist_by_length" : twist_by_length, "symmetry" : symmetry, "last_normal" : last_normal, "tangent" : tangent, "uniform" : uniform, "relaxed" : relaxed, "caps" : caps, "style" : style, "convexity" : convexity, "anchor" : anchor, "cp" : cp, "spin" : spin, "orient" : orient, "atype" : atype, "profiles" : profiles, "width" : width, **kwargs})

class path_sweep2d(_Bosl2Base):
    def __init__(self, profile=None, path=None, closed=None, caps=None, quality=None, style=None, convexity=None, anchor=None, cp=None, spin=None, orient=None, atype=None, **kwargs):
       super().__init__("path_sweep2d", {"profile" : profile, "path" : path, "closed" : closed, "caps" : caps, "quality" : quality, "style" : style, "convexity" : convexity, "anchor" : anchor, "cp" : cp, "spin" : spin, "orient" : orient, "atype" : atype, **kwargs})

class sweep(_Bosl2Base):
    def __init__(self, shape=None, transforms=None, closed=None, caps=None, style=None, convexity=None, anchor=None, cp=None, spin=None, orient=None, atype=None, **kwargs):
       super().__init__("sweep", {"shape" : shape, "transforms" : transforms, "closed" : closed, "caps" : caps, "style" : style, "convexity" : convexity, "anchor" : anchor, "cp" : cp, "spin" : spin, "orient" : orient, "atype" : atype, **kwargs})

class textured_linear_sweep(_Bosl2Base):
    def __init__(self, path=None, texture=None, tex_size=None, h=None, inset=None, rot=None, tscale=None, twist=None, scale=None, shift=None, style=None, reverse=None, l=None, counts=None, anchor=None, spin=None, orient=None, convexity=None, **kwargs):
       super().__init__("textured_linear_sweep", {"path" : path, "texture" : texture, "tex_size" : tex_size, "h" : h, "inset" : inset, "rot" : rot, "tscale" : tscale, "twist" : twist, "scale" : scale, "shift" : shift, "style" : style, "reverse" : reverse, "l" : l, "counts" : counts, "anchor" : anchor, "spin" : spin, "orient" : orient, "convexity" : convexity, **kwargs})

class textured_revolution(_Bosl2Base):
    def __init__(self, path=None, texture=None, tex_size=None, tscale=None, inset=None, rot=None, caps=None, wrap=None, shift=None, style=None, reverse=None, atype=None, convexity=None, counts=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("textured_revolution", {"path" : path, "texture" : texture, "tex_size" : tex_size, "tscale" : tscale, "inset" : inset, "rot" : rot, "caps" : caps, "wrap" : wrap, "shift" : shift, "style" : style, "reverse" : reverse, "atype" : atype, "convexity" : convexity, "counts" : counts, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

class textured_cylinder(_Bosl2Base):
    def __init__(self, h=None, r=None, texture=None, tex_size=None, counts=None, tscale=None, inset=None, rot=None, style=None, reverse=None, shift=None, l=None, r1=None, r2=None, d=None, d1=None, d2=None, chamfer=None, chamfer1=None, chamfer2=None, rounding=None, rounding1=None, rounding2=None, convexity=None, anchor=None, spin=None, orient=None, **kwargs):
       super().__init__("textured_cylinder", {"h" : h, "r" : r, "texture" : texture, "tex_size" : tex_size, "counts" : counts, "tscale" : tscale, "inset" : inset, "rot" : rot, "style" : style, "reverse" : reverse, "shift" : shift, "l" : l, "r1" : r1, "r2" : r2, "d" : d, "d1" : d1, "d2" : d2, "chamfer" : chamfer, "chamfer1" : chamfer1, "chamfer2" : chamfer2, "rounding" : rounding, "rounding1" : rounding1, "rounding2" : rounding2, "convexity" : convexity, "anchor" : anchor, "spin" : spin, "orient" : orient, **kwargs})

