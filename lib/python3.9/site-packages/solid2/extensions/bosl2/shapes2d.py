from solid2.core.object_base import OpenSCADConstant as _OpenSCADConstant
from solid2.core.scad_import import extra_scad_include as _extra_scad_include
from pathlib import Path as _Path

from .bosl2_base import Bosl2Base as _Bosl2Base

_extra_scad_include(f"{_Path(__file__).parent.parent / '../libs/BOSL2/shapes2d.scad'}", use_not_include=False)

class square(_Bosl2Base):
    def __init__(self, size=None, center=None, anchor=None, spin=None, **kwargs):
       super().__init__("square", {"size" : size, "center" : center, "anchor" : anchor, "spin" : spin, **kwargs})

class rect(_Bosl2Base):
    def __init__(self, size=None, rounding=None, chamfer=None, atype=None, anchor=None, spin=None, **kwargs):
       super().__init__("rect", {"size" : size, "rounding" : rounding, "chamfer" : chamfer, "atype" : atype, "anchor" : anchor, "spin" : spin, **kwargs})

class circle(_Bosl2Base):
    def __init__(self, r=None, d=None, points=None, corner=None, anchor=None, spin=None, **kwargs):
       super().__init__("circle", {"r" : r, "d" : d, "points" : points, "corner" : corner, "anchor" : anchor, "spin" : spin, **kwargs})

class _ellipse_refine(_Bosl2Base):
    def __init__(self, a=None, b=None, N=None, _theta=None, **kwargs):
       super().__init__("_ellipse_refine", {"a" : a, "b" : b, "N" : N, "_theta" : _theta, **kwargs})

class _ellipse_refine_realign(_Bosl2Base):
    def __init__(self, a=None, b=None, N=None, _theta=None, i=None, **kwargs):
       super().__init__("_ellipse_refine_realign", {"a" : a, "b" : b, "N" : N, "_theta" : _theta, "i" : i, **kwargs})

class ellipse(_Bosl2Base):
    def __init__(self, r=None, d=None, realign=None, circum=None, uniform=None, anchor=None, spin=None, **kwargs):
       super().__init__("ellipse", {"r" : r, "d" : d, "realign" : realign, "circum" : circum, "uniform" : uniform, "anchor" : anchor, "spin" : spin, **kwargs})

class regular_ngon(_Bosl2Base):
    def __init__(self, n=None, r=None, d=None, _or=None, od=None, ir=None, id=None, side=None, rounding=None, realign=None, align_tip=None, align_side=None, anchor=None, spin=None, _mat=None, _anchs=None, **kwargs):
       super().__init__("regular_ngon", {"n" : n, "r" : r, "d" : d, "_or" : _or, "od" : od, "ir" : ir, "id" : id, "side" : side, "rounding" : rounding, "realign" : realign, "align_tip" : align_tip, "align_side" : align_side, "anchor" : anchor, "spin" : spin, "_mat" : _mat, "_anchs" : _anchs, **kwargs})

class pentagon(_Bosl2Base):
    def __init__(self, r=None, d=None, _or=None, od=None, ir=None, id=None, side=None, rounding=None, realign=None, align_tip=None, align_side=None, anchor=None, spin=None, **kwargs):
       super().__init__("pentagon", {"r" : r, "d" : d, "_or" : _or, "od" : od, "ir" : ir, "id" : id, "side" : side, "rounding" : rounding, "realign" : realign, "align_tip" : align_tip, "align_side" : align_side, "anchor" : anchor, "spin" : spin, **kwargs})

class hexagon(_Bosl2Base):
    def __init__(self, r=None, d=None, _or=None, od=None, ir=None, id=None, side=None, rounding=None, realign=None, align_tip=None, align_side=None, anchor=None, spin=None, **kwargs):
       super().__init__("hexagon", {"r" : r, "d" : d, "_or" : _or, "od" : od, "ir" : ir, "id" : id, "side" : side, "rounding" : rounding, "realign" : realign, "align_tip" : align_tip, "align_side" : align_side, "anchor" : anchor, "spin" : spin, **kwargs})

class octagon(_Bosl2Base):
    def __init__(self, r=None, d=None, _or=None, od=None, ir=None, id=None, side=None, rounding=None, realign=None, align_tip=None, align_side=None, anchor=None, spin=None, **kwargs):
       super().__init__("octagon", {"r" : r, "d" : d, "_or" : _or, "od" : od, "ir" : ir, "id" : id, "side" : side, "rounding" : rounding, "realign" : realign, "align_tip" : align_tip, "align_side" : align_side, "anchor" : anchor, "spin" : spin, **kwargs})

class right_triangle(_Bosl2Base):
    def __init__(self, size=None, center=None, anchor=None, spin=None, **kwargs):
       super().__init__("right_triangle", {"size" : size, "center" : center, "anchor" : anchor, "spin" : spin, **kwargs})

class trapezoid(_Bosl2Base):
    def __init__(self, h=None, w1=None, w2=None, angle=None, shift=None, chamfer=None, rounding=None, flip=None, anchor=None, spin=None, **kwargs):
       super().__init__("trapezoid", {"h" : h, "w1" : w1, "w2" : w2, "angle" : angle, "shift" : shift, "chamfer" : chamfer, "rounding" : rounding, "flip" : flip, "anchor" : anchor, "spin" : spin, **kwargs})

class star(_Bosl2Base):
    def __init__(self, n=None, r=None, ir=None, d=None, _or=None, od=None, id=None, step=None, realign=None, align_tip=None, align_pit=None, anchor=None, spin=None, atype=None, _mat=None, _anchs=None, **kwargs):
       super().__init__("star", {"n" : n, "r" : r, "ir" : ir, "d" : d, "_or" : _or, "od" : od, "id" : id, "step" : step, "realign" : realign, "align_tip" : align_tip, "align_pit" : align_pit, "anchor" : anchor, "spin" : spin, "atype" : atype, "_mat" : _mat, "_anchs" : _anchs, **kwargs})

class _path_add_jitter(_Bosl2Base):
    def __init__(self, path=None, dist=None, closed=None, **kwargs):
       super().__init__("_path_add_jitter", {"path" : path, "dist" : dist, "closed" : closed, **kwargs})

class teardrop2d(_Bosl2Base):
    def __init__(self, r=None, ang=None, cap_h=None, d=None, anchor=None, spin=None, **kwargs):
       super().__init__("teardrop2d", {"r" : r, "ang" : ang, "cap_h" : cap_h, "d" : d, "anchor" : anchor, "spin" : spin, **kwargs})

class egg(_Bosl2Base):
    def __init__(self, length=None, r1=None, r2=None, R=None, d1=None, d2=None, D=None, anchor=None, spin=None, **kwargs):
       super().__init__("egg", {"length" : length, "r1" : r1, "r2" : r2, "R" : R, "d1" : d1, "d2" : d2, "D" : D, "anchor" : anchor, "spin" : spin, **kwargs})

class glued_circles(_Bosl2Base):
    def __init__(self, r=None, spread=None, tangent=None, d=None, anchor=None, spin=None, **kwargs):
       super().__init__("glued_circles", {"r" : r, "spread" : spread, "tangent" : tangent, "d" : d, "anchor" : anchor, "spin" : spin, **kwargs})

class _superformula(_Bosl2Base):
    def __init__(self, theta=None, m1=None, m2=None, n1=None, n2=None, n3=None, a=None, b=None, **kwargs):
       super().__init__("_superformula", {"theta" : theta, "m1" : m1, "m2" : m2, "n1" : n1, "n2" : n2, "n3" : n3, "a" : a, "b" : b, **kwargs})

class supershape(_Bosl2Base):
    def __init__(self, step=None, m1=None, m2=None, n1=None, n2=None, n3=None, a=None, b=None, r=None, d=None, anchor=None, spin=None, atype=None, **kwargs):
       super().__init__("supershape", {"step" : step, "m1" : m1, "m2" : m2, "n1" : n1, "n2" : n2, "n3" : n3, "a" : a, "b" : b, "r" : r, "d" : d, "anchor" : anchor, "spin" : spin, "atype" : atype, **kwargs})

class reuleaux_polygon(_Bosl2Base):
    def __init__(self, n=None, r=None, d=None, anchor=None, spin=None, **kwargs):
       super().__init__("reuleaux_polygon", {"n" : n, "r" : r, "d" : d, "anchor" : anchor, "spin" : spin, **kwargs})

class square(_Bosl2Base):
    def __init__(self, size=None, center=None, anchor=None, spin=None, **kwargs):
       super().__init__("square", {"size" : size, "center" : center, "anchor" : anchor, "spin" : spin, **kwargs})

class rect(_Bosl2Base):
    def __init__(self, size=None, rounding=None, atype=None, chamfer=None, anchor=None, spin=None, **kwargs):
       super().__init__("rect", {"size" : size, "rounding" : rounding, "atype" : atype, "chamfer" : chamfer, "anchor" : anchor, "spin" : spin, **kwargs})

class circle(_Bosl2Base):
    def __init__(self, r=None, d=None, points=None, corner=None, anchor=None, spin=None, **kwargs):
       super().__init__("circle", {"r" : r, "d" : d, "points" : points, "corner" : corner, "anchor" : anchor, "spin" : spin, **kwargs})

class ellipse(_Bosl2Base):
    def __init__(self, r=None, d=None, realign=None, circum=None, uniform=None, anchor=None, spin=None, **kwargs):
       super().__init__("ellipse", {"r" : r, "d" : d, "realign" : realign, "circum" : circum, "uniform" : uniform, "anchor" : anchor, "spin" : spin, **kwargs})

class regular_ngon(_Bosl2Base):
    def __init__(self, n=None, r=None, d=None, _or=None, od=None, ir=None, id=None, side=None, rounding=None, realign=None, align_tip=None, align_side=None, anchor=None, spin=None, **kwargs):
       super().__init__("regular_ngon", {"n" : n, "r" : r, "d" : d, "_or" : _or, "od" : od, "ir" : ir, "id" : id, "side" : side, "rounding" : rounding, "realign" : realign, "align_tip" : align_tip, "align_side" : align_side, "anchor" : anchor, "spin" : spin, **kwargs})

class pentagon(_Bosl2Base):
    def __init__(self, r=None, d=None, _or=None, od=None, ir=None, id=None, side=None, rounding=None, realign=None, align_tip=None, align_side=None, anchor=None, spin=None, **kwargs):
       super().__init__("pentagon", {"r" : r, "d" : d, "_or" : _or, "od" : od, "ir" : ir, "id" : id, "side" : side, "rounding" : rounding, "realign" : realign, "align_tip" : align_tip, "align_side" : align_side, "anchor" : anchor, "spin" : spin, **kwargs})

class hexagon(_Bosl2Base):
    def __init__(self, r=None, d=None, _or=None, od=None, ir=None, id=None, side=None, rounding=None, realign=None, align_tip=None, align_side=None, anchor=None, spin=None, **kwargs):
       super().__init__("hexagon", {"r" : r, "d" : d, "_or" : _or, "od" : od, "ir" : ir, "id" : id, "side" : side, "rounding" : rounding, "realign" : realign, "align_tip" : align_tip, "align_side" : align_side, "anchor" : anchor, "spin" : spin, **kwargs})

class octagon(_Bosl2Base):
    def __init__(self, r=None, d=None, _or=None, od=None, ir=None, id=None, side=None, rounding=None, realign=None, align_tip=None, align_side=None, anchor=None, spin=None, **kwargs):
       super().__init__("octagon", {"r" : r, "d" : d, "_or" : _or, "od" : od, "ir" : ir, "id" : id, "side" : side, "rounding" : rounding, "realign" : realign, "align_tip" : align_tip, "align_side" : align_side, "anchor" : anchor, "spin" : spin, **kwargs})

class right_triangle(_Bosl2Base):
    def __init__(self, size=None, center=None, anchor=None, spin=None, **kwargs):
       super().__init__("right_triangle", {"size" : size, "center" : center, "anchor" : anchor, "spin" : spin, **kwargs})

class trapezoid(_Bosl2Base):
    def __init__(self, h=None, w1=None, w2=None, angle=None, shift=None, chamfer=None, rounding=None, flip=None, anchor=None, spin=None, **kwargs):
       super().__init__("trapezoid", {"h" : h, "w1" : w1, "w2" : w2, "angle" : angle, "shift" : shift, "chamfer" : chamfer, "rounding" : rounding, "flip" : flip, "anchor" : anchor, "spin" : spin, **kwargs})

class star(_Bosl2Base):
    def __init__(self, n=None, r=None, ir=None, d=None, _or=None, od=None, id=None, step=None, realign=None, align_tip=None, align_pit=None, anchor=None, spin=None, atype=None, **kwargs):
       super().__init__("star", {"n" : n, "r" : r, "ir" : ir, "d" : d, "_or" : _or, "od" : od, "id" : id, "step" : step, "realign" : realign, "align_tip" : align_tip, "align_pit" : align_pit, "anchor" : anchor, "spin" : spin, "atype" : atype, **kwargs})

class jittered_poly(_Bosl2Base):
    def __init__(self, path=None, dist=None, **kwargs):
       super().__init__("jittered_poly", {"path" : path, "dist" : dist, **kwargs})

class teardrop2d(_Bosl2Base):
    def __init__(self, r=None, ang=None, cap_h=None, d=None, anchor=None, spin=None, **kwargs):
       super().__init__("teardrop2d", {"r" : r, "ang" : ang, "cap_h" : cap_h, "d" : d, "anchor" : anchor, "spin" : spin, **kwargs})

class egg(_Bosl2Base):
    def __init__(self, length=None, r1=None, r2=None, R=None, d1=None, d2=None, D=None, anchor=None, spin=None, **kwargs):
       super().__init__("egg", {"length" : length, "r1" : r1, "r2" : r2, "R" : R, "d1" : d1, "d2" : d2, "D" : D, "anchor" : anchor, "spin" : spin, **kwargs})

class glued_circles(_Bosl2Base):
    def __init__(self, r=None, spread=None, tangent=None, d=None, anchor=None, spin=None, **kwargs):
       super().__init__("glued_circles", {"r" : r, "spread" : spread, "tangent" : tangent, "d" : d, "anchor" : anchor, "spin" : spin, **kwargs})

class supershape(_Bosl2Base):
    def __init__(self, step=None, m1=None, m2=None, n1=None, n2=None, n3=None, a=None, b=None, r=None, d=None, anchor=None, spin=None, atype=None, **kwargs):
       super().__init__("supershape", {"step" : step, "m1" : m1, "m2" : m2, "n1" : n1, "n2" : n2, "n3" : n3, "a" : a, "b" : b, "r" : r, "d" : d, "anchor" : anchor, "spin" : spin, "atype" : atype, **kwargs})

class reuleaux_polygon(_Bosl2Base):
    def __init__(self, n=None, r=None, d=None, anchor=None, spin=None, **kwargs):
       super().__init__("reuleaux_polygon", {"n" : n, "r" : r, "d" : d, "anchor" : anchor, "spin" : spin, **kwargs})

class text(_Bosl2Base):
    def __init__(self, text=None, size=None, font=None, halign=None, valign=None, spacing=None, direction=None, language=None, script=None, anchor=None, spin=None, **kwargs):
       super().__init__("text", {"text" : text, "size" : size, "font" : font, "halign" : halign, "valign" : valign, "spacing" : spacing, "direction" : direction, "language" : language, "script" : script, "anchor" : anchor, "spin" : spin, **kwargs})

class round2d(_Bosl2Base):
    def __init__(self, r=None, _or=None, ir=None, **kwargs):
       super().__init__("round2d", {"r" : r, "_or" : _or, "ir" : ir, **kwargs})

class shell2d(_Bosl2Base):
    def __init__(self, thickness=None, _or=None, ir=None, **kwargs):
       super().__init__("shell2d", {"thickness" : thickness, "_or" : _or, "ir" : ir, **kwargs})

