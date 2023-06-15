include </home/jeff/code/solid_stuff/SolidPython/solid2/extensions/../libs/BOSL2/version.scad>;
include </home/jeff/code/solid_stuff/SolidPython/solid2/extensions/../libs/BOSL2/constants.scad>;
include </home/jeff/code/solid_stuff/SolidPython/solid2/extensions/../libs/BOSL2/transforms.scad>;
include </home/jeff/code/solid_stuff/SolidPython/solid2/extensions/../libs/BOSL2/distributors.scad>;
include </home/jeff/code/solid_stuff/SolidPython/solid2/extensions/../libs/BOSL2/mutators.scad>;
include </home/jeff/code/solid_stuff/SolidPython/solid2/extensions/../libs/BOSL2/color.scad>;
include </home/jeff/code/solid_stuff/SolidPython/solid2/extensions/../libs/BOSL2/attachments.scad>;
include </home/jeff/code/solid_stuff/SolidPython/solid2/extensions/../libs/BOSL2/shapes3d.scad>;
include </home/jeff/code/solid_stuff/SolidPython/solid2/extensions/../libs/BOSL2/shapes2d.scad>;
include </home/jeff/code/solid_stuff/SolidPython/solid2/extensions/../libs/BOSL2/drawing.scad>;
include </home/jeff/code/solid_stuff/SolidPython/solid2/extensions/../libs/BOSL2/masks3d.scad>;
include </home/jeff/code/solid_stuff/SolidPython/solid2/extensions/../libs/BOSL2/masks2d.scad>;
include </home/jeff/code/solid_stuff/SolidPython/solid2/extensions/../libs/BOSL2/math.scad>;
include </home/jeff/code/solid_stuff/SolidPython/solid2/extensions/../libs/BOSL2/paths.scad>;
include </home/jeff/code/solid_stuff/SolidPython/solid2/extensions/../libs/BOSL2/lists.scad>;
include </home/jeff/code/solid_stuff/SolidPython/solid2/extensions/../libs/BOSL2/comparisons.scad>;
include </home/jeff/code/solid_stuff/SolidPython/solid2/extensions/../libs/BOSL2/linalg.scad>;
include </home/jeff/code/solid_stuff/SolidPython/solid2/extensions/../libs/BOSL2/trigonometry.scad>;
include </home/jeff/code/solid_stuff/SolidPython/solid2/extensions/../libs/BOSL2/vectors.scad>;
include </home/jeff/code/solid_stuff/SolidPython/solid2/extensions/../libs/BOSL2/affine.scad>;
include </home/jeff/code/solid_stuff/SolidPython/solid2/extensions/../libs/BOSL2/coords.scad>;
include </home/jeff/code/solid_stuff/SolidPython/solid2/extensions/../libs/BOSL2/geometry.scad>;
include </home/jeff/code/solid_stuff/SolidPython/solid2/extensions/../libs/BOSL2/regions.scad>;
include </home/jeff/code/solid_stuff/SolidPython/solid2/extensions/../libs/BOSL2/strings.scad>;
include </home/jeff/code/solid_stuff/SolidPython/solid2/extensions/../libs/BOSL2/skin.scad>;
include </home/jeff/code/solid_stuff/SolidPython/solid2/extensions/../libs/BOSL2/vnf.scad>;
include </home/jeff/code/solid_stuff/SolidPython/solid2/extensions/../libs/BOSL2/utility.scad>;
include </home/jeff/code/solid_stuff/SolidPython/solid2/extensions/../libs/BOSL2/partitions.scad>;
include </home/jeff/code/solid_stuff/SolidPython/solid2/extensions/../libs/BOSL2/gears.scad>;
include </home/jeff/code/solid_stuff/SolidPython/solid2/extensions/../libs/BOSL2/beziers.scad>;
include </home/jeff/code/solid_stuff/SolidPython/solid2/extensions/../libs/BOSL2/screws.scad>;
include </home/jeff/code/solid_stuff/SolidPython/solid2/extensions/../libs/BOSL2/cubetruss.scad>;



$fa = 1;
$fs = 1;

xdistribute(spacing = 50) {
	recolor(c = "#f77") {
		diff(remove = "hole") {
			cuboid(anchor = FRONT, chamfer = 10, edges = [(RIGHT + BACK), (RIGHT + FRONT)], size = [45, 45, 10]) {
				tag(tag = "hole") {
					cuboid(chamfer = 5, edges = [(RIGHT + BACK), (RIGHT + FRONT)], size = [30, 30, 11]);
				}
				attach(from = FRONT, overlap = 5, to = BACK) {
					cuboid(edges = [(RIGHT + BACK), (RIGHT + FRONT)], rounding = 15, size = [45, 45, 10]) {
						tag(tag = "hole") {
							cuboid(edges = [(RIGHT + BACK), (RIGHT + FRONT)], rounding = 10, size = [30, 30, 11]);
						}
					}
				}
			}
		}
	}
	recolor(c = "#7f7") {
		bevel_gear(face_width = 12, pitch = 8, pitch_angle = 45, shaft_diam = 25, slices = 12, spiral_angle = 30, teeth = 20);
	}
	recolor(c = "#99f") {
		path_sweep(path = bezpath_curve(bezpath = [[-18, -20], [-18, -45], [18, -45], [18, -20], [18, 0], [-18, 0], [-18, 20], [-18, 45], [18, 45], [18, 20]]), shape = regular_ngon(d = 10, n = 3, spin = 90));
	}
	recolor(c = "#0bf") {
		move(v = [-15, -35, 0]) {
			cubetruss_corner(bracing = false, clipthick = 0, extents = [3, 8, 0, 0, 0], h = 1, size = 10, strut = 1);
		}
	}
	recolor(c = "#777") {
		xdistribute(spacing = 24) {
			screw(anchor = "origin", head = "hex", name = "M12,70", orient = BACK) {
				attach(from = BOT, to = CENTER) {
					nut(diameter = 20, name = "M12", thickness = 10);
				}
			}
			screw(anchor = "origin", head = "hex", name = "M12,70", orient = BACK) {
				attach(from = BOT, to = CENTER) {
					nut(diameter = 20, name = "M12", thickness = 10);
				}
			}
		}
	}
}
