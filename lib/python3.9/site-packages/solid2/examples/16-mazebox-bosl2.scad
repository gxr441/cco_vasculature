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





difference() {
	union() {
		tube(or = 30.3577490736, center = true, h = 88, ir = 28.3577490736);
		down(z = 40) {
			zcyl(h = 10, r = 35.5577490736);
		}
	}
	union() {
		cylindrical_extrude(or = 30.5577490736, ir = 30.3577490736) {
			offset(r = 0.0000000000) {
				projection(cut = true) {
					translate(v = [0, 0, 10]) {
						surface(center = true, file = "maze7.png", invert = true);
					}
				}
			}
		}
		cylindrical_extrude(or = 30.3577490736, ir = 30.1577490736) {
			offset(r = -0.0123116594) {
				projection(cut = true) {
					translate(v = [0, 0, 10]) {
						surface(center = true, file = "maze7.png", invert = true);
					}
				}
			}
		}
		cylindrical_extrude(or = 30.1577490736, ir = 29.9577490736) {
			offset(r = -0.0489434837) {
				projection(cut = true) {
					translate(v = [0, 0, 10]) {
						surface(center = true, file = "maze7.png", invert = true);
					}
				}
			}
		}
		cylindrical_extrude(or = 29.9577490736, ir = 29.7577490736) {
			offset(r = -0.1089934758) {
				projection(cut = true) {
					translate(v = [0, 0, 10]) {
						surface(center = true, file = "maze7.png", invert = true);
					}
				}
			}
		}
		cylindrical_extrude(or = 29.7577490736, ir = 29.5577490736) {
			offset(r = -0.1909830056) {
				projection(cut = true) {
					translate(v = [0, 0, 10]) {
						surface(center = true, file = "maze7.png", invert = true);
					}
				}
			}
		}
		cylindrical_extrude(or = 29.5577490736, ir = 29.3577490736) {
			offset(r = -0.2928932188) {
				projection(cut = true) {
					translate(v = [0, 0, 10]) {
						surface(center = true, file = "maze7.png", invert = true);
					}
				}
			}
		}
		cylindrical_extrude(or = 29.3577490736, ir = 29.1577490736) {
			offset(r = -0.4122147477) {
				projection(cut = true) {
					translate(v = [0, 0, 10]) {
						surface(center = true, file = "maze7.png", invert = true);
					}
				}
			}
		}
		cylindrical_extrude(or = 29.1577490736, ir = 28.9577490736) {
			offset(r = -0.5460095003) {
				projection(cut = true) {
					translate(v = [0, 0, 10]) {
						surface(center = true, file = "maze7.png", invert = true);
					}
				}
			}
		}
		cylindrical_extrude(or = 28.9577490736, ir = 28.7577490736) {
			offset(r = -0.6909830056) {
				projection(cut = true) {
					translate(v = [0, 0, 10]) {
						surface(center = true, file = "maze7.png", invert = true);
					}
				}
			}
		}
		cylindrical_extrude(or = 28.7577490736, ir = 28.5577490736) {
			offset(r = -0.8435655350) {
				projection(cut = true) {
					translate(v = [0, 0, 10]) {
						surface(center = true, file = "maze7.png", invert = true);
					}
				}
			}
		}
		zrot(a = 120.0000000000) {
			union() {
				cylindrical_extrude(or = 30.5577490736, ir = 30.3577490736) {
					offset(r = 0.0000000000) {
						projection(cut = true) {
							translate(v = [0, 0, 10]) {
								surface(center = true, file = "maze7.png", invert = true);
							}
						}
					}
				}
				cylindrical_extrude(or = 30.3577490736, ir = 30.1577490736) {
					offset(r = -0.0123116594) {
						projection(cut = true) {
							translate(v = [0, 0, 10]) {
								surface(center = true, file = "maze7.png", invert = true);
							}
						}
					}
				}
				cylindrical_extrude(or = 30.1577490736, ir = 29.9577490736) {
					offset(r = -0.0489434837) {
						projection(cut = true) {
							translate(v = [0, 0, 10]) {
								surface(center = true, file = "maze7.png", invert = true);
							}
						}
					}
				}
				cylindrical_extrude(or = 29.9577490736, ir = 29.7577490736) {
					offset(r = -0.1089934758) {
						projection(cut = true) {
							translate(v = [0, 0, 10]) {
								surface(center = true, file = "maze7.png", invert = true);
							}
						}
					}
				}
				cylindrical_extrude(or = 29.7577490736, ir = 29.5577490736) {
					offset(r = -0.1909830056) {
						projection(cut = true) {
							translate(v = [0, 0, 10]) {
								surface(center = true, file = "maze7.png", invert = true);
							}
						}
					}
				}
				cylindrical_extrude(or = 29.5577490736, ir = 29.3577490736) {
					offset(r = -0.2928932188) {
						projection(cut = true) {
							translate(v = [0, 0, 10]) {
								surface(center = true, file = "maze7.png", invert = true);
							}
						}
					}
				}
				cylindrical_extrude(or = 29.3577490736, ir = 29.1577490736) {
					offset(r = -0.4122147477) {
						projection(cut = true) {
							translate(v = [0, 0, 10]) {
								surface(center = true, file = "maze7.png", invert = true);
							}
						}
					}
				}
				cylindrical_extrude(or = 29.1577490736, ir = 28.9577490736) {
					offset(r = -0.5460095003) {
						projection(cut = true) {
							translate(v = [0, 0, 10]) {
								surface(center = true, file = "maze7.png", invert = true);
							}
						}
					}
				}
				cylindrical_extrude(or = 28.9577490736, ir = 28.7577490736) {
					offset(r = -0.6909830056) {
						projection(cut = true) {
							translate(v = [0, 0, 10]) {
								surface(center = true, file = "maze7.png", invert = true);
							}
						}
					}
				}
				cylindrical_extrude(or = 28.7577490736, ir = 28.5577490736) {
					offset(r = -0.8435655350) {
						projection(cut = true) {
							translate(v = [0, 0, 10]) {
								surface(center = true, file = "maze7.png", invert = true);
							}
						}
					}
				}
			}
		}
		zrot(a = 240.0000000000) {
			union() {
				cylindrical_extrude(or = 30.5577490736, ir = 30.3577490736) {
					offset(r = 0.0000000000) {
						projection(cut = true) {
							translate(v = [0, 0, 10]) {
								surface(center = true, file = "maze7.png", invert = true);
							}
						}
					}
				}
				cylindrical_extrude(or = 30.3577490736, ir = 30.1577490736) {
					offset(r = -0.0123116594) {
						projection(cut = true) {
							translate(v = [0, 0, 10]) {
								surface(center = true, file = "maze7.png", invert = true);
							}
						}
					}
				}
				cylindrical_extrude(or = 30.1577490736, ir = 29.9577490736) {
					offset(r = -0.0489434837) {
						projection(cut = true) {
							translate(v = [0, 0, 10]) {
								surface(center = true, file = "maze7.png", invert = true);
							}
						}
					}
				}
				cylindrical_extrude(or = 29.9577490736, ir = 29.7577490736) {
					offset(r = -0.1089934758) {
						projection(cut = true) {
							translate(v = [0, 0, 10]) {
								surface(center = true, file = "maze7.png", invert = true);
							}
						}
					}
				}
				cylindrical_extrude(or = 29.7577490736, ir = 29.5577490736) {
					offset(r = -0.1909830056) {
						projection(cut = true) {
							translate(v = [0, 0, 10]) {
								surface(center = true, file = "maze7.png", invert = true);
							}
						}
					}
				}
				cylindrical_extrude(or = 29.5577490736, ir = 29.3577490736) {
					offset(r = -0.2928932188) {
						projection(cut = true) {
							translate(v = [0, 0, 10]) {
								surface(center = true, file = "maze7.png", invert = true);
							}
						}
					}
				}
				cylindrical_extrude(or = 29.3577490736, ir = 29.1577490736) {
					offset(r = -0.4122147477) {
						projection(cut = true) {
							translate(v = [0, 0, 10]) {
								surface(center = true, file = "maze7.png", invert = true);
							}
						}
					}
				}
				cylindrical_extrude(or = 29.1577490736, ir = 28.9577490736) {
					offset(r = -0.5460095003) {
						projection(cut = true) {
							translate(v = [0, 0, 10]) {
								surface(center = true, file = "maze7.png", invert = true);
							}
						}
					}
				}
				cylindrical_extrude(or = 28.9577490736, ir = 28.7577490736) {
					offset(r = -0.6909830056) {
						projection(cut = true) {
							translate(v = [0, 0, 10]) {
								surface(center = true, file = "maze7.png", invert = true);
							}
						}
					}
				}
				cylindrical_extrude(or = 28.7577490736, ir = 28.5577490736) {
					offset(r = -0.8435655350) {
						projection(cut = true) {
							translate(v = [0, 0, 10]) {
								surface(center = true, file = "maze7.png", invert = true);
							}
						}
					}
				}
			}
		}
	}
}
