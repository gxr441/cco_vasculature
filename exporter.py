# cco_vasculature
#exporter.py
#     Copyright (C) 2023  Anirudh P K
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.
#

import csv
import math
import binarytree
from lib import Point
import solid2
import lib
import euclid

import constants

def export_tree_to_csv(v_root):
    with open('v_output.csv', 'w') as csvfile:
        csv_writer = csv.writer(csvfile, delimiter=',', quoting = csv.QUOTE_MINIMAL)
        tpl = []
        binarytree.get_tree_preorder_list(v_root, tpl)

        for node in tpl:
            segment_p1 = None
            if node.parent is None:
                segment_p1 = Point(0.0,0.0)
            else:
                segment_p1 = node.parent.data.distal_point

            segment_p2 = node.data.distal_point
            csv_writer.writerow([segment_p2.x, segment_p2.y])
            csv_writer.writerow([segment_p1.x, segment_p1.y])
def export_to_openscad_csg(v_root):
    tpl = []
    binarytree.get_tree_preorder_list(v_root, tpl)

    solid2_segments_list = []

    for node in tpl:
        segment_p1 = None
        r1_l = None
        if node.parent is None:
            segment_p1 = Point(0.0,0.0)
            r1_l = node.data.r
        else:
            segment_p1 = node.parent.data.distal_point
            r1_l = node.parent.data.r

        segment_p2 = node.data.distal_point

        lj = lib.calculate_distance_bw_points(segment_p1, segment_p2)

        segment_local_cartesian_quadrant = lib.get_segment_quadrant_in_local_cartesian(node)
        segment_angle = math.degrees(math.atan((segment_p2.y - segment_p1.y)/(segment_p2.x - segment_p1.x)))
        segment_angle_absolute = None

        if segment_local_cartesian_quadrant == 1:
            segment_angle_absolute = segment_angle + 90
        else:
            segment_angle_absolute = segment_angle - 90




        d = solid2.cylinder(r1 = r1_l*constants.OPENSCAD_EXPORT_SCALING_FACTOR, r2=node.data.r*constants.OPENSCAD_EXPORT_SCALING_FACTOR, h = lj*constants.OPENSCAD_EXPORT_SCALING_FACTOR)
        #d = solid2.cylinder(r=node.data.r*1000, h = lj*1000)
        #id = solid2.cylinder(r1 = r1_l*(1-constants.WALL_THICKNESS_LUMEN_RATIO)*1000, r2 = node.data.r*(1-constants.WALL_THICKNESS_LUMEN_RATIO)*1000, h = lj*1000)
        #d = d-id
        d = d.rotate([90,0,segment_angle_absolute]).translate([segment_p1.x*constants.OPENSCAD_EXPORT_SCALING_FACTOR, segment_p1.y*constants.OPENSCAD_EXPORT_SCALING_FACTOR, 0])



        solid2_segments_list.append(d)

    unioned_body = solid2_segments_list[0]

    for segment in solid2_segments_list[1:]:
        unioned_body = unioned_body + segment


    bifurcation_geos = []

    for node in tpl:
        parent_segment_p1 = None
        if node.parent is None:
            parent_segment_p1 = Point(0.0,0.0)
        else:
            parent_segment_p1 = node.parent.data.distal_point

        parent_segment_p2 = node.data.distal_point

        if node.left is not None:
            left_daughter_segment_p1 = node.data.distal_point
            left_daughter_segment_p2 = node.left.data.distal_point
            lleft = lib.calculate_distance_bw_points(left_daughter_segment_p1, left_daughter_segment_p2)
            left_daughter_segment_angle = math.degrees(math.atan((left_daughter_segment_p2.y - left_daughter_segment_p1.y)/(left_daughter_segment_p2.x - left_daughter_segment_p1.x)))

            if lib.get_segment_quadrant_in_local_cartesian(node.left) == 1:
                left_daughter_segment_angle_absolute = left_daughter_segment_angle + 90
            else:
                left_daughter_segment_angle_absolute = left_daughter_segment_angle - 90
        if node.right is not None:
            right_daughter_segment_p1 = node.data.distal_point
            right_daughter_segment_p2 = node.right.data.distal_point
            lright = lib.calculate_distance_bw_points(right_daughter_segment_p1, right_daughter_segment_p2)
            right_daughter_segment_angle = math.degrees(math.atan((right_daughter_segment_p2.y - right_daughter_segment_p1.y)/(right_daughter_segment_p2.x - right_daughter_segment_p1.x)))

            if lib.get_segment_quadrant_in_local_cartesian(node.right) == 1:
                right_daughter_segment_angle_absolute = right_daughter_segment_angle + 90
            else:
                right_daughter_segment_angle_absolute = right_daughter_segment_angle - 90





        parent_segment_angle = math.degrees(math.atan((parent_segment_p2.y - parent_segment_p1.y)/(parent_segment_p2.x - parent_segment_p1.x)))

        if lib.get_segment_quadrant_in_local_cartesian(node) == 1:
            parent_segment_angle_absolute = parent_segment_angle + 90
        else:
            parent_segment_angle_absolute = parent_segment_angle - 90

        lparent = lib.calculate_distance_bw_points(parent_segment_p1, parent_segment_p2)

        parent_segment_end_cylinder_solid2 = solid2.cylinder(r=node.data.r*constants.OPENSCAD_EXPORT_SCALING_FACTOR, h=constants.EPSILON * lparent*constants.OPENSCAD_EXPORT_SCALING_FACTOR)
        parent_segment_end_cylinder_solid2 = parent_segment_end_cylinder_solid2.rotate([0, parent_segment_angle_absolute, 0])

        parent_segment_end_cylinder_mask_solid2 = solid2.cylinder(r=node.data.r*constants.OPENSCAD_EXPORT_SCALING_FACTOR*(1-constants.SEGMENT_ASPECT_RATIO), h=constants.EPSILON * lparent * constants.OPENSCAD_EXPORT_SCALING_FACTOR)
        parent_segment_end_cylinder_mask_solid2 = parent_segment_end_cylinder_mask_solid2.rotate([0, parent_segment_angle_absolute, 0])

        right_daughter_start_cylinder_solid2 = solid2.cylinder(r=node.data.r*constants.OPENSCAD_EXPORT_SCALING_FACTOR, h=constants.EPSILON * lright * constants.OPENSCAD_EXPORT_SCALING_FACTOR).rotate([0, right_daughter_segment_angle_absolute, 0])
        right_daughter_start_cylinder_mask_solid2 = solid2.cylinder(r=node.data.r*constants.OPENSCAD_EXPORT_SCALING_FACTOR*(1-constants.SEGMENT_ASPECT_RATIO), h=constants.EPSILON * lright * constants.OPENSCAD_EXPORT_SCALING_FACTOR)
        left_daughter_start_cylinder_solid2 = solid2.cylinder(r=node.data.r*constants.OPENSCAD_EXPORT_SCALING_FACTOR, h=constants.EPSILON * lleft * constants.OPENSCAD_EXPORT_SCALING_FACTOR).rotate([0, left_daughter_segment_angle_absolute, 0])
        left_daughter_start_cylinder_mask_solid2 = solid2.cylinder(r=node.data.r*constants.OPENSCAD_EXPORT_SCALING_FACTOR*(1-constants.SEGMENT_ASPECT_RATIO), h=constants.EPSILON * lleft * constants.OPENSCAD_EXPORT_SCALING_FACTOR)


        if not (node.left is None or node.right is None):
            unioned_geo = parent_segment_end_cylinder_solid2 + left_daughter_start_cylinder_solid2 + right_daughter_start_cylinder_solid2
            unioned_mask_geo = parent_segment_end_cylinder_mask_solid2 + left_daughter_start_cylinder_mask_solid2 + right_daughter_start_cylinder_mask_solid2
            hulled_geo = unioned_geo.hull().rotate([90,0,0]).translate([node.data.distal_point.x*constants.OPENSCAD_EXPORT_SCALING_FACTOR, node.data.distal_point.y*constants.OPENSCAD_EXPORT_SCALING_FACTOR, 0])
            hulled_mask_geo = unioned_mask_geo.hull().rotate([90,0,0]).translate([node.data.distal_point.x*constants.OPENSCAD_EXPORT_SCALING_FACTOR, node.data.distal_point.y*constants.OPENSCAD_EXPORT_SCALING_FACTOR, 0])
            hulled_geo = hulled_geo #- hulled_mask_geo
            bifurcation_geos.append(hulled_geo)

    for geo in bifurcation_geos:
        unioned_body += geo
    #unioned_body = bifurcation_geos[0]

    #Create mask to hollow out the Vasculature
    solid2_segment_masks_list = []
    for node in tpl:
        segment_p1 = None
        r1_l = None
        if node.parent is None:
            segment_p1 = Point(0.0,0.0)
            r1_l = node.data.r
        else:
            segment_p1 = node.parent.data.distal_point
            r1_l = node.parent.data.r

        segment_p2 = node.data.distal_point

        lj = lib.calculate_distance_bw_points(segment_p1, segment_p2)

        segment_local_cartesian_quadrant = lib.get_segment_quadrant_in_local_cartesian(node)
        segment_angle = math.degrees(math.atan((segment_p2.y - segment_p1.y)/(segment_p2.x - segment_p1.x)))
        segment_angle_absolute = None

        if segment_local_cartesian_quadrant == 1:
            segment_angle_absolute = segment_angle + 90
        else:
            segment_angle_absolute = segment_angle - 90


        line_segment_dy = (segment_p2.y-segment_p1.y)
        line_segment_dx = (segment_p2.x - segment_p1.x)

        line_vector = euclid.Vector2(line_segment_dx, line_segment_dy)
        line_vector_normalized = line_vector.normalize()


        d = solid2.cylinder(r1 = r1_l*constants.OPENSCAD_EXPORT_SCALING_FACTOR*(1-constants.WALL_THICKNESS_LUMEN_RATIO), r2=node.data.r*constants.OPENSCAD_EXPORT_SCALING_FACTOR*(1-constants.WALL_THICKNESS_LUMEN_RATIO), h = lj*1000)
        #d = solid2.cylinder(r=node.data.r*1000, h = lj*1000)
        #id = solid2.cylinder(r1 = r1_l*(1-constants.WALL_THICKNESS_LUMEN_RATIO)*1000, r2 = node.data.r*(1-constants.WALL_THICKNESS_LUMEN_RATIO)*1000, h = lj*1000)
        #d = d-id
        d = d.rotate([90,0,segment_angle_absolute]).translate([segment_p1.x*constants.OPENSCAD_EXPORT_SCALING_FACTOR, segment_p1.y*constants.OPENSCAD_EXPORT_SCALING_FACTOR, 0])

        solid2_segment_start_cap = solid2.cylinder(r = node.data.r*(1-constants.WALL_THICKNESS_LUMEN_RATIO)*constants.OPENSCAD_EXPORT_SCALING_FACTOR, h = lj * constants.EPSILON * 2  * constants.OPENSCAD_EXPORT_SCALING_FACTOR)
        solid2_segment_start_cap = solid2_segment_start_cap.rotate([90,0,segment_angle_absolute]).translate([segment_p1.x*constants.OPENSCAD_EXPORT_SCALING_FACTOR-constants.EPSILON*line_vector_normalized.x, segment_p1.y*constants.OPENSCAD_EXPORT_SCALING_FACTOR-0.01*line_vector_normalized.y, 0])
        d = d + solid2_segment_start_cap

        solid2_segment_end_cap = solid2.cylinder(r = node.data.r*(1-constants.WALL_THICKNESS_LUMEN_RATIO)*constants.OPENSCAD_EXPORT_SCALING_FACTOR, h = lj*constants.EPSILON*2*constants.OPENSCAD_EXPORT_SCALING_FACTOR)
        solid2_segment_end_cap = solid2_segment_end_cap.rotate([90, 0, segment_angle_absolute]).translate([segment_p2.x*constants.OPENSCAD_EXPORT_SCALING_FACTOR - constants.EPSILON*line_vector_normalized.x, segment_p2.y*constants.OPENSCAD_EXPORT_SCALING_FACTOR - 0.01*line_vector_normalized.y, 0])
        d = d + solid2_segment_end_cap

        #Correct OpenSCAD Z-fighting artefact


        solid2_segment_masks_list.append(d)


    unioned_mask = solid2_segment_masks_list[0]
    for segment_mask in solid2_segment_masks_list:
        unioned_mask += segment_mask

    unioned_body = unioned_body# - unioned_mask

    #Correct OpenSCAD Z-fighting artefact using a thing cylinder at end of each node





    unioned_body.save_as_scad("test.scad")
