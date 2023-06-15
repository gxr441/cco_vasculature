# cco_vasculature
#cli.py
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

from lib import Node, Segment, Circle, Segment, Point
import lib
import math
from copy import deepcopy
import binarytree
import constants
import exporter
import visualizer
import time

import random



def get_new_terminal_point(circle_rsupp, v_root):
    spacing = lib.compute_dthresh(r_supp, k_term)
    prospective_distal_point = circle_rsupp.get_random_coordinate(spacing) #used len-1 because node_list should never be 0

    tree_inorder_list = []
    binarytree.get_tree_inorder_list(v_root, tree_inorder_list)


    done = False
    lap = 0
    while not done:
        for k, node in enumerate(tree_inorder_list):
            segment_data = node.data
            if node.parent is None:
                segment_start = Point(0.0,0.0)
            else:
                segment_start = node.parent.data.distal_point
            segment_stop = segment_data.distal_point


            dproj = lib.compute_dproj(prospective_distal_point, segment_start, segment_stop)
            dcrit = lib.compute_dcrit(prospective_distal_point, segment_start, segment_stop, dproj)
            dthresh = lib.compute_dthresh(r_supp, k_term, lap)

            if dcrit is not None:
                if dcrit > dthresh:
                    #print("New terminal point found: ", prospective_distal_point)
                    done = True
                    return {'new_point': prospective_distal_point, 'parent': k}
        lap = lap + 1


visualizer.init_pygame()

print("Initializing root...")

current_segment_index = 0

k_term = 1
k_tot = 2*k_term - 1
N_TOT = 2*constants.N_TERM - 1

r_supp = ((constants.R_PERF**2) / constants.N_TERM) ** (1/2)
spacing = lib.compute_dthresh(r_supp, k_term) / 1.1
perfusion_circle = Circle(r_supp, Point(0.0, -r_supp))
plant_root_coordinate = perfusion_circle.get_random_coordinate(spacing)
Q_root = constants.Q_PERF / constants.N_TERM
l_root = lib.calculate_distance_bw_points(Point(0.0,0.0), plant_root_coordinate)
R_red_root = (8*constants.ETA*l_root) / math.pi


R_root = (constants.P_PERF - constants.P_TERM) / Q_root
r_root = (((8*constants.ETA)/math.pi) * (l_root/R_root))**(1/4)

root_segment = Segment(index = current_segment_index, distal_point = plant_root_coordinate, b_l = None, b_r = None, Q = Q_root, R_red = R_red_root)
root_segment.r = r_root



v_root = Node(root_segment)
ibif_node_c = None
print("Tree initialized")

wl_n_term = 0
#TODO: check for intersections!!
while wl_n_term < (constants.N_TERM - 1):
    #visualizer.draw_tree(v_root, r_supp)
    local_min_table = {}

    o_r_supp = r_supp

    k_tot = 2*binarytree.getLeafCount(v_root) - 1
    r_supp = ((k_tot + 1) * ((constants.R_PERF**2)/N_TOT)) ** (1/2)

    new_virtual_circle = Circle(r_supp, Point(0.0, -r_supp))

    scale_factor_inflation = r_supp / o_r_supp
    if ibif_node_c is not None:
        tpl = []
        binarytree.get_tree_preorder_list(v_root, tpl)

        for node in tpl:
            node.data.distal_point.x = node.data.distal_point.x * scale_factor_inflation
            node.data.distal_point.y = node.data.distal_point.y * scale_factor_inflation
            #node.data.r = node.data.r * scale_factor_inflation
        #v_root = lib.recompute_tree_characteristics_after_bifurcation(v_root, ibif_node_c)
    wl_done = False


    while len(local_min_table) <= 0:
        non_intersecting_indices_list = []

        while len(non_intersecting_indices_list) <= 0:
            new_terminal_point_with_parent = get_new_terminal_point(new_virtual_circle, v_root)

            non_intersecting_indices_list = lib.get_non_intersecting_segment_indices(v_root, new_terminal_point_with_parent['new_point'])

        for non_intersecting_index in non_intersecting_indices_list:
        ##########Segment addition round
            v_root_virt = deepcopy(v_root)
            visualizer.draw_tree(v_root_virt, r_supp)
            #time.sleep(2)
            new_tree_and_bifurcating_segment = lib.add_segment_to_node(v_root_virt, new_terminal_point_with_parent, current_segment_index, non_intersecting_index)

            v_root_virt = new_tree_and_bifurcating_segment['v_root']
            ibif_node = new_tree_and_bifurcating_segment['ibif']

            ibif_node_c = ibif_node

            #visualizer.draw_tree(v_root_virt)

            v_root_virt = lib.recompute_tree_characteristics_after_bifurcation(v_root_virt, ibif_node)

            visualizer.draw_tree(v_root_virt, r_supp)

            #Optimize ibif distal_point
            c_cost_normalized_grad_vector = lib.get_cost_function_normalized_gradient_vector(v_root_virt, ibif_node, new_virtual_circle.radius)
            dx = lib.compute_dthresh(new_virtual_circle.radius, binarytree.getLeafCount(v_root_virt))
            dPx = dx * c_cost_normalized_grad_vector.x
            dPy = dx * c_cost_normalized_grad_vector.y
            v_root_virt_geo_opt = deepcopy(v_root_virt)


            wl_s_c_cost = 0.0
            wl_s_n_c_cost = 0.0
            wl_count = 0
            v_final_geo = None
            #print("init cost", lib.compute_cost_function(v_root_virt_geo_opt))
            while not wl_s_n_c_cost > wl_s_c_cost:
                v_final_geo = deepcopy(v_root_virt_geo_opt)
                visualizer.draw_tree(v_final_geo, r_supp)
                wl_s_c_cost = lib.compute_cost_function(v_root_virt_geo_opt)
                wl_ibif_node = binarytree.get_node_with_index(v_root_virt_geo_opt, ibif_node.data.index)

                wl_ibif_node.data.distal_point.x = wl_ibif_node.data.distal_point.x - dPx
                wl_ibif_node.data.distal_point.y = wl_ibif_node.data.distal_point.y - dPy

                #if not lib.is_bifurcation_sane(v_root_virt_geo_opt, wl_ibif_node, current_segment_index):
                #    break

                v_root_virt_geo_opt = lib.recompute_tree_characteristics_after_bifurcation(v_root_virt_geo_opt, wl_ibif_node)

                #Optimization constraintds
                if not new_virtual_circle._is_point_in_circle(wl_ibif_node.data.distal_point):
                    print('degeneracy1')
                    break
                l_ibif = None
                if wl_ibif_node.parent is None:
                    l_ibif = lib.calculate_distance_bw_points(Point(0.0,0.0), wl_ibif_node.data.distal_point)
                else:
                    l_ibif = lib.calculate_distance_bw_points(wl_ibif_node.parent.data.distal_point, wl_ibif_node.data.distal_point)

                if l_ibif <= constants.SEGMENT_ASPECT_RATIO*wl_ibif_node.data.r:
                    print('degeneracy2')
                    break

                l_ibif_l = lib.calculate_distance_bw_points(wl_ibif_node.data.distal_point, wl_ibif_node.left.data.distal_point)
                l_ibif_r = lib.calculate_distance_bw_points(wl_ibif_node.data.distal_point, wl_ibif_node.right.data.distal_point)

                if l_ibif_l <= constants.SEGMENT_ASPECT_RATIO*wl_ibif_node.left.data.r or l_ibif_r <= constants.SEGMENT_ASPECT_RATIO*wl_ibif_node.right.data.r:
                    print('degeneracy3')
                    break





                wl_s_n_c_cost = lib.compute_cost_function(v_root_virt_geo_opt)
                wl_count = wl_count+1
                #print('cost', wl_s_n_c_cost)

                #print(lib.is_segments_intersecting(wl_ibif_node, v_root_virt_geo_opt))


            tpl = []
            binarytree.get_tree_preorder_list(v_final_geo, tpl)
            print(tpl)
            tpl = [node for node in tpl if not (node.data.index==wl_ibif_node.data.index or node.data.index==wl_ibif_node.left.data.index or node.data.index == wl_ibif_node.right.data.index)]
            print(tpl)
            #time.sleep(2)
            is_intersecting = False
            for node in tpl:
                #print(lib.is_segments_intersecting(wl_ibif_node, node))
                #time.sleep(10)
                if lib.is_segments_intersecting(wl_ibif_node.left, node) or lib.is_segments_intersecting(wl_ibif_node.right, node):
                    is_intersecting = True
                    break
            if (not is_intersecting is False and wl_ibif_node.parent is not None):
                print(tpl)
                tpl = [node for node in tpl if not node.data.index == wl_ibif_node.parent.data.index]
                for node in tpl:
                    if lib.is_segments_intersecting(wl_ibif_node, node):
                        is_intersecting = True
                        break

            if is_intersecting is False:
                print('Segment accept')
                #time.sleep(2)
                local_min_table[deepcopy(v_final_geo)] = wl_s_c_cost
            #print(local_min_table)

    local_min_tree_list = []
    local_min_cost_list = []

    for local_min in local_min_table.items():
        local_min_tree_list.append(local_min[0])
        local_min_cost_list.append(local_min[1])

    global_min_cost_index = local_min_cost_list.index(min(local_min_cost_list))
    print(global_min_cost_index)
    #global_min_cost_index = local_min_cost_list[global_min_cost_index:].index(min(local_min_cost_list[global_min_cost_index:]))
    #global_min_tree = local_min_tree_list[random.randint(0, len(local_min_table) - 1)] #TODO: Do it right!
    global_min_tree = local_min_tree_list[global_min_cost_index]

    v_root = deepcopy(global_min_tree)

    current_segment_index = current_segment_index + 2
    wl_n_term = wl_n_term + 1

tpl = []
binarytree.get_tree_preorder_list(v_root, tpl)
print([((x.data.distal_point.x, x.data.distal_point.y), (x.parent.data.distal_point.x, x.parent.data.distal_point.y), x.data.r, x.data.index, x.parent.data.index) for x in tpl if x.parent is not None])
exporter.export_tree_to_csv(v_root)
exporter.export_to_openscad_csg(v_root)

while True:
    time.sleep(0.04)
    visualizer.draw_tree(v_root, r_supp)
