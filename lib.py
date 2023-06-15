# cco_vasculature
#lib.py
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


from copy import copy, deepcopy
import random
import math
import numpy
import binarytree
import constants
import euclid
import shapely

import time

from multiprocessing import Process, Queue
import itertools


def get_segment_quadrant_in_local_cartesian(node):
    segment_p1 = None
    if node.parent is None:
        segment_p1 = Point(0.0,0.0)
    else:
        segment_p1 = node.parent.data.distal_point

    segment_p2 = node.data.distal_point
    print(segment_p2.x, segment_p1.x, segment_p2.y, segment_p1.y)
    if (segment_p2.x >= segment_p1.x):
        return 1
    else:
        return 2

def is_segments_intersecting(node1, node2):
    segment1_corner_points = get_bounding_rectangle_for_segment(node1)
    segment2_corner_points = get_bounding_rectangle_for_segment(node2)

    rect1_coords_tuple = None
    rect2_coords_tuple = None
    coords_list = []
    for point in segment1_corner_points:
        coords_list.append((point.x, point.y))
    coords_list.append(copy(coords_list[0]))

    rect1_coords_tuple = tuple(coords_list)
    coords_list.clear()
    for point in segment2_corner_points:
        coords_list.append((point.x, point.y))
    coords_list.append(copy(coords_list[0]))
    rect2_coords_tuple = tuple(coords_list)


    rect1_polygon = shapely.Polygon(rect1_coords_tuple)
    rect2_polygon = shapely.Polygon(rect2_coords_tuple)


    if (rect1_polygon.is_valid or rect2_polygon.is_valid) is False:
        print("Invalid bounding rectangle, stopping")
        while True:
            time.sleep(60)

    return rect1_polygon.intersects(rect2_polygon)
def get_gross_segment_radius(node):
    return node.data.r + (constants.WALL_THICKNESS_LUMEN_RATIO) * node.data.r
def get_bounding_rectangle_for_segment(node):
    segment_p1 = None

    if node.parent is None:
        segment_p1 = Point(0.0,0.0)
    else:
        segment_p1 = node.parent.data.distal_point
    segment_p2 = node.data.distal_point

    corner_points = []

    x3 = [0.0,0.0]
    y3 = [0.0,0.0]

    deltax = (segment_p2.x - segment_p1.x)
    deltay = (segment_p2.y - segment_p1.y)

    gross_segment_radius = get_gross_segment_radius(node)

    if deltax == 0:
        x3[0] = segment_p1.x - gross_segment_radius
        y3[0] = segment_p1.y

        x3[1] = segment_p1.x + gross_segment_radius
        y3[1] = segment_p1.y

        corner_points.append(Point(x3[0], y3[0]))
        corner_points.append(Point(x3[1], y3[1]))
        corner_points.append(Point(x3[1], y3[1] + deltay))
        corner_points.append(Point(x3[0], y3[0] + deltay))
        return corner_points


    if deltay == 0:
        x3[0] = segment_p1.x
        y3[0] = segment_p1.y + gross_segment_radius

        x3[1] = segment_p1.x
        y3[1] = segment_p1.y - gross_segment_radius

        corner_points.append(Point(x3[0], y3[0]))
        corner_points.append(Point(x3[1], y3[1]))
        corner_points.append(Point(x3[1] + deltax, y3[1]))
        corner_points.append(Point(x3[0] + deltax, y3[0]))
        return corner_points


    m = (segment_p2.y - segment_p1.y) / (segment_p2.x - segment_p1.x)
    mp = (-1/m)

    k = mp**2 + 1

    perpendicular_dist_quadeq = QuadraticEquation(k, -2*k*segment_p1.x, (k*segment_p1.x**2) - gross_segment_radius**2)

    x3_sol_list = perpendicular_dist_quadeq.get_solution_list()
    #print('x3sol', x3_sol_list)
    x3[0] = x3_sol_list[0]
    x3[1] = x3_sol_list[1]

    y3[0] = mp*(x3[0] - segment_p1.x) + segment_p1.y

    y3[1] = mp*(x3[1] - segment_p1.x) + segment_p1.y

    corner_points.append(Point(x3[0], y3[0]))
    corner_points.append(Point(x3[1], y3[1]))
    corner_points.append(Point(x3[1]+deltax, y3[1]+deltay))
    corner_points.append(Point(x3[0] + deltax, y3[0] + deltay))

    return corner_points








def recompute_tree_characteristics_after_bifurcation(v_root, ibif_node):
    #balance bifurcation ratios at ibif

    ibif_node = balance_bifurcation_ratio_for_node(ibif_node)



    #balance bifurcations upto root
    v_root = balance_bifurcations_upto_root(v_root, ibif_node)
    #rescale root
    v_root.data.r = get_root_radius_rescaled(v_root)




    #Update radii for the tree

    v_root = update_radii_for_tree(v_root)

    return v_root


def update_radii_for_tree(v_root):
    tree_preorder_list = []
    binarytree.get_tree_preorder_list(v_root, tree_preorder_list)

    for node in tree_preorder_list:
        if node.left is not None:
            node.left.data.r = node.data.b_l * node.data.r
        if node.right is not None:
            node.right.data.r = node.data.b_r * node.data.r
    return v_root


def balance_bifurcations_upto_root(v_root, ibif_node):
    if ibif_node.parent is not None:
        current_parent = ibif_node.parent
        last_parent = current_parent
        while current_parent is not None:
            current_parent = balance_bifurcation_ratio_for_node(current_parent)
            last_parent = current_parent
            current_parent = current_parent.parent
        v_root = last_parent
    else:
        v_root = balance_bifurcation_ratio_for_node(v_root)
    return v_root

def get_cost_function_normalized_gradient_vector(v_root, ibif_node, r_supp):
    cost_list_x = [compute_cost_function(v_root)]
    cost_list_y = deepcopy(cost_list_x)

    dx = compute_dthresh(r_supp, binarytree.getLeafCount(v_root)) / 10
    dy = deepcopy(dx)

    v_root_copy = deepcopy(v_root)
    ibif_node_copy = binarytree.get_node_with_index(v_root_copy, ibif_node.data.index)

    for i in range(3):


        ibif_node_copy.data.distal_point.x = ibif_node_copy.data.distal_point.x + dx

        v_root_copy = recompute_tree_characteristics_after_bifurcation(v_root_copy, ibif_node_copy)
        cost_list_x.append(compute_cost_function(v_root_copy))

    v_root_copy = deepcopy(v_root)
    ibif_node_copy = binarytree.get_node_with_index(v_root_copy, ibif_node.data.index)

    for i in range(3):
        ibif_node_copy.data.distal_point.y = ibif_node_copy.data.distal_point.y + dy

        v_root_copy = recompute_tree_characteristics_after_bifurcation(v_root_copy, ibif_node_copy)
        cost_list_y.append(compute_cost_function(v_root_copy))



    pgx_list = [cost_list_x[i+1] - cost_list_x[i] for i in range(len(cost_list_x) - 1)]
    pgy_list = [cost_list_y[i+1] - cost_list_y[i] for i in range(len(cost_list_y) - 1)]


    pgx_avg = sum(pgx_list)/len(pgx_list)
    pgy_avg = sum(pgy_list)/len(pgy_list)


    gradient_vector = euclid.Vector2(pgx_avg, pgy_avg)
    normalized_gradient_vector = gradient_vector.normalize()
    return normalized_gradient_vector

def compute_cost_function(v_root):
    tree_inorder_list = []
    binarytree.get_tree_inorder_list(v_root, tree_inorder_list)

    c_cost = 0.0
    for node in tree_inorder_list:
        segment_p1 = None
        if node.parent is None:
            segment_p1 = Point(0.0,0.0)
        else:
            segment_p1 = node.parent.data.distal_point
        segment_p2 = node.data.distal_point

        lj = calculate_distance_bw_points(segment_p1, segment_p2)

        rj = node.data.r

        c_cost = c_cost + (lj * (rj**2))
    return c_cost


def _get_non_intersecting_segment_indices_task(v_root, new_terminal_point, q, chunk):
    non_intersecting_indices_list = []

    for node in chunk:
        is_intersecting = False

        til_copy = deepcopy(chunk)

        til_copy = [lc_node for lc_node in til_copy if not lc_node.data.index == node.data.index]

        prospective_connecting_segment_p1 = None
        if node.parent is None:
            prospective_connecting_segment_p1 = shapely.Point(0.0,0.0)
        else:
            prospective_connecting_segment_p1 = shapely.Point(node.parent.data.distal_point.x, node.parent.data.distal_point.y)
        prospective_connecting_segment_p2 = shapely.Point(node.data.distal_point.x, node.data.distal_point.y)

        prospective_connecting_segment_midpoint = shapely.Point((prospective_connecting_segment_p1.x+prospective_connecting_segment_p2.x)/2, (prospective_connecting_segment_p1.y+prospective_connecting_segment_p2.y)/2)
        new_terminal_point = shapely.Point(new_terminal_point.x, new_terminal_point.y)
        prospective_segment_linestring = shapely.LineString([prospective_connecting_segment_midpoint, new_terminal_point])

        for c_node in til_copy:
            test_segment_p1 = None
            if c_node.parent is None:
                test_segment_p1 = shapely.Point(0.0,0.0)
            else:
                test_segment_p1 = shapely.Point(c_node.parent.data.distal_point.x, c_node.parent.data.distal_point.y)
            test_segment_p2 = shapely.Point(c_node.data.distal_point.x, c_node.data.distal_point.y)

            test_segment_linestring = shapely.LineString([test_segment_p1, test_segment_p2])

            if test_segment_linestring.intersects(prospective_segment_linestring):
                is_intersecting = True
        if not is_intersecting:
            non_intersecting_indices_list.append(node.data.index)
    print(non_intersecting_indices_list)
    q.put(non_intersecting_indices_list)


def get_non_intersecting_segment_indices(v_root, new_terminal_point):
    til = []
    binarytree.get_tree_inorder_list(v_root, til)
    non_intersecting_indices_list = []


    chunk_size = math.ceil(len(til)/constants.CORE_COUNT)

    chunks_list = []
    for i in range(0, len(til), chunk_size):
        chunks_list.append(til[i:i+chunk_size])
    print(chunks_list)



    task_p_list=[]
    task_q = Queue()
    for chunk in chunks_list:
        fl_task_p = Process(target = _get_non_intersecting_segment_indices_task, args=(v_root, new_terminal_point, task_q, chunk))
        fl_task_p.start()
        task_p_list.append(fl_task_p)
    for task in task_p_list:
        task.join()
    for i in range(len(chunks_list)):
        non_intersecting_indices_list = non_intersecting_indices_list + task_q.get()
    return non_intersecting_indices_list


def add_segment_to_node(v_root, new_terminal_point_with_parent, current_segment_index, parent_index):
    #print('pi', parent_index, current_segment_index)
    tpl = []
    binarytree.get_tree_preorder_list(v_root, tpl)
    #print([x.data.index for x in tpl])
    node_iconn = binarytree.get_node_with_index(v_root, parent_index)
    Q_iconn = binarytree.getLeafCount(node_iconn) * constants.Q_TERM
    Q_inew = constants.Q_TERM

    new_point = new_terminal_point_with_parent['new_point']
    segment_start = None
    if node_iconn.parent is None:
        segment_start = Point(0.0,0.0)
    else:
        segment_start = node_iconn.parent.data.distal_point

    segment_stop = node_iconn.data.distal_point

    virt_current_segment_index = current_segment_index + 1

    mid_point = Point((segment_start.x + segment_stop.x)/2, (segment_start.y + segment_stop.y)/2)
    segment_ibif = Segment(index = virt_current_segment_index, distal_point = mid_point, b_l = None, b_r = None, Q = Q_iconn + Q_inew, R_red = None)
    virt_current_segment_index = virt_current_segment_index + 1
    l_inew = calculate_distance_bw_points(mid_point, new_terminal_point_with_parent['new_point'])
    segment_inew = Segment(index = virt_current_segment_index, distal_point = new_terminal_point_with_parent['new_point'], b_l = None, b_r = None, Q = Q_inew, R_red = (8*constants.ETA*l_inew)/math.pi)

    iconn_dc = deepcopy(node_iconn)

    ibif_node = Node(segment_ibif)
    ibif_node.left = iconn_dc
    iconn_dc.parent = ibif_node

    inew_node = Node(segment_inew)
    inew_node.parent = ibif_node
    ibif_node.right = inew_node

    if node_iconn.parent is None:
        v_root = ibif_node
        ibif_node.parent = None
    else:
        ibif_node.parent = node_iconn.parent
        if node_iconn.isLeftChild():
            node_iconn.parent.left = ibif_node
        elif node_iconn.isRightChild():
            node_iconn.parent.right = ibif_node

    return {'ibif': ibif_node, 'v_root': v_root} #return the bifurcating segment node to update radii

def get_root_radius_rescaled(v_root):
    segment_p1 = Point(0.0,0.0)
    segment_p2 = v_root.data.distal_point
    l_root = calculate_distance_bw_points(segment_p1, segment_p2)
    k_term = binarytree.getLeafCount(v_root)
    R_root = (constants.P_PERF - constants.P_TERM) / ((k_term + 1)* constants.Q_TERM)
    r_root = (((8*constants.ETA)/math.pi) * (l_root/R_root))**(1/4)
    return r_root

def balance_bifurcation_ratio_for_node(node): # WARNING: This function is also responsible for keeping the hydrodynamic resistances updated, bad design, refactor needed
    R_red_updated_subtree_right = update_R_red_for_subtree(node.right)
    node.right = R_red_updated_subtree_right

    R_red_updated_subtree_left = update_R_red_for_subtree(node.left)
    node.left = R_red_updated_subtree_left
    r_icon_inew_ratio = ((node.left.data.Q / node.right.data.Q) * (node.left.data.R_red/node.right.data.R_red)) ** (1/4)

    b_l = (1 + ((r_icon_inew_ratio) ** (-constants.BE))) ** (-1/constants.BE)
    b_r = (1 + ((r_icon_inew_ratio) ** (constants.BE))) ** (-1/constants.BE)

    node.data.b_l = b_l
    node.data.b_r = b_r

    R_red_update_subtree = update_R_red_for_subtree(node)



    return node

def calculate_distance_bw_points(point1, point2):
    return ((point2.x - point1.x)**2 + (point2.y - point1.y)**2)**(1/2)

def update_R_red_for_subtree(subtree_node):
    leaves = []
    binarytree.getLeafNodes(subtree_node, leaves)
    tree_postorder = []
    binarytree.get_tree_postorder_list(subtree_node, tree_postorder)
    tree_postorder_wo_leaves = [node for node in tree_postorder if node not in leaves]


    #Update reduced hydrodynamic resistances for leaves
    for leaf_node in leaves:
        segment_stop = leaf_node.data.distal_point
        segment_start = leaf_node.parent.data.distal_point
        lj = calculate_distance_bw_points(segment_start, segment_stop)
        leaf_node.data.R_red = (8*constants.ETA*lj) / math.pi
    for node in tree_postorder_wo_leaves:
        segment_stop = node.data.distal_point
        segment_start = None
        if node.parent is None:
            segment_start = Point(0.0,0.0)
        else:
            segment_start = node.parent.data.distal_point

        lj = calculate_distance_bw_points(segment_start, segment_stop)
        node.data.R_red = ((8*constants.ETA*lj) / math.pi) + (((node.data.b_l**4) / node.left.data.R_red) + ((node.data.b_r**4) / node.right.data.R_red)) ** (-1)
    return subtree_node


def compute_dcrit(prospective_distal_point, segment_start, segment_stop, dproj):
    dcrit = None
    if dproj is None:
        return
    if(dproj >= 0 and dproj <= 1):
        a = numpy.array([
            -segment_start.y + segment_stop.y,
            segment_start.x - segment_stop.y
        ])

        b = numpy.array([
            prospective_distal_point.x - segment_stop.x,
            prospective_distal_point.y - segment_stop.y
        ])

        length = calculate_distance_bw_points(segment_start, segment_stop)

        dcrit = abs(numpy.dot(a, b)) * (1/length)
    else:
        a = ((prospective_distal_point.x - segment_stop.x)**2 + (prospective_distal_point.y - segment_stop.y)**2)**(1/2)
        b = ((prospective_distal_point.x - segment_start.x)**2 + (prospective_distal_point.y - segment_start.y)**2)**(1/2)
        dcrit = min(a,b)
    #print("dcrit: ", dcrit)
    return dcrit

def compute_dproj(prospective_distal_point, segment_start, segment_stop):

    print(segment_start, segment_stop)

    a = numpy.array([
        segment_start.x - segment_stop.x,
        segment_start.y - segment_stop.y
    ])

    b = numpy.array([
        prospective_distal_point.x - segment_stop.x,
        prospective_distal_point.y - segment_stop.y
    ])

    length = calculate_distance_bw_points(segment_start, segment_stop)
    if round(length, 7) == 0.0000000:
        print("Segment with length zero encountered!!")
        return
    dproj = numpy.dot(a, b) * (1/(length**2))

    print(dproj)

    return dproj







def compute_dthresh(r_supp, k_term, lap=0):
    dthresh = ((math.pi*r_supp**2) / k_term**(1/2))
    #print("dthresh: ", (0.9**lap)*dthresh)

    return (0.9**lap)*dthresh

class Point:
    x: float = None
    y: float = None
    def __init__(self, x: float, y: float):
        self.x = x
        self.y = y

class Circle:
    radius: float = None
    center: Point = None
    def __init__(self, radius: float, center: Point):
        self.radius = radius
        self.center = center
        print(self.center.x, self.center.y, self.radius)
    def _is_point_in_circle(self, point: Point):
        rsquared_curr = (point.x - self.center.x)**2 + (point.y - self.center.y)**2
        if(rsquared_curr**(1/2) < self.radius):
            return True
        else:
            return False
    def get_2d_mesh_points(self, spacing: float):
        node_list=[]
        x = self.center.x - self.radius
        y = self.center.y - self.radius
        while y < self.center.y + self.radius:
            while x < self.center.x + self.radius:
                node = Point(x, y)
                if self._is_point_in_circle(node):
                    node_list.append(node)
                x = round(x+spacing, 7)
            y = round(y + spacing, 7)
            x = self.center.x - self.radius
        return node_list

    def get_random_coordinate(self, spacing):
        r = random.uniform(0,self.radius)
        theta = random.uniform(0, 2*math.pi)

        return Point(self.center.x + r*math.cos(theta), self.center.y + r*math.sin(theta))

class Segment:
    distal_point: Point = None
    b_l: float = None
    b_r: float = None
    Q: float = None
    R_red: float = None
    parent = None
    r = None

    def __init__(self, index = None, distal_point = None, b_l = None, b_r = None, Q = None, R_red = None):
        self.index = index
        self.distal_point = distal_point
        self.b_l = b_l
        self.b_r = b_r
        self.Q = Q
        self.R_red = R_red



class Node:
    parent = None
    left = None
    right = None
    data: Segment = None

    def __init__(self, segment: Segment):
        self.data = segment
    def isLeftChild(self):
        return self is self.parent.left
    def isRightChild(self):
        return self is self.parent.right

class QuadraticEquation:
    a = None
    b = None
    c = None

    def __init__(self, a, b, c):
        self.a = a
        self.b = b
        self.c = c

    def get_solution_list(self):
        abc = [self.a, self.b, self.c]
        if None in abc:
            return None
        else:
            sol_list = []
            sol_list.append((-self.b + math.sqrt( (self.b**2) - (4*self.a*self.c))) / (2*self.a))
            sol_list.append((-self.b - math.sqrt( (self.b**2) - (4*self.a*self.c))) / (2*self.a))

            return sol_list
