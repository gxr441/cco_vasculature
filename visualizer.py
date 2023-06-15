# cco_vasculature
#visualizer.py
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


import pygame
import lib

import binarytree
import math


screen = None
ZF = 10e2

SCROLL_SPEED = 0.5

def init_pygame():
    pygame.init()
    global screen
    screen = pygame.display.set_mode((1024,1028))
    screen.fill((255,255,255))

def quit_pygame():
    pygame.quit()

def draw_tree(v_root, r_supp):
    global ZF
    screen.fill((255,255,255))
    for event in pygame.event.get():
        #print(event)
        if event.type == pygame.QUIT:
            quit_pygame()
        elif event.type == pygame.MOUSEWHEEL:
            direction = event.y
            ZF = ZF + direction * SCROLL_SPEED * 0.2 * ZF

    tpl = []
    binarytree.get_tree_preorder_list(v_root, tpl)

    for node in tpl:
        segment_p1 = None
        if node.parent is None:
            segment_p1 = lib.Point(0.0,0.0)
        else:
            segment_p1 = node.parent.data.distal_point

        segment_p2 = node.data.distal_point
        pygame.draw.circle(screen, (0,0,0), (1024/2, r_supp*ZF), r_supp*ZF, width=1)
        line_width = math.floor(node.data.r*ZF)
        if line_width <= 1:
            line_width = 1
        pygame.draw.line(screen, (0,0,0), ((segment_p1.x*ZF+(1024/2)), -segment_p1.y*ZF), ((segment_p2.x*ZF+(1024/2)), -segment_p2.y*ZF), width=line_width)

    pygame.display.flip()
