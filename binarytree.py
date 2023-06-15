# cco_vasculature
#binarytree.py
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

def hasPath(root, arr, index):

    # if root is None there is no path
    if (not root):
        return False

    # push the node's value in 'arr'
    arr.append(root.data.index)

    # if it is the required node
    # return true
    if (root.data.index == index):
        return True

    # else check whether the required node
    # lies in the left subtree or right
    # subtree of the current node
    if (hasPath(root.left, arr, index) or
        hasPath(root.right, arr, index)):
        return True

    # required node does not lie either in
    # the left or right subtree of the current
    # node. Thus, remove current node's value
    # from 'arr'and then return false
    arr.pop(-1)
    return False
def get_tree_preorder_list(root, list):
    if root:

        list.append(root)
        get_tree_preorder_list(root.left, list)
        get_tree_preorder_list(root.right, list)

def get_tree_inorder_list(root, list):
    if root:
        get_tree_inorder_list(root.left, list)
        list.append(root)
        get_tree_inorder_list(root.right, list)

def get_node_with_index(root, index):
    tree_inorder_list = []
    get_tree_inorder_list(root, tree_inorder_list)

    for node in tree_inorder_list:
        if node.data.index == index:
            return node

def getLeafCount(node):
    if node is None:
        return 0
    if(node.left is None and node.right is None):
        return 1
    else:
        return getLeafCount(node.left) + getLeafCount(node.right)

def get_tree_postorder_list(root, list):
    if root:
        get_tree_inorder_list(root.left, list)
        get_tree_inorder_list(root.right, list)
        list.append(root)

def getLeafNodes(root, out_list):

    # If node is null, return
    if (not root):
        return

    # If node is leaf node,
    # print its data
    if (not root.left and
        not root.right):
        out_list.append(root)
        return

    # If left child exists,
    # check for leaf recursively
    if root.left:
        getLeafNodes(root.left, out_list)

    # If right child exists,
    # check for leaf recursively
    if root.right:
        getLeafNodes(root.right, out_list)
