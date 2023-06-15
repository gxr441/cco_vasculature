__nothing__ = None

from .std import union, attachable

#============ attachable add =============
#enhance the add function of the attachable OpenSCADObject so it can be used
#properly: cf. 07-libs-bosl2-attachable.py
attachable_default_add = attachable.add

def attachable_add(self, c):
    if len(self._children) == 0:
        attachable_default_add(self, c)
    elif len(self._children) == 1:
        attachable_default_add(self, union()(c))
    else:
        assert(len(self._children) == 2)
        self._children[1].add(c)

attachable.add = attachable_add
#============ attachable add end =============

