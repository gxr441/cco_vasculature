from .customizer_widgets import *
from .scad_variable import *
from .scad_interface import *

# register this extension
from ...core.extension_manager import default_extension_manager
default_extension_manager.register_pre_render(lambda root : get_scad_header())

