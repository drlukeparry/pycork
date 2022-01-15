# -*- coding: utf-8 -*-

import sys
import os
import pyslm
import platform
import pyslm

#sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


# are we on linux
is_linux = 'linux' in platform.system().lower()

PY_VER = (sys.version_info.major, sys.version_info.minor)
