from mecode import G
import math
import os
from collections import defaultdict

import numpy as np




def spiral(radius, over, x_center, y_center, direction = 'CW'):
    g.abs_move(x=x_center, y=(y_center-over*0.5))
    repeats = 2*(radius/over)
    count = 0
    sign=-1
    for i in range(repeats):
        count=count+1
        sign=(sign)^(count+1)
        g.arc(direction=direction, radius = 0.5*(over*count), x=0, y=(over*count*sign))
        