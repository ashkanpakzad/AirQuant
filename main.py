import numpy as np

class Road:

    def __init__(self, generation, origin_node, terminal_node, voxels):
        self.gen = generation
        self.onode = origin_node
        self.tnode = terminal_node
        self.vox = voxels
    # airway generation, node that airway segment originates from, node that airway segment terminates at.

    class Spline:
        # TODO: need to reformat to attributes of spline.
        def __init__(self, interval):
            self.interval = interval

        def tangent(self, instance):
            road_arclenth = self.interval * instance
            CT_loc = spline(road_arclenth)
            CT_dir =