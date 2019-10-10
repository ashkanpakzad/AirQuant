import numpy as np
import scipy

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

        # TODO: Reformat once spline understood.
        def tangent(self, instance):
            # find the tangential direction along the spline at a given point (instance).
            road_arclenth = self.interval * instance
            # tangent position in CT
            CT_pos = spline(road_arclenth)
            # B-spline rep of ND curve.
            tck = scipy.interpolate.splprep
            # tangent direction in CT
            ct_dir = scipy.interpolate.splev(x, tck, der=1)
            return ct_dir

        def CT_perp(self, instance, ct_dir)
            # compute plane perpendicular to tangent given direction and instance.
            # TODO: create optimisation so that the size of this plane can be changed.