"""Unit tests for the force calculation module"""

import sys
sys.path.append('../')

import unittest
import random
import math
import numpy as np

import calibration_generate_data as calgen
import data_generation as datagen
import force_calc as forcecalc
import plotting as plt

class TestForce(unittest.TestCase):
    """Tests scripts which calculate forces on trapped particle and check if they correspond to expected values.
    """
    
    def setUp(self):
        random.seed(123)
        num_points = 10**3

        pos_1 = (-60.,12.)
        pos_2 = (5.,30.)
        self.k = (2.5e-6,0.5e-6)
        self.phi = 0.

        self.expected_coeffs = (0.0082505622,0.0082126888)
        self.expected_means_calc = (0.00217595, 0.00207052)        
        self.expected_means_axis = (0.82491509, 0.2130912, 0.19448143, -0.322048 )
        self.expected_means_displacements = (0.29525034, 0.85102993, 0.10935019, -0.51692861)
        self.expected_distance = 67.446274915669
        self.expected_variances = (0.01061657, 0.00394103, 0.01069194, 0.00372818)

        data_1 = calgen.generate(self.k,273,random.uniform(0,2*math.pi),(pos_1[0]+random.uniform(0,1),pos_1[1]+random.uniform(0,1)),num_points)
        data_2 = calgen.generate(self.k,273,random.uniform(0,2*math.pi),(pos_2[0]+random.uniform(-1,0),pos_2[1]+random.uniform(-1,0)),num_points)
        #random displacements represent inter-particle interaction

        self.times = calgen.generate_time(10**3)
        self.trajectories = np.hstack((data_1,data_2))
        position = np.hstack((pos_1,pos_2))
        self.positions = np.tile(position,(num_points,1))

    def test_force_single(self):
        random.seed(123)
        
        time,traps,trajectories = datagen.generate(1000,0.005, self.k,(2,1),(1e-6,1e-6),0.5e-6,9.7e-4,self.phi,(0.,0.),300, 1,False)
        
        _,means = forcecalc.calculate(time, trajectories[:, 0:2], traps[:, 0:2], self.k,self.phi)

        self.assertTrue(np.allclose(means, self.expected_means_calc, rtol=1e-05, atol=1e-08))

    def test_force_axial(self):
        random.seed(123)
        _,means,_,mean_distance = forcecalc.calculate_axial(self.times, self.trajectories, self.positions, self.k, self.k, self.phi, self.phi)

        np.testing.assert_approx_equal(mean_distance,self.expected_distance,6)
        self.assertTrue(np.allclose(means, self.expected_means_axis, rtol=1e-05, atol=1e-08))

    def test_force_displacements(self):
        random.seed(123)
        _,means,variances = forcecalc.displacement_calculation(self.trajectories, self.positions)

        self.assertTrue(np.allclose(means, self.expected_means_displacements, rtol=1e-05, atol=1e-08))
        self.assertTrue(np.allclose(variances, self.expected_variances, rtol=1e-05, atol=1e-08))

if __name__ == "__main__":
    unittest.main()