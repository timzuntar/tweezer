"""Unit tests for the data simulation functions"""

import sys
sys.path.append('../')

import unittest
import random
import math
import numpy as np

import calibration_generate_data as calgen
import data_generation as datagen

class TestSim(unittest.TestCase):
    """Tests scripts which simulate the Brownian motion of particle in trap.
    """

    def setUp(self):

        self.num_points = 10**3
        self.position = (3.,3.)
        self.k = (1e-6,2e-6)
        self.phi = 0.5

        self.expected_mean_data_prob = (3.209274, 3.4459371)
        self.last_time = 4.995
        self.expected_mean_traj_prob = (6.0879, 12.7741)

    def test_probabilistic_gen(self):

        random.seed(123)
        data = calgen.generate(self.k,293,self.phi,self.position,self.num_points,True)
        mean_data = np.mean(np.fabs(data), axis=0)

        self.assertTrue(np.allclose(mean_data, self.expected_mean_data_prob, rtol=1e-05, atol=1e-08))

    def test_simulation_gen(self):
        
        random.seed(123)
        time,trap,trajectory = datagen.generate(self.num_points,0.005,self.k,(1.,1.),(1e-5,2e-5),1e-6,8.9e-4,self.phi,self.position,293,1,True)
        mean_traj = np.mean(np.fabs(trajectory), axis=0)

        self.assertTrue(np.isclose(time[self.num_points - 1], self.last_time, rtol=1e-05, atol=1e-08))
        self.assertTrue(np.allclose(mean_traj, self.expected_mean_traj_prob, rtol=1e-05, atol=1e-08))


if __name__ == "__main__":
    unittest.main()