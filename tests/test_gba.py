#!/usr/bin/env python
"""
Created on Jun 23, 2014

@author: behry
"""
import os
import unittest

import numpy as np
import scipy.stats as stats

import gba


class GbATestCase(unittest.TestCase):

    def setUp(self):
        # load the test filterbank data
        self.path = os.path.dirname(__file__)
        self.fbdata = np.loadtxt(os.path.join(self.path, 'data', 'az_obs.txt'),
                                 unpack=True)
        training_data = os.path.join(self.path, '..', 'data', 'az_training.nc')
        self.g = gba.GbA()
        self.g.init(training_data)

    def test_likelihood(self):
        """
        Test that the likelihood pdf has the correct mean and
        covariance matrix.
        """
        nsim = 30
        mean = np.zeros((2))
        cov = np.zeros((2, 2))
        self.g.compute_likelihood(self.fbdata, 0.5, 'z', nsim, mean, cov)
        np.testing.assert_almost_equal(mean[0], 1.574, decimal=3)
        np.testing.assert_almost_equal(mean[1], 5.726, decimal=3)
        np.testing.assert_almost_equal(cov[0, 0], 0.039, decimal=3)
        np.testing.assert_almost_equal(cov[0, 1], 0.059, decimal=3)
        np.testing.assert_almost_equal(cov[1, 1], 0.152, decimal=3)
        self.assertEqual(cov[0, 1], cov[1, 0])

    def test_marginal_likelihood(self):
        """
        Test that the maximum likelihood estimates for the marginal likelihood
        pdfs are correct.
        """
        m = np.linspace(2, 8, 31)
        r = np.linspace(0, 2.0, 21)
        M, R = np.meshgrid(m, r)
        pos = np.empty(M.shape + (2,))
        pos[:, :, 0] = R; pos[:, :, 1] = M
        nsim = 30
        mean = np.zeros((2))
        cov = np.zeros((2, 2))
        self.g.compute_likelihood(self.fbdata, 0.5, 'z', nsim, mean, cov)
        rv = stats.multivariate_normal(mean, cov)
        p = rv.pdf(pos)
        mp = np.trapz(p, x=r, axis=0)
        rp = np.trapz(p, x=m, axis=1)
        mp_normed = mp / np.trapz(mp, m)
        rp_normed = rp / np.trapz(rp, r)
        mhat = m[np.argmax(mp_normed)]
        rhat = r[np.argmax(rp_normed)]
        np.testing.assert_almost_equal(mhat, 5.8, decimal=1)
        np.testing.assert_almost_equal(rhat, 1.6, decimal=1)

def suite():
    return unittest.makeSuite(GbATestCase, 'test')

if __name__ == '__main__':
    unittest.main(defaultTest='suite')
