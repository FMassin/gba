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
        self.nbands = 9
        self.nsim = 30
        self.g = gba.GbA(self.nbands, self.nsim)
        self.g.init(training_data)

    def test_likelihood(self):
        """
        Test that the likelihood pdf has the correct mean and
        covariance matrix.
        """
        mags = np.zeros(2 * self.nsim)
        r = np.zeros(2 * self.nsim)
        mean = np.zeros((2))
        cov = np.zeros((2, 2))
        self.g.process(self.fbdata, 0.5, 0)
        self.g.get_m_r(mags, r)
        self.g.get_mean_cov(mean, cov)
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
        pos[:, :, 0] = R
        pos[:, :, 1] = M
        mean = np.zeros((2))
        cov = np.zeros((2, 2))
        self.g.process(self.fbdata, 0.5, 0)
        self.g.get_mean_cov(mean, cov)
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

    def test_pdf(self):
        """
        Test the marginal pdfs computed by the C++ code.
        """
        m = np.linspace(2, 8, 31)
        r = np.linspace(0, 2.0, 21)
        mcp = m.copy()
        rcp = r.copy()
        pdf = np.zeros((m.size, r.size))
        self.g.process(self.fbdata, 0.5, 0)
        self.g.get_pdf(pdf, mcp, rcp)
        mhat = m[np.argmax(mcp)]
        rhat = r[np.argmax(rcp)]
        # Test the MAP values
        np.testing.assert_almost_equal(mhat, 5.8, decimal=1)
        np.testing.assert_almost_equal(rhat, 1.6, decimal=1)
        # Test that the pdfs have been properly normalized
        np.testing.assert_almost_equal(np.trapz(mcp, x=m), 1.0, decimal=1)
        np.testing.assert_almost_equal(np.trapz(rcp, x=r), 1.0, decimal=1)


def suite():
    return unittest.makeSuite(GbATestCase, 'test')

if __name__ == '__main__':
    unittest.main(defaultTest='suite')
