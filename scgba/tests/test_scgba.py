#!/usr/bin/env python
"""
Created on Apr 14, 2016

@author: behry
"""

import os
import re
import unittest

import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from obspy import UTCDateTime, read
import scipy.stats as stats
import ipdb
import pyproj

from run_test import TestSCGbA

class SCGbATestCase(unittest.TestCase):

    def setUp(self):
        self.path = os.path.dirname(__file__)
        station = 'ABK'
        channel = 'HGZ'
        self.gba_data = os.path.join(self.path, 'data',
                                     'BO.%s..%s_gba.txt' % (station, channel))
        wf = os.path.join(self.path, 'data', 'test2.mseed.sorted')

        if not os.path.isfile(self.gba_data):
            faketimepath = \
            '/usr/lib/x86_64-linux-gnu/faketime/libfaketime.so.1'
            inventory = os.path.join(self.path, 'data', 'Inventory.xml')
            scgba_bin = os.path.join(self.path, '..', './scgba')
            msrtsimul_bin = os.path.join(self.path, './msrtsimul.py')
            seiscomp_bin = "/home/behry/eewamps/bin/seiscomp"
            tg = TestSCGbA(faketimepath, 'BO', station, '', channel,
                           inventory, self.gba_data, scgba_bin, msrtsimul_bin,
                           seiscomp_bin)
            tg.run(wf)
        abk_loc = (140.451, 39.0384)
        self.ot = UTCDateTime('1996-08-10T18:12:17.3Z')
        # load the sgba result
        st = read(wf)
        self.tr_abk = st.select(channel='HGZ', station='ABK')[0]
        m_data = os.path.join(self.path, 'data', 'AKT_test.nc')
        nc = Dataset(m_data)
        evlat = nc.groups['meta'].groups['event'].variables['latitude'][:]
        evlon = nc.groups['meta'].groups['event'].variables['longitude'][:]
        self.ev_mag = nc.groups['meta'].groups['event'].variables['magnitude'][:]
        pid_abk = nc.groups['filterbank'].groups['AKT019'].variables['ppxIdx'][:]
        self.pt_abk = self.tr_abk.stats.starttime + self.tr_abk.stats.delta * (pid_abk - 1)
        # self.pt_abk = pt_abk - self.ot
        self.toffset_abk = self.tr_abk.stats.starttime - self.pt_abk
        self.abk_mhat = nc.groups['filterbank'].groups['AKT019'].variables['mHat'][:]
        self.abk_rhat = nc.groups['filterbank'].groups['AKT019'].variables['rHat'][:]
        self.abk_t1 = nc.groups['filterbank'].groups['AKT019'].variables['time1'][:]

        g = pyproj.Geod(ellps='WGS84')
        az, baz, dist = g.inv(abk_loc[0], abk_loc[1], evlon, evlat)
        self.evdist = dist / 1000.

        fh = open(self.gba_data)
        pfloat = r'(\d+\.?\d*)'
        p = r'(\S+); (\S+); Mhat: ' + pfloat
        p += r'; Mbar: ' + pfloat
        p += r'; Sigma\^2: ' + pfloat
        p += r'; Rhat: ' + pfloat
        p += r'; Rbar: ' + pfloat
        p += r'; Sigma\^2: ' + pfloat
        self.mag = []
        self.magHat = []
        self.mag_unc = []
        self.dist = []
        self.dist_unc = []
        self.times = []
        for _l in fh.readlines():
            match = re.match(p, _l)
            if match:
                try:
                    ts = UTCDateTime(match.group(1))
                    ct = UTCDateTime(match.group(2))
                except:
                    print match.groups(1)
                    raise
                mhat, mbar, msig, rhat, rbar, rsig = map(float, match.groups()[2::])
                self.mag.append(mbar)
                self.magHat.append(mhat)
                self.mag_unc.append(np.sqrt(msig))
                self.dist.append(rbar)
                self.dist_unc.append(np.sqrt(rsig))
                self.times.append(ct - self.pt_abk)




    def test_ABK(self):
        """
        Test that the likelihood pdf has the correct mean and
        covariance matrix.
        """
        fig = plt.figure()
        xmin = self.ot - self.pt_abk
        xmax = 20
        ax1 = fig.add_axes([0.07, 0.1, 0.85, 0.5])
        ax1.errorbar(self.times, self.mag, yerr=self.mag_unc, color='g')
        # ax1.plot(self.abk_t1, self.abk_mhat, 'g:')
        ax1.vlines(0, 0, 7, 'r')
        ax1.hlines(self.ev_mag, xmin, xmax, color='g', linewidths=2,
                   linestyle='--')
        ax1.set_yticks([2, 3, 4, 5, 6])
        ax1.set_xticks(range(-4, 22, 2))
        ax1.set_xlim(xmin, xmax)
        ax1.set_ylim(2, 7)
        ax1.set_ylabel('Magnitude', color='g')
        for tl in ax1.get_yticklabels():
            tl.set_color('g')
        ax1.set_xlabel('Time since Pick [s]')
        ax11 = ax1.twinx()
        ax11.set_ylabel('Distance [km]', color='b')
        ax11.errorbar(self.times, self.dist, yerr=self.dist_unc, color='b')
        ax11.hlines(np.log10(self.evdist), xmin, xmax, color='b', linewidths=2,
                   linestyle='--')
        # ax11.plot(self.abk_t1, self.abk_rhat, 'b:')
        for tl in ax11.get_yticklabels():
            tl.set_color('b')
        ax11.set_yticks([np.log10(x) for x in [20, 30, 40, 50, 100, 200]])
        ax11.set_yticklabels([20, 30, 40, 50, 100, 200])
        ax11.set_ylim(np.log10(10), np.log10(200))
        ax11.set_xlim(xmin, xmax)

        ax2 = fig.add_axes([0.07, 0.6, 0.85, 0.3])
        ax2.plot(self.tr_abk.times() + self.toffset_abk, self.tr_abk.data, 'k')
        ax2.set_xticks([])
        ax2.set_yticks([])
        ymin, ymax = ax2.get_ylim()
        ax2.vlines(0, ymin, ymax, 'r')
        ax2.set_xlim(xmin, xmax)
        fig.savefig('scgba_abk.png', bbox_inches='tight', dpi=300)
        plt.show()

def suite():
    return unittest.makeSuite(SCGbATestCase, 'test')

if __name__ == '__main__':
    unittest.main(defaultTest='suite')
