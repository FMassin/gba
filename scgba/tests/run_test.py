#!/usr/bin/env python
"""
Created on Apr 13, 2016

@author: behry
"""

import datetime
import os
import subprocess as sp
import sys
import time

from obspy import UTCDateTime

from seiscomp3 import Config, System
import seiscomp3.Kernel
import seiscomp3.IO


class PBError(Exception): pass


class TestSCGbA:

    def __init__(self, libfaketime, network, station, location, channel,
                 inventory, fout, scgba_bin, msrtsimul_bin, seiscomp_bin):
        self.libfaketime = libfaketime
        self.network = network
        self.station = station
        self.location = location
        self.channel = channel
        self.seiscomp_bin = seiscomp_bin
        self.msrtsimul_bin = msrtsimul_bin
        self.scgba_bin = scgba_bin
        self.stream_allow = ".".join((network, station, location, channel))
        self.inventory = inventory
        self.fout = fout
        # "/home/behry/eewamps/bin/seiscomp"
        # "./tests/msrtsimul.py"
        # "./scgba"
        # "./tests/data/Inventory.xml"
        # "BO.ABK..HGZ"
        #

    def get_start_time(self, fn):
        """
        Find the earliest end-time of all records in the waveform file.
        """
        stream = seiscomp3.IO.RecordStream.Open('file://%s' % fn)
        input = seiscomp3.IO.RecordInput(stream, seiscomp3.Core.Array.INT,
                                         seiscomp3.Core.Record.SAVE_RAW)
        tmin = datetime.datetime.utcnow()
        while True:
            try:
                rec = input.next()
            except:
                break
            if not rec:
                break
            te = rec.endTime().toString("%FT%T.%4fZ")
            ts = rec.startTime().toString("%FT%T.%4fZ")
            dts = datetime.datetime.strptime(ts, "%Y-%m-%dT%H:%M:%S.%fZ")
            dte = datetime.datetime.strptime(te, "%Y-%m-%dT%H:%M:%S.%fZ")
            if dte < tmin:
                tmin = dte
                Id = rec.streamID()
        return tmin


    def run(self, wf, delays=None, startupdelay=5):
        """
        Start msrtsimul and the waveform playback.
        """
        if not os.path.isfile(wf):
            raise PBError('%s does not exist.' % wf)

        # construct msrtsimul command
        msrt_cmnd = [self.seiscomp_bin, "exec", self.msrtsimul_bin, "-c"]
        if delays is not None:
            msrt_cmnd += ["-d", delays]
        msrt_cmnd += ['-m', 'historic']
        t0 = self.get_start_time(wf)
        msrt_cmnd += ['-t', str(t0)]
        msrt_cmnd += [wf]
        print msrt_cmnd
        t0 -= datetime.timedelta(seconds=startupdelay)
        print "Start time %s" % t0
        os.environ['LD_PRELOAD'] = self.libfaketime
        ts = time.time()
        # Set system time in seconds relative to UTC now
        os.environ['FAKETIME'] = "%f" % (float(UTCDateTime(t0)) - ts)
        gba_cmnd = [self.seiscomp_bin, "exec", self.scgba_bin,
                    "--streams-allow", self.stream_allow, "--inventory",
                    self.inventory, "--debug", "-I", "-"]
        msrt_proc = sp.Popen(msrt_cmnd, shell=False, env=os.environ,
                             stdout=sp.PIPE, stderr=sys.stderr)
        scgba_proc = sp.Popen(gba_cmnd, shell=False, env=os.environ,
                              stdin=msrt_proc.stdout, stdout=sp.PIPE)
        output = scgba_proc.communicate()[0]
        fh = open(self.fout, 'w')
        fh.write(output)
        fh.close()

# run("./tests/data/test2.mseed.sorted")

