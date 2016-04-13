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


def get_start_time(fn):
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


def run(wf, delays=None, startupdelay=5):
    """
    Start msrtsimul and the waveform playback.
    """
    if not os.path.isfile(wf):
        raise PBError('%s does not exist.' % wf)

    # construct msrtsimul command
    msrt_cmnd = ["/home/behry/eewamps/bin/seiscomp", "exec",
                 "./tests/msrtsimul.py", "-c"]
    if delays is not None:
        msrt_cmnd += ["-d", delays]
    msrt_cmnd += ['-m', 'historic']
    t0 = get_start_time(wf)
    msrt_cmnd += ['-t', str(t0)]
    msrt_cmnd += [wf]
    t0 -= datetime.timedelta(seconds=startupdelay)
    print "Start time %s" % t0
    os.environ['LD_PRELOAD'] = '/usr/lib/x86_64-linux-gnu/faketime/libfaketime.so.1'
    ts = time.time()
    # Set system time in seconds relative to UTC now
    os.environ['FAKETIME'] = "%f" % (float(UTCDateTime(t0)) - ts)
    gba_cmnd = ["/home/behry/eewamps/bin/seiscomp", "exec", "./scgba", "--streams-allow",
                "BO.ABK..HGZ", "--inventory", "./tests/data/Inventory.xml",
                "--debug", "-I", "-"]
    msrt_proc = sp.Popen(msrt_cmnd, shell=False, env=os.environ,
                         stdout=sp.PIPE, stderr=sys.stderr)
    scgba_proc = sp.Popen(gba_cmnd, shell=False, env=os.environ,
                          stdin=msrt_proc.stdout)
    scgba_proc.wait()


run("./tests/data/test2.mseed.sorted")

