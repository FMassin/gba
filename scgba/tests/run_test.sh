#!/bin/bash

faketime -f "@1996-08-10 18:12:22.82" ~/eewamps/bin/seiscomp exec ./scgba \
--streams-allow "BO.*..HG?" \
--inventory-db ./tests/data/Inventory.xml \
--debug -I ./tests/data/test2.mseed.sorted
