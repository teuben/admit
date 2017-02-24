#!/usr/bin/env casarun
#
# This ADMIT script was auto-generated. Edit at your own risk.
# It expects to run from /data2/data1/teuben/ADMIT/N6503/.
# REMOVE existing project directory before running script for the first time!
#
import admit
import os

loglevel = 10            # 10=DEBUG, 15=TIMING 20=INFO 30=WARNING 40=ERROR 50=FATAL


# always start from clean slate
os.system('rm -rf test6503.admit')

# Master project.
p = admit.Project('test6503.admit', loglevel=loglevel)


# Flow tasks.
t0  = p.addtask(admit.Ingest_AT(alias='x', file='test6503.fits'))    
t1  = p.addtask(admit.CubeStats_AT(ppp=True), ['x'])
t2  = p.addtask(admit.CubeSum_AT(numsigma=4.0, sigma=99.0), ['x', t1])

#   can't run this, since moment maps are 2D
#t3  = p.addtask(admit.CubeSpectrum_AT(), ['x', t2])

if False:


    t4  = p.addtask(admit.PVSlice_AT(clip=0.3, smooth=[10, 10], width=5), ['x', t2])
    t5  = p.addtask(admit.PVCorr_AT(), [t4, t1])
    t6  = p.addtask(admit.LineID_AT(csub=[0, 0], minchan=3, numsigma=4.0, references='/home/teuben/ADMIT/admit/etc/co_lines.list', vlsr=25.0), [t3, t1, t5])
    t7  = p.addtask(admit.LineCube_AT(grow=10), ['x', t6])
    t8  = p.addtask(admit.Moment_AT(mom0clip=2.0, moments=[0, 1, 2]), [t7, t1])
    t9  = p.addtask(admit.CubeSpectrum_AT(), [t7, t8])

p.run()
