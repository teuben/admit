""" .. _ContinuumSub-at-api:

   **ContinuumSub_AT** --- Subtracts continuum emission from a cube.
   -----------------------------------------------------------------

   This module defines the ContinuumSub_AT class.
"""
from admit.AT import AT
from admit.Summary import SummaryEntry
import admit.util.bdp_types as bt
import admit.util.casautil as casautil
import admit.util.Image as Image
import admit.util.Line as Line
import admit.util.ImPlot as ImPlot
import admit.util.Segments as Segments
from admit.bdp.SpwCube_BDP import SpwCube_BDP
from admit.bdp.Image_BDP import Image_BDP
from admit.bdp.LineList_BDP import LineList_BDP
from admit.bdp.LineSegment_BDP import LineSegment_BDP
import admit.util.utils as utils
import admit.util.filter.Filter1D as Filter1D
from admit.util.AdmitLogging import AdmitLogging as logging
import numpy as np
import numpy.ma as ma
from copy import deepcopy

import types
import os

class ContinuumSub_AT(AT):
    """Continuum subtraction from a cube. Produces a line cube and continuum map.

    Based on line segments found (usually from LineSegments_AT from a CubeStats_BDP)
    this AT will fit the continuum in channels not covered by the line segments.
    The continuum segments can also be explicitly given, as can also be done in
    Ingest_AT if they are known beforehand.  This AT is meant for the automated
    continuum subtraction via LineSegments_AT.

    Although both are optional, you need to given either a LineSegment list, or
    explicitly define the **contsub** continuum segments.
    
    **Keywords**

        **contsub**: list of tuples 
            List a set of channel segments, 0 based and edges included,
            where the continuum is fitted. For example [(100,200),(800,900)]
            See also Ingest_AT for an alternate method.
            **Default**: []

        **pad**: integer
            Widen the line segments from a LineList_BDP if that was given.
            For insane reasons negative numbers are also allowed to narrow
            the segments. It will apply pad channels on either side of the segments.
            **Default**: 5
        
        **fitorder**: integer
            Order of continuum fit polynomial.
            **Default**: 0

    **Input BDPs**

        **SpwCube_BDP**: count: 1
            Input spectral window cube. 

        **LineList_BDP**: count: 1 (optional)
            Optional linelist, usually derived from LineSegments_AT, although
            LineID_AT should also work. If given, the contsub= is ignored.
 
    **Output BDPs**

        **SpwCube_BDP**: 1
            Output Line Cube which should now be continuum free.
            New extension will be ".lim"

        **Image_BDP**: 1
            Output Continuum Map. 
            New extension will be ".cim"

    **Graphics Produced**
        TBD

    Parameters
    ----------
        keyval : dictionary, optional

    Attributes
    ----------
        _version : string
    """

    def __init__(self, **keyval):
        keys = {
            "contsub"    : [],      # see also Ingest_AT (although deprecated there now)
            "pad"        : 5,       # see also LineCube_AT
            "fitorder"   : 0,       # polynomial order
        }
        AT.__init__(self,keys,keyval)
        self._version = "1.0.2"
        self.set_bdp_in([(SpwCube_BDP,      1, bt.REQUIRED),        # input spw cube 
                         (LineList_BDP,     1, bt.OPTIONAL),        # will catch SegmentList as well
                        ])
        self.set_bdp_out([(SpwCube_BDP,  1),                        # output line cube (.lim)
                          (Image_BDP,    1)],                       # output cont map  (.cim)
                        )

    def run(self):
        """ The run method creates the BDP.

            Parameters
            ----------
            None

            Returns
            -------
            None
        """
        dt = utils.Dtime("ContinuumSub")         # tagging time
        self._summary = {}                       # an ADMIT summary will be created here

        contsub = self.getkey("contsub")
        pad = self.getkey("pad")
        fitorder = self.getkey("fitorder")

        # x.im -> x.cim + x.lim

        # b1  = input spw BDP
        # b1a = optional input {Segment,Line}List
        # b1b = optional input Cont Map (now deprecated)
        # b2  = output line cube
        # b3  = output cont map
        b1 = self._bdp_in[0]
        f1 = b1.getimagefile(bt.FITS)

        b1a = self._bdp_in[1]
        # b1b = self._bdp_in[2]      
        b1b = None                   # do not allow continuum maps to be input

        f2 = self.mkext(f1,'lim')
        f3 = self.mkext(f1,'cim')
        f3a = self.mkext(f1,'cim3d')      # temporary cube name, map is needed
        b2 = SpwCube_BDP(f2)
        b3 = Image_BDP(f3)

        self.addoutput(b2)
        self.addoutput(b3)

          
        dt.tag("done")
        dt.end()

    def summary(self):
        """Returns the summary dictionary from the AT, for merging
           into the ADMIT Summary object.

           ContinuumSub_AT adds the following to ADMIT summary:

           .. table::
              :class: borderless

              +-------------+--------+---------------------------------------------+
              |   Key       | type   |    Description                              |
              +=============+========+=============================================+
              |continuumsub | list   |    Info about ContinuumSub produced         |
              +-------------+--------+---------------------------------------------+

           Parameters
           ----------
           None

           Returns
           -------
           dict
               Dictionary of SummaryEntry
        """
        if hasattr(self,"_summary"):
            return self._summary
        else:
            return {}
