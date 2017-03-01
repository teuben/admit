""" .. _Ingest-at-api:

   **Ingest_AT** --- Ingests a (FITS) data cube.
   ---------------------------------------------

   This module defines the Ingest_AT class.
"""
import os, sys

import admit
from admit.AT import AT
from admit.Summary import SummaryEntry
import admit.util.bdp_types as bt
from admit.util.Image import Image
import admit.util.utils as utils
import admit.util.ratutil as casautil
import admit.util.ImPlot  as ImPlot
from admit.bdp.SpwCube_BDP import SpwCube_BDP
from admit.bdp.Image_BDP   import Image_BDP
from admit.util.AdmitLogging import AdmitLogging as logging

import numpy as np
import math

# @todo 
# - rel path should be to ../basename.fits
# = edge/box are the keywords to cut a cube
#   len(box) = 0       edge allowed
#              2       no edge allowed, they are z1,z2
#              4       edge allowed, they are blc,trc
#              6       no edge allowed
# - vlsr needs to be stored, in a veldef?   For this we need access to ASDM/Source.xml
#   and will need a small tool in util.py to parse the xml and return what we need
# - resolve the confusion between symlink= and copy=
# - allow LineCube instead of the default SpwCube
# - autobox?
# - if not a fits file, there is no smoothing
#   (obviously if it remains a symlink, of course we cannot change the data)
# - check if spectral reference is LSRK (e.g. some ALMA is TOPO, that's bad)
#   although some may not have it, e.g. test6503
#   SPECSYS = 'TOPO' or 'TOPOCENT' (casa 3.3.0)
#             'BARY'  (helio)
#   imhead->['reffreqtype']
#  - smooth and decimate option?
#  - vlsr=0.0 cannot be given?


class Ingest_AT(AT):
    """Ingest an image (cube) into ADMIT, normally to bootstrap a flow.

    See also :ref:`Ingest-AT-Design` for the design document.

    Ingest an image cube (usually FITS, but a CASA image or MIRIAD
    image are also natively supported by CASA) into CASA.  Expect I/O
    penalties if you use FITS or MIRIAD because CASA images are tiled.

    A number of selections and corrections to the cube can be made, as
    specified by the keywords. Notably, a sub-cube can be taken out of
    the cube by trimming off spatial and spectral edges, a primary beam
    correction can be made as well.

    Internally ADMIT will store images as 4D CASA images, with any missing
    3rd or 4th axis created redundantly (FREQ as axis 3, and POL as axis 4)

    **Keywords**

      **file**: string
               Input filename.

               Usually 'basename.fits' where 'basename' can be long, and can
               involve a directory hierarchy.  The admit directory will be
               'basename.admit', which will include the whole directory path
               if one was given.

               A symbolic link within ADMIT will then be used to resolve
               this as a local file.
    
               An absolute address is advised to be used if your working directory
               can change in a flow.

               If a CASA image is given, a symlink is used when no modifications
               (e.g. box=) are used, but this will cause uncertain integrity of the
               input BDP, as other clients could be modifying it.  A MIRIAD image
               can also be given, but I/O (like with FITS) will be slower.

               [no default]

      **basename**: string
               New short basename (like an alias) of the output file.
    
               This is meant to override a possibly long basename in the input
               file. The **alias** keyword in the baseclass can still be used
               to create **basename-alias** type filenames.
               Warning: if you use a dash in your basename, there is a good risk
               this will be confused with the alias separator in the basename
               in a flow, as these are "basename-alias.extension".
               Default: empty string, basename inherited from the input file name.

      **pb**:  string
               If given, this is the filename of the Primary Beam by which the
               input file needs to be multiplied to get a noise flat image for
               processing the ADMIT flow.
               The ALMA pipeline product of the Primary Beam should be in
               'longname.pb.fits' where the flux flat cube is
               'longname.pbcor.fits'  Note: PB correction is slow.
               Default:  empty string, no PB correction is done.

      **usepb**: boolean
               If True, the PB is actually used in an assumed flux flat input file to
               create a noise flat input BDP.  If False, it is assumed the user has
               given a noise flat file, and the PB can be used downstream to compute
               proper PB corrected fluxes.  Note that the PB is always stored as an 2D image,
               not as a cube.
               Default: True

      **box**:  blc,tlc (a list of 2, 4 or 6 integers)
               Select a box region from the cube.
               For example box=[xmin,ymin,xmax,ymax], which takes all channels, or
               [xmin,ymin,zmin,xmax,ymax,zmax], which also selects a range in channels.
               You can also select just some channels, with box=[zmin,zmax].
               As always, pixels and channels are 0 based in CASA.
               Arbitrary CASA regions are not implemented here, we only support
               a box/edge selection.

      **edge**:  Z_start,Z_end (a list of 1 or 2 integers)
               You can use edge= to remove edge channels, e.g. if box= was not specified,
               or when only an XY box was given. If box contains 2 or 6 numbers, any edge
               specification would be ignored. If one number is given, the edge rejection
               is the same at the upper and lower end. Default: not used.

      **smooth**: [nx,ny,[nz]] (a list of 2 or 3 integers)
               You can convolve your cube spatially, and optionally spectrally as well,
               by supplying the number of pixels by which it is convolved. Spatially the
               FWHM is used. 
               If nz is 1, a Hanning smooth is applied, else a boxcar of size nz is used.
               See also :ref:`Smooth-AT-api`, where a common beam can be computed or supplied.
               A future version should contain a decimation option.
               By default no smoothing is applied.
    
      **mask**: boolean
               If True, as mask needs to be created where the 
               cube has 0's. This option is automatically bypassed if the input
               CASA image had a mask. 
               [False]

      **vlsr**: float (km/s)
               VLSR of the source (km/s).  If not set, the frequency at the center of the
               band is compared with the rest frequency from the header, which we call VLSRc.
               If the input file is not FITS, or header items are missing if the input
               file is CASA or MIRIAD already, unexpected things may happen.
               This VLSR (or VLSRc) is added to the ADMIT summary, which will be used
               downstream in the flow by other AT's (e.g. LineID)
               Default: -999999.99 (not set).

      **restfreq**: float (GHz)
               An alternative method providing the source VLSR would be to specify the true
               restfreq (f0) where the fits header has a 'fake' restfreq (f). This technique
               is sometimes used by the PI to avoid complex high-z doppler calculations and
               supply the redshifted line directly.
               In this case VLSR = c * (1-f/f0), in the radio definition, with z in the optical
               convention of course. We call this VLSRf.
               NOTE: clarify/check if the "1+z" velocity scale of the high-z object is correct.
               Default: -1.0 (method not used). Units must be GHz!
      

    **Input BDPs**
      None. The input is specific via the file= keyword.

    **Output BDPs**

      **SpwCube**: count: 1
        The output spectral window cube. The convention is that the name of the BDP
        inherits from the basename of the input fits cube,and adding an extension
        "im". Each AT has a hidden keyword called alias=, use this keyword if
        you want to modify the cube name to "alias.im", and hence the BDP to
        "alias.im.bdp".   Note this is an exception from the usual rule, where
        alias= is used to create a dashed-prefix to an extension, e.g. "x.alias-y".

    Parameters
    ----------

    keyval : dictionary, optional
      Keyword-value pairs, directly passed to the contructor for ease of
      assignment.

    Attributes
    ----------

    _version : string
        Version ID for some future TBD use.  Also should not be documented here,
        as underscore attributes are for internal usage only.

    """

    #### DEPRECATED KEYWORDS BUT STILL ACTIVE IN CODE BY THEIR DEFAULT
    """

          **symlink** : True/False:
               If True, A symlink is kept to the input file without
               any conversion (if that was needed) This is used in
               those cases where your whole flow can work with the
               fits file, without need to convert to a CASA image
               Setting to True, also disabled all other processing
               (mask/region/pbcor) Use with caution!  [False]
               DEPRECATION

    """

    def __init__(self,**keyval):
        keys = {
            'file'    : "",        # fitsfile cube or map (or casa/miriad)
            'basename': "",        # override basename (useful for shorter names)
            'pb'      : "",        # PB cube or map
            'usepb'   : True,      # use PB, or was it just given for downstream
            'mask'    : True,      # define a mask where data==0.0 if no mask present
            'box'     : [],        # [] or z1,z2 or x1,y1,x2,y2  or x1,y1,z1,x2,y2,z2 
            'edge'    : [],        # [] or zl,zr - number of edge channels
            'smooth'  : [],        # pixel smoothing size applied to data (can be slow) - see also Smooth_AT
            'vlsr'    : -999999.0, # force a VLSR (see also LineID)
            'restfreq': -1.0,      # alternate VLSRf specification
            # 'symlink' : False,   # 
            # 'autobox' : False,   # automatically cut away spatial and spectral slices that are masked
            # 'cbeam'   : 0.5,     # channel beam variation allowed in terms of pixel size to use median beam
        }
        AT.__init__(self,keys,keyval)
        self._version = "1.2.2"
        self.set_bdp_in()                            # no input BDP
        self.set_bdp_out([(SpwCube_BDP, 1),          # one or two output BDPs
                        ])

    def summary(self):
        """Returns the summary dictionary from the AT, for merging
           into the ADMIT Summary object.

        Ingest_AT adds the following to ADMIT summary::

           Key      type        Description
         --------  ------       -----------
         fitsname  string      Pathless filename of FITS cube
         casaname  string      Pathless filename of CASA cube
         object    string      Object (or field) name
         naxis     integer     Number of axes
         naxisn    integer     size of axis n (n=1 to naxis0)
         crpix1    float       Reference pixel axis 1 (n=1 to naxis)
         crvaln    float       axis value at CRPIX1 (n=1 to naxis)
         ctypen    string      axis type 1 (n=1 to naxis0
         cdeltn    float       axis increment 1 (n=1 to naxis)
         cunitn    string      axis unit 1 (n=1 to naxis)
         equinox   string      equinox
         restfreq  float       rest frequency, Hz
         bmaj      float       beam major axis, radians
         bmin      float       beam minor axis, radians
         bpa       float       beam position angle, deg
         bunit     string      units of pixel values
         telescop  string      telescope name
         observer  string      observer name
         date-obs  string      date of observation
         datamax   float       maximum data value
         datamin   float       minimum data value
         badpixel  float       fraction of invalid pixels in the cube (a number between 0 and 1)
         vlsr      float       Object line-of-sight velocity (km/s)
        """
        # @todo the master list is in $ADMIT/admit/summary_defs.tab
        #       1) we duplicate that here....
        #       2) should it be with code?  or in $ADMIT/etc ?
        if hasattr(self,"_summary"):
            return self._summary
        else:
            return {}

    def run(self):
        # 
        self._summary = {}                  # prepare to make a summary here
        dt = utils.Dtime("Ingest")          # timer for debugging

        do_cbeam = True                     # enforce a common beam
        #
        pb = self.getkey('pb')
        do_pb = len(pb) > 0
        use_pb = self.getkey("usepb")
        # 
        create_mask = self.getkey('mask')   # create a new mask ?
        box   = self.getkey("box")          # corners in Z, XY or XYZ
        edge  = self.getkey("edge")         # number of edge channels to remove
        restfreq = self.getkey("restfreq")  # < 0 means not activated

        # smooth=  could become deprecated, and/or include a decimation option to make it useful
        #          again, Smooth_AT() does this also , at the cost of an extra cube to store
        smooth = self.getkey("smooth")      # 
        #
        vlsr = self.getkey("vlsr")          # see also LineID, where this could be given again

        # first place a fits file in the admit project directory (symlink)
        # this is a bit involved, depending on if an absolute or relative path was
        # give to Ingest_AT(file=)
        fitsfile = self.getkey('file')
        if fitsfile[0] != os.sep:
            fitsfile = os.path.abspath(os.getcwd() + os.sep + fitsfile)
        logging.debug('FILE=%s' % fitsfile)
        if fitsfile[0] != os.sep:
            raise Exception,"Bad file=%s, expected absolute name",fitsfile

        loc = fitsfile.rfind(os.sep)               # find the last '/'
        ffile0 = fitsfile[loc+1:]                  # basename.fits
        basename = self.getkey('basename')         # (new) basename allowed (allow no dots?)
        if len(basename) == 0:
            basename = ffile0[:ffile0.rfind('.')]  # basename
        logging.info("basename=%s" % basename)
        target = self.dir(ffile0)
        bdpfile = self.mkext(basename,"fits")
        #print 'target',target
        #print 'bdpfile',bdpfile

        if not os.path.exists(target) :
            #cmd = 'ln -s "%s" "%s"' % (fitsfile, target)
            cmd = 'ln -s "%s" "%s"' % (fitsfile, self.dir(bdpfile))
            logging.debug("CMD: %s" % cmd)
            os.system(cmd)


        if bdpfile == basename:
            raise Exception,"basename and bdpfile are the same, Ingest_AT needs a fix for this"
        b1  = SpwCube_BDP(bdpfile)
        b1.setkey("image", Image(images={bt.FITS:bdpfile}))
        self.addoutput(b1)
        if do_pb:
            raise Exception,"PB correction not implemented"

        dt.tag("start")

        # self._summarize(fitsfile, bdpfile, h, shape, taskargs)
        #self._summary[k].setTaskArgs(taskargs)
        self._summary['fitsname'] = SummaryEntry(fitsfile)
        taskargs = "bugus=True"
        for k in self._summary:
            self._summary[k].setTaskname("Ingest_AT")
            self._summary[k].setTaskID(self.id(True))
            self._summary[k].setTaskArgs(taskargs)
        
        

        dt.tag("done")
        dt.end()

    def _summarize(self, fitsname, casaname, header, shape, taskargs):
        """Convenience function to populate dictionary for
           items to add to the ADMIT Summary. The contract
           function is self.summary(), called by AT()

        """

        self._summary = {}
        self._summary['fitsname'] = SummaryEntry(fitsname)
        self._summary['casaname'] = SummaryEntry(casaname)

        # these are one-to-one match keywords 
        easy = [ 'object'  , 'equinox', 
                 'observer', 'date-obs',   'datamax', 
                 'datamin' , 'badpixel',   'vlsr',
               ]

        naxis = len(shape)
        self._summary['naxis'] = SummaryEntry(naxis)
        for i in range(naxis):
            j = i+1
            jay = str(j)
            easy.append('crpix'+jay)
            easy.append('ctype'+jay)
            easy.append('crval'+jay)
            easy.append('cdelt'+jay)
            easy.append('cunit'+jay)
            self._summary['naxis'+jay] = SummaryEntry(int(shape[i]))

        # FITS is only 8 chars.
        if 'telescope' in header:
            self._summary['telescop'] = SummaryEntry(header['telescope'])

        if 'imtype' in header:
            self._summary['bunit'] = SummaryEntry(header['bunit'])

        for k in easy:
            if k in header:
                self._summary[k] = SummaryEntry(header[k])

        if 'restfreq' in header:
            self._summary['restfreq']  = SummaryEntry(header['restfreq'][0])
            
        # These are in imhead returned as dictionaries {'unit','value'} 
        # so we have to munge them

        # convert beam parameters
        if 'beampa' in header:
            self._summary['bpa']  = SummaryEntry(taskinit.qa.convert(header['beampa'],'deg')['value'])
        if 'beammajor' in header:
            self._summary['bmaj'] = SummaryEntry(taskinit.qa.convert(header['beammajor'],'rad')['value'])
        if 'beamminor' in header:
            self._summary['bmin'] = SummaryEntry(taskinit.qa.convert(header['beamminor'],'rad')['value'])
        
        # Now tag all summary items with task name and task ID.

        for k in self._summary:
            self._summary[k].setTaskname("Ingest_AT")
            self._summary[k].setTaskID(self.id(True))
            self._summary[k].setTaskArgs(taskargs)
