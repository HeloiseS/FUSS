del Zero*  ver-
del Master_Flat*  ver-
del Flat*  ver-
del ARC.fit*  ver-

noao
imred
ccdred

ls BIAS* > list_bias
 
zerocombine ("@list_bias",
output="Zero", combine="average", reject="minmax", ccdtype="", process=no,
delete=no, clobber=no, scale="none", statsec="", nlow=0, nhigh=1, nkeep=1,
mclip=yes, lsigma=3., hsigma=3., rdnoise="0.", gain="1.", snoise="0.",
pclip=-0.5, blank=0.)

ls ARC* > list_arc
ls SCIENCE* > list_obj
ls FLAT* > list_flat

ccdproc ("@list_flat",
output="", ccdtype="", max_cache=0, noproc=no, fixpix=no,
overscan=no, trim=no, zerocor=yes, darkcor=no, flatcor=no, illumcor=no,
fringecor=no, readcor=no, scancor=no, readaxis="line", fixfile="", biassec="",
trimsec="", zero="Zero", dark="", flat="", illum="", fringe="", minreplace=1.,
scantype="shortscan", nscan=1, interactive=no, function="legendre", order=1,
sample="*", naverage=1, niterate=1, low_reject=3., high_reject=3., grow=0.)

flatcombine ("@list_flat",
output="Flat", combine="average", reject="avsigclip", ccdtype="", process=no,
subsets=yes, delete=no, clobber=no, scale="mode", statsec="", nlow=1, nhigh=1,
nkeep=1, mclip=yes, lsigma=3., hsigma=3., rdnoise="0.", gain="1.",
snoise="0.", pclip=-0.5, blank=1.)

specred

apnormalize ("Flat ",
"Master_Flat", apertures="", references="", interactive=yes, find=yes,
recenter=yes, resize=yes, edit=yes, trace=yes, fittrace=yes, normalize=yes,
fitspec=yes, line=INDEF, nsum=10, cennorm=no, threshold=10.,
background="none", weights="none", pfit="fit1d", clean=no, skybox=1,
saturation=INDEF, readnoise="0.", gain="1.", lsigma=4., usigma=4.,
function="spline3", order=35, sample="*", naverage=1, niterate=0,
low_reject=3., high_reject=3., grow=0.)

ccdred

ccdproc ("@list_obj",
output="", ccdtype="", max_cache=0, noproc=no, fixpix=no, overscan=no,
trim=no, zerocor=yes, darkcor=no, flatcor=yes, illumcor=no, fringecor=no,
readcor=no, scancor=no, readaxis="line", fixfile="", biassec="", trimsec="",
zero="Zero", dark="", flat="Master_Flat", illum="", fringe="", minreplace=1.,
scantype="shortscan", nscan=1, interactive=no, function="legendre", order=1,
sample="*", naverage=1, niterate=1, low_reject=3., high_reject=3., grow=0.)

ccdproc ("@list_arc",
output="", ccdtype="", max_cache=0, noproc=no, fixpix=no,
overscan=no, trim=no, zerocor=yes, darkcor=no, flatcor=yes, illumcor=no,
fringecor=no, readcor=no, scancor=no, readaxis="line", fixfile="", biassec="",
trimsec="", zero="Zero", dark="", flat="Master_Flat", illum="", fringe="", minreplace=1.,
scantype="shortscan", nscan=1, interactive=no, function="legendre", order=1,
sample="*", naverage=1, niterate=1, low_reject=3., high_reject=3., grow=0.)

ccdproc ("@list_std",
output="", ccdtype="", max_cache=0, noproc=no, fixpix=no,
overscan=no, trim=no, zerocor=yes, darkcor=no, flatcor=yes, illumcor=no,
fringecor=no, readcor=no, scancor=no, readaxis="line", fixfile="", biassec="",
trimsec="", zero="Zero", dark="", flat="Master_Flat", illum="", fringe="", minreplace=1.,
scantype="shortscan", nscan=1, interactive=no, function="legendre", order=1,
sample="*", naverage=1, niterate=1, low_reject=3., high_reject=3., grow=0.)


imcombine ("@list_arc",
"ARC", headers="", bpmasks="", rejmasks="", nrejmasks="", expmasks="",
sigmas="", imcmb="$I", logfile="STDOUT", combine="average", reject="none",
project=no, outtype="real", outlimits="", offsets="none", masktype="none",
maskvalue="0", blank=0., scale="none", zero="none", weight="none", statsec="",
expname="", lthreshold=INDEF, hthreshold=INDEF, nlow=1, nhigh=1, nkeep=1,
mclip=yes, lsigma=3., hsigma=3., rdnoise="0.", gain="1.", snoise="0.",
sigscale=0.1, pclip=-0.5, grow=0.)

