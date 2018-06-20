# My Cheat Sheet
---

## Data Reduction
---

Open terminal in directory containing the data.
```
>>> import FUSS.datred as r
>>> r.sort_red()
>>> metadataframe = r.Meta().data
File Created
>>> metadataframe = r.Meta(target_flag = 'TOO_SN').data
```

In IRAF
```
> display SCIENCE_##.fits
> imexam 
	→ ‘l’ on the spectrum to check it’s the right object and levels. Make note of loc on detector
 	→ ‘q’ to exit imexam 
> firstcal 
>Find apertures for Flat?  (yes): 
>Number of apertures to be found automatically: 0
>Edit apertures for Flat?  (yes): 
	→’n’ to create new apertures
	→:lower -35 (or other number)
	→:upper 35 (or other number just make sure aperture wide enough) 
	→ ‘q’
		* Trace apertures for Flat? (yes)
		* Fit traced positions for Flat interactively? (yes)
		* Fit curve to aperture 1 of Flat interactively? (yes)
			→ func=legendre, order=1, ‘f’ to refit if needed
			→ ‘q’
		* [Do the same for the other apertures]
		* Normalise blabla? (yes)
		* Fit spectra from Flat interactively? (yes)
		* Fit spectrum for aperture 1 for Flat.fits interactively? (yes)
			→ Use a splien function. Usually order = 35 is fine
			→ ‘d’ to delete points, ‘u’ to undo, ‘f’ ro refit
			→ ‘q’ to exit when happy with the fit 
		* Repeat fit for other apertures
> noao
> imred
> specred
> epar apall
	→  make sure input = @list_obj
	→ change the RON and 1/gain according to what is in metadata 
		(for other options see “how_to_iraf”)
	→ :q to save and exit, :go to save and run
		*Number of apertures to be found automatically (1): 
		*Resize apertures for SCIENCE_05?  (yes): 
		*Edit apertures for SCIENCE_05?  (yes): 
		* ‘m’ to create new apertures
		* ‘wj’, ‘wk’ to crop left and right
		* :upper and :lower to change bounds of apertures
		* ‘b’ to fit the background
			→ sample -low1:-low2, up1:up2 to change the background sample ranges
			→  :order 7 (ord other number) to change order of fit
			→ ‘f’ to fit again, ‘q’ to exit
		* ‘q’ to exit
			→ Trace apertures for SCIENCE_05? (yes)
			→ Fit traced positions for SCIENCE_05 interactively? (yes)
			→ Fit curve to aperture 1 of SCIENCE_05 interactively? (yes)
				* ‘d’ to delete discrepant points, ‘f’ to refit, ‘q’ when happy
			→ Repeat for other apertures
			→ Write apertures for SCIENCE_05 to database? (yes)
			→ Extract apertures spectra for SCIENCE_05? (yes)
			→ Review extracted spectrum for aperture 1 from SCIENCE_05?(yes)
			→ Check the extraction is sensible, then ‘q’
			→ Repeat visual check for other apertures
			→ Find apertures for SCIENCE_06? (yes) #next image
		* Repeat all images
 
> arc_aps
> epar identify
	→input = @list_cal
	→ :go
> epar dispcor
	→input   =   @list_ms  
	→output  =  d//@list_ms 
(First time around)
	→  (w1     =                INDEF) Starting wavelength
	→ (w2     =                INDEF) Ending wavelength
	→ (dw     =                INDEF) Wavelength interval per pixel

> disp2
(Can use other numbers than the ones below)
	→ (w1     =                2780.) Starting wavelength
	→ (w2     =                9330.) Ending wavelength
	→ (dw     =                  3.3) Wavelength interval per pixel
> ascii_pol
```

Back to the  terminal
```
>>> pol = r.LinearSpecPol(metadata = metadataframe, bin_size=15)
>>> poldataframe = pol.calculate()
4 Files per images... All good here
4 Files per images... All good here
4 Files per images... All good here
4 Files per images... All good here
Binning to  15 Angstrom
Index Error at  9320.6
Index Error at  9327.2
Index Error at  9320.6
Index Error at  9327.2

======== BEFORE BINNING ======
MEDIAN SNR 
94.3581085739
CENTRAL SNR at ( 6205.4  A)
108.462402142
======== AFTER BINNING ======
MEDIAN SNR / EXPECTED 
197.889145246 201.172163419
CENTRAL SNR / EXPECTED (at  6205.4  A)
232.628315614 231.242618344


Binning to  15 Angstrom
Index Error at  9320.6
Index Error at  9327.2
Index Error at  9320.6
Index Error at  9327.2
Binning to  15 Angstrom
Index Error at  9320.6
Index Error at  9327.2
Index Error at  9320.6
Index Error at  9327.2
Binning to  15 Angstrom
Index Error at  9320.6
Index Error at  9327.2
Index Error at  9320.6
Index Error at  9327.2
What do you want to name the polarisation file? 10as.pol
```

## PolData Object and Simple QU plot
---

Example in the python interpreter in a terminal:
```
>>> import FUSS.polmisc as F
>>> data = F.PolData(filename)
==== PolData - instance:   ====
Polarisation data initialised. If you want to add Stokes I use add_flux_data(). To find ISP use find_isp(). 

```

Note that the file must either be a text file with only numeric values and 9 columns in the following order (in which case it will be read with np.loadtxt):
wl	p	p\_r	 q	q\_r	u u\_r theta	theta\_r

Or it can have these as header and be tab separated (in which case it will be read with pandas).

The first option is for backwards compatibility with rpevious versions of FUSS, and the second option is the format that comes out of FUSS.datred currently. 

A third option is to give one array containing all the quantities we want in the right order: `[wl, p, p_r, q, q_r, u, u_r, theta, theta_r]`

You can then make a qu plot simply calling the `qu_plt` method:
```
>>> import matplotlib.pyplot as plt
>>> data.qu_plt()
>>> plt.show()
```

More details on the options of the `qu_plt` method, and more methods available with `PolData` can be found in the docstrings.


