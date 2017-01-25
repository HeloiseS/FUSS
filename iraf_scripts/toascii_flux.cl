ls c_dS* > list_scopy_flux 

scopy ("@list_scopy_flux",
"1D_//@list_scopy_flux", w1=INDEF, w2=INDEF, apertures="", bands="", beams="", apmodulus=0,
format="onedspec", renumber=no, offset=0, clobber=no, merge=no, rebin=yes,
verbose=no)

mv *100* 1D_spectra
mv *200* 1D_spectra

ls 1D_*.fits > list_txt_flux

wspectext ("@list_txt_flux",
"@list_txt_flux//.txt", header=no, wformat="")

mv *00?.fits 1D_spectra

ls *000* > list_dat_flux
ls *300* > list_err_flux

!python /home/heloise/iraf/my_scripts/rename_flux_txt.py



