ls dS* > list_scopy 

scopy ("@list_scopy",
"1D_//@list_scopy", w1=INDEF, w2=INDEF, apertures="", bands="", beams="", apmodulus=0,
format="onedspec", renumber=no, offset=0, clobber=no, merge=no, rebin=yes,
verbose=no)

mv *100* 1D_spectra
mv *200* 1D_spectra

ls 1D_*.fits > list_txt

wspectext ("@list_txt",
"@list_txt//.txt", header=no, wformat="")

mv *00?.fits 1D_spectra

ls *000* > list_dat
ls *300* > list_err

!python /home/heloise/iraf/my_scripts/rename_txt.py

