!python /home/heloise/iraf/my_scripts/airmass.py

list="list_calibrate"
s1="c_"
while(fscan (list, s2, x, y) != EOF) {
calibrate (input=s2, output='c_'//s2, airmass=x, exptime=y, extinct=yes, flux=yes, extinction="home$paranal.txt", observatory="paranal", ignoreaps=no, sensitivity = "sens", fnu=no) 
}




