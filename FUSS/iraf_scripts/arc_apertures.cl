noao
imred
specred

list="list_obj"
s2="cal_"
while(fscan (list, s1) != EOF) {
apall("ARC.fits", 2, output=s2//s1, references=s1, interactive=no, recenter=no, trace=no, background="none")
}

ls cal* > list_cal

ls SCIENCE*ms* > list_ms

!python /home/heloise/iraf/my_scripts/edit_list.py

list="list_hedit"
while(fscan (list, s1, s2) != EOF) {
hedit(images=s1, fields="REFSPEC1", value=s2)
}

