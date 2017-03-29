import os

f1 = open('list_dat_flux', 'r')
f2 = open('list_err_flux', 'r')

for line in f1:
    name = line[:-23]
    ap = line [-11:-10]
    new_name = name+"_ap"+ap+".txt"
    os.system("mv "+line[:-1]+" "+new_name)
    #print new_name

for line in f2:
    name = line[:-23]
    ap = line [-11:-10]
    new_name = name+"_ap"+ap+"_err.txt"
    #print new_name
    os.system("mv "+line[:-1]+" "+new_name)


