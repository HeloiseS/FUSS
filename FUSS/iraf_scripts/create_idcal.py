import os
import fileinput

template = raw_input("What is the name of (or path to) the tempalte file: ")
template_num = raw_input("NUMBER in name of the template file (e.g 02): ")


with open('list_obj') as f:
    filenames = f.readlines() 
for filename in filenames:
    os.system('cp '+template+' '+'idcal_'+filename[:-6])

#template_num='01'#just for the test
keyword1 = "SCIENCE_"+str(template_num)

for filename in filenames:
    replacement1 = filename[:-6]
    for line in fileinput.input('idcal_'+filename[:-6], inplace=True):
        line = line.rstrip()
        if keyword1 in line:
            line = line.replace(keyword1, replacement1)
        print line

os.system('mv idcal_* database')
