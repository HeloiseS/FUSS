import os
import numpy as np

f = open('list_ms', 'r')
i = 0
try:
    os.remove("list_hedit")
except:
    print "kittens"

for line in f:
    cal = "cal_"+line
    with open("list_hedit", 'a') as f2:
        f2.write(line[:-1]+" "+cal[:-9]+"\n") 
f.close()
#f2.close
