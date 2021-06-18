#!/usr/bin/env python2
import os
i=0
while i<30:
	i+=1
	os.remove("readingpressure.txt")
	out_file = open("readingpressure.txt","w")
	var = i*0.1*0.156
	out_file.write("%lf" % var)
	out_file.close()
	tuamamma = "nohup mosrun ./DNA_duplex > results.",i ," & "
	os.system(tuamamma)
