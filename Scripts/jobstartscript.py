import subprocess
import glob
import os
koverk2="2.0"
for doverp in ["1.0"]:
	# make the file
	outputdir="doverp_"+doverp+"_"+"koverk2_"+koverk2
	# open the submission script template and modify it
	subscript = open("myscript_template.pbs").read()
	subscript = subscript.replace("INSERT_DOVERP",doverp)
	subscript = subscript.replace("INSERT_KOVERK2",koverk2)
	subscript = subscript.replace("INSERT_OUTPUTDIR",outputdir)
	scriptname = "myscript_"+outputdir+".pbs"
	f = open(scriptname,"w")
	f.write(subscript)
	f.close()
	#launch it
	subprocess.call(["sbatch", "-p", "devel",scriptname])

