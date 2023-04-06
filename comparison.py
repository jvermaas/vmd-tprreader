import subprocess
import glob
import numpy as np
def loaddumpcoordinates(dumpfile):
	with open(dumpfile) as fin:
		output = subprocess.check_output("grep -a \"\\sx\\[.*\\]\" %s" % dumpfile, shell=True).decode("utf-8")
		outdata = []
		for line in output.split("\n")[:-1]:
			data = np.array(line[line.index("{")+1:line.index("}")].split(","), dtype=float)
			outdata.append(data)
		return np.vstack(outdata)
def loadtprtestcoordinates(dumpfile):
	with open(dumpfile) as fin:
		output = subprocess.check_output("grep -a \"x\\[.*\\]\" %s" % dumpfile, shell=True).decode("utf-8")
		outdata = []
		for line in output.split("\n")[:-1]:
			data = np.array(line[line.index(" ")+1:].split(" "), dtype=float)
			outdata.append(data)
		return np.vstack(outdata)
for tprfile in glob.glob("testtprs/*tpr"):
	subprocess.call(f"gmx dump -s {tprfile} > dump.log", shell=True, stderr=subprocess.DEVNULL)
	gmxcoordinates = 10 * loaddumpcoordinates("dump.log") # Multiply by 10 to convert positions to Angstroms
	subprocess.call(f"./tprtest {tprfile} > tprtest.log", shell=True)
	tprtestcoordinates = loadtprtestcoordinates("tprtest.log")
	print(tprfile, np.amax(np.abs(gmxcoordinates-tprtestcoordinates)))
exit()
