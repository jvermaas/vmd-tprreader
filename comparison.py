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
maxerror = None
identity = None
for tprfile in glob.glob("testtprs/*tpr"):
	subprocess.call(f"gmx dump -s {tprfile} > dump.log", shell=True)
	gmxcoordinates = 10 * loaddumpcoordinates("dump.log") # Multiply by 10 to convert positions to Angstroms
	subprocess.call(f"./tprtest {tprfile} > tprtest.log", shell=True)
	tprtestcoordinates = loadtprtestcoordinates("tprtest.log")
	error = np.amax(np.abs(gmxcoordinates-tprtestcoordinates))
	print(tprfile, error)
	if maxerror is None or error > maxerror:
		maxerror = error
		identity = tprfile
print(f"The overall largest coordinate error relative to gmx dump was {maxerror}, from file {tprfile}")