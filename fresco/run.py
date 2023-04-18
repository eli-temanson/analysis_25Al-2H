import os
from sys import platform

def run_fresco(filename):
	if platform == "linux" or platform == "linux2" or platform == "darwin":
		os.system("fresco < {0}.in > {0}.out".format(filename))
	elif platform == "win32":
		os.system("fresco < {0}.in > {0}.out".format(filename))
	print("Fresco run completed!")


def run_sfresco(filename):
	if platform == "linux" or platform == "linux2" or platform == "darwin":
		os.system("sfresco < {0}.min".format(filename))
	elif platform == "win32":
		os.system("sfresco < {0}.min".format(filename))


if __name__ == "__main__":

	filename='2-25Al_dn-error'
	# filename='1-25Al_dn-eli'
	run_fresco(filename)
	# run_sfresco(filename)
	os.system("python ../scripts/parse_dfcs.py 5")