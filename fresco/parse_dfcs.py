import csv  
import sys

def Filter(stateNum):
	dataFlag=False

	outfile = open('data/state'+str(stateNum)+'.csv', 'w', encoding='UTF8', newline='') 
	writer = csv.writer(outfile)
	writer.writerow(['theta','xsec'])

	#lines is a list with each item representing a line of the file
	for l in lines: 
		if l[0:8]=='#  Theta' in l:
			dataFlag=True
			# outfile = open('data/state'+str(stateNum)+'.csv', 'w') 
			continue
		if l[0:4]==' END' in l:
			dataFlag=False
			outfile.close()
			continue
		if dataFlag==True:
			x,y = l.split()
			# outfile.write(l)
			writer.writerow([x,y])
	infile.close()


if __name__ == "__main__":

	name_of_script = sys.argv[0]
	number_of_files = int(sys.argv[1])
	index = range(1, number_of_files+1) # start with file 201

	for i in index:
		infile = open('fort.20{0}'.format(i), 'r') 
		lines = infile.readlines() 
		Filter(i)

	# infile = open('fort.203', 'r') 
	# lines = infile.readlines() 
	# Filter(2)

	# infile = open('fort.204', 'r') 
	# lines = infile.readlines() 
	# Filter(3)

	# infile = open('fort.205', 'r') 
	# lines = infile.readlines() 
	# Filter(4)

	# infile = open('fort.206', 'r') 
	# lines = infile.readlines() 
	# Filter(5)