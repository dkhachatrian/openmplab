import re
import sys #recognize command line args
import __future__ #in case using Python 2

def main():
	""" Takes in input string of file in current directory.
	Calculates average numbers and appends to end of file.
	"""
	if len(sys.argv) != 2:
		print("Not the correct number of arguments!")
		return

	f_string = sys.argv[1]

	with open(f_string, 'r+') as inf:
		f_sum = 0
		t_sum = 0
		i = 0
		for line in inf:
			line = re.sub(r'\s', r'', line)
			parts = line.split(':') # how it's formatted
			num = float(parts[-1])
			if parts[0] == 'FUNCTIME':
				i+=1
				f_sum += num
			elif parts[0] == 'TOTALTIME':
				t_sum += num

		f_ave, t_ave = f_sum/i, t_sum/i

		inf.seek(0,2) #be at end of file
		inf.write('\n')
		inf.write("Number of iterations: {0}\n".format(i))
		inf.write("Average for FUNC TIME: {0}\n".format(f_ave))
		inf.write("Average for TOTAL TIME: {0}".format(t_ave))

	print("Statistics appended to {0}".format(f_string))
	return

main()