import sys

#print "Nbr args : %d" % (len(sys.argv))

if len(sys.argv)>=3:
	readfrom=sys.argv[1]
	writeto=sys.argv[2]
else:
	print("Syntax: python %s <file_to_read> <file_to_write> [<id_label> <rest_of_cols> ['<char_sep>'] ]" % sys.argv[0])
	print(" <rest_of_cols> : columns to include in output file (ex: 1,2,3), starting in 1. Use - for including all.")
	exit(0)

if len(sys.argv)>=5:
	label=int(sys.argv[3])
	rest=sys.argv[4]
	if len(sys.argv)>=6:
		char_sep = sys.argv[5]
	else:
		char_sep = ' '
else:
	label=1
	rest="-"
	char_sep = ' '

if rest!="-" :
	some_cols=True
	id_cols=rest.rstrip().split(",")
else:
	some_cols=False

g = open(writeto, 'w')

print("Character separation: [%c]"%char_sep)

with open(readfrom, "r") as f:
	for line in f:
		cols = line.rstrip().split(char_sep)
		if label>len(cols) :
			print("ERROR, indicated column for label does not exist (%d columns detected in file)"% len(cols))
		if cols[0].find("#") != -1 :
			continue
		cad='%s ' % (cols[label-1])
		#raw_input("[%s]"% cad)
		index=1
		ind_w=1
		for item in cols:

			if (some_cols==True and (str(index) in id_cols)==False):
				add = ""
			elif (index==label):
				add = ""
			elif index!=label:
				add = "%d:%s " % (ind_w,item)
				ind_w=ind_w+1				
			
			index=index+1

			if index-1 == 1 :	# 1st col is nbr line
				continue
			if '\n' in item:
				raw_input("ERROR: There should be no \n in item")
			try:
				if float(item)==0.:
					continue
			except ValueError:
				continue
			cad = cad + add 
		#g.write(cad+'\r')
		g.write(cad+"\n")
		print(cad)
	#g.write('\r')
	f.close()
g.close()