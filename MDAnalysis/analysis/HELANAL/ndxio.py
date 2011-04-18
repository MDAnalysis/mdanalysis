def read_ndx(filename):
	ndx = open(filename,"r")
	group_dictionary = {}
	current_group = ""
	for line in ndx:
		if line[0] == "[":
			current_group = line[2:-3]
			group_dictionary[current_group] = []
		else:
			contents = line.split()
			atom_numbers = [(int(item)-1) for item in contents]
			group_dictionary[current_group] += atom_numbers
	return group_dictionary