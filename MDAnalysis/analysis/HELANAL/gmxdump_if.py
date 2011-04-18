#!/usr/bin/python 

#gmxdump interface
#Modeled on xtcio

try:
	import psyco
	psyco.full()
except:
	pass	

import sys,os,subprocess

class BinaryError(Exception):
	pass


#Test for gmxdump in the path
location = subprocess.Popen(["which", "gmxdump"], stdout=subprocess.PIPE).communicate()[0]

if len(location)>0:
	pass
else:
	raise BinaryError, "Cannot find gmxdump in path"
	

class read_xtc:
	def __init__(self,filename):
		#self.source = open(filename, "rb")
		command = "gmxdump -f " +  filename
		p = subprocess.Popen(command, shell=True, bufsize=1,
          stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)

		self.filename = filename
		self.source = p.stdout
		self.contents = ''
		self.time = 0
		self.step = 0
		self.frame = [0]
		self.natoms = 0
		self.box = [0,0,0]
		self.cartesian = [0]
		self.precision = None
		
	def next_frame(self):
		cartesian = []
		frame = []
		current_atom = 0
		read_to_end = False
		for line in self.source:
			if line[3:9]== "natoms":
				contents=line.split()
				self.natoms = int(contents[1])
				self.time   = float(contents[4][5:])
				self.precision = float(contents[6])
			elif line[6:10] == "box[":
				contents = line.split()
				selection = int(contents[1][0])
				self.box[selection] == float(contents[2+selection][:-1])
			elif line[6:8] == "x[":
				#print line
				contents = [float(item[:-1]) for item in line[16:].replace("{"," ").split()]
				#print contents
				frame.append(contents[0])
				frame.append(contents[1])
				frame.append(contents[2])
				cartesian.append([ contents[0],contents[1],contents[2] ])
				#print int(contents[1][:-3]), self.natoms
				if int( line[:16].split("[")[1][:-3] ) == (self.natoms-1):
					#print "Final Atom"
					read_to_end = True
					break
		if read_to_end:
			self.cartesian = cartesian
			self.frame = frame
		return read_to_end

	def convert_coordinates(self):
		pass
		
	def next(self):
		correctly_read = self.next_frame()
		if not correctly_read: raise StopIteration
		return self.cartesian
	
	def __iter__(self):
		return self


if __name__ == "__main__":
		print "Test"
		if len(sys.argv) >1:
			print sys.argv[1]
			if os.path.exists(sys.argv[1]):
				print "First Frame"
				trajectory = read_xtc(sys.argv[1])
				trajectory.next_frame()
				for coordinate in trajectory.cartesian:
					print coordinate
				print "Last frame"
				while(trajectory.next_frame()):
					pass
				for coordinate in trajectory.cartesian:
					print coordinate
		