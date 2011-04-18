#!/usr/bin/python
#
#Note to self. In nanometers
import xdrlib
import struct
import copy
import itertools
try:
	import psyco
	psyco.full()
except:
	pass	

class read_xtc:
	def __init__(self,filename):
		self.source = open(filename, "rb")
		self.contents = ''
		self.unpacker = xdrlib.Unpacker(self.contents)
		self.time = 0
		self.step = 0
		self.frame = [0]
		self.natoms = 0
		self.box = 0
		self.cartesian = [0]

	def sizeofint(self,intsize):
		num = 1
		num_of_bits = 0
		while(intsize >= num and num_of_bits < 32):
			num_of_bits = num_of_bits + 1
			num <<=  1
		return num_of_bits

	def sizeofints(self,length, sizearray):
		num_of_bytes = 1
		num_of_bits = 0
		bytes = []
		bytes.append(1)
		i = 0
		for i in range(length):
			tmp = 0
			for bytecnt in range(num_of_bytes):
				tmp = bytes[bytecnt] * sizearray[i] + tmp
				bytes[bytecnt] = tmp & 0xff
				tmp >>= 8
			bytecnt = num_of_bytes
			while(tmp != 0):
				bytecnt += 1
				bytes.append(tmp & 0xff)
				tmp >>= 8
			num_of_bytes = bytecnt
		num = 1
		num_of_bytes -= 1
		while (bytes[num_of_bytes] >= num):
			num_of_bits = num_of_bits + 1
			num *= 2
		return (num_of_bits + num_of_bytes * 8)

	def decodebits(self,stream, num_of_bits):
		mask = (1 << num_of_bits) - 1
		num = 0
		cnt = stream[0]
		lastbits = stream[1]
		lastbyte = stream[2]
		for item in range(num_of_bits,7,-8):
		#while(num_of_bits >= 8):
			lastbyte = (lastbyte << 8) | stream[cnt+3]
			#lastbyte &= ((1 << 32) - 1)
			cnt += 1
			num |= (lastbyte >> lastbits) << (num_of_bits - 8)
			#num_of_bits -= 8
		num_of_bits -= num_of_bits/8 * 8
		if (num_of_bits > 0):
			if(lastbits < num_of_bits):
				lastbits += 8
				lastbyte = (lastbyte << 8) | stream[cnt+3]
				#lastbyte &= ((1 << 32) - 1)
				cnt+=1
			lastbits -= num_of_bits
			num |= (lastbyte >> lastbits) & ((1 << num_of_bits) -1)
		lastbyte &= ((1 << 32) - 1)
		num &= mask
		stream[0] = cnt
		stream[1] = lastbits
		stream[2] = lastbyte
		return num, stream

	def decodeints(self,stream, num_of_ints, num_of_bits, sizes):
		lrange = copy.copy(range)
		di_bytes = [0] * 32
		num_of_bytes = 0
		coordinates = [0] * 3
		for item in lrange(num_of_bits,8,-8):
		#while(num_of_bits >8):
			#[bytes[num_of_bytes],stream] = self.decodebits(stream, 8)
			'''decodebits moved inline for speed'''
			#copy shared variables
			db_num_of_bits = 8
			mask = (1 << db_num_of_bits) - 1
			num = 0
			cnt = stream[0]
			lastbits = stream[1]
			lastbyte = stream[2]
			for item in lrange(db_num_of_bits,7,-8):
			#while(db_num_of_bits >= 8):
				lastbyte = (lastbyte << 8) | stream[cnt+3]
				#lastbyte &= ((1 << 32) - 1)
				cnt += 1
				num |= (lastbyte >> lastbits) << (db_num_of_bits - 8)
				#db_num_of_bits -= 8
			
			db_num_of_bits -= db_num_of_bits/8 * 8
			if (db_num_of_bits > 0):
				if(lastbits < db_num_of_bits):
					lastbits += 8
					lastbyte = (lastbyte << 8) | stream[cnt+3]
					#lastbyte &= ((1 << 32) - 1)
					cnt+=1
				lastbits -= db_num_of_bits
				num |= (lastbyte >> lastbits) & ((1 << db_num_of_bits) -1)
			lastbyte &= ((1 << 32) - 1)
			num &= mask
			stream[0] = cnt
			stream[1] = lastbits
			stream[2] = lastbyte
			di_bytes[num_of_bytes] = num
			
			num_of_bytes += 1
			#num_of_bits -= 8
		num_of_bits -= num_of_bytes*8
		if (num_of_bits > 0):
			#[bytes[num_of_bytes],stream] = self.decodebits(stream, num_of_bits)

			'''decodebits moved inline for speed'''
			#copy shared variables
			
			db_num_of_bits = num_of_bits
			mask = (1 << db_num_of_bits) - 1
			num = 0
			cnt = stream[0]
			lastbits = stream[1]
			lastbyte = stream[2]
			for item in lrange(db_num_of_bits,7,-8):
			#while(db_num_of_bits >= 8):
				lastbyte = (lastbyte << 8) | stream[cnt+3]
				#lastbyte &= ((1 << 32) - 1)
				cnt += 1
				num |= (lastbyte >> lastbits) << (db_num_of_bits - 8)
				#db_num_of_bits -= 8
			db_num_of_bits -= db_num_of_bits/8 * 8
			if (db_num_of_bits > 0):
				if(lastbits < db_num_of_bits):
					lastbits += 8
					lastbyte = (lastbyte << 8) | stream[cnt+3]
					#lastbyte &= ((1 << 32) - 1)
					cnt+=1
				lastbits -= db_num_of_bits
				num |= (lastbyte >> lastbits) & ((1 << db_num_of_bits) -1)
			lastbyte &= ((1 << 32) - 1)
			num &= mask
			stream[0] = cnt
			stream[1] = lastbits
			stream[2] = lastbyte
			di_bytes[num_of_bytes] = num


			num_of_bytes += 1
		for i in lrange( (num_of_ints-1), 0, -1):
			num = 0
			for j in lrange( (num_of_bytes - 1), -1, -1):
				num = (num << 8) | di_bytes[j]
				p = num / sizes[i]
				di_bytes[j] = p
				num = num - p * sizes[i]
			coordinates[i] = num
		coordinates[0] = di_bytes[0] | (di_bytes[1] << 8) | (di_bytes[2] << 16) | (di_bytes[3] << 24)
		return stream,coordinates

	def next_frame(self):
		#This is for the first 52 bytes
		self.content = self.source.read(52)
		self.unpacker.reset(self.content)
		
		#first must be equal to 1995
		try:
		   flag = self.unpacker.unpack_int()
		except:
		   return 0
		if ( flag != 1995): return 0; print "Magic Number Error"
		
		#number of atoms is the next int
		self.natoms = self.unpacker.unpack_int()
		#print natoms
		
		#this is the timestep in ps
		self.step = self.unpacker.unpack_int()
		#print step
		
		#this is the frame time in ps
		self.time = self.unpacker.unpack_float()
		#print time
		
		#reading 9 floats representing box vectors
		self.box = [0.0] * 9
		for i in range(9):
			self.box[i] = self.unpacker.unpack_float()

		#Now read the coordinates
		#self.frame = self.decompress_coordinates()
		'''print "decompression"'''
		try:
			self.frame = self.decompress_coordinates()
		except:
			return 0
		
		if (len(self.frame) == 1): return 0
		
		return 1
		
	def decompress_coordinates(self):
		###############################################################
		#Compressed coordinates; no internal stuff for this so we'll do it manually...
		###############################################################
		#A tuple of magic ints
		magicints = (0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 10, 12, 16, 20, 25, 32, 40, 50, 64,80, 101, 128, 161, 203, 256, 322, 406, 512, 645, 812, 1024, 1290, 1625, 2048, 2580, 3250, 4096, 5060, 6501, 8192, 10321, 13003, 16384, 20642, 26007, 32768, 41285, 52015, 65536,82570, 104031, 131072, 165140, 208063, 262144, 330280, 416127, 524287, 660561, 832255, 1048576, 1321122, 1664510, 2097152, 2642245, 3329021, 4194304, 5284491, 6658042, 8388607, 10568983, 13316085, 16777216)
		FIRSTIDX = 9
		LASTIDX = len(magicints)
		
		frame_coordinates = []
		
		self.content = self.source.read(4)
		self.unpacker.reset(self.content)
		
		#int- number of atoms. If its not equal to natoms we have problems
		lsize = self.unpacker.unpack_int()
		if (lsize != self.natoms and lsize != 0):
			#Don't know why its good to be equal to zero, but hey
			return 0
		#Number of coordinates is 3*lsize - xyz
		size3 = lsize * 3
		#comments in xdrfile.c state that files with fewer than 3 atoms 
		##should not have compression. However the test is against size, not
		##size3 and therefore it is testing for fewer or equal to 9
		if (lsize <= 9):
			self.content = self.source.read(4*lsize)
			self.unpacker.reset(self.content)
			for i in range(size3):
				#probably unpacker.unpack_float()
				frame_coordinates.append(self.unpacker.unpack_float())
			return frame_coordinates
		
		self.content = self.source.read(36)
		self.unpacker.reset(self.content)
		
		#Ok now, lets start unpacking the compression
		#Find the precision
		precision = self.unpacker.unpack_float()
				
		#Here we work out the size of the integers
		
		minint = []
		maxint = []
		sizeint = []
		
		for i in range(3):
			minint.append(self.unpacker.unpack_int())
		for i in range(3):
			maxint.append(self.unpacker.unpack_int())
		for i in range(3):
			sizeint.append(maxint[i] - minint[i] +1)
		
		#Now determine if we need a big int size. 

		if ( (sizeint[0] > 0xffffff) or (sizeint[1] > 0xffffff) or (sizeint[2] > 0xffffff) ):
			bitsizeint = []
			bitsizeint.append(self.sizeofint(sizeint[0]))
			bitsizeint.append(self.sizeofint(sizeint[1]))
			bitsizeint.append(self.sizeofint(sizeint[2]))
			bitsize = 0 #this indicates that we should use the large sizes
		else:
			bitsize = self.sizeofints(len(sizeint), sizeint)
			
		
		smallidx = self.unpacker.unpack_int()
		if (smallidx == 0):
			return 0;

		#Not clear what happens here. Straight code translation
		tmp = smallidx + 8
		if (LASTIDX < tmp):
			maxidx = LASTIDX
		else:
			maxidx = tmp
		minidx = maxidx - 8
		tmp = smallidx - 1
		if (FIRSTIDX > tmp):
			tmp = FIRSTIDX
		smaller = magicints[tmp] / 2
		smallnum = magicints[smallidx] / 2
		sizesmall = [magicints[smallidx]] * 3
		larger = magicints[maxidx]

		#length_of_stream holds the length in bytes
		length_of_stream = self.unpacker.unpack_int()
		
		self.content = self.source.read((length_of_stream+3)//4*4)
		self.unpacker.reset(self.content)
		
		if (length_of_stream == 0):
			return 0;
		compressed_coordinates = self.unpacker.unpack_fopaque(length_of_stream)
		
		#by this point we have read all of the data from this frame
		#Now we translate this to floats
		
		#We multiply the resulting ints by the inverse of the prec to
		#find the true value	
		inv_precision = 1.0 / precision
		
		#now we add a buffer of 3 zeros to the compressed_coorinates

		#compressed_coordinates = compressed_coordinates[0] * 3 + compressed_coordinates
		#Converting to bytes is necessary to convert the results appropriately
		bytes = list(struct.unpack('B'*len(compressed_coordinates), compressed_coordinates))
		bytes = [0] * 3 + bytes
		
		lrange = copy.copy(range)
		
		#Loop over every atom to get the coordinates
		i=0
		run=0
		modifier = 0
		#while(i < self.natoms):
		for i in lrange(0,self.natoms):
			#print i
			i += modifier
			if i >= self.natoms:
				break
			this_coord = [0] * 3
			prev_coord = [0] * 3
			if (bitsize == 0):
				for j in lrange(3):
					#[num,bytes] = self.decodebits(bytes, bitsizeint[j])
					
					'''decodebits moved inline for speed'''
					#copy shared variables
					
					stream = bytes
					db_num_of_bits = bitsizeint[j]
					mask = (1 << db_num_of_bits) - 1
					num = 0
					cnt = stream[0]
					lastbits = stream[1]
					lastbyte = stream[2]
					for item in lrange(db_num_of_bits,7,-8):
					#while(db_num_of_bits >= 8):
						lastbyte = (lastbyte << 8) | stream[cnt+3]
						#lastbyte &= ((1 << 32) - 1)
						cnt += 1
						num |= (lastbyte >> lastbits) << (db_num_of_bits - 8)
						#db_num_of_bits -= 8
					db_num_of_bits -= db_num_of_bits/8 * 8
					if (db_num_of_bits > 0):
						if(lastbits < db_num_of_bits):
							lastbits += 8
							lastbyte = (lastbyte << 8) | stream[cnt+3]
							#lastbyte &= ((1 << 32) - 1)
							cnt+=1
						lastbits -= db_num_of_bits
						num |= (lastbyte >> lastbits) & ((1 << db_num_of_bits) -1)
					lastbyte &= ((1 << 32) - 1)
					num &= mask
					stream[0] = cnt
					stream[1] = lastbits
					stream[2] = lastbyte
					
					
					this_coord[j] = (num)
			else:
				#[bytes,this_coord] = self.decodeints(bytes, 3, bitsize, sizeint)
				########
				'''Inline code here def decodeints(self,stream, num_of_ints, num_of_bits, sizes)'''
				num_of_ints = 3
				num_of_bits = bitsize
				sizes = sizeint
				stream = bytes
				
				di_bytes = [0] * 32
				num_of_bytes = 0
				coordinates = [0] * 3
				for item in lrange(num_of_bits,8,-8):
				#while(num_of_bits >8):
					#[bytes[num_of_bytes],stream] = self.decodebits(stream, 8)
					#decodebits moved inline for speed
					#copy shared variables
					db_num_of_bits = 8
					mask = (1 << db_num_of_bits) - 1
					num = 0
					cnt = stream[0]
					lastbits = stream[1]
					lastbyte = stream[2]
					for item in lrange(db_num_of_bits,7,-8):
					#while(db_num_of_bits >= 8):
						lastbyte = (lastbyte << 8) | stream[cnt+3]
						#lastbyte &= ((1 << 32) - 1)
						cnt += 1
						num |= (lastbyte >> lastbits) << (db_num_of_bits - 8)
						#db_num_of_bits -= 8
					
					db_num_of_bits -= db_num_of_bits/8 * 8
					if (db_num_of_bits > 0):
						if(lastbits < db_num_of_bits):
							lastbits += 8
							lastbyte = (lastbyte << 8) | stream[cnt+3]
							#lastbyte &= ((1 << 32) - 1)
							cnt+=1
						lastbits -= db_num_of_bits
						num |= (lastbyte >> lastbits) & ((1 << db_num_of_bits) -1)
					lastbyte &= ((1 << 32) - 1)
					num &= mask
					stream[0] = cnt
					stream[1] = lastbits
					stream[2] = lastbyte
					di_bytes[num_of_bytes] = num
					
					num_of_bytes += 1
					#num_of_bits -= 8
				num_of_bits -= num_of_bytes*8
				if (num_of_bits > 0):
					#[bytes[num_of_bytes],stream] = self.decodebits(stream, num_of_bits)
		
					#decodebits moved inline for speed
					#copy shared variables
					
					db_num_of_bits = num_of_bits
					mask = (1 << db_num_of_bits) - 1
					num = 0
					cnt = stream[0]
					lastbits = stream[1]
					lastbyte = stream[2]
					for item in lrange(db_num_of_bits,7,-8):
					#while(db_num_of_bits >= 8):
						lastbyte = (lastbyte << 8) | stream[cnt+3]
						#lastbyte &= ((1 << 32) - 1)
						cnt += 1
						num |= (lastbyte >> lastbits) << (db_num_of_bits - 8)
						#db_num_of_bits -= 8
					db_num_of_bits -= db_num_of_bits/8 * 8
					if (db_num_of_bits > 0):
						if(lastbits < db_num_of_bits):
							lastbits += 8
							lastbyte = (lastbyte << 8) | stream[cnt+3]
							#lastbyte &= ((1 << 32) - 1)
							cnt+=1
						lastbits -= db_num_of_bits
						num |= (lastbyte >> lastbits) & ((1 << db_num_of_bits) -1)
					lastbyte &= ((1 << 32) - 1)
					num &= mask
					stream[0] = cnt
					stream[1] = lastbits
					stream[2] = lastbyte
					di_bytes[num_of_bytes] = num
		
		
					num_of_bytes += 1
				for di_k in lrange( (num_of_ints-1), 0, -1):
					num = 0
					for di_j in lrange( (num_of_bytes - 1), -1, -1):
						num = (num << 8) | di_bytes[di_j]
						p = num / sizes[di_k]
						di_bytes[di_j] = p
						num = num - p * sizes[di_k]
					coordinates[di_k] = num
				coordinates[0] = di_bytes[0] | (di_bytes[1] << 8) | (di_bytes[2] << 16) | (di_bytes[3] << 24)
				this_coord = coordinates

				#return stream,coordinates
				''''''
				########
			
			i +=1
			#modifier += 1
			for j in lrange(3):
				this_coord[j] += minint[j]
				prev_coord[j] = this_coord[j]
			[flag,bytes] = self.decodebits(bytes, 1)
			is_smaller = 0
			if (flag == 1):
				[run,bytes] = self.decodebits(bytes,5)
				is_smaller = run % 3
				run -= is_smaller
				is_smaller -= 1
			if (run > 0):
				for k in lrange(0,run,3):
					#[bytes,this_coord] = self.decodeints(bytes, 3 , smallidx, sizesmall)
					
					num_of_ints = 3
					num_of_bits = smallidx
					sizes = sizesmall
					stream = bytes
					
					di_bytes = [0] * 32
					num_of_bytes = 0
					coordinates = [0] * 3
					for item in lrange(num_of_bits,8,-8):
					#while(num_of_bits >8):
						#[bytes[num_of_bytes],stream] = self.decodebits(stream, 8)
						'''decodebits moved inline for speed'''
						#copy shared variables
						db_num_of_bits = 8
						mask = (1 << db_num_of_bits) - 1
						num = 0
						cnt = stream[0]
						lastbits = stream[1]
						lastbyte = stream[2]
						for item in lrange(db_num_of_bits,7,-8):
						#while(db_num_of_bits >= 8):
							lastbyte = (lastbyte << 8) | stream[cnt+3]
							#lastbyte &= ((1 << 32) - 1)
							cnt += 1
							num |= (lastbyte >> lastbits) << (db_num_of_bits - 8)
							#db_num_of_bits -= 8
						
						db_num_of_bits -= db_num_of_bits/8 * 8
						if (db_num_of_bits > 0):
							if(lastbits < db_num_of_bits):
								lastbits += 8
								lastbyte = (lastbyte << 8) | stream[cnt+3]
								#lastbyte &= ((1 << 32) - 1)
								cnt+=1
							lastbits -= db_num_of_bits
							num |= (lastbyte >> lastbits) & ((1 << db_num_of_bits) -1)
						lastbyte &= ((1 << 32) - 1)
						num &= mask
						stream[0] = cnt
						stream[1] = lastbits
						stream[2] = lastbyte
						di_bytes[num_of_bytes] = num
						
						num_of_bytes += 1
						#num_of_bits -= 8
					num_of_bits -= num_of_bytes*8
					if (num_of_bits > 0):
						#[bytes[num_of_bytes],stream] = self.decodebits(stream, num_of_bits)
			
						'''decodebits moved inline for speed'''
						#copy shared variables
						
						db_num_of_bits = num_of_bits
						mask = (1 << db_num_of_bits) - 1
						num = 0
						cnt = stream[0]
						lastbits = stream[1]
						lastbyte = stream[2]
						for item in lrange(db_num_of_bits,7,-8):
						#while(db_num_of_bits >= 8):
							lastbyte = (lastbyte << 8) | stream[cnt+3]
							#lastbyte &= ((1 << 32) - 1)
							cnt += 1
							num |= (lastbyte >> lastbits) << (db_num_of_bits - 8)
							#db_num_of_bits -= 8
						db_num_of_bits -= db_num_of_bits/8 * 8
						if (db_num_of_bits > 0):
							if(lastbits < db_num_of_bits):
								lastbits += 8
								lastbyte = (lastbyte << 8) | stream[cnt+3]
								#lastbyte &= ((1 << 32) - 1)
								cnt+=1
							lastbits -= db_num_of_bits
							num |= (lastbyte >> lastbits) & ((1 << db_num_of_bits) -1)
						lastbyte &= ((1 << 32) - 1)
						num &= mask
						stream[0] = cnt
						stream[1] = lastbits
						stream[2] = lastbyte
						di_bytes[num_of_bytes] = num
			
			
						num_of_bytes += 1
					for di_k in lrange( (num_of_ints-1), 0, -1):
						num = 0
						for di_j in lrange( (num_of_bytes - 1), -1, -1):
							num = (num << 8) | di_bytes[di_j]
							p = num / sizes[di_k]
							di_bytes[di_j] = p
							num = num - p * sizes[di_k]
						coordinates[di_k] = num
					coordinates[0] = di_bytes[0] | (di_bytes[1] << 8) | (di_bytes[2] << 16) | (di_bytes[3] << 24)
					this_coord = coordinates					
					
					i += 1
					modifier += 1
					for l in lrange(3):
						this_coord[l] += prev_coord[l] - smallnum
					if (k==0):
						#interchange first with second atom for better compression of water molecules
						for l in lrange(3):
							tmp = this_coord[l]
							this_coord[l] = prev_coord[l]
							prev_coord[l] = tmp
							frame_coordinates.append(prev_coord[l] * inv_precision)
					else:
						for l in lrange(3):
							prev_coord[l] = this_coord[l]
					for l in lrange(3):
						frame_coordinates.append(this_coord[l] * inv_precision)
			else:
				for k in lrange(3):
					frame_coordinates.append(this_coord[k] * inv_precision)
			smallidx += is_smaller
			if (is_smaller < 0):
				smallnum = smaller
				if (smallidx > FIRSTIDX):
					smaller = magicints[smallidx - 1] /2
				else:
					smaller = 0
			elif (is_smaller > 0):
				smaller = smallnum
				smallnum = magicints[smallidx] / 2
			sizesmall = [magicints[smallidx]] * 3
		
		return frame_coordinates
	def convert_coordinates(self):
		self.cartesian =  [ [self.frame[i*3], self.frame[i*3+1], self.frame[i*3+2]] for i in range(len(self.frame)/3) ]
	
	def next(self):
		correctly_read = self.next_frame()
		if not correctly_read: raise StopIteration
		self.convert_coordinates()
		return self.cartesian
	
	def __iter__(self):
		return self
	
