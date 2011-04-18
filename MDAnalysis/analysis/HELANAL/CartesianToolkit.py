#!/usr/bin/python

import math

def center(list_of_coordinates):
	transposed_list = zip(*list_of_coordinates)
	return [sum(cart)/float(len(cart)) for  cart in transposed_list]
	
def vecsub(cart_list_a, cart_list_b):
	"""a - b"""
	zipped_carts = zip(cart_list_a, cart_list_b)
	return [(cart[0] - cart[1]) for cart in zipped_carts]

def vecadd(cart_list_a, cart_list_b):
	"""a + b"""
	zipped_carts = zip(cart_list_a, cart_list_b)
	return [(cart[0] + cart[1]) for cart in zipped_carts]

def veclength(cart_list):
	"""|a|"""
	squares = [cart ** 2 for cart in cart_list]
	return float(sum(squares))**0.5
	
def vecnorm(cart_list):
	"""a/|a|"""
	length = veclength(cart_list)
	return [cart/length for cart in cart_list]
	
def vecscaler(cart_list_a,cart_list_b):
	"""a.b"""
	return sum([cart[0] * cart[1] for cart in zip(cart_list_a,cart_list_b)])
	
def vecscale(cart_list,value):
	"""value * a"""
	return [(item * value) for item in cart_list]

def vecangle(cart_list_a,cart_list_b):
	"""angle(a,b)"""
	length_product = veclength(cart_list_a) * veclength(cart_list_b)
	return math.acos(round(vecscaler(cart_list_a,cart_list_b)/length_product,8))

def vecdist(cart_list_a,cart_list_b):
	"""|a-b|"""
	relating_vector = vecsub(cart_list_a, cart_list_b)
	return veclength(relating_vector)
	
def veccross(cart_list_a,cart_list_b):
	return [ ( cart_list_a[1]*cart_list_b[2] - cart_list_a[2]*cart_list_b[1] ), ( cart_list_a[2]*cart_list_b[0] - cart_list_a[0]*cart_list_b[2] ), ( cart_list_a[0]*cart_list_b[1] - cart_list_a[1]*cart_list_b[0] ) ]
	
def rotation_matrix_to_x(cart_list):
	angle = vecangle(cart_list,[1,0,0])
	[x,y,z] = vecnorm(veccross(cart_list,[1,0,0]))
	rotation_matrix = [	[ (1 + (1-math.cos(angle))*(x*x-1)),	(-z*math.sin(angle)+(1-math.cos(angle))*x*y),	(y*math.sin(angle)+(1-math.cos(angle))*x*z) ],
						[ (z*math.sin(angle)+(1-math.cos(angle))*x*y), 	(1 + (1-math.cos(angle))*(y*y-1)), 	(-x*math.sin(angle)+(1-math.cos(angle))*y*z) ],
						[ (-y*math.sin(angle)+(1-math.cos(angle))*x*z),	(x*math.sin(angle)+(1-math.cos(angle))*y*z),	(1 + (1-math.cos(angle))*(z*z-1)) ] ]
	return rotation_matrix

def rotate_by_matrix(cart_list,rotation_matrix):
	rotated_vector = [	rotation_matrix[0][0]*cart_list[0] + rotation_matrix[0][1]*cart_list[1] + rotation_matrix[0][2]*cart_list[2],
						rotation_matrix[1][0]*cart_list[0] + rotation_matrix[1][1]*cart_list[1] + rotation_matrix[1][2]*cart_list[2],
						rotation_matrix[2][0]*cart_list[0] + rotation_matrix[2][1]*cart_list[1] + rotation_matrix[2][2]*cart_list[2] ]
	return rotated_vector

def transform_by_rot_trans_matrix(cart_list,rotation_matrix):
	rotated_vector = [	rotation_matrix[0][0]*cart_list[0] + rotation_matrix[0][1]*cart_list[1] + rotation_matrix[0][2]*cart_list[2] + rotation_matrix[0][3], 
						rotation_matrix[1][0]*cart_list[0] + rotation_matrix[1][1]*cart_list[1] + rotation_matrix[1][2]*cart_list[2] + rotation_matrix[1][3],
						rotation_matrix[2][0]*cart_list[0] + rotation_matrix[2][1]*cart_list[1] + rotation_matrix[2][2]*cart_list[2] + rotation_matrix[2][3] ]
	return rotated_vector

def vecdot(cart_list_a,cart_list_b):
	return vecscaler(cart_list_a,cart_list_b)

def wrapangle(angle):
	if angle > math.pi:
		angle -= 2*math.pi
	elif angle < -math.pi:
		angle += 2*math.pi
	return angle

def rotation_matrix_relation(cart_list_a,cart_list_b):
	pass

##### Numpy versions #####

import numpy

def numpy_center(array_of_coordinates):
	return [sum(cart[:,0])/float(len(cart)) , sum(cart[:,1])/float(len(cart)), sum(cart[:,2])/float(len(cart)) ]
