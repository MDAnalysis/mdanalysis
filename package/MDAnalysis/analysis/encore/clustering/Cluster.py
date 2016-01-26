# Cluster.py --- classes to handle results of clustering runs
# Copyright (C) 2014 Wouter Boomsma, Matteo Tiberti
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
Ensemble representation --- :mod:`MDAnalysis.analysis.ensemble.ensemble`
=====================================================================

The module contains the Cluster and ClusterCollection classes which are designed
to store results from clustering algorithms. 
"""

import numpy as np

class Cluster:
    """
    Generic Cluster class for clusters with centroids.

**Attributes**:

   `id` : int
	Cluster ID number. Useful for the ClustersCollection class

   `metadata` : iterable
	dict of lists, containing metadata for the cluster elements. The iterable must return the same number of elements as those that belong to the cluster.

   `size` : int
	number of elements.

   `centroid` : element object
	cluster centroid.

   `elements` : numpy.array
	array containing the cluster elements.
   """

    def __init__(self, elem_list=None, centroid=None, idn=None, metadata=None):
	"""Class constructor. If elem_list is None, an empty cluster is created and the remaning arguments ignored.

	**Arguments:**

	`elem_list` : numpy.array or None
	    numpy array of cluster elements. if None, the cluster will be initialized as empty.
	
	`centroid` : None or element object
	    centroid object

        `idn` : int
	    cluster ID

	`metadata` : {str:iterable, ...}
	    metadata, one value for each cluster element. The iterable must have the same length as the elements array.

"""
	

        self.id = idn
        
        if elem_list == None:
            self.size = 0
            self.elements = np.array([])
            self.centroid = None
            self.metadata = {}
            return
        
        self.metadata = {}
        self.elements = elem_list
        if not centroid in self.elements:
            raise LookupError
        
        self.centroid = centroid
        self.size = self.elements.shape[0]
        if metadata:
            for k,v in metadata.iteritems():
                if len(v) != self.size:
                    raise TypeError
            self.metadata[k] = np.array(v)
        
    def __iter__(self):
        return iter(self.elements)
        
    def add_metadata(self, name, data):
        if len(data) != self.size:
            raise TypeError
        self.metadata[name] = np.array(data)
        
class ClustersCollection():
    """Clusters collection class; this class represents the results of a full clustering run. It stores a group of clusters defined as encore.clustering.Cluster objects.

    **Attributes:**
    
    `clusters` : list of Cluster objects
	clusters object which are part of the Cluster collection

"""

       
    def __init__(self, elements=None, metadata=None):
	"""Class constructor. If elements is None, an empty cluster collection will be created. Otherwise, the constructor takes as input an iterable of ints with the following format:

	[ a, a, a, a, b, b, b, c, c, ... , z, z ]

	the variables a,b,c,...,z are clusters centroids, represented as cluster element numbers (i.e. 3 means the 4th element of the ordered input data for clustering).
 The array has the same number of elements as input data. The array maps a correspondence between cluster elements (which are implicitly associated with the position in the array) with centroids, i. e. defines clusters. For instance:

	[ 1, 1, 1, 4, 4, 5 ]
	
	means that elements 0, 1, 2 form a cluster which has 1 as centroids, elements 3 and 4 form a cluster which has 4 as centroid, and element 5 has its own cluster.

**arguments**:

`elements` : iterable of ints or None
	clustering results. See the previous description for details

`metadata` : {str:list, str:list,...} or None
	metadata for the data elements. The list must be of the same size as the elements array, with one value per element.
	
"""
        idn = 0
        if elements == None:
            self.clusters=None
            return
        
        if not len(set(map(type,elements))) == 1:
            raise TypeError
        self.clusters = []            
        elements_array = np.array(elements)
        centroids = np.unique(elements_array)
        for i in centroids:
            if elements[i] != i:
                raise AssertionError
        for c in centroids:
            this_metadata = {}
            this_array = np.where(elements_array == c)
            if metadata:
                for k,v in metadata.iteritems():
                    this_metadata[k] = np.array(v)[this_array]
            self.clusters.append(Cluster(elem_list=this_array[0], idn=idn, centroid=c, metadata=this_metadata))
            
            idn += 1
            
    def get_ids(self):
	"""Get the ID numbers of the clusters

	**Returns:**

	`ids` : list of int 
	    list of cluster ids """
        return [v.idn for v in self.clusters]
        
    def get_centroids(self):
        """Get the centroids of the clusters

        **Returns:**

        `centroids` : list of cluster element objects
            list of cluster centroids """

        return [v.centroid for v in self.clusters]

    def __setitiem__(self, name, val):
        if type(val) != Cluster:
            raise TypeError
        self.clusters[name] = val
        
    def __getitem__(self, name):
        return self.clusters[name]
        
    def __iter__(self):
        return iter(self.clusters)

