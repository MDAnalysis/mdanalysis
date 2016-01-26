# similarity.py --- Simularity measures between protein ensembles
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
Ensemble similarity calculations --- :mod:`encore.similarity`
=====================================================================

The module contains implementations of similary measures between
protein ensembles described in:

     Similarity Measures for Protein Ensembles. Lindorff-Larsen, K.; 
     Ferkinghoff-Borg, J. PLoS ONE 2009, 4, e4203.

"""
import optparse
import numpy
import warnings
import logging
from time import sleep
from MDAnalysis import Universe
from .Ensemble import Ensemble
from .clustering.Cluster import ClustersCollection
from .clustering.affinityprop import AffinityPropagation
from .dimensionality_reduction.stochasticproxembed import StochasticProximityEmbedding, kNNStochasticProximityEmbedding
from .confdistmatrix import MinusRMSDMatrixGenerator, RMSDMatrixGenerator
from .covariance import covariance_matrix, EstimatorShrinkage, EstimatorML
from multiprocessing import cpu_count
from .utils import *
from scipy.stats import gaussian_kde
from random import randint

# Silence deprecation warnings - scipy problem
warnings.filterwarnings("ignore", category=DeprecationWarning) 
warnings.filterwarnings("ignore", category=RuntimeWarning) 
warnings.filterwarnings("ignore", category=FutureWarning) 



# Low boundary value for log() argument - ensure no nans 
EPSILON=1E-15

# x*log(y) with the assumption that 0*(log(0)) = 0
xlogy = numpy.vectorize(lambda x,y : 0.0 if (x<=EPSILON and y<=EPSILON) else x*numpy.log(y))     

# discrete dKL
def discrete_kullback_leibler_divergence(pA, pB):
    """Kullback-Leibler divergence between discrete probability distribution. Notice that since this measure is not symmetric  :math:`d_{KL}(p_A,p_B) != d_{KL}(p_B,p_A)`

    **Arguments:**
	
	`pA` : iterable of floats
		First discrete probability density function
	
	`pB` : iterable of floats
		Second discrete probability density function

    **Returns:**
	
	`dkl` : float
		Discrete Kullback-Liebler divergence
	"""

    return numpy.sum( xlogy(pA, pA/pB) )

# discrete dJS
def discrete_jensen_shannon_divergence(pA, pB):
    """Jensen-Shannon divergence between discrete probability distributions.

    **Arguments:**
        
	`pA` : iterable of floats
                First discrete probability density function
        
	`pB` : iterable of floats
                Second discrete probability density function

    **Returns:**
        
	`djs` : float
                Discrete Jensen-Shannon divergence
"""
    return 0.5*( discrete_kullback_leibler_divergence(pA, (pA+pB)*0.5) + 
                 discrete_kullback_leibler_divergence(pB, (pA+pB)*0.5) )

# calculate harmonic similarity
def harmonic_ensemble_similarity(ensemble1=None,
                                 ensemble2=None,
                                 sigma1=None,
                                 sigma2=None,
                                 x1=None,
                                 x2=None,
                                 mass_weighted=True,
                                 covariance_estimator = EstimatorShrinkage()):
    ''' 
    Calculate the harmonic ensemble similarity measure
    as defined in 

	    Similarity Measures for Protein Ensembles. Lindorff-Larsen, K.; 
	    Ferkinghoff-Borg, J. PLoS ONE 2009, 4, e4203.

    **Arguments:**

	`ensemble1` : encore.Ensemble or None
		First ensemble to be compared. If this is None, sigma1 and x1 must be provided.

	`ensemble2` : encore.Ensemble or None
		Second ensemble to be compared. If this is None, sigma2 and x2 must be provided.

	`sigma1` : numpy.array
		Covariance matrix for the first ensemble. If this None, calculate it from ensemble1 using covariance_estimator

	`sigma2` : numpy.array
		Covariance matrix for the second ensemble. If this None, calculate it from ensemble1 using covariance_estimator

	`x1`: numpy.array 
		Mean for the estimated normal multivariate distribution of the first ensemble. If this is None, calculate it from ensemble1

	`x2`: numpy.array
                Mean for the estimated normal multivariate distribution of the first ensemble.. If this is None, calculate it from ensemble2

	`mass_weighted` : bool
		Whether to perform mass-weighted covariance matrix estimation

	`covariance_estimator` : either EstimatorShrinkage or EstimatorML objects
		Which covariance estimator to use
	
    **Returns:**
	
	`dhes` : float
		harmonic similarity measure
    '''

    # If matrices and means are specified, use them
    if x1 == None or x2 == None or sigma1 == None or sigma2 == None:
        if ensemble1 == None or ensemble2 == None:
            raise RuntimeError

        # Extract coordinates from ensembles
        coordinates_system1 = ensemble1.coordinates
        coordinates_system2 = ensemble2.coordinates
    
        # Average coordinates in the two systems
        x1 = numpy.average(coordinates_system1, axis=0).flatten()
        x2 = numpy.average(coordinates_system2, axis=0).flatten()

        # Covariance matrices in the two systems
        sigma1 = covariance_matrix(ensemble1, 
                               mass_weighted=mass_weighted,
                               estimator = covariance_estimator)
        sigma2 = covariance_matrix(ensemble2, 
                               mass_weighted=mass_weighted,
                               estimator = covariance_estimator)

    # Inverse covariance matrices
    sigma1_inv = numpy.linalg.pinv(sigma1)
    sigma2_inv = numpy.linalg.pinv(sigma2)

    # Difference between average vectors
    d_avg = x1 - x2

    # Sigma
    sigma = sigma1_inv + sigma2_inv

    # Distance measure
    trace = numpy.trace(numpy.dot(sigma1, sigma2_inv) + 
                        numpy.dot(sigma2, sigma1_inv)
                        - 2*numpy.identity(sigma1.shape[0]))

    d_hes = 0.25*(numpy.dot(numpy.transpose(d_avg), 
                            numpy.dot(sigma1_inv + sigma2_inv,
                                      d_avg)) + trace)
    return d_hes

def clustering_ensemble_similarity(cc, ens1, ens1_id, ens2, ens2_id):
    """Clustering ensemble similarity: calculate the probability densities from the clusters and calculate discrete Jensen-Shannon divergence.
	
	**Arguments:**

	`cc` : encore.ClustersCollection 
		Collection from cluster calculated by a clustering algorithm (e.g. Affinity propagation)
	
	`ens1` : encore.Ensemble
		First ensemble to be used in comparison
	
        `ens2` : encore.Ensemble
                Second ensemble to be used in comparison
		
	`ens1_id` : int
		First ensemble id as detailed in the ClustersCollection metadata

	`ens2_id` : int
		Second ensemble id as detailed in the ClustersCollection metadata

	**Returns:**

	`djs` : float
		Jensen-Shannon divergence between the two ensembles, as calculated by the clustering ensemble similarity method
	"""
    tmpA = numpy.array( [ numpy.where(c.metadata['ensemble'] == ens1_id)[0].shape[0]/float(ens1.coordinates.shape[0]) for c in cc ] )
    tmpB = numpy.array( [ numpy.where(c.metadata['ensemble'] == ens2_id)[0].shape[0]/float(ens2.coordinates.shape[0]) for c in cc ] )
                    
    # Exclude clusters which have 0 elements in both ensembles    
    pA=tmpA[tmpA+tmpB > EPSILON]
    pB=tmpB[tmpA+tmpB > EPSILON]

    return discrete_jensen_shannon_divergence(pA, pB)

def cumulative_clustering_ensemble_similarity(cc, ens1, ens1_id, ens2, ens2_id, ens1_id_min=1, ens2_id_min=1):
    """ Calculate clustering ensemble similarity between joined ensembles. This means that, after clustering has been performed, some ensembles are merged and the dJS is calculated between the probability distributions of the two clusters groups. In particular, the two ensemble groups are defined by their ensembles id: one of the two joined ensembles will comprise all the ensembles with id [ens1_id_min, ens1_id], and the other ensembles will comprise all the ensembles with id [ens2_id_min, ens2_id].

**Arguments:**

        `cc` : encore.ClustersCollection
                Collection from cluster calculated by a clustering algorithm (e.g. Affinity propagation)

        `ens1` : encore.Ensemble
                First ensemble to be used in comparison

        `ens2` : encore.Ensemble
                Second ensemble to be used in comparison

        `ens1_id` : int
                First ensemble id as detailed in the ClustersCollection metadata

        `ens2_id` : int
                Second ensemble id as detailed in the ClustersCollection metadata

        **Returns:**

        `djs` : float
                Jensen-Shannon divergence between the two ensembles, as calculated by the clustering ensemble similarity method

"""

    ensA = [ numpy.where( numpy.logical_and(c.metadata['ensemble'] <= ens1_id, c.metadata['ensemble']) >= ens1_id_min)[0].shape[0] for c in cc ]
    ensB = [ numpy.where( numpy.logical_and(c.metadata['ensemble'] <= ens2_id, c.metadata['ensemble']) >= ens2_id_min)[0].shape[0] for c in cc ]
    sizeA = float(numpy.sum(ensA))
    sizeB = float(numpy.sum(ensB))
    #sizeA = float( numpy.sum( [numpy.where( numpy.logical_and(c.metadata['ensemble'] <= ens1_id, c.metadata['ensemble']) >= ens1_id_min)[0].shape[0] for c in cc])
    #sizeB = float(numpy.sum( [numpy.where( numpy.logical_and(c.metadata['ensemble'] <= ens2_id, c.metadata['ensemble']) >= ens2_id_min)[0].shape[0] for c in cc])

    tmpA = numpy.array( ensA )/sizeA
    tmpB = numpy.array( ensB  )/sizeB

    # Exclude clusters which have 0 elements in both ensembles
    pA=tmpA[tmpA+tmpB > EPSILON]
    pB=tmpB[tmpA+tmpB > EPSILON]

    return discrete_jensen_shannon_divergence(pA, pB)

def gen_kde_pdfs(embedded_space, ensemble_assignment, nensembles,  nsamples=None, **kwargs):
    """ 
    Generate Kernel Density Estimates (KDE) from embedded spaces and elaborate the coordinates for later use.

**Arguments:**

`embedded_space` : numpy.array
	Array containing the coordinates of the embedded space

`ensemble_assignment` : numpy.array
	Array containing one int per ensemble conformation. These allow to distinguish, in the complete embedded space, which conformations belong to each ensemble. For instance if ensemble_assignment is [1,1,1,1,2,2], it means that the first four conformations belong to ensemble 1 and the last two to ensemble 2

`nesensembles` : int
	Number of ensembles

`nsamples` : int samples to be drawn from the ensembles. Will be required in a later stage in order to calculate dJS.`

**Returns:**

`kdes` : scipy.stats.gaussian_kde
	KDEs calculated from ensembles

`resamples` : list of numpy.array
	For each KDE, draw samples according to the probability distribution of the KDE mixture model

`embedded_ensembles` : list of numpy.array
	List of numpy.array containing, each one, the elements of the embedded space belonging to a certain ensemble
"""
    kdes = []
    embedded_ensembles = []
    resamples = []
    
    for i in range(1,nensembles+1):
        this_embedded = embedded_space.transpose()[numpy.where(ensemble_assignment == i)].transpose()
        embedded_ensembles.append(this_embedded)
        kdes.append(gaussian_kde(this_embedded)) # XXX support different bandwidth values

    # Set number of samples
    if not nsamples:
        nsamples = this_embedded.shape[1]*10

    # Resample according to probability distributions
    for this_kde in kdes:
        resamples.append(this_kde.resample(nsamples))

    return (kdes, resamples, embedded_ensembles)
    
def dimred_ensemble_similarity(kde1, resamples1, kde2, resamples2, ln_P1_exp_P1=None, ln_P2_exp_P2=None, ln_P1P2_exp_P1=None, ln_P1P2_exp_P2=None):
    """ Calculate the Jensen-Shannon divergence according the the Dimensionality reduction method. In this case, we have continuous probability densities we have to integrate over the measureable space. Our target is calculating Kullback-Liebler, which is defined as:

.. math::
	D_{KL}(P(x) || Q(x)) = \\int_{-\\infty}^{\\infty}P(x_i) ln(P(x_i)/Q(x_i)) = \\langle{}ln(P(x))\\rangle{}_P - \\langle{}ln(Q(x))\\rangle{}_P

where the :math:`\\langle{}.\\rangle{}_P` denotes an expectation calculated under the 
distribution P. We can, thus, just estimate the expectation values of the components to get an estimate of dKL.
Since the Jensen-Shannon distance is actually  more complex, we need to estimate four expectation values:

.. math::	
     \\langle{}log(P(x))\\rangle{}_P

     \\langle{}log(Q(x))\\rangle{}_Q

     \\langle{}log(0.5*(P(x)+Q(x)))\\rangle{}_P

     \\langle{}log(0.5*(P(x)+Q(x)))\\rangle{}_Q

**Arguments:**

`kde1` : scipy.stats.gaussian_kde
	Kernel density estimation for ensemble 1

`resamples1` : numpy.array
	Samples drawn according do kde1. Will be used as samples to calculate the expected values according to 'P' as detailed before.

`kde2` : scipy.stats.gaussian_kde
        Kernel density estimation for ensemble 2

`resamples2` : numpy.array
        Samples drawn according do kde2. Will be used as sample to calculate the expected values according to 'Q' as detailed before.	

`ln_P1_exp_P1` : float or None
	Use this value for :math:`\\langle{}log(P(x))\\rangle{}_P`; if None, calculate it instead

`ln_P2_exp_P2` : float or None
        Use this value for :math:`\\langle{}log(Q(x))\\rangle{}_Q`; if None, calculate it instead

`ln_P1P2_exp_P1` : float or None
        Use this value for :math:`\\langle{}log(0.5*(P(x)+Q(x)))\\rangle{}_P`;  if None, calculate it instead

`ln_P1P2_exp_P1` : float or None
        Use this value for :math:`\\langle{}log(0.5*(P(x)+Q(x)))\\rangle{}_Q`; if None, calculate it instead	

**Returns:**

`djs` : float
	Jensen-Shannon divergence calculated according to the dimensionality reduction method 
"""

    if not ln_P1_exp_P1 and not ln_P2_exp_P2 and not ln_P1P2_exp_P1 and not ln_P1P2_exp_P2:
        ln_P1_exp_P1 = numpy.average(numpy.log(kde1.evaluate(resamples1)))
        ln_P2_exp_P2 = numpy.average(numpy.log(kde2.evaluate(resamples2)))
        ln_P1P2_exp_P1 = numpy.average(numpy.log(0.5*(kde1.evaluate(resamples1)+kde2.evaluate(resamples1))))
        ln_P1P2_exp_P2 = numpy.average(numpy.log(0.5*(kde1.evaluate(resamples2)+kde2.evaluate(resamples2))))

    return 0.5 * (ln_P1_exp_P1 - ln_P1P2_exp_P1 + ln_P2_exp_P2 - ln_P1P2_exp_P2)

def cumulative_gen_kde_pdfs(embedded_space, ensemble_assignment, nensembles,  nsamples=None, ens_id_min=1, ens_id_max=None):
    """
    Generate Kernel Density Estimates (KDE) from embedded spaces and elaborate the coordinates for later use. However, consider more than one ensemble as the space on which the KDE will be generated. In particular, will use ensembles with ID [ens_id_min, ens_id_max]. 

**Arguments:**

`embedded_space` : numpy.array
        Array containing the coordinates of the embedded space

`ensemble_assignment` : numpy.array
        array containing one int per ensemble conformation. These allow to distinguish, in the complete embedded space, which conformations belong to each ensemble. For instance if ensemble_assignment is [1,1,1,1,2,2], it means that the first four conformations belong to ensemble 1 and the last two to ensemble 2

`nesensembles` : int
        Number of ensembles

`nsamples : int 
	Samples to be drawn from the ensembles. Will be required in a later stage in order to calculate dJS.`

`ens_id_min` : int 
	Minimum ID of the ensemble to be considered; see description

`ens_id_max` : int
	Maximum ID of the ensemble to be considered; see description

**Returns:**

`kdes` : scipy.stats.gaussian_kde
        KDEs calculated from ensembles

`resamples` : list of numpy.array
        For each KDE, draw samples according to the probability distribution of the kde mixture model

`embedded_ensembles` : list of numpy.array
        List of numpy.array containing, each one, the elements of the embedded space belonging to a certain ensemble
    """

    kdes = []
    embedded_ensembles = []
    resamples = []
    if not ens_id_max:
        ens_id_max = nensembles+1
    for i in range(ens_id_min, ens_id_max+1):
        this_embedded = embedded_space.transpose()[numpy.where(numpy.logical_and(ensemble_assignment >= ens_id_min, ensemble_assignment <= i))].transpose()
        embedded_ensembles.append(this_embedded)
        kdes.append(gaussian_kde(this_embedded)) # XXX support different bandwidth values

    # Set number of samples
    if not nsamples:
        nsamples = this_embedded.shape[1]*10

    # Resample according to probability distributions
    for this_kde in kdes:
        resamples.append(this_kde.resample(nsamples))

    return (kdes, resamples, embedded_ensembles)

def write_output(matrix, base_fname=None, header="", suffix="", extension="dat"):
    """
    Write output matrix with a nice format, to stdout and optionally a file.  

**Arguments:**

`matrix` : encore.utils.TriangularMatrix
        Matrix containing the values to be printed

`base_fname` : str
	Basic filename for output. If None, no files will be written, and the matrix will be just printed on screen

`header` : str
        Line to be written just before the matrix

`suffix` : str 
	String to be concatenated to basename, in order to get the final file name        

`extension` : str 
	Extension for the output file       

    """

    if base_fname != None:
        fname = base_fname+"-"+suffix+"."+extension
    else:
        fname = None
    matrix.square_print(header=header, fname=fname)
        
def write_output_line(value, fhandler=None, suffix="", label="win.", number=0, rawline=None):
    """
    Write a line of data with a fixed format to standard output and optionally file. The line will be appended or written to a file object.
The format is (in the Python str.format specification language): '{:s}{:d}\t{:.3f}', with the first element being the label, the second being
a number that identifies the data point, and the third being the number itself. For instance:

win.3	0.278

**Arguments:**

`value` : float
        Value to be printed.

`fhandler` : file object
	File object in which the line will be written. if None, nothing will be written to file, and the value will be just printed on screen

`label` : str
        Label to be written before the data 

`number` : int 
	Number that identifies the data being written in this line.        

`rawline` : str
	If rawline is not None, write rawline to fhandler instead of the formatted number line. rawline can be any arbitrary string.       
    """

    if fhandler == None:
        fh = Tee(sys.stdout)
    else:
        fh = Tee(sys.stdout, fhandler)

    if rawline != None:
        print >>fh, rawline
        return

    print >>fh, "{:s}{:d}\t{:.3f}".format(label, number, value)

def bootstrap_coordinates(coords, times):
    """
    Bootstrap conformations in a encore.Ensemble. This means drawing from the encore.Ensemble.coordinates numpy array with replacement "times" times and returning the outcome. 

**Arguments:**

`coords` : numpy.array
        3-dimensional coordinates array

`times` : int
        number of times the coordinates will be bootstrapped

**Returns:**

`out` : list
        Bootstrapped coordinates list. len(out) = times.
    """
    out = []
    for t in range(times):
        this_coords = numpy.zeros(coords.shape)
        for c in range(this_coords.shape[0]):
            this_coords[c,:,:] = coords[numpy.random.randint(low=0, high=this_coords.shape[0]),:,:]
        out.append(this_coords)
    return out

def bootstrap_matrix(matrix):
    """
    Bootstrap an input square matrix. The resulting matrix will have the same shape as the original one, but the order of its elements will be drawn (with repetition). Separately bootstraps each ensemble.

**Arguments:**

`matrix` : encore.utils.TriangularMatrix
        similarity/dissimilarity matrix

**Returns:**

`this_m` : encore.utils.TriangularMatrix
        bootstrapped similarity/dissimilarity matrix
    """
    ensemble_identifiers = numpy.unique(ensemble_assignment)
    this_m = TriangularMatrix(size = matrix.size)
    indexes = []
    for ens in ensemble_identifiers:
        old_indexes = numpy.where(ensemble_assignment == ens)[0]
        indexes.append( numpy.random.randint(low=numpy.min(old_indexes), high=numpy.max(old_indexes)+1, size=old_indexes.shape[0] ) )

    indexes = numpy.hstack(indexes)
    for j in range(this_m.size):
        for k in range(j):
            this_m[j, k] = matrix[indexes[j], indexes[k]]
        
    logging.info("Matrix bootstrapped.")
    return this_m





def get_similarity_matrix(  ensembles,
                            similarity_mode="minusrmsd",
                            load_matrix = None,
                            change_sign = None,
                            save_matrix = None,
                            superimpose = True,
                            superimposition_subset = "name CA",
                            mass_weighted = True,
                            bootstrap_matrix = False,
                            bootstrapping_samples = 100,
                            np = 1):

    trajlist = []
    ensemble_assignment = []

    nensembles = len(ensembles)

    # Define ensemble assignments as required on the joined ensemble
    for i in range(1, nensembles+1):
        ensemble_assignment += [i for j in ensembles[i-1].coordinates]
    ensemble_assignment = numpy.array(ensemble_assignment)

    # Joined ensemble
    joined_ensemble = Ensemble(topology = ensembles[0].topology_filename,
                               trajectory = [ensembles[0].topology_filename],
                               atom_selection_string = "all",
                               superimposition_selection_string = ensembles[0].superimposition_selection_string)

    # Joined ensemble coordinates as a concatenation of single ensembles - faster this way
    joined_ensemble.coordinates = numpy.concatenate(tuple([ e.coordinates for e in ensembles ]) )
    joined_ensemble.superimposition_coordinates = numpy.concatenate(tuple([ e.superimposition_coordinates for e in ensembles ]) )
   
    # Define metadata dictionary
    metadata = {'ensemble': ensemble_assignment}
    
    # Choose distance metric
    if similarity_mode == "minusrmsd":
        logging.info("    Similarity matrix: -RMSD matrix")
        matrix_builder = MinusRMSDMatrixGenerator()
    elif similarity_mode == "rmsd":
        logging.info("    Similarity matrix: RMSD matrix")
        matrix_builder = RMSDMatrixGenerator()
    else: 
        logging.error("Supported conformational distance measures are rmsd and minusrmsd")
        return None

    # Load the matrix if required
    if load_matrix:
        logging.info("        Loading similarity matrix from: %s"%load_matrix)
        confdistmatrix = TriangularMatrix(size=joined_ensemble.coordinates.shape[0], loadfile=load_matrix)
        logging.info("        Done!")
        for key in confdistmatrix.metadata.dtype.names:
            logging.info("        %s : %s" % (key, str(confdistmatrix.metadata[key][0])) )

        # Change matrix sign if required. Useful to switch between similarity/distance matrix.
        if change_sign:
            logging.info("        The matrix sign will be changed.")
            for k,v in enumerate(confdistmatrix._elements):
                confdistmatrix._elements[k] = -v

        # Check matrix size for consistency
        if not confdistmatrix.size == joined_ensemble.coordinates.shape[0]:
            logging.error("ERROR: The size of the loaded matrix and of the ensemble do not match")
            return None

    # Calculate the matrix  
    else:
        logging.info("        Perform pairwise alignment: %s"       % str(superimpose))
        logging.info("        Mass-weighted alignment and RMSD: %s" % str(mass_weighted))
        if superimpose:
            logging.info("        Atoms subset for alignment: %s" % superimposition_subset ) 
        logging.info("    Calculating similarity matrix . . .")

        # Use superimposition subset, if necessary. If the pairwise alignment is not required, it will not be performed anyway.
        if superimposition_subset:
            confdistmatrix = matrix_builder(joined_ensemble, 
                                pairwise_align = superimpose, 
                                align_subset_coordinates = joined_ensemble.superimposition_coordinates,
                                mass_weighted = mass_weighted,
                                ncores = np)

        else:
            confdistmatrix = matrix_builder(joined_ensemble, 
                                pairwise_align = superimpose,
                                mass_weighted = mass_weighted,
                                ncores = np)                            
        
        logging.info("    Done!")

        if save_matrix:
            logging.info("    Similarity matrix will be saved in %s.%s"%(parser_phase3.options.save_matrix, "" if parser_phase3.options.save_matrix[-3:] == "npz" else "npz"))
            confdistmatrix.savez(parser_phase3.options.save_matrix)

    if bootstrap_matrix:
        logging.info("Error estimation mode: Bootstrapping")
        logging.info("the similarity matrix will be bootstrapped %d times." % parser_phase3.options.bootstrapping_runs)

        bs_args = [tuple([confdistmatrix]) for i in range(bootstrapping_samples)]

        pc = ParallelCalculation(parser_phase3.options.coresn, bootstrap_matrix, bs_args)
        
        pc_results = pc.run()

        bootstrap_matrices = zip(*pc_results)[1]

        return bootstrap_matrices

    return confdistmatrix





def prepare_ensembles_for_convergence_increasing_window(ensembles, window_size):

    ens_size = ensembles[0].coordinates.shape[0]

    rest_slices = ens_size / window_size
    residuals = ens_size % window_size
    slices_n = [0]

    for rs in range(rest_slices-1):
        slices_n.append(slices_n[-1] + window_size)
    if residuals != 0:
        slices_n.append(slices_n[-1] + residuals + window_size)
        logging.warning("the last window will be shorter than the prescribed window size (%s frames)"%residuals)
    else:
        slices_n.append(slices_n[-1] + window_size)
            
    for s in range(len(slices_n)-1):
        tmp_ensembles.append( Ensemble(topology = ensembles[0].topology,
                                       trajectory = [ensembles[0].topology],
                                       atom_selection_string = ensembles[0].atom_selection_string,
                                       superimposition_selection_string = ensembles[0].superimposition_subset))
        #print slices_n
        tmp_ensembles[-1].coordinates = ensembles[0].coordinates[slices_n[s]:slices_n[s+1],:,:]

    return tmp_ensembles





def hes(ensembles, 
        cov_estimator = "shrinkage",
        mass_weighted = True,
        details = None,
        estimate_error = False,
        error_estimation_mode = "bootstrapping",
        bootstrapping_runs = 100,):

    logging.info("Chosen metric: Harmonic similarity")
    if cov_estimator == "shrinkage":
        covariance_estimator = EstimatorShrinkage()
        logging.info("    Covariance matrix estimator: Shrinkage")
    elif cov_estimator == "ml":
        covariance_estimator = EstimatorML()
        logging.info("    Covariance matrix estimator: Maximum Likelihood")
    else:
        logging.error("Covariance estimator %s is not supported. Choose between 'shrinkage' and 'ml'." % cov_estimator)
        return None

    xs = []
    sigmas = []

    if estimate_error:
        if error_estimation_mode == "bootstrapping":
            data = []
            for t in range(parser_phase3.options.bootstrapping_runs):
                logging.info("The coordinates will be bootstrapped.")
                xs = []
                sigmas = []
                values = numpy.zeros((out_matrix_eln,out_matrix_eln))
                for e in ensembles:
                    this_coords = bootstrap_coordinates(e.coordinates, 1)[0]
                    xs.append(numpy.average(this_coords, axis=0).flatten())
                    sigmas.append( covariance_matrix(e,
                                                     mass_weighted=True,
                                                     estimator = covariance_estimator) )
                for i,j in pairs_indeces:
                    value = harmonic_ensemble_similarity(x1 = xs[i],
                                                         x2 = xs[j],
                                                         sigma1 = sigmas[i],
                                                         sigma2 = sigmas[j])
                    values[i,j] = value
                    values[j,i] = value
                data.append(values)
            outs = numpy.array(data)
            avgs = np.average(data, axis=0)
            stds = np.std(data, axis=0)

            return (avgs, stds)

        else: 
            logging.error("Only bootstrapping mode is supported so far.")
            return None



    # Calculate the parameters for the multivariate normal distribution of each ensemble
    values = numpy.zeros((len(ensembles), len(ensembles)))
    pairs_indeces = list( trm_indeces_nodiag(len(ensembles)) )


    for e in ensembles:
        print e    
        # Extract coordinates from each ensemble
        coordinates_system = e.coordinates

        # Average coordinates in each system
        xs.append(numpy.average(coordinates_system, axis=0).flatten())

        # Covariance matrices in each system
        sigmas.append( covariance_matrix(e, 
                                        mass_weighted = mass_weighted,
                                        estimator = covariance_estimator) )
                    
    for i,j in pairs_indeces:
        value = harmonic_ensemble_similarity(x1 = xs[i],
                                             x2 = xs[j],
                                             sigma1 = sigmas[i],
                                             sigma2 = sigmas[j])
        values[i,j] = value
        values[j,i] = value

    # Save details as required
    if details:
        kwds = {}
        for i in range(len(ensembles)):
            kwds['ensemble%d_mean'%(i+1)] = xs[i]
            kwds['ensemble%d_covariance_matrix'%(i+1)] = sigmas[i]
        numpy.savez(details, **kwds)

    return values


def ces(ensembles,
                            preference_values=[-1.0],
                            max_iterations = 500,
                            convergence = 50,
                            damping = 0.9,
                            noise = True,
                            mode = "ap",
                            similarity_matrix = None,
                            cluster_collections = None,
                            estimate_error = False,
                            error_estimation_mode = "bootstrapping",
                            boostrapped_matrices = None,
                            details = False,
                            np = 1,
                            **kwargs):




    ensemble_assignment = []
    for i in range(1, len(ensembles)+1):
        ensemble_assignment += [i for j in ensembles[i-1].coordinates]
    ensemble_assignment = numpy.array(ensemble_assignment)

    metadata = {'ensemble': ensemble_assignment}

    pairs_indeces = list( trm_indeces_nodiag(len(ensembles)) )

    if not conf_dist_matrix:
        confdistmatrix = get_similarity_matrix(  ensembles, **kwargs)
    else:
        confdistmatrix = similarity_matrix

    print confdistmatrix, "CDM"

    if mode == "ap":
        
        preferences = map(float, preference_values)
                                
        logging.info("    Clustering algorithm: Affinity Propagation")
        logging.info("        Preference values: %s" % ", ".join(map(lambda x: "%3.2f"%x ,preferences)))
        logging.info("        Maximum iterations: %d" % max_iterations)
        logging.info("        Convergence: %d" % convergence)
        logging.info("        Damping: %1.2f"%  damping)
        logging.info("        Apply noise to similarity matrix: %s" % str(noise))

        # Choose clustering algorithm
        clustalgo = AffinityPropagation()

        # Prepare input for parallel calculation
        if estimate_error:
            if error_estimation_mode == "bootstrapping":
                confdistmatrixs = []
                lams = []
                max_iterationss = []
                convergences = []
                noises = []
                real_prefs = []
                nmat = len(bootstrap_matrices)
                for p in preferences: 
                    confdistmatrixs.extend(bootstrap_matrices)
                    lams.extend([damping]*nmat)
                    max_iterationss.extend([max_iterations]*nmat)
                    noises.extend([noise]*nmat)
                    convergences.extend([convergence]*nmat)
                    real_prefs.extend([p]*nmat)
                old_prefs = preferences
                preferences = real_prefs
        else:
            confdistmatrixs = [ confdistmatrix for i in preferences ]
            lams = [ damping for i in preferences ]
            max_iterationss = [ max_iterations for i in preferences ]
            convergences = [ convergence for i in preferences ]
            noises = [ int(noise) for i in preferences ]

        args = zip(confdistmatrixs, preferences, lams, max_iterationss, convergences, noises)
        logging.info("    Starting affinity propagation runs . . .")

        # Do it
        pc = ParallelCalculation(np, clustalgo, args)

        results = pc.run()
        
        logging.info("\n    Done!")

        # Create clusters collections from clustering results, one for each cluster. None if clustering didn't work.
        ccs = [ ClustersCollection(clusters[1], metadata=metadata) for clusters in results ]
        
        if estimate_error:
            if error_estimation_mode == "bootstrapping":
                preferences = old_prefs
                k = 0
                for i,p in enumerate(preferences):
                    failed_runs = 0
                    values = []
                    for j in range(parser_phase3.options.bootstrapping_runs):
                        if ccs[k].clusters == None:
                            failed_runs += 1
                            k += 1
                            continue
                        values.append(numpy.zeros((out_matrix_eln,out_matrix_eln)))
        
                        for pair in pairs_indeces:
                            # Calculate dJS
                            this_djs = clustering_ensemble_similarity( ccs[k], ensembles[pair[0]], pair[0]+1, ensembles[pair[1]], pair[1]+1 )
                            values[-1][pair[0],pair[1]] = this_djs
                            values[-1][pair[1],pair[0]] = this_djs
                        k += 1
                    outs = numpy.array(values)
                    avgs = numpy.average(outs, axis=0)
                    stds = numpy.std(outs, axis=0)

            return (avgs, stds)

        values = {}
        for i,p in enumerate(preferences):
            if ccs[i].clusters == None:
                continue
            else:
                values[p] = numpy.zeros((len(ensembles),len(ensembles)))

                for pair in pairs_indeces:
                # Calculate dJS
                    this_val = clustering_ensemble_similarity( ccs[i], ensembles[pair[0]], pair[0]+1, ensembles[pair[1]], pair[1]+1)
                    values[p][pair[0],pair[1]] = this_val
                    values[p][pair[1],pair[0]] = this_val

            if details:
                kwds = {}
                kwds['centroids'] = numpy.array([c.centroid for c in ccs[i]])
                kwds['ensemble_sizes'] = numpy.array([e.coordinates.shape[0] for e in ensembles])
                for cln,cluster in enumerate(ccs[i]):
                    kwds["cluster%d"%(cln+1)] = numpy.array(cluster.elements)
                details_array = np.array(kwds)

                return values, details

    return values
        

def dres(   ensembles,
            conf_dist_matrix = None,
            mode='vanilla',
            dimensions = [3],
            maxlam = 2.0,
            minlam = 0.1,
            ncycle = 100,
            nstep = 10000,
            neighborhood_cutoff = 1.5,
            kn = 100,
            estimate_error = False,
            boostrapped_matrices = None,
            nsamples=1000,
            details = None,
            np=1,
            **kwargs):

    stressfreq = -1

    ensemble_assignment = []
    for i in range(1, len(ensembles)+1):
        ensemble_assignment += [i for j in ensembles[i-1].coordinates]
    ensemble_assignment = numpy.array(ensemble_assignment)

    metadata = {'ensemble': ensemble_assignment}

    pairs_indeces = list( trm_indeces_nodiag(len(ensembles)) )

    if not conf_dist_matrix:
        confdistmatrix = get_similarity_matrix(  ensembles, **kwargs)
    else:
        confdistmatrix = conf_dist_matrix

    dimensions = map(int, dimensions)

    # prepare runs. (e.g.: runs = [1,2,3,1,2,3,1,2,3, ...])
    if estimate_error:        
        runs = []
        for d in dimensions: 
            runs.extend([d]*len(bootstrapped_matrices))
        matrices = bootstrap_matrices*len(bootstrapped_matrices)
    else:
        runs = dimensions
        matrices = [confdistmatrix for i in runs]
    for d in dimensions:
        if d > confdistmatrix.size:
            logging.error("ERROR: The embedded space must have a number of dimensions inferior to the original space.")
            exit(1)
    
    # Choose algorithm and prepare options
    embedding_options = []
    if mode == 'vanilla':
        embedder = StochasticProximityEmbedding()
        for r in range(len(runs)):
            embedding_options += [(matrices[r], 
                              neighborhood_cutoff, 
                              runs[r],
                              maxlam,
                              minlam,
                              ncycle,
                              nstep,
                              stressfreq)]

    if mode == 'rn':
        embedder = RandomNeighborhoodStochasticProximityEmbedding()
        for r in range(len(runs)):
            embedding_options += [(matrices[r],
                              neighborhood_cutoff,
                              kn,
                              runs[r],
                              maxlam,
                              minlam,
                              ncycle,
                              stressfreq)]

    if mode == 'knn':
        embedder = kNNStochasticProximityEmbedding()
        for r in range(len(runs)):
            embedding_options += [(matrices[r],
                              kn,
                              runs[r],
                              maxlam,
                              minlam,
                              ncycle,
                              nstep,
                              stressfreq)]

    pc = ParallelCalculation(np, embedder, embedding_options)
    

    # Run parallel calculation
    results = pc.run()
    sleep(1)

    embedded_spaces_perdim = {}
    stresses_perdim = {}

    # Sort out obtained spaces and their residual stress values

    if estimate_error: # if bootstrap
        k = 0
        values = {}
        for ndim in dimensions:
            values[ndim] = []
            for i in range(len(bootstrapped_matrices)):

                header = "# ==== Number of dimensions: %d ==="%ndim
                values.append(numpy.zeros((len(ensembles),len(ensembles))))

                embedded_stress  = results[k][1][0]
                embedded_space = results[k][1][1]
                
                kdes, resamples, embedded_ensembles = gen_kde_pdfs(embedded_space, ensemble_assignment, parser_phase3.options.nensembles, nsamples = parser_phase3.options.samples)

                for pair in pairs_indeces:
                    this_value = dimred_ensemble_similarity(kdes[pair[0]], resamples[pair[0]], kdes[pair[1]],resamples[pair[1]])
                    values[-1][pair[0],pair[1]] = this_value
                    values[-1][pair[1],pair[0]] = this_value
                
                k += 1
                outs = numpy.array(values)
                avgs = numpy.average(outs, axis=0)
                stds = numpy.std(outs, axis=0)

        return (avgs, stds)                    

    values = {}

    for i in range(len(dimensions)):
        stresses_perdim[dimensions[i]] = []
        embedded_spaces_perdim[dimensions[i]] = []
        for j in range(1):
            stresses_perdim[dimensions[i]].append(results[j*len(dimensions)+i][1][0])
            embedded_spaces_perdim[dimensions[i]].append(results[j*len(dimensions)+i][1][1])

    for ndim in dimensions:

        values[ndim] = numpy.zeros((len(ensembles),len(ensembles)))

        embedded_spaces = embedded_spaces_perdim[ndim]
        embedded_stresses = stresses_perdim[ndim]

        embedded_stress = embedded_stresses[numpy.argmin(embedded_stresses)]
        embedded_space  = embedded_spaces[numpy.argmin(embedded_stresses)]

        kdes, resamples, embedded_ensembles = gen_kde_pdfs(embedded_space, ensemble_assignment, len(ensembles), nsamples = nsamples)
                             
        for pair in pairs_indeces:
            this_value = dimred_ensemble_similarity(kdes[pair[0]], resamples[pair[0]], kdes[pair[1]],resamples[pair[1]])
            values[ndim][pair[0],pair[1]] = this_value
            values[ndim][pair[1],pair[0]] = this_value

    return values

    if parser_phase3.options.details:
        kwds = {}
        kwds["stress"] = numpy.array([embedded_stress])
        for en,e in enumerate(embedded_ensembles):
            kwds[("ensemble%d"%en)] = e
        details_array = np.array(kwds)

    return values, details_array



def ces_ensemble_convergence(   original_ensembles, 
                                window_size,
                                preferences = [1.0],                            
                                max_iterations = 500,
                                convergence = 50,
                                damping = 0.9,
                                noise = True,   
                                save_matrix = None,
                                load_matrix = None):

    ensembles = prepare_ensembles_for_convergence_increasing_window(original_ensembles, window_size)

    if not similarity_matrix:
        similarity_matrix = get_similarity_matrix(arguments)

    preferences = preference_values
                            
    logging.info("    Clustering algorithm: Affinity Propagation")
    logging.info("        Preference values: %s" % ", ".join(map(lambda x: "%3.2f"%x ,preferences)))
    logging.info("        Maximum iterations: %d" % parser_phase3.options.max_iterations)
    logging.info("        Convergence: %d" % parser_phase3.options.convergence)
    logging.info("        Damping: %1.2f"%  parser_phase3.options.lam)
    logging.info("        Apply noise to similarity matrix: %s" % str(parser_phase3.options.noise))

    if len(preferences) % parser_phase3.options.coresn != 0:
        logging.warning("WARNING: for optimal performance, the number of cores should be a factor of the number of preference values.")

    confdistmatrixs = [ confdistmatrix for i in preferences ]
    lams = [ parser_phase3.options.lam for i in preferences ]
    max_iterationss = [ parser_phase3.options.max_iterations for i in preferences ]
    convergences = [ parser_phase3.options.convergence for i in preferences ]
    noises = [ int(parser_phase3.options.noise) for i in preferences ]

    clustalgo = AffinityPropagation()

    args = zip(confdistmatrixs, preferences, lams, max_iterationss, convergences, noises)

    logging.info("    Starting affinity propagation runs . . .")

    pc = ParallelCalculation(parser_phase3.options.coresn, clustalgo, args)

    results = pc.run()
    
    logging.info("\n    Done!")
    ccs = [ ClustersCollection(clusters[1], metadata=metadata) for clusters in results ]

    out = {}

    for i,p in enumerate(preferences):
        if ccs[i].clusters == None:
            continue
        out[p] = np.zeros(len(ensembles))
        for j in range(0,len(ensembles)):
            out[p][j] = cumulative_clustering_ensemble_similarity( ccs[i], 
                                                                        ensembles[-1], 
                                                                        len(ensembles)+1,
                                                                        ensembles[j], j+1)

    return out




def dres_ensemble_convergence():
    embedding_options = []
    if parser_phase3.options.spe_mode == 'vanilla':
        embedder = StochasticProximityEmbedding()
        for r in range(len(runs)):
            embedding_options += [(matrices[r], 
                              parser_phase3.options.neighborhood_cutoff, 
                              runs[r],
                              parser_phase3.options.maxlam,
                              parser_phase3.options.minlam,
                              parser_phase3.options.ncycle,
                              parser_phase3.options.nstep,
                              parser_phase3.options.stressfreq)]

    if parser_phase3.options.spe_mode == 'rn':
        embedder = RandomNeighborhoodStochasticProximityEmbedding()
        for r in range(len(runs)):
            embedding_options += [(matrices[r],
                              parser_phase3.options.neighborhood_cutoff,
                              parser_phase3.options.kn,
                              runs[r],
                              parser_phase3.options.maxlam,
                              parser_phase3.options.minlam,
                              parser_phase3.options.ncycle,
                              parser_phase3.options.stressfreq)]

    if parser_phase3.options.spe_mode == 'knn':
        embedder = kNNStochasticProximityEmbedding()
        for r in range(len(runs)):
            embedding_options += [(matrices[r],
                              parser_phase3.options.kn,
                              runs[r],
                              parser_phase3.options.maxlam,
                              parser_phase3.options.minlam,
                              parser_phase3.options.ncycle,
                              parser_phase3.options.nstep,
                              parser_phase3.options.stressfreq)]

    pc = ParallelCalculation(parser_phase3.options.coresn, embedder, embedding_options)

    for i in range(len(dimensions)):
        stresses_perdim[dimensions[i]] = []
        embedded_spaces_perdim[dimensions[i]] = []
        for j in range(1):
            stresses_perdim[dimensions[i]].append(results[j*len(dimensions)+i][1][0])
            embedded_spaces_perdim[dimensions[i]].append(results[j*len(dimensions)+i][1][1])

    out = {}

    for ndim in dimensions:

        embedded_spaces = embedded_spaces_perdim[ndim]
        embedded_stresses = stresses_perdim[ndim]

        embedded_stress = embedded_stresses[numpy.argmin(embedded_stresses)]
        embedded_space  = embedded_spaces[numpy.argmin(embedded_stresses)]

        kdes, resamples, embedded_ensembles = gen_kde_pdfs(embedded_space, ensemble_assignment, parser_phase3.options.nensembles, nsamples = parser_phase3.options.samples)


    # Run parallel calculation
    results = pc.run()
    sleep(1)

    embedded_spaces_perdim = {}
    stresses_perdim = {}
    out = {}

    for ndim in dimensions:

        out[ndim] = np.zeros(len(ensembles))

        embedded_spaces = embedded_spaces_perdim[ndim]
        embedded_stresses = stresses_perdim[ndim]

        embedded_stress = embedded_stresses[numpy.argmin(embedded_stresses)]
        embedded_space  = embedded_spaces[numpy.argmin(embedded_stresses)]


    # For every chosen dimension value:

        kdes, resamples, embedded_ensembles = cumulative_gen_kde_pdfs(embedded_space, ensemble_assignment, parser_phase3.options.nensembles-1, nsamples = parser_phase3.options.samples)
                
        for j in range(0,len(ensembles)):
            out[ndim][j] = dimred_ensemble_similarity(kdes[-1], 
                                                        resamples[-1], 
                                                        kdes[j], 
                                                        resamples[j])

    return out