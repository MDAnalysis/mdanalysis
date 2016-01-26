# covariance.py --- Covariance matrix calculations
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
Covariance calculation --- :mod:`encore.covariance`
=====================================================================

The module contains functions to estimate the covariance matrix of
an ensemble of structures.

"""

import sys
import MDAnalysis
import MDAnalysis.analysis
import MDAnalysis.analysis.align
import numpy

from Ensemble import Ensemble


class EstimatorML:
    ''' 
    Standard maximum likelihood estimator of the covariance matrix. 
    The generated object acts as a functor.
    '''
    def calculate(self, coordinates, reference_coordinates=None):
        '''
	**Arguments:**

        `coordinates` : numpy.array
            Flattened array of coordiantes
        
	`reference_coordinates` : numpy.array 
            Optional reference to use instead of mean

	**Returns:**

	`cov_mat` : numpy.array
            Estimate of  covariance matrix
	
        '''

        if reference_coordinates != None:

            # Offset from reference (for a normal covariance calculation
            # this would be the distance to the average)
            coordinates_offset = coordinates - reference_coordinates

            # Calculate covariance manually
            coordinates_cov = numpy.zeros((coordinates.shape[1], coordinates.shape[1]))
            for frame in coordinates_offset:
                coordinates_cov += numpy.outer(frame, frame)
            coordinates_cov /= coordinates.shape[0]

            return coordinates_cov

        else:
            return numpy.cov(coordinates, rowvar=0)

    __call__ = calculate
        
class EstimatorShrinkage:
    ''' 
    Shrinkage estimator of the covariance matrix using the method described in

    Improved Estimation of the Covariance Matrix of Stock Returns With an Application 
    to Portfolio Selection. Ledoit, O.; Wolf, M., Journal of Empirical Finance, 10, 5, 2003

    This implementation is based on the matlab code made available by Olivier Ledoit on
    his website:

    http://www.ledoit.net/ole2_abstract.htm
        
    The generated object acts as a functor.

    '''

    def __init__(self, shrinkage_parameter=None):
        '''
        Constructor.
        
        **Arguments:**

        `shrinkage_parameter` : float
            Makes it possible to set the shrinkage parameter explicitly, rather than having it estimated automatically.
        '''
        self.shrinkage_parameter = shrinkage_parameter


    def calculate(self, coordinates, reference_coordinates=None):
        '''
        
        **Arguments:**
        
        `coordinates` : numpy.array
            Flattened array of coordiantes
        `reference_coordinates`: numpy.array
            Optional reference to use instead of mean
            
        **Returns**:

        `cov_mat` : nump.array
            Covariance matrix
        '''

        x = coordinates
        t = x.shape[0]
        n = x.shape[1]

        mean_x = numpy.average(x, axis=0)

        # Use provided coordinates as "mean" if provided
        if reference_coordinates != None:
            mean_x = reference_coordinates

        x = x - mean_x
        xmkt = numpy.average(x, axis=1)

        # Call maximum likelihood estimator (note the additional column)
        sample = EstimatorML()(numpy.hstack([x,xmkt[:,numpy.newaxis]]), 0) * (t-1)/float(t)

        # Split covariance matrix into components
        covmkt = sample[0:n,n]
        varmkt = sample[n,n]
        sample = sample[:n,:n]

        # Prior
        prior = numpy.outer(covmkt,covmkt)/varmkt
        prior[numpy.ma.make_mask(numpy.eye(n))] = numpy.diag(sample)

        # If shrinkage parameter is not set, estimate it
        if self.shrinkage_parameter == None:

            # Frobenius norm
            c = numpy.linalg.norm(sample - prior, ord='fro')**2

            y = x**2
            p=1/float(t)*numpy.sum(numpy.dot(numpy.transpose(y),y))-numpy.sum(numpy.sum(sample**2))
            rdiag=1/float(t)*numpy.sum(numpy.sum(y**2))-numpy.sum(numpy.diag(sample)**2)
            z = x * numpy.repeat(xmkt[:,numpy.newaxis], n, axis=1)
            v1 = 1/float(t) * numpy.dot(numpy.transpose(y),z) - numpy.repeat(covmkt[:,numpy.newaxis],n, axis=1)*sample
            roff1 = (numpy.sum(v1*numpy.transpose(numpy.repeat(covmkt[:,numpy.newaxis],n, axis=1)))/varmkt -
                     numpy.sum(numpy.diag(v1)*covmkt)/varmkt)
            v3 = 1/float(t)*numpy.dot(numpy.transpose(z),z) - varmkt*sample
            roff3 = (numpy.sum(v3*numpy.outer(covmkt, covmkt))/varmkt**2 - 
                     numpy.sum(numpy.diag(v3)*covmkt**2)/varmkt**2)
            roff=2*roff1-roff3;
            r=rdiag+roff;

            # Shrinkage constant
            k=(p-r)/c;
            self.shrinkage_parameter=max(0,min(1,k/float(t)))

        # calculate covariance matrix
        sigma=self.shrinkage_parameter*prior+(1-self.shrinkage_parameter)*sample;

        return sigma
    
    __call__ = calculate

def covariance_matrix(ensemble,
                      estimator = EstimatorShrinkage(),
                      mass_weighted=True,
                      reference = None,
                      start=0,
                      end=None):

    '''
    Calculates (optionally mass weighted) covariance matrix

    **Arguments:**
    
    `ensemble` : Ensemble object
        The structural ensemble
    
    `estimator` : MLEstimator or ShrinkageEstimator object
        Which estimator type to use (maximum likelihood, shrinkage). This object is required to have a __call__ function defined.
    
    `mass_weighted` : bool
        Whether to do a mass-weighted analysis
    
    `reference` : MDAnalysis.Universe object
        Use the distances to a specific reference structure rather than the distance to the mean.
        
    **Returns:**

    cov_mat : numpy.array
        Covariance matrix

    '''

    # Extract coordinates from ensemble
    # coordinates = ensemble.get_coordinates(start=start, end=end)
    coordinates = ensemble.coordinates


    # Flatten coordinate matrix into n_frame x n_coordinates
    coordinates = numpy.reshape(coordinates, (coordinates.shape[0], -1))

    # Extract coordinates from reference structure, if specified
    reference_coordinates = None
    if reference:

        # Select the same atoms in reference structure
        reference_atom_selection = reference.select_atoms(ensemble.get_atom_selection_string())
        reference_coordinates = reference_atom_selection.atoms.coordinates()

        # Flatten reference coordinates 
        reference_coordinates = reference_coordinates.flatten()

    sigma = estimator(coordinates, reference_coordinates)

    # Optionally correct with mass-weighting
    if mass_weighted:
        # Calculate mass-weighted covariance matrix
        masses = numpy.repeat(ensemble.atom_selection.masses, 3)
        mass_matrix = numpy.sqrt(numpy.identity(len(masses))*masses)
        sigma = numpy.dot(mass_matrix, numpy.dot(sigma, mass_matrix))

    return sigma
    

    
if __name__ == "__main__":

    import optparse
    import utils

    usage = "usage %prog [options]"
    parser = optparse.OptionParser(usage=usage)
    
    parser.add_option("--topology", dest="topology_filename", default="",
                      help="Path to a topology file (supported formats: PDB,PSF,CRD,GRO)")
    parser.add_option("--trajectory", dest="trajectory_filename", action="callback", type="string", nargs=0, default=[],
                      callback=utils.vararg_callback,
                      metavar="TRAJECTORY_FILENAME(S)",
                      help="Add trajectory filenames")
    parser.add_option("--atom-selection", dest="atom_selection_string", default="(name CA)",
                      help="CHARMM-style atom selection")
    parser.add_option("--mass-weighted-analysis", dest="mass_weighted_analysis", action="store_true", default=False,
                      help="Use mass-weighting in the calculation of the covariance matrix.")
    parser.add_option("--no-align", dest="no_align", action="store_true", default=False,
                      help="Do not superimpose structures before calculating covariance.")
    parser.add_option("--frame-interval", dest="frame_interval", default=1,
                      help="The interval between frames used for the analysis")
    parser.add_option("--use-distance-to-reference", dest="use_distance_to_reference", action="store_true", default=False,
                      help="Whether to use the distance to the reference structure rather than the distance to the average structure when calculating covariance matrix.")
    parser.add_option("--output", dest="output_filename", default="", 
                      help="Output file for covariance matrix")
    parser.add_option("--covariance-estimator", type="choice", dest="covariance_estimator", default="shrinkage",
                      choices=["ml","shrinkage"],
                      help="Type of estimator (maximum likelihood (ml) or shrinkage")
    (options, args) = parser.parse_args()

    if not options.trajectory_filename or not options.topology_filename:
        parser.error("--trajectory and --topology options must be specified")


    # Construct reference if available
    try:
        reference = MDAnalysis.Universe(options.topology_filename, 
                                        options.topology_filename)
    except:
        reference = None

    # Construct ensemble
    ensemble = Ensemble(topology=options.topology_filename,
                        trajectory=options.trajectory_filename, 
                        atom_selection_string=options.atom_selection_string,
                        frame_interval=options.frame_interval)

    # Align ensemble to reference
    if not options.no_align:
        ensemble.align(reference)

    # Select covariance estimator
    estimator = EstimatorML()
    if options.covariance_estimator == "shrinkage":
        estimator = EstimatorShrinkage()

    # Disable reference unless use_distance_to_reference is set
    if not options.use_distance_to_reference:
        reference = None
    
    # Calculate covariance matrix
    cov_matrix = covariance_matrix(ensemble,
                                   estimator = estimator,
                                   mass_weighted=options.mass_weighted_analysis,
                                   reference=reference)
    # Output results
    print >> sys.stderr, "Dimensions: (%s,%s)" %  cov_matrix.shape
    print >> sys.stderr, "Trace : %s" % numpy.trace(cov_matrix)
 
    if options.output_filename:
        print >> sys.stderr, cov_matrix
        output_file = open(options.output_filename, 'w')
        numpy.savetxt(output_file, cov_matrix, delimiter=" ")
        output_file.close()
    else:
        numpy.savetxt(sys.stdout, cov_matrix, delimiter=" ")

