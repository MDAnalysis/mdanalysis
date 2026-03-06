"""
This module contains all of the qmmm namelist variables for the
amber programs and automatically loads those dictionaries with the
default values found in that program (sander only).
"""
from typing import Dict

from .typing import InputDataType
class qmmm:

    def __init__(self):
        self.sander: Dict[str, InputDataType] = {
            'qmcut' : -1, 'iqmatoms' : '', 'qmmask' : '', 'qmgb' : 2,
            'qm_theory' : '', 'qmcharge' : 0, 'qmqmdx' : 1,
            'verbosity' : 0, 'tight_p_conv' : 0, 'scfconv' : 1.0e-8,
            'errconv' : 1e-1, 'ndiis_matrices' :6, 'ndiis_attempts' : 0,
            'printcharges' : 0, 'printdipole' : 0, 'peptide_corr' : 0,
            'itrmax' : 1000, 'qmshake' : 1, 'qmqm_erep_incore' : 1,
            'qmmmrij_incore' : 1, 'lnk_dis' : 1.09, 'lnk_atomic_no' : 1,
            'lnk_method' : 1, 'spin' : 1, 'pseudo_diag' : 1,
            'pseudo_diag_criteria' : 0.05, 'qm_ewald' : 1, 'qm_pme' : 1,
            'kmaxqx' : 5, 'kmaxqy' : 5, 'kmaxqz' : 5, 'ksqmaxq' : 27,
            'writepdb' : 0, 'qmmm_int' : 1, 'adjust_q' : 2,
            'diag_routine' : 1, 'density_predict' : 0,
            'fock_predict' : 0, 'fockp_d1' : 2.4, 'fockp_d2' : -1.2,
            'fockp_d3' : -0.8, 'fockp_d4' : 0.6, 'idc' : 0, 'divpb' : 0,
            'dftb_maxiter' : 70, 'dftb_disper' : 0,
            'dftb_3rd_order' : 'NONE', 'dftb_chg' : 0,
            'dftb_telec' : 0, 'dftb_telec_step' : 0,
            'qmmm_omp_max_threads' : 1, 'chg_lambda' : 1,
            'nearest_qm_solvent' : 0, 'nearest_qm_solvent_fq' : 1,
            'nearest_qm_solvent_resname' : 'WAT ',
        }
