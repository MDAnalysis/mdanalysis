from libc.stdint cimport uint64_t, UINT64_MAX


cdef extern from "string.h":
    void* memcpy(void* dst, void* src, int len)

cdef extern from "calc_distances.h":
    ctypedef float coordinate[3]
    cdef bint USED_OPENMP
    void _calc_distance_array(coordinate* ref, uint64_t numref, coordinate* conf, uint64_t numconf, double* distances)
    void _calc_distance_array_ortho(coordinate* ref, uint64_t numref, coordinate* conf, uint64_t numconf, float* box, double* distances)
    void _calc_distance_array_triclinic(coordinate* ref, uint64_t numref, coordinate* conf, uint64_t numconf, float* box, double* distances)
    void _calc_self_distance_array(coordinate* ref, uint64_t numref, double* distances)
    void _calc_self_distance_array_ortho(coordinate* ref, uint64_t numref, float* box, double* distances)
    void _calc_self_distance_array_triclinic(coordinate* ref, uint64_t numref, float* box, double* distances)
    void _coord_transform(coordinate* coords, uint64_t numCoords, double* box)
    void _calc_bond_distance(coordinate* atom1, coordinate* atom2, uint64_t numatom, double* distances)
    void _calc_bond_distance_ortho(coordinate* atom1, coordinate* atom2, uint64_t numatom, float* box, double* distances)
    void _calc_bond_distance_triclinic(coordinate* atom1, coordinate* atom2, uint64_t numatom, float* box, double* distances)
    void _calc_angle(coordinate* atom1, coordinate* atom2, coordinate* atom3, uint64_t numatom, double* angles)
    void _calc_angle_ortho(coordinate* atom1, coordinate* atom2, coordinate* atom3, uint64_t numatom, float* box, double* angles)
    void _calc_angle_triclinic(coordinate* atom1, coordinate* atom2, coordinate* atom3, uint64_t numatom, float* box, double* angles)
    void _calc_dihedral(coordinate* atom1, coordinate* atom2, coordinate* atom3, coordinate* atom4, uint64_t numatom, double* angles)
    void _calc_dihedral_ortho(coordinate* atom1, coordinate* atom2, coordinate* atom3, coordinate* atom4, uint64_t numatom, float* box, double* angles)
    void _calc_dihedral_triclinic(coordinate* atom1, coordinate* atom2, coordinate* atom3, coordinate* atom4, uint64_t numatom, float* box, double* angles)
    void _ortho_pbc(coordinate* coords, uint64_t numcoords, float* box)
    void _triclinic_pbc(coordinate* coords, uint64_t numcoords, float* box)
    void minimum_image(double* x, float* box, float* inverse_box)
    void minimum_image_triclinic(float* x, float* box, float* inverse_box)