# Copyright (c) 2016-2025 The Regents of the University of Michigan
# Part of GSD, released under the BSD 2-Clause License.

from libc.stdint cimport uint8_t, int8_t, uint16_t, int16_t, uint32_t, int32_t,\
    uint64_t, int64_t

cdef extern from "gsd.h" nogil:
    cdef enum gsd_type:
        GSD_TYPE_UINT8=1
        GSD_TYPE_UINT16
        GSD_TYPE_UINT32
        GSD_TYPE_UINT64
        GSD_TYPE_INT8
        GSD_TYPE_INT16
        GSD_TYPE_INT32
        GSD_TYPE_INT64
        GSD_TYPE_FLOAT
        GSD_TYPE_DOUBLE
        GSD_TYPE_CHARACTER

    cdef enum gsd_open_flag:
        GSD_OPEN_READWRITE=1
        GSD_OPEN_READONLY
        GSD_OPEN_APPEND

    cdef enum gsd_error:
        GSD_SUCCESS = 0
        GSD_ERROR_IO = -1
        GSD_ERROR_INVALID_ARGUMENT = -2
        GSD_ERROR_NOT_A_GSD_FILE = -3
        GSD_ERROR_INVALID_GSD_FILE_VERSION = -4
        GSD_ERROR_FILE_CORRUPT = -5
        GSD_ERROR_MEMORY_ALLOCATION_FAILED = -6
        GSD_ERROR_NAMELIST_FULL = -7
        GSD_ERROR_FILE_MUST_BE_WRITABLE = -8
        GSD_ERROR_FILE_MUST_BE_READABLE = -9

    cdef struct gsd_header:
        uint64_t magic
        uint32_t gsd_version
        char application[64]
        char schema[64]
        uint32_t schema_version
        uint64_t index_location
        uint64_t index_allocated_entries
        uint64_t namelist_location
        uint64_t namelist_allocated_entries
        char reserved[80]

    cdef struct gsd_index_entry:
        uint64_t frame
        uint64_t N
        int64_t location
        uint32_t M
        uint16_t id
        uint8_t type
        uint8_t flags

    cdef struct gsd_namelist_entry:
        char name[64]

    cdef struct gsd_index_buffer:
        gsd_index_entry *data
        size_t size
        size_t reserved
        void *mapped_data
        size_t mapped_len

    cdef struct gsd_name_id_map:
        void *v
        size_t size

    cdef struct gsd_write_buffer:
        char *data
        size_t size
        size_t reserved

    cdef struct gsd_handle:
        int fd
        gsd_header header
        gsd_index_buffer file_index
        gsd_index_buffer frame_index
        gsd_index_buffer buffer_index
        gsd_write_buffer write_buffer
        gsd_namelist_entry *namelist
        uint64_t namelist_num_entries
        uint64_t cur_frame
        int64_t file_size
        gsd_open_flag open_flags
        gsd_name_id_map name_map
        uint64_t namelist_written_entries

    uint32_t gsd_make_version(unsigned int major, unsigned int minor)
    int gsd_create(const char *fname,
                   const char *application,
                   const char *schema,
                   uint32_t schema_version)
    int gsd_create_and_open(gsd_handle* handle,
                            const char *fname,
                            const char *application,
                            const char *schema,
                            uint32_t schema_version,
                            const gsd_open_flag flags,
                            int exclusive_create)
    int gsd_open(gsd_handle* handle, const char *fname,
                 const gsd_open_flag flags)
    int gsd_truncate(gsd_handle* handle)
    int gsd_close(gsd_handle* handle)
    int gsd_end_frame(gsd_handle* handle)
    int gsd_flush(gsd_handle* handle)
    int gsd_write_chunk(gsd_handle* handle,
                        const char *name,
                        gsd_type type,
                        uint64_t N,
                        uint8_t M,
                        uint8_t flags,
                        const void *data)
    const gsd_index_entry* gsd_find_chunk(gsd_handle* handle,
                                          uint64_t frame,
                                          const char *name)
    int gsd_read_chunk(gsd_handle* handle, void* data,
                       const gsd_index_entry* chunk)
    uint64_t gsd_get_nframes(gsd_handle* handle)
    size_t gsd_sizeof_type(gsd_type type)
    const char *gsd_find_matching_chunk_name(gsd_handle* handle,
                                             const char *match,
                                             const char *prev)
    int gsd_upgrade(gsd_handle *handle)
    uint64_t gsd_get_maximum_write_buffer_size(gsd_handle* handle)
    int gsd_set_maximum_write_buffer_size(gsd_handle* handle, uint64_t size)
