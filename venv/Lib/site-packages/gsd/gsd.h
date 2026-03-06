// Copyright (c) 2016-2025 The Regents of the University of Michigan
// Part of GSD, released under the BSD 2-Clause License.

#ifndef GSD_H
#define GSD_H

#include <stdbool.h>
#include <stdint.h>
#include <string.h>

#ifdef __cplusplus
extern "C"
    {
#endif

    /*! \file gsd.h
        \brief Declare GSD data types and C API
    */

    /// Identifiers for the gsd data chunk element types
    enum gsd_type
        {
        /// Unsigned 8-bit integer.
        GSD_TYPE_UINT8 = 1,

        /// Unsigned 16-bit integer.
        GSD_TYPE_UINT16,

        /// Unsigned 32-bit integer.
        GSD_TYPE_UINT32,

        /// Unsigned 53-bit integer.
        GSD_TYPE_UINT64,

        /// Signed 8-bit integer.
        GSD_TYPE_INT8,

        /// Signed 16-bit integer.
        GSD_TYPE_INT16,

        /// Signed 32-bit integer.
        GSD_TYPE_INT32,

        /// Signed 64-bit integer.
        GSD_TYPE_INT64,

        /// 32-bit floating point number.
        GSD_TYPE_FLOAT,

        /// 64-bit floating point number.
        GSD_TYPE_DOUBLE,

        /// 8-bit character.
        GSD_TYPE_CHARACTER
        };

    /// Flag for GSD file open options
    enum gsd_open_flag
        {
        /// Open for both reading and writing
        GSD_OPEN_READWRITE = 1,

        /// Open only for reading
        GSD_OPEN_READONLY,

        /// Open only for writing
        GSD_OPEN_APPEND
        };

    /// Error return values
    enum gsd_error
        {
        /// Success.
        GSD_SUCCESS = 0,

        /// IO error. Check ``errno`` for details
        GSD_ERROR_IO = -1,

        /// Invalid argument passed to function.
        GSD_ERROR_INVALID_ARGUMENT = -2,

        /// The file is not a GSD file.
        GSD_ERROR_NOT_A_GSD_FILE = -3,

        /// The GSD file version cannot be read.
        GSD_ERROR_INVALID_GSD_FILE_VERSION = -4,

        /// The GSD file is corrupt.
        GSD_ERROR_FILE_CORRUPT = -5,

        /// GSD failed to allocated memory.
        GSD_ERROR_MEMORY_ALLOCATION_FAILED = -6,

        /// The GSD file cannot store any additional unique data chunk names.
        GSD_ERROR_NAMELIST_FULL = -7,

        /** This API call requires that the GSD file opened in with the mode GSD_OPEN_APPEND or
            GSD_OPEN_READWRITE.
        */
        GSD_ERROR_FILE_MUST_BE_WRITABLE = -8,

        /** This API call requires that the GSD file opened the mode GSD_OPEN_READ or
            GSD_OPEN_READWRITE.
        */
        GSD_ERROR_FILE_MUST_BE_READABLE = -9,
        };

    enum
        {
        /** v1 file: Size of a GSD name in memory. v2 file: The name buffer size is a multiple of
            GSD_NAME_SIZE.
        */
        GSD_NAME_SIZE = 64
        };

    enum
        {
        /// Reserved bytes in the header structure
        GSD_RESERVED_BYTES = 80
        };

    /** GSD file header

        The in-memory and on-disk storage of the GSD file header. Stored in the first 256 bytes of
        the file.

        @warning All members are **read-only** to the caller.
    */
    struct gsd_header
        {
        /// Magic number marking that this is a GSD file.
        uint64_t magic;

        /// Location of the chunk index in the file.
        uint64_t index_location;

        /// Number of index entries that will fit in the space allocated.
        uint64_t index_allocated_entries;

        /// Location of the name list in the file.
        uint64_t namelist_location;

        /// Number of bytes in the namelist divided by GSD_NAME_SIZE.
        uint64_t namelist_allocated_entries;

        /// Schema version: from gsd_make_version().
        uint32_t schema_version;

        /// GSD file format version from gsd_make_version().
        uint32_t gsd_version;

        /// Name of the application that generated this file.
        char application[GSD_NAME_SIZE];

        /// Name of data schema.
        char schema[GSD_NAME_SIZE];

        /// Reserved for future use.
        char reserved[GSD_RESERVED_BYTES];
        };

    /** Index entry

        An index entry for a single chunk of data.

        @warning All members are **read-only** to the caller.
    */
    struct gsd_index_entry
        {
        /// Frame index of the chunk.
        uint64_t frame;

        /// Number of rows in the chunk.
        uint64_t N;

        /// Location of the chunk in the file.
        int64_t location;

        /// Number of columns in the chunk.
        uint32_t M;

        /// Index of the chunk name in the name list.
        uint16_t id;

        /// Data type of the chunk: one of gsd_type.
        uint8_t type;

        /// Flags (for internal use).
        uint8_t flags;
        };

    /** Name/id mapping

        A string name paired with an ID. Used for storing sorted name/id mappings in a hash map.
    */
    struct gsd_name_id_pair
        {
        /// Pointer to name (actual name storage is allocated in gsd_handle)
        char* name;

        /// Next name/id pair with the same hash
        struct gsd_name_id_pair* next;

        /// Entry id
        uint16_t id;
        };

    /** Name/id hash map

        A hash map of string names to integer identifiers.
    */
    struct gsd_name_id_map
        {
        /// Name/id mappings
        struct gsd_name_id_pair* v;

        /// Number of entries in the mapping
        size_t size;
        };

    /** Array of index entries

        May point to a mapped location of index entries in the file or an in-memory buffer.
    */
    struct gsd_index_buffer
        {
        /// Indices in the buffer
        struct gsd_index_entry* data;

        /// Number of entries in the buffer
        size_t size;

        /// Number of entries available in the buffer
        size_t reserved;

        /// Pointer to mapped data (NULL if not mapped)
        void* mapped_data;

        /// Number of bytes mapped
        size_t mapped_len;
        };

    /** Byte buffer

        Used to buffer of small data chunks held for a buffered write at the end of a frame. Also
        used to hold the names.
    */
    struct gsd_byte_buffer
        {
        /// Data
        char* data;

        /// Number of bytes in the buffer
        size_t size;

        /// Number of bytes available in the buffer
        size_t reserved;
        };

    /** Name buffer

        Holds a list of string names in order separated by NULL terminators. In v1 files, each name
        is 64 bytes. In v2 files, only one NULL terminator is placed between each name.
    */
    struct gsd_name_buffer
        {
        /// Data
        struct gsd_byte_buffer data;

        /// Number of names in the list
        size_t n_names;
        };

    /** File handle

        A handle to an open GSD file.

        This handle is obtained when opening a GSD file and is passed into every method that
        operates on the file.

        @warning All members are **read-only** to the caller.
    */
    struct gsd_handle
        {
        /// File descriptor
        int fd;

        /// The file header
        struct gsd_header header;

        /// Mapped data chunk index
        struct gsd_index_buffer file_index;

        /// Buffered index entries to append to the current frame
        struct gsd_index_buffer buffer_index;

        /// Buffered write data
        struct gsd_byte_buffer write_buffer;

        /// List of names stored in the file
        struct gsd_name_buffer file_names;

        /// List of names added in the current frame
        struct gsd_name_buffer frame_names;

        /// The index of the last frame in the buffer
        uint64_t buffer_frame;

        /// The index of the last frame comitted to the file
        uint64_t file_frame;

        /// Size of the file (in bytes)
        int64_t file_size;

        /// Flags passed to gsd_open() when opening this handle
        enum gsd_open_flag open_flags;

        /// Access the names in the namelist
        struct gsd_name_id_map name_map;

        /// Number of index entries pending in the current frame.
        uint64_t pending_index_entries;

        /// Maximum write buffer size (bytes).
        uint64_t maximum_write_buffer_size;
        };

    /** Specify a version.

        @param major major version
        @param minor minor version

        @return a packed version number aaaa.bbbb suitable for storing in a gsd file version entry.
    */
    uint32_t gsd_make_version(unsigned int major, unsigned int minor);

    /** Create a GSD file.

        @param fname File name (UTF-8 encoded).
        @param application Generating application name (truncated to 63 chars).
        @param schema Schema name for data to be written in this GSD file (truncated to 63 chars).
        @param schema_version Version of the scheme data to be written (make with
        gsd_make_version()).

        @post Create an empty gsd file in a file of the given name. Overwrite any existing file at
        that location.

        The generated gsd file is not opened. Call gsd_open() to open it for writing.

        @return
          - GSD_SUCCESS (0) on success. Negative value on failure:
          - GSD_ERROR_IO: IO error (check errno).
    */
    int gsd_create(const char* fname,
                   const char* application,
                   const char* schema,
                   uint32_t schema_version);

    /** Create and open a GSD file.

        @param handle Handle to open.
        @param fname File name (UTF-8 encoded).
        @param application Generating application name (truncated to 63 chars).
        @param schema Schema name for data to be written in this GSD file (truncated to 63 chars).
        @param schema_version Version of the scheme data to be written (make with
            gsd_make_version()).
        @param flags Either GSD_OPEN_READWRITE, or GSD_OPEN_APPEND.
        @param exclusive_create Set to non-zero to force exclusive creation of the file.

        @post Create an empty gsd file with the given name. Overwrite any existing file at that
        location.

        Open the generated gsd file in *handle*.

        The file descriptor is closed if there when an error opening the file.

        @return
          - GSD_SUCCESS (0) on success. Negative value on failure:
          - GSD_ERROR_IO: IO error (check errno).
          - GSD_ERROR_NOT_A_GSD_FILE: Not a GSD file.
          - GSD_ERROR_INVALID_GSD_FILE_VERSION: Invalid GSD file version.
          - GSD_ERROR_FILE_CORRUPT: Corrupt file.
          - GSD_ERROR_MEMORY_ALLOCATION_FAILED: Unable to allocate memory.
    */
    int gsd_create_and_open(struct gsd_handle* handle,
                            const char* fname,
                            const char* application,
                            const char* schema,
                            uint32_t schema_version,
                            enum gsd_open_flag flags,
                            int exclusive_create);

    /** Open a GSD file.

        @param handle Handle to open.
        @param fname File name to open (UTF-8 encoded).
        @param flags Either GSD_OPEN_READWRITE, GSD_OPEN_READONLY, or GSD_OPEN_APPEND.

        @pre The file name *fname* is a GSD file.

        @post Open a GSD file and populates the handle for use by API calls.

        The file descriptor is closed if there is an error opening the file.

        @return
          - GSD_SUCCESS (0) on success. Negative value on failure:
          - GSD_ERROR_IO: IO error (check errno).
          - GSD_ERROR_NOT_A_GSD_FILE: Not a GSD file.
          - GSD_ERROR_INVALID_GSD_FILE_VERSION: Invalid GSD file version.
          - GSD_ERROR_FILE_CORRUPT: Corrupt file.
          - GSD_ERROR_MEMORY_ALLOCATION_FAILED: Unable to allocate memory.
    */
    int gsd_open(struct gsd_handle* handle, const char* fname, enum gsd_open_flag flags);

    /** Truncate a GSD file.

        @param handle Open GSD file to truncate.

        After truncating, a file will have no frames and no data chunks. The file size will be that
        of a newly created gsd file. The application, schema, and schema version metadata will be
        kept. Truncate does not close and reopen the file, so it is suitable for writing restart
        files on Lustre file systems without any metadata access.

        @return
          - GSD_SUCCESS (0) on success. Negative value on failure:
          - GSD_ERROR_IO: IO error (check errno).
          - GSD_ERROR_NOT_A_GSD_FILE: Not a GSD file.
          - GSD_ERROR_INVALID_GSD_FILE_VERSION: Invalid GSD file version.
          - GSD_ERROR_FILE_CORRUPT: Corrupt file.
          - GSD_ERROR_MEMORY_ALLOCATION_FAILED: Unable to allocate memory.
    */
    int gsd_truncate(struct gsd_handle* handle);

    /** Close a GSD file.

        @param handle GSD file to close.

        @pre *handle* was opened by gsd_open().

        @post Writable files: All data and index entries buffered before the previous call to
              gsd_end_frame() is written to the file by implicitly calling gsd_flush().
        @post The file is closed.
        @post *handle* is freed and can no longer be used.

        @note There is no need to manually call gsd_flush() before gsd_close().

        @warning Ensure that all gsd_write_chunk() calls are completed with gsd_end_frame() before
        closing the file.

        @return
          - GSD_SUCCESS (0) on success. Negative value on failure:
          - GSD_ERROR_IO: IO error (check errno).
          - GSD_ERROR_INVALID_ARGUMENT: *handle* is NULL.
    */
    int gsd_close(struct gsd_handle* handle);

    /** Complete the current frame.

        @param handle Handle to an open GSD file

        @pre *handle* was opened by gsd_open().

        @post The current frame counter is increased by 1.

        @note Starting with GSD 4.0, gsd_end_frame() does NOT automatically flush buffered frames
        to the filesystem. Callers must manually call gsd_flush() or gsd_close() to commit the
        buffered frames to the filesystem.

        @return
          - GSD_SUCCESS (0) on success. Negative value on failure:
          - GSD_ERROR_IO: IO error (check errno).
          - GSD_ERROR_INVALID_ARGUMENT: *handle* is NULL.
          - GSD_ERROR_FILE_MUST_BE_WRITABLE: The file was opened read-only.
          - GSD_ERROR_MEMORY_ALLOCATION_FAILED: Unable to allocate memory.
    */
    int gsd_end_frame(struct gsd_handle* handle);

    /** Flush the write buffer.

        @param handle Handle to an open GSD file

        @pre *handle* was opened by gsd_open().

        @post All data buffered by gsd_write_chunk() are present in the file.
        @post All index entries buffered by gsd_write_chunk() prior to the last call to
              gsd_end_frame() are present in the file.

        @return
          - GSD_SUCCESS (0) on success. Negative value on failure:
          - GSD_ERROR_IO: IO error (check errno).
          - GSD_ERROR_INVALID_ARGUMENT: *handle* is NULL.
          - GSD_ERROR_FILE_MUST_BE_WRITABLE: The file was opened read-only.
          - GSD_ERROR_MEMORY_ALLOCATION_FAILED: Unable to allocate memory.
    */
    int gsd_flush(struct gsd_handle* handle);

    /** Add a data chunk to the current frame.

        @param handle Handle to an open GSD file.
        @param name Name of the data chunk.
        @param type type ID that identifies the type of data in *data*.
        @param N Number of rows in the data.
        @param M Number of columns in the data.
        @param flags set to 0, non-zero values reserved for future use.
        @param data Data buffer.

        @pre *handle* was opened by gsd_open().
        @pre *name* is a unique name for data chunks in the given frame.
        @pre data is allocated and contains at least `N * M * gsd_sizeof_type(type)` bytes.

        @post When there is space in the buffer: The given data is present in the write buffer.
              Otherwise, the data is present at the end of the file.
        @post The index is present in the buffer.

        @note If the GSD file is version 1.0, the chunk name is truncated to 63 bytes. GSD version
        2.0 files support arbitrarily long names.

        @note *N* == 0 is allowed. When *N* is 0, *data* may be NULL.

        @return
          - GSD_SUCCESS (0) on success. Negative value on failure:
          - GSD_ERROR_IO: IO error (check errno).
          - GSD_ERROR_INVALID_ARGUMENT: *handle* is NULL, *N* == 0, *M* == 0, *type* is invalid, or
            *flags* != 0.
          - GSD_ERROR_FILE_MUST_BE_WRITABLE: The file was opened read-only.
          - GSD_ERROR_NAMELIST_FULL: The file cannot store any additional unique chunk names.
          - GSD_ERROR_MEMORY_ALLOCATION_FAILED: failed to allocate memory.
    */
    int gsd_write_chunk(struct gsd_handle* handle,
                        const char* name,
                        enum gsd_type type,
                        uint64_t N,
                        uint32_t M,
                        uint8_t flags,
                        const void* data);

    /** Find a chunk in the GSD file.

        @param handle Handle to an open GSD file
        @param frame Frame to look for chunk
        @param name Name of the chunk to find

        @pre *handle* was opened by gsd_open() in read or readwrite mode.

        The found entry contains size and type metadata and can be passed to gsd_read_chunk() to
        read the data.

        @return A pointer to the found chunk, or NULL if not found.

        @note In read/write files gsd_find_chunk() can only find chunks that have been committed
        to the file with gsd_flush().
    */
    const struct gsd_index_entry*
    gsd_find_chunk(struct gsd_handle* handle, uint64_t frame, const char* name);

    /** Read a chunk from the GSD file.

        @param handle Handle to an open GSD file.
        @param data Data buffer to read into.
        @param chunk Chunk to read.

        @pre *handle* was opened in read or readwrite mode.
        @pre *chunk* was found by gsd_find_chunk().
        @pre *data* points to an allocated buffer with at least `N * M * gsd_sizeof_type(type)`
       bytes.

        @return
          - GSD_SUCCESS (0) on success. Negative value on failure:
          - GSD_ERROR_IO: IO error (check errno).
          - GSD_ERROR_INVALID_ARGUMENT: *handle* is NULL, *data* is NULL, or *chunk* is NULL.
          - GSD_ERROR_FILE_MUST_BE_READABLE: The file was opened in append mode.
          - GSD_ERROR_FILE_CORRUPT: The GSD file is corrupt.

        @note In read/write files gsd_read_chunk() can only read chunks that have been committed
        to the file with gsd_flush().
    */
    int gsd_read_chunk(struct gsd_handle* handle, void* data, const struct gsd_index_entry* chunk);

    /** Get the number of frames in the GSD file.

        @param handle Handle to an open GSD file

        @pre *handle* was opened by gsd_open().

        @return The number of frames in the file, or 0 on error.
    */
    uint64_t gsd_get_nframes(struct gsd_handle* handle);

    /** Query size of a GSD type ID.

        @param type Type ID to query.

        @return Size of the given type in bytes, or 0 for an unknown type ID.
    */
    size_t gsd_sizeof_type(enum gsd_type type);

    /** Search for chunk names in a gsd file.

        @param handle Handle to an open GSD file.
        @param match String to match.
        @param prev Search starting point.

        @pre *handle* was opened by gsd_open()
        @pre *prev* was returned by a previous call to gsd_find_matching_chunk_name()

        To find the first matching chunk name, pass NULL for prev. Pass in the previous found string
        to find the next after that, and so on. Chunk names match if they begin with the string in
        *match*. Chunk names returned by this function may be present in at least one frame.

        @return Pointer to a string, NULL if no more matching chunks are found found, or NULL if
        *prev* is invalid

        @note In read/write files gsd_find_matching_chunk_names() can only find names that have
        been committed to the file with gsd_flush().
    */
    const char*
    gsd_find_matching_chunk_name(struct gsd_handle* handle, const char* match, const char* prev);

    /** Upgrade a GSD file to the latest specification.

        @param handle Handle to an open GSD file

        @pre *handle* was opened by gsd_open() with a writable mode.
        @pre There are no pending data to write to the file in gsd_end_frame()

        @return
          - GSD_SUCCESS (0) on success. Negative value on failure:
          - GSD_ERROR_IO: IO error (check errno).
          - GSD_ERROR_INVALID_ARGUMENT: *handle* is NULL
          - GSD_ERROR_FILE_MUST_BE_WRITABLE: The file was opened in read-only mode.
    */
    int gsd_upgrade(struct gsd_handle* handle);

    /** Get the maximum write buffer size.

        @param handle Handle to an open GSD file

        @pre *handle* was opened by gsd_open().

        @return The maximum write buffer size in bytes, or 0 on error.
    */
    uint64_t gsd_get_maximum_write_buffer_size(struct gsd_handle* handle);

    /** Set the maximum write buffer size.

        @param handle Handle to an open GSD file
        @param size Maximum number of bytes to allocate in the write buffer (must be greater than
        0).

        @pre *handle* was opened by gsd_open().

        @return
          - GSD_SUCCESS (0) on success. Negative value on failure:
          - GSD_ERROR_INVALID_ARGUMENT: *handle* is NULL
          - GSD_ERROR_INVALID_ARGUMENT: size == 0
    */
    int gsd_set_maximum_write_buffer_size(struct gsd_handle* handle, uint64_t size);

#ifdef __cplusplus
    }
#endif

#endif // #ifndef GSD_H
