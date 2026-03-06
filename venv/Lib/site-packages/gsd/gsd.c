// Copyright (c) 2016-2025 The Regents of the University of Michigan
// Part of GSD, released under the BSD 2-Clause License.

#include <sys/stat.h>
#ifdef _WIN32

#pragma warning(push)
#pragma warning(disable : 4996)

#define GSD_USE_MMAP 0
#define WIN32_LEAN_AND_MEAN
#include <io.h>
#include <windows.h>

#else // linux / mac

#define _XOPEN_SOURCE 500
#include <sys/mman.h>
#include <unistd.h>
#define GSD_USE_MMAP 1

#endif

#ifdef __APPLE__
#include <limits.h>
#endif

#include <errno.h>
#include <fcntl.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "gsd.h"

/** @file gsd.c
    @brief Implements the GSD C API
*/

/// Magic value identifying a GSD file
const uint64_t GSD_MAGIC_ID = 0x65DF65DF65DF65DF;

/// Initial index size
enum
    {
    GSD_INITIAL_INDEX_SIZE = 128
    };

/// Initial namelist size
enum
    {
    GSD_INITIAL_NAME_BUFFER_SIZE = 1024
    };

/// Size of initial frame index
enum
    {
    GSD_INITIAL_FRAME_INDEX_SIZE = 16
    };

/// Initial size of write buffer
enum
    {
    GSD_INITIAL_WRITE_BUFFER_SIZE = 1024
    };

/// Default maximum size of write buffer (in bytes)
enum
    {
    GSD_DEFAULT_MAXIMUM_WRITE_BUFFER_BYTES = 1024 * 1024
    };

/// Size of hash map
enum
    {
    GSD_NAME_MAP_SIZE = 57557
    };

/// Current GSD file specification
enum
    {
    GSD_CURRENT_FILE_VERSION_MAJOR = 2
    };

enum
    {
    GSD_CURRENT_FILE_VERSION_MINOR = 1
    };

// define windows wrapper functions
#ifdef _WIN32
#define lseek _lseeki64
#define ftruncate _chsize
typedef int64_t ssize_t;

int S_IRUSR = _S_IREAD;
int S_IWUSR = _S_IWRITE;
int S_IRGRP = _S_IREAD;
int S_IWGRP = _S_IWRITE;

int gsd_fsync(int fd)
    {
    return _commit(fd);
    }

inline ssize_t pread(int fd, void* buf, size_t count, int64_t offset)
    {
    // Note: _read only accepts unsigned int values
    if (count > UINT_MAX)
        return GSD_ERROR_IO;

    int64_t oldpos = _telli64(fd);
    _lseeki64(fd, offset, SEEK_SET);
    ssize_t result = _read(fd, buf, (unsigned int)count);
    _lseeki64(fd, oldpos, SEEK_SET);
    return result;
    }

inline ssize_t pwrite(int fd, const void* buf, size_t count, int64_t offset)
    {
    // Note: _write only accepts unsigned int values
    if (count > UINT_MAX)
        return GSD_ERROR_IO;

    int64_t oldpos = _telli64(fd);
    _lseeki64(fd, offset, SEEK_SET);
    ssize_t result = _write(fd, buf, (unsigned int)count);
    _lseeki64(fd, oldpos, SEEK_SET);
    return result;
    }

#elif defined(__APPLE__)
int gsd_fsync(int fd)
    {
    return fcntl(fd, F_FULLFSYNC);
    }
#else
int gsd_fsync(int fd)
    {
    return fsync(fd);
    }
#endif

/** Zero memory

    @param d pointer to memory region
    @param size_to_zero size of the area to zero in bytes
*/
inline static void gsd_util_zero_memory(void* d, size_t size_to_zero)
    {
    memset(d, 0, size_to_zero);
    }

/** @internal
    @brief Write large data buffer to file

    The system call pwrite() fails to write very large data buffers. This method calls pwrite() as
    many times as necessary to completely write a large buffer.

    @param fd File descriptor.
    @param buf Data buffer.
    @param count Number of bytes to write.
    @param offset Location in the file to start writing.

    @returns The total number of bytes written or a negative value on error.
*/
inline static ssize_t gsd_io_pwrite_retry(int fd, const void* buf, size_t count, int64_t offset)
    {
    size_t total_bytes_written = 0;
    const char* ptr = (char*)buf;

    // perform multiple pwrite calls to complete a large write successfully
    while (total_bytes_written < count)
        {
        size_t to_write = count - total_bytes_written;
#if defined(_WIN32) || defined(__APPLE__)
        // win32 and apple raise an error for writes greater than INT_MAX
        if (to_write > INT_MAX / 2)
            {
            to_write = INT_MAX / 2;
            }
#endif

        errno = 0;
        ssize_t bytes_written
            = pwrite(fd, ptr + total_bytes_written, to_write, offset + total_bytes_written);

        if (bytes_written == -1 || (bytes_written == 0 && errno != 0))
            {
            return GSD_ERROR_IO;
            }

        total_bytes_written += bytes_written;
        }

    return total_bytes_written;
    }

/** @internal
    @brief Read large data buffer to file

    The system call pread() fails to read very large data buffers. This method calls pread() as many
    times as necessary to completely read a large buffer.

    @param fd File descriptor.
    @param buf Data buffer.
    @param count Number of bytes to read.
    @param offset Location in the file to start reading.

    @returns The total number of bytes read or a negative value on error.
*/
inline static ssize_t gsd_io_pread_retry(int fd, void* buf, size_t count, int64_t offset)
    {
    size_t total_bytes_read = 0;
    char* ptr = (char*)buf;

    // perform multiple pread calls to complete a large write successfully
    while (total_bytes_read < count)
        {
        size_t to_read = count - total_bytes_read;
#if defined(_WIN32) || defined(__APPLE__)
        // win32 and apple raise errors for reads greater than INT_MAX
        if (to_read > INT_MAX / 2)
            {
            to_read = INT_MAX / 2;
            }
#endif

        errno = 0;
        ssize_t bytes_read = pread(fd, ptr + total_bytes_read, to_read, offset + total_bytes_read);
        if (bytes_read == -1 || (bytes_read == 0 && errno != 0))
            {
            return GSD_ERROR_IO;
            }

        total_bytes_read += bytes_read;

        // handle end of file
        if (bytes_read == 0)
            {
            return total_bytes_read;
            }
        }

    return total_bytes_read;
    }

/** @internal
    @brief Allocate a name/id map

    @param map Map to allocate.
    @param size Number of entries in the map.

    @returns GSD_SUCCESS on success, GSD_* error codes on error.
*/
inline static int gsd_name_id_map_allocate(struct gsd_name_id_map* map, size_t size)
    {
    if (map == NULL || map->v || size == 0 || map->size != 0)
        {
        return GSD_ERROR_INVALID_ARGUMENT;
        }

    map->v = calloc(size, sizeof(struct gsd_name_id_pair));
    if (map->v == NULL)
        {
        return GSD_ERROR_MEMORY_ALLOCATION_FAILED;
        }

    map->size = size;

    return GSD_SUCCESS;
    }

/** @internal
    @brief Free a name/id map

    @param map Map to free.

    @returns GSD_SUCCESS on success, GSD_* error codes on error.
*/
inline static int gsd_name_id_map_free(struct gsd_name_id_map* map)
    {
    if (map == NULL || map->v == NULL || map->size == 0)
        {
        return GSD_ERROR_INVALID_ARGUMENT;
        }

    // free all of the linked lists
    size_t i;
    for (i = 0; i < map->size; i++)
        {
        free(map->v[i].name);

        struct gsd_name_id_pair* cur = map->v[i].next;
        while (cur != NULL)
            {
            struct gsd_name_id_pair* prev = cur;
            cur = cur->next;
            free(prev->name);
            free(prev);
            }
        }

    // free the main map
    free(map->v);

    map->v = 0;
    map->size = 0;

    return GSD_SUCCESS;
    }

/** @internal
    @brief Hash a string

    @param str String to hash

    @returns Hashed value of the string.
*/
inline static unsigned long gsd_hash_str(const unsigned char* str)
    {
    unsigned long hash = 5381; // NOLINT
    int c;

    while ((c = *str++))
        {
        hash = ((hash << 5) + hash) + c; /* hash * 33 + c NOLINT */
        }

    return hash;
    }

/** @internal
    @brief Insert a string into a name/id map

    @param map Map to insert into.
    @param str String to insert.
    @param id ID to associate with the string.

    @returns GSD_SUCCESS on success, GSD_* error codes on error.
*/
inline static int gsd_name_id_map_insert(struct gsd_name_id_map* map, const char* str, uint16_t id)
    {
    if (map == NULL || map->v == NULL || map->size == 0)
        {
        return GSD_ERROR_INVALID_ARGUMENT;
        }

    size_t hash = gsd_hash_str((const unsigned char*)str) % map->size;

    // base case: no conflict
    if (map->v[hash].name == NULL)
        {
        map->v[hash].name = calloc(strlen(str) + 1, sizeof(char));
        if (map->v[hash].name == NULL)
            {
            return GSD_ERROR_MEMORY_ALLOCATION_FAILED;
            }
        memcpy(map->v[hash].name, str, strlen(str) + 1);
        map->v[hash].id = id;
        map->v[hash].next = NULL;
        }
    else
        {
        // go to the end of the conflict list
        struct gsd_name_id_pair* insert_point = map->v + hash;

        while (insert_point->next != NULL)
            {
            insert_point = insert_point->next;
            }

        // allocate and insert a new entry
        insert_point->next = malloc(sizeof(struct gsd_name_id_pair));
        if (insert_point->next == NULL)
            {
            return GSD_ERROR_MEMORY_ALLOCATION_FAILED;
            }

        insert_point->next->name = calloc(strlen(str) + 1, sizeof(char));
        if (insert_point->next->name == NULL)
            {
            return GSD_ERROR_MEMORY_ALLOCATION_FAILED;
            }
        memcpy(insert_point->next->name, str, strlen(str) + 1);
        insert_point->next->id = id;
        insert_point->next->next = NULL;
        }

    return GSD_SUCCESS;
    }

/** @internal
    @brief Find an ID in a name/id mapping

    @param map Map to search.
    @param str String to search.

    @returns The ID if found, or UINT16_MAX if not found.
*/
inline static uint16_t gsd_name_id_map_find(struct gsd_name_id_map* map, const char* str)
    {
    if (map == NULL || map->v == NULL || map->size == 0)
        {
        return UINT16_MAX;
        }

    size_t hash = gsd_hash_str((const unsigned char*)str) % map->size;

    struct gsd_name_id_pair* cur = map->v + hash;

    while (cur != NULL)
        {
        if (cur->name == NULL)
            {
            // not found
            return UINT16_MAX;
            }

        if (strcmp(str, cur->name) == 0)
            {
            // found
            return cur->id;
            }

        // keep looking
        cur = cur->next;
        }

    // not found in any conflict
    return UINT16_MAX;
    }

/** @internal
    @brief Utility function to validate index entry
    @param handle handle to the open gsd file
    @param idx index of entry to validate

    @returns 1 if the entry is valid, 0 if it is not
*/
inline static int gsd_is_entry_valid(struct gsd_handle* handle, size_t idx)
    {
    const struct gsd_index_entry entry = handle->file_index.data[idx];

    // check for valid type
    if (gsd_sizeof_type((enum gsd_type)entry.type) == 0)
        {
        return 0;
        }

    // validate that we don't read past the end of the file
    size_t size = entry.N * entry.M * gsd_sizeof_type((enum gsd_type)entry.type);
    if ((entry.location + size) > (uint64_t)handle->file_size)
        {
        return 0;
        }

    // check for valid frame (frame cannot be more than the number of index entries)
    if (entry.frame >= handle->header.index_allocated_entries)
        {
        return 0;
        }

    // check for valid id
    if (entry.id >= (handle->file_names.n_names + handle->frame_names.n_names))
        {
        return 0;
        }

    // check for valid flags
    if (entry.flags != 0)
        {
        return 0;
        }

    return 1;
    }

/** @internal
    @brief Allocate a write buffer

    @param buf Buffer to allocate.
    @param reserve Number of bytes to allocate.

    @returns GSD_SUCCESS on success, GSD_* error codes on error.
*/
inline static int gsd_byte_buffer_allocate(struct gsd_byte_buffer* buf, size_t reserve)
    {
    if (buf == NULL || buf->data || reserve == 0 || buf->reserved != 0 || buf->size != 0)
        {
        return GSD_ERROR_INVALID_ARGUMENT;
        }

    buf->data = calloc(reserve, sizeof(char));
    if (buf->data == NULL)
        {
        return GSD_ERROR_MEMORY_ALLOCATION_FAILED;
        }

    buf->size = 0;
    buf->reserved = reserve;

    return GSD_SUCCESS;
    }

/** @internal
    @brief Append bytes to a byte buffer

    @param buf Buffer to append to.
    @param data Data to append.
    @param size Number of bytes in *data*.

    @returns GSD_SUCCESS on success, GSD_* error codes on error.
*/
inline static int gsd_byte_buffer_append(struct gsd_byte_buffer* buf, const char* data, size_t size)
    {
    if (buf == NULL || buf->data == NULL || size == 0 || buf->reserved == 0)
        {
        return GSD_ERROR_INVALID_ARGUMENT;
        }

    if (buf->size + size > buf->reserved)
        {
        // reallocate by doubling
        size_t new_reserved = buf->reserved * 2;
        while (buf->size + size >= new_reserved)
            {
            new_reserved = new_reserved * 2;
            }

        char* old_data = buf->data;
        // NOLINTNEXTLINE(bugprone-suspicious-realloc-usage): realloc is used correctly
        buf->data = realloc(buf->data, sizeof(char) * new_reserved);
        if (buf->data == NULL)
            {
            // this free should not be necessary, but clang-tidy disagrees
            free(old_data);
            return GSD_ERROR_MEMORY_ALLOCATION_FAILED;
            }

        // zero the new memory, but only the portion after the end of the new section to be appended
        gsd_util_zero_memory(buf->data + (buf->size + size),
                             sizeof(char) * (new_reserved - (buf->size + size)));
        buf->reserved = new_reserved;
        }

    memcpy(buf->data + buf->size, data, size);
    buf->size += size;

    return GSD_SUCCESS;
    }

/** @internal
    @brief Free the memory allocated by the write buffer or unmap the mapped memory.

    @param buf Buffer to free.

    @returns GSD_SUCCESS on success, GSD_* error codes on error.
*/
inline static int gsd_byte_buffer_free(struct gsd_byte_buffer* buf)
    {
    if (buf == NULL || buf->data == NULL)
        {
        return GSD_ERROR_INVALID_ARGUMENT;
        }

    free(buf->data);

    gsd_util_zero_memory(buf, sizeof(struct gsd_byte_buffer));
    return GSD_SUCCESS;
    }

/** @internal
    @brief Allocate a buffer of index entries

    @param buf Buffer to allocate.
    @param reserve Number of entries to allocate.

    @post The buffer's data element has *reserve* elements allocated in memory.

    @returns GSD_SUCCESS on success, GSD_* error codes on error.
*/
inline static int gsd_index_buffer_allocate(struct gsd_index_buffer* buf, size_t reserve)
    {
    if (buf == NULL || buf->mapped_data || buf->data || reserve == 0 || buf->reserved != 0
        || buf->size != 0)
        {
        return GSD_ERROR_INVALID_ARGUMENT;
        }

    buf->data = calloc(reserve, sizeof(struct gsd_index_entry));
    if (buf->data == NULL)
        {
        return GSD_ERROR_MEMORY_ALLOCATION_FAILED;
        }

    buf->size = 0;
    buf->reserved = reserve;
    buf->mapped_data = NULL;
    buf->mapped_len = 0;

    return GSD_SUCCESS;
    }

/** @internal
    @brief Map index entries from the file

    @param buf Buffer to map.
    @param handle GSD file handle to map.

    @post The buffer's data element contains the index data from the file.

    On some systems, this will use mmap to efficiently access the file. On others, it may result in
    an allocation and read of the entire index from the file.

    @returns GSD_SUCCESS on success, GSD_* error codes on error.
*/
inline static int gsd_index_buffer_map(struct gsd_index_buffer* buf, struct gsd_handle* handle)
    {
    if (buf == NULL || buf->mapped_data || buf->data || buf->reserved != 0 || buf->size != 0)
        {
        return GSD_ERROR_INVALID_ARGUMENT;
        }

    // validate that the index block exists inside the file
    if (handle->header.index_location
            + sizeof(struct gsd_index_entry) * handle->header.index_allocated_entries
        > (uint64_t)handle->file_size)
        {
        return GSD_ERROR_FILE_CORRUPT;
        }

#if GSD_USE_MMAP
    // map the index in read only mode
    size_t index_size = sizeof(struct gsd_index_entry) * handle->header.index_allocated_entries;
    buf->mapped_data = mmap(NULL,
                            index_size + handle->header.index_location,
                            PROT_READ,
                            MAP_SHARED,
                            handle->fd,
                            0);

    if (buf->mapped_data == MAP_FAILED)
        {
        return GSD_ERROR_IO;
        }

    buf->data
        = (struct gsd_index_entry*)(((char*)buf->mapped_data) + handle->header.index_location);

    buf->mapped_len = index_size + handle->header.index_location;
    buf->reserved = handle->header.index_allocated_entries;
#else
    // mmap not supported, read the data from the disk
    int retval = gsd_index_buffer_allocate(buf, handle->header.index_allocated_entries);
    if (retval != GSD_SUCCESS)
        {
        return retval;
        }

    ssize_t bytes_read = gsd_io_pread_retry(handle->fd,
                                            buf->data,
                                            sizeof(struct gsd_index_entry)
                                                * handle->header.index_allocated_entries,
                                            handle->header.index_location);

    if (bytes_read == -1
        || bytes_read != sizeof(struct gsd_index_entry) * handle->header.index_allocated_entries)
        {
        return GSD_ERROR_IO;
        }
#endif

    // determine the number of index entries in the list
    // file is corrupt if first index entry is invalid
    if (buf->data[0].location != 0 && !gsd_is_entry_valid(handle, 0))
        {
        return GSD_ERROR_FILE_CORRUPT;
        }

    if (buf->data[0].location == 0)
        {
        buf->size = 0;
        }
    else
        {
        // determine the number of index entries (marked by location = 0)
        // binary search for the first index entry with location 0
        size_t L = 0;
        size_t R = buf->reserved;

        // progressively narrow the search window by halves
        do
            {
            size_t m = (L + R) / 2;

            // file is corrupt if any index entry is invalid or frame does not increase
            // monotonically
            if (buf->data[m].location != 0
                && (!gsd_is_entry_valid(handle, m) || buf->data[m].frame < buf->data[L].frame))
                {
                return GSD_ERROR_FILE_CORRUPT;
                }

            if (buf->data[m].location != 0)
                {
                L = m;
                }
            else
                {
                R = m;
                }
            } while ((R - L) > 1);

        // this finds R = the first index entry with location = 0
        buf->size = R;
        }

    return GSD_SUCCESS;
    }

/** @internal
    @brief Free the memory allocated by the index buffer or unmap the mapped memory.

    @param buf Buffer to free.

    @returns GSD_SUCCESS on success, GSD_* error codes on error.
*/
inline static int gsd_index_buffer_free(struct gsd_index_buffer* buf)
    {
    if (buf == NULL || buf->data == NULL)
        {
        return GSD_ERROR_INVALID_ARGUMENT;
        }

#if GSD_USE_MMAP
    if (buf->mapped_data)
        {
        int retval = munmap(buf->mapped_data, buf->mapped_len);

        if (retval != 0)
            {
            return GSD_ERROR_IO;
            }
        }
    else
#endif
        {
        free(buf->data);
        }

    gsd_util_zero_memory(buf, sizeof(struct gsd_index_buffer));
    return GSD_SUCCESS;
    }

/** @internal
    @brief Add a new index entry and provide a pointer to it.

    @param buf Buffer to add too.
    @param entry [out] Pointer to set to the new entry.

    Double the size of the reserved space as needed to hold the new entry. Does not accept mapped
    indices.

    @returns GSD_SUCCESS on success, GSD_* error codes on error.
*/
inline static int gsd_index_buffer_add(struct gsd_index_buffer* buf, struct gsd_index_entry** entry)
    {
    if (buf == NULL || buf->mapped_data || entry == NULL || buf->reserved == 0)
        {
        return GSD_ERROR_INVALID_ARGUMENT;
        }

    if (buf->size == buf->reserved)
        {
        // grow the array
        size_t new_reserved = buf->reserved * 2;
        // NOLINTNEXTLINE(bugprone-suspicious-realloc-usage): realloc is used correctly
        buf->data = realloc(buf->data, sizeof(struct gsd_index_entry) * new_reserved);
        if (buf->data == NULL)
            {
            return GSD_ERROR_MEMORY_ALLOCATION_FAILED;
            }

        // zero the new memory
        gsd_util_zero_memory(buf->data + buf->reserved,
                             sizeof(struct gsd_index_entry) * (new_reserved - buf->reserved));
        buf->reserved = new_reserved;
        }

    size_t insert_pos = buf->size;
    buf->size++;
    *entry = buf->data + insert_pos;

    return GSD_SUCCESS;
    }

inline static int gsd_cmp_index_entry(const struct gsd_index_entry* a,
                                      const struct gsd_index_entry* b)
    {
    int result = 0;

    if (a->frame < b->frame)
        {
        result = -1;
        }

    if (a->frame > b->frame)
        {
        result = 1;
        }

    if (a->frame == b->frame)
        {
        if (a->id < b->id)
            {
            result = -1;
            }

        if (a->id > b->id)
            {
            result = 1;
            }

        if (a->id == b->id)
            {
            result = 0;
            }
        }

    return result;
    }

/** @internal
    @brief Compute heap parent node.
    @param i Node index.
*/
inline static size_t gsd_heap_parent(size_t i)
    {
    return (i - 1) / 2;
    }

/** @internal
    @brief Compute heap left child.
    @param i Node index.
*/
inline static size_t gsd_heap_left_child(size_t i)
    {
    return 2 * i + 1;
    }

/** @internal
    @brief Swap the nodes *a* and *b* in the buffer
    @param buf Buffer.
    @param a First index to swap.
    @param b Second index to swap.
*/
inline static void gsd_heap_swap(struct gsd_index_buffer* buf, size_t a, size_t b)
    {
    struct gsd_index_entry tmp = buf->data[a];
    buf->data[a] = buf->data[b];
    buf->data[b] = tmp;
    }

/** @internal
    @brief Shift heap node downward
    @param buf Buffer.
    @param start First index of the valid heap in *buf*.
    @param end Last index of the valid hep in *buf*.
*/
inline static void gsd_heap_shift_down(struct gsd_index_buffer* buf, size_t start, size_t end)
    {
    size_t root = start;

    while (gsd_heap_left_child(root) <= end)
        {
        size_t child = gsd_heap_left_child(root);
        size_t swap = root;

        if (gsd_cmp_index_entry(buf->data + swap, buf->data + child) < 0)
            {
            swap = child;
            }
        if (child + 1 <= end && gsd_cmp_index_entry(buf->data + swap, buf->data + child + 1) < 0)
            {
            swap = child + 1;
            }

        if (swap == root)
            {
            return;
            }

        gsd_heap_swap(buf, root, swap);
        root = swap;
        }
    }

/** @internal
    @brief Convert unordered index buffer to a heap
    @param buf Buffer.
*/
inline static void gsd_heapify(struct gsd_index_buffer* buf)
    {
    ssize_t start = gsd_heap_parent(buf->size - 1);

    while (start >= 0)
        {
        gsd_heap_shift_down(buf, start, buf->size - 1);
        start--;
        }
    }

/** @internal
    @brief Sort the index buffer.

    @param buf Buffer to sort.

    Sorts an in-memory index buffer. Does not accept mapped indices.

    @returns GSD_SUCCESS on success, GSD_* error codes on error.
*/
inline static int gsd_index_buffer_sort(struct gsd_index_buffer* buf)
    {
    if (buf == NULL || buf->mapped_data || buf->reserved == 0)
        {
        return GSD_ERROR_INVALID_ARGUMENT;
        }

    // arrays of size 0 or 1 are already sorted
    if (buf->size <= 1)
        {
        return GSD_SUCCESS;
        }

    gsd_heapify(buf);

    size_t end = buf->size - 1;
    while (end > 0)
        {
        gsd_heap_swap(buf, end, 0);
        end = end - 1;
        gsd_heap_shift_down(buf, 0, end);
        }

    return GSD_SUCCESS;
    }

/** @internal
    @brief Utility function to expand the memory space for the index block in the file.

    @param handle Handle to the open gsd file.
    @param size_required The new index must be able to hold at least this many elements.

    @returns GSD_SUCCESS on success, GSD_* error codes on error.
*/
inline static int gsd_expand_file_index(struct gsd_handle* handle, size_t size_required)
    {
    if (handle->open_flags == GSD_OPEN_READONLY)
        {
        return GSD_ERROR_FILE_MUST_BE_WRITABLE;
        }

    // multiply the index size each time it grows
    // this allows the index to grow rapidly to accommodate new frames
    const int multiplication_factor = 2;

    // save the old size and update the new size
    size_t size_old = handle->header.index_allocated_entries;
    size_t size_new = size_old * multiplication_factor;

    while (size_new <= size_required)
        {
        size_new *= multiplication_factor;
        }

    // Mac systems deadlock when writing from a mapped region into the tail end of that same region
    // unmap the index first and copy it over by chunks
    int retval = gsd_index_buffer_free(&handle->file_index);
    if (retval != 0)
        {
        return retval;
        }

    int64_t new_index_location = lseek(handle->fd, 0, SEEK_END);
    int64_t old_index_location = handle->header.index_location;

    // fill the new index space with 0s
    size_t new_index_bytes = size_new * sizeof(struct gsd_index_entry);
    retval = ftruncate(handle->fd, new_index_location + new_index_bytes);
    if (retval != 0)
        {
        return GSD_ERROR_IO;
        }

    // write the current index to the end of the file
    uint64_t copy_buffer_size = GSD_DEFAULT_MAXIMUM_WRITE_BUFFER_BYTES;
    if (copy_buffer_size > size_old * sizeof(struct gsd_index_entry))
        {
        copy_buffer_size = size_old * sizeof(struct gsd_index_entry);
        }
    char* buf = malloc(copy_buffer_size);

    size_t total_bytes_written = 0;
    size_t old_index_bytes = size_old * sizeof(struct gsd_index_entry);
    while (total_bytes_written < old_index_bytes)
        {
        size_t bytes_to_copy = copy_buffer_size;
        if (old_index_bytes - total_bytes_written < copy_buffer_size)
            {
            bytes_to_copy = old_index_bytes - total_bytes_written;
            }

        ssize_t bytes_read = gsd_io_pread_retry(handle->fd,
                                                buf,
                                                bytes_to_copy,
                                                old_index_location + total_bytes_written);

        if (bytes_read == -1 || bytes_read != bytes_to_copy)
            {
            free(buf);
            return GSD_ERROR_IO;
            }

        ssize_t bytes_written = gsd_io_pwrite_retry(handle->fd,
                                                    buf,
                                                    bytes_to_copy,
                                                    new_index_location + total_bytes_written);

        if (bytes_written == -1 || bytes_written != bytes_to_copy)
            {
            free(buf);
            return GSD_ERROR_IO;
            }

        total_bytes_written += bytes_written;
        }

    // sync the expanded index
    retval = gsd_fsync(handle->fd);
    if (retval != 0)
        {
        free(buf);
        return GSD_ERROR_IO;
        }

    // free the copy buffer
    free(buf);

    // update the header
    handle->header.index_location = new_index_location;
    handle->file_size = new_index_location + new_index_bytes;
    handle->header.index_allocated_entries = size_new;

    // write the new header out
    ssize_t bytes_written
        = gsd_io_pwrite_retry(handle->fd, &(handle->header), sizeof(struct gsd_header), 0);
    if (bytes_written != sizeof(struct gsd_header))
        {
        return GSD_ERROR_IO;
        }

    // sync the updated header
    retval = gsd_fsync(handle->fd);
    if (retval != 0)
        {
        return GSD_ERROR_IO;
        }

    // remap the file index
    retval = gsd_index_buffer_map(&handle->file_index, handle);
    if (retval != 0)
        {
        return retval;
        }

    return GSD_SUCCESS;
    }

/** @internal
    @brief Flush the write buffer.

    gsd_write_frame() writes small data chunks into the write buffer. It adds index entries for
    these chunks to gsd_handle::buffer_index with locations offset from the current end of file.
    gsd_flush_write_buffer() writes the buffer to the end of the file. The entries in the buffered
    index have valid locations because nothing else writes past the end of the file.

    @param handle Handle to flush the write buffer.
    @returns GSD_SUCCESS on success or GSD_* error codes on error
*/
inline static int gsd_flush_write_buffer(struct gsd_handle* handle)
    {
    if (handle == NULL)
        {
        return GSD_ERROR_INVALID_ARGUMENT;
        }

    if (handle->write_buffer.size == 0)
        {
        // nothing to do
        return GSD_SUCCESS;
        }

    // write the buffer to the end of the file
    uint64_t offset = handle->file_size;
    ssize_t bytes_written = gsd_io_pwrite_retry(handle->fd,
                                                handle->write_buffer.data,
                                                handle->write_buffer.size,
                                                offset);

    if (bytes_written == -1 || bytes_written != handle->write_buffer.size)
        {
        return GSD_ERROR_IO;
        }

    handle->file_size += handle->write_buffer.size;

    // reset write_buffer for new data
    handle->write_buffer.size = 0;

    return GSD_SUCCESS;
    }

/** @internal
    @brief Flush the name buffer.

    gsd_write_frame() adds new names to the frame_names buffer. gsd_flush_name_buffer() flushes
    this buffer at the end of a frame write and commits the new names to the file. If necessary,
    the namelist is written to a new location in the file.

    @param handle Handle to flush the write buffer.
    @returns GSD_SUCCESS on success or GSD_* error codes on error
*/
inline static int gsd_flush_name_buffer(struct gsd_handle* handle)
    {
    if (handle == NULL)
        {
        return GSD_ERROR_INVALID_ARGUMENT;
        }

    if (handle->frame_names.n_names == 0)
        {
        // nothing to do
        return GSD_SUCCESS;
        }

    if (handle->frame_names.data.size == 0)
        {
        // error: bytes in buffer, but no names for them
        return GSD_ERROR_INVALID_ARGUMENT;
        }

    size_t old_reserved = handle->file_names.data.reserved;
    size_t old_size = handle->file_names.data.size;

    // add the new names to the file name list and zero the frame list
    int retval = gsd_byte_buffer_append(&handle->file_names.data,
                                        handle->frame_names.data.data,
                                        handle->frame_names.data.size);
    if (retval != GSD_SUCCESS)
        {
        return retval;
        }

    handle->file_names.n_names += handle->frame_names.n_names;
    handle->frame_names.n_names = 0;
    handle->frame_names.data.size = 0;
    gsd_util_zero_memory(handle->frame_names.data.data, handle->frame_names.data.reserved);

    // reserved space must be a multiple of the GSD name size
    if (handle->file_names.data.reserved % GSD_NAME_SIZE != 0)
        {
        return GSD_ERROR_INVALID_ARGUMENT;
        }

    if (handle->file_names.data.reserved > old_reserved)
        {
        // write the new name list to the end of the file
        uint64_t offset = handle->file_size;
        ssize_t bytes_written = gsd_io_pwrite_retry(handle->fd,
                                                    handle->file_names.data.data,
                                                    handle->file_names.data.reserved,
                                                    offset);

        if (bytes_written == -1 || bytes_written != handle->file_names.data.reserved)
            {
            return GSD_ERROR_IO;
            }

        // sync the updated name list
        retval = gsd_fsync(handle->fd);
        if (retval != 0)
            {
            return GSD_ERROR_IO;
            }

        handle->file_size += handle->file_names.data.reserved;
        handle->header.namelist_location = offset;
        handle->header.namelist_allocated_entries
            = handle->file_names.data.reserved / GSD_NAME_SIZE;

        // write the new header out
        bytes_written
            = gsd_io_pwrite_retry(handle->fd, &(handle->header), sizeof(struct gsd_header), 0);
        if (bytes_written != sizeof(struct gsd_header))
            {
            return GSD_ERROR_IO;
            }
        }
    else
        {
        // write the new name list to the old index location
        uint64_t offset = handle->header.namelist_location;
        ssize_t bytes_written = gsd_io_pwrite_retry(handle->fd,
                                                    handle->file_names.data.data + old_size,
                                                    handle->file_names.data.reserved - old_size,
                                                    offset + old_size);
        if (bytes_written != (handle->file_names.data.reserved - old_size))
            {
            return GSD_ERROR_IO;
            }
        }

    // sync the updated name list or header
    retval = gsd_fsync(handle->fd);
    if (retval != 0)
        {
        return GSD_ERROR_IO;
        }

    return GSD_SUCCESS;
    }

/** @internal
    @brief utility function to append a name to the namelist

    @param id [out] ID of the new name
    @param handle handle to the open gsd file
    @param name string name

    Append a name to the names in the current frame. gsd_end_frame() will add this list to the
    file names.

    @return
      - GSD_SUCCESS (0) on success. Negative value on failure:
      - GSD_ERROR_IO: IO error (check errno).
      - GSD_ERROR_MEMORY_ALLOCATION_FAILED: Unable to allocate memory.
      - GSD_ERROR_FILE_MUST_BE_WRITABLE: File must not be read only.
*/
inline static int gsd_append_name(uint16_t* id, struct gsd_handle* handle, const char* name)
    {
    if (handle->open_flags == GSD_OPEN_READONLY)
        {
        return GSD_ERROR_FILE_MUST_BE_WRITABLE;
        }

    if (handle->file_names.n_names + handle->frame_names.n_names == UINT16_MAX)
        {
        // no more names may be added
        return GSD_ERROR_NAMELIST_FULL;
        }

    // Provide the ID of the new name
    *id = (uint16_t)(handle->file_names.n_names + handle->frame_names.n_names);

    if (handle->header.gsd_version < gsd_make_version(2, 0))
        {
        // v1 files always allocate GSD_NAME_SIZE bytes for each name and put a NULL terminator
        // at address 63
        char name_v1[GSD_NAME_SIZE];
        strncpy(name_v1, name, GSD_NAME_SIZE - 1);
        name_v1[GSD_NAME_SIZE - 1] = 0;
        gsd_byte_buffer_append(&handle->frame_names.data, name_v1, GSD_NAME_SIZE);
        handle->frame_names.n_names++;

        // update the name/id mapping with the truncated name
        int retval = gsd_name_id_map_insert(&handle->name_map, name_v1, *id);
        if (retval != GSD_SUCCESS)
            {
            return retval;
            }
        }
    else
        {
        gsd_byte_buffer_append(&handle->frame_names.data, name, strlen(name) + 1);
        handle->frame_names.n_names++;

        // update the name/id mapping
        int retval = gsd_name_id_map_insert(&handle->name_map, name, *id);
        if (retval != GSD_SUCCESS)
            {
            return retval;
            }
        }

    return GSD_SUCCESS;
    }

/** @internal
    @brief Cross-platform wrapper for the POSIX open() system function.
    @param pathname file path using UTF-8 encoding on all platforms
    @return file descriptor
*/
inline static int gsd_open_file(const char* pathname, int flags, int mode)
    {
#ifndef _WIN32
    return open(pathname, flags, mode);
#else
    // On Windows, we call the _wopen() function, which requires converting the UTF-8 input path to
    // UTF-16 wide-character encoding.
    int count_wchars;
    wchar_t* wpathname;
    int fd;

    // First, determine the number of wide characters needed to represent the input string.
    count_wchars = MultiByteToWideChar(CP_UTF8, 0, pathname, -1, NULL, 0);
    // Then allocate temporary wchar_t buffer and perform the string conversion.
    wpathname = malloc(sizeof(wchar_t) * count_wchars);
    MultiByteToWideChar(CP_UTF8, 0, pathname, -1, wpathname, count_wchars);
    fd = _wopen(wpathname, flags, mode);
    free(wpathname);
    return fd;
#endif
    }

/** @internal
    @brief Truncate the file and write a new gsd header.

    @param fd file descriptor to initialize
    @param application Generating application name (truncated to 63 chars)
    @param schema Schema name for data to be written in this GSD file (truncated to 63 chars)
    @param schema_version Version of the scheme data to be written (make with gsd_make_version())
*/
inline static int
gsd_initialize_file(int fd, const char* application, const char* schema, uint32_t schema_version)
    {
    // check if the file was created
    if (fd == -1)
        {
        return GSD_ERROR_IO;
        }

    int retval = ftruncate(fd, 0);
    if (retval != 0)
        {
        return GSD_ERROR_IO;
        }

    // populate header fields
    struct gsd_header header;
    gsd_util_zero_memory(&header, sizeof(header));

    header.magic = GSD_MAGIC_ID;
    header.gsd_version
        = gsd_make_version(GSD_CURRENT_FILE_VERSION_MAJOR, GSD_CURRENT_FILE_VERSION_MINOR);
    strncpy(header.application, application, sizeof(header.application) - 1);
    header.application[sizeof(header.application) - 1] = 0;
    strncpy(header.schema, schema, sizeof(header.schema) - 1);
    header.schema[sizeof(header.schema) - 1] = 0;
    header.schema_version = schema_version;
    header.index_location = sizeof(header);
    header.index_allocated_entries = GSD_INITIAL_INDEX_SIZE;
    header.namelist_location
        = header.index_location + sizeof(struct gsd_index_entry) * header.index_allocated_entries;
    header.namelist_allocated_entries = GSD_INITIAL_NAME_BUFFER_SIZE / GSD_NAME_SIZE;
    gsd_util_zero_memory(header.reserved, sizeof(header.reserved));

    // write the header out
    ssize_t bytes_written = gsd_io_pwrite_retry(fd, &header, sizeof(header), 0);
    if (bytes_written != sizeof(header))
        {
        return GSD_ERROR_IO;
        }

    // allocate and zero default index memory
    struct gsd_index_entry index[GSD_INITIAL_INDEX_SIZE];
    gsd_util_zero_memory(index, sizeof(index));

    // write the empty index out
    bytes_written = gsd_io_pwrite_retry(fd, index, sizeof(index), sizeof(header));
    if (bytes_written != sizeof(index))
        {
        return GSD_ERROR_IO;
        }

    // allocate and zero the namelist memory
    char names[GSD_INITIAL_NAME_BUFFER_SIZE];
    gsd_util_zero_memory(names, sizeof(char) * GSD_INITIAL_NAME_BUFFER_SIZE);

    // write the namelist out
    bytes_written = gsd_io_pwrite_retry(fd, names, sizeof(names), sizeof(header) + sizeof(index));
    if (bytes_written != sizeof(names))
        {
        return GSD_ERROR_IO;
        }

    // sync file
    retval = gsd_fsync(fd);
    if (retval != 0)
        {
        return GSD_ERROR_IO;
        }

    return GSD_SUCCESS;
    }

/** @internal
    @brief Read in the file index and initialize the handle.

    @param handle Handle to read the header

    @pre handle->fd is an open file.
    @pre handle->open_flags is set.
*/
inline static int gsd_initialize_handle(struct gsd_handle* handle)
    {
    // check if the file was created
    if (handle->fd == -1)
        {
        return GSD_ERROR_IO;
        }

    // read the header
    ssize_t bytes_read
        = gsd_io_pread_retry(handle->fd, &handle->header, sizeof(struct gsd_header), 0);
    if (bytes_read == -1)
        {
        return GSD_ERROR_IO;
        }
    if (bytes_read != sizeof(struct gsd_header))
        {
        return GSD_ERROR_NOT_A_GSD_FILE;
        }

    // validate the header
    if (handle->header.magic != GSD_MAGIC_ID)
        {
        return GSD_ERROR_NOT_A_GSD_FILE;
        }

    if (handle->header.gsd_version < gsd_make_version(1, 0)
        && handle->header.gsd_version != gsd_make_version(0, 3))
        {
        return GSD_ERROR_INVALID_GSD_FILE_VERSION;
        }

    if (handle->header.gsd_version >= gsd_make_version(3, 0))
        {
        return GSD_ERROR_INVALID_GSD_FILE_VERSION;
        }

    // determine the file size
    handle->file_size = lseek(handle->fd, 0, SEEK_END);

    // validate that the namelist block exists inside the file
    if (handle->header.namelist_location
            + (GSD_NAME_SIZE * handle->header.namelist_allocated_entries)
        > (uint64_t)handle->file_size)
        {
        return GSD_ERROR_FILE_CORRUPT;
        }

    // allocate the hash map
    int retval = gsd_name_id_map_allocate(&handle->name_map, GSD_NAME_MAP_SIZE);
    if (retval != GSD_SUCCESS)
        {
        return retval;
        }

    // read the namelist block
    size_t namelist_n_bytes = GSD_NAME_SIZE * handle->header.namelist_allocated_entries;
    retval = gsd_byte_buffer_allocate(&handle->file_names.data, namelist_n_bytes);
    if (retval != GSD_SUCCESS)
        {
        return retval;
        }
    bytes_read = gsd_io_pread_retry(handle->fd,
                                    handle->file_names.data.data,
                                    namelist_n_bytes,
                                    handle->header.namelist_location);

    if (bytes_read == -1 || bytes_read != namelist_n_bytes)
        {
        return GSD_ERROR_IO;
        }

    // The name buffer must end in a NULL terminator or else the file is corrupt
    if (handle->file_names.data.data[handle->file_names.data.reserved - 1] != 0)
        {
        return GSD_ERROR_FILE_CORRUPT;
        }

    // Add the names to the hash map. Also determine the number of used bytes in the namelist.
    size_t name_start = 0;
    handle->file_names.n_names = 0;
    while (name_start < handle->file_names.data.reserved)
        {
        char* name = handle->file_names.data.data + name_start;

        // an empty name notes the end of the list
        if (name[0] == 0)
            {
            break;
            }

        retval
            = gsd_name_id_map_insert(&handle->name_map, name, (uint16_t)handle->file_names.n_names);
        if (retval != GSD_SUCCESS)
            {
            return retval;
            }
        handle->file_names.n_names++;

        if (handle->header.gsd_version < gsd_make_version(2, 0))
            {
            // gsd v1 stores names in fixed 64 byte segments
            name_start += GSD_NAME_SIZE;
            }
        else
            {
            size_t len = strnlen(name, handle->file_names.data.reserved - name_start);
            name_start += len + 1;
            }
        }

    handle->file_names.data.size = name_start;

    // read in the file index
    retval = gsd_index_buffer_map(&handle->file_index, handle);
    if (retval != GSD_SUCCESS)
        {
        return retval;
        }

    // determine the current frame counter
    if (handle->file_index.size == 0)
        {
        handle->file_frame = 0;
        }
    else
        {
        handle->file_frame = handle->file_index.data[handle->file_index.size - 1].frame + 1;
        }
    handle->buffer_frame = handle->file_frame;

    // if this is a write mode, allocate the initial frame index and the name buffer
    if (handle->open_flags != GSD_OPEN_READONLY)
        {
        retval = gsd_index_buffer_allocate(&handle->buffer_index, GSD_INITIAL_FRAME_INDEX_SIZE);
        if (retval != GSD_SUCCESS)
            {
            return retval;
            }

        retval = gsd_byte_buffer_allocate(&handle->write_buffer, GSD_INITIAL_WRITE_BUFFER_SIZE);
        if (retval != GSD_SUCCESS)
            {
            return retval;
            }

        handle->frame_names.n_names = 0;
        retval = gsd_byte_buffer_allocate(&handle->frame_names.data, GSD_NAME_SIZE);
        if (retval != GSD_SUCCESS)
            {
            return retval;
            }
        }

    handle->pending_index_entries = 0;
    handle->maximum_write_buffer_size = GSD_DEFAULT_MAXIMUM_WRITE_BUFFER_BYTES;

    // Silently upgrade writable files from a previous matching major version to the latest
    // minor version.
    if ((handle->open_flags == GSD_OPEN_READWRITE || handle->open_flags == GSD_OPEN_APPEND)
        && (handle->header.gsd_version
            <= gsd_make_version(GSD_CURRENT_FILE_VERSION_MAJOR, GSD_CURRENT_FILE_VERSION_MINOR))
        && (handle->header.gsd_version >> (sizeof(uint32_t) * 4) == GSD_CURRENT_FILE_VERSION_MAJOR))
        {
        handle->header.gsd_version
            = gsd_make_version(GSD_CURRENT_FILE_VERSION_MAJOR, GSD_CURRENT_FILE_VERSION_MINOR);
        size_t bytes_written
            = gsd_io_pwrite_retry(handle->fd, &(handle->header), sizeof(struct gsd_header), 0);

        if (bytes_written != sizeof(struct gsd_header))
            {
            return GSD_ERROR_IO;
            }
        }

    return GSD_SUCCESS;
    }

uint32_t gsd_make_version(unsigned int major, unsigned int minor)
    {
    return major << (sizeof(uint32_t) * 4) | minor;
    }

int gsd_create(const char* fname,
               const char* application,
               const char* schema,
               uint32_t schema_version)
    {
    int extra_flags = 0;
#ifdef _WIN32
    extra_flags = _O_BINARY;
#endif

    // create the file
    int fd = gsd_open_file(fname,
                           O_RDWR | O_CREAT | O_TRUNC | extra_flags,
                           S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP);
    int retval = gsd_initialize_file(fd, application, schema, schema_version);
    if (fd != -1)
        {
        close(fd);
        }
    return retval;
    }

int gsd_create_and_open(struct gsd_handle* handle,
                        const char* fname,
                        const char* application,
                        const char* schema,
                        uint32_t schema_version,
                        const enum gsd_open_flag flags,
                        int exclusive_create)
    {
    // zero the handle
    gsd_util_zero_memory(handle, sizeof(struct gsd_handle));

    int extra_flags = 0;
#ifdef _WIN32
    extra_flags = _O_BINARY;
#endif

    // set the open flags in the handle
    if (flags == GSD_OPEN_READWRITE)
        {
        handle->open_flags = GSD_OPEN_READWRITE;
        }
    else if (flags == GSD_OPEN_READONLY)
        {
        return GSD_ERROR_FILE_MUST_BE_WRITABLE;
        }
    else if (flags == GSD_OPEN_APPEND)
        {
        handle->open_flags = GSD_OPEN_APPEND;
        }

    // set the exclusive create bit
    if (exclusive_create)
        {
        extra_flags |= O_EXCL;
        }

    // create the file
    handle->fd = gsd_open_file(fname,
                               O_RDWR | O_CREAT | O_TRUNC | extra_flags,
                               S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP);
    int retval = gsd_initialize_file(handle->fd, application, schema, schema_version);
    if (retval != 0)
        {
        if (handle->fd != -1)
            {
            close(handle->fd);
            }
        return retval;
        }

    retval = gsd_initialize_handle(handle);
    if (retval != 0)
        {
        if (handle->fd != -1)
            {
            close(handle->fd);
            }
        }
    return retval;
    }

int gsd_open(struct gsd_handle* handle, const char* fname, const enum gsd_open_flag flags)
    {
    // zero the handle
    gsd_util_zero_memory(handle, sizeof(struct gsd_handle));

    int extra_flags = 0;
#ifdef _WIN32
    extra_flags = _O_BINARY;
#endif

    // open the file
    if (flags == GSD_OPEN_READWRITE)
        {
        handle->fd = gsd_open_file(fname, O_RDWR | extra_flags, 0);
        handle->open_flags = GSD_OPEN_READWRITE;
        }
    else if (flags == GSD_OPEN_READONLY)
        {
        handle->fd = gsd_open_file(fname, O_RDONLY | extra_flags, 0);
        handle->open_flags = GSD_OPEN_READONLY;
        }
    else if (flags == GSD_OPEN_APPEND)
        {
        handle->fd = gsd_open_file(fname, O_RDWR | extra_flags, 0);
        handle->open_flags = GSD_OPEN_APPEND;
        }

    int retval = gsd_initialize_handle(handle);
    if (retval != 0)
        {
        if (handle->fd != -1)
            {
            close(handle->fd);
            }
        }
    return retval;
    }

int gsd_truncate(struct gsd_handle* handle)
    {
    if (handle == NULL)
        {
        return GSD_ERROR_INVALID_ARGUMENT;
        }
    if (handle->open_flags == GSD_OPEN_READONLY)
        {
        return GSD_ERROR_FILE_MUST_BE_WRITABLE;
        }

    int retval = 0;

    // deallocate indices
    if (handle->frame_names.data.reserved > 0)
        {
        retval = gsd_byte_buffer_free(&handle->frame_names.data);
        if (retval != GSD_SUCCESS)
            {
            return retval;
            }
        }

    if (handle->file_names.data.reserved > 0)
        {
        retval = gsd_byte_buffer_free(&handle->file_names.data);
        if (retval != GSD_SUCCESS)
            {
            return retval;
            }
        }

    retval = gsd_name_id_map_free(&handle->name_map);
    if (retval != GSD_SUCCESS)
        {
        return retval;
        }

    retval = gsd_index_buffer_free(&handle->file_index);
    if (retval != GSD_SUCCESS)
        {
        return retval;
        }

    if (handle->buffer_index.reserved > 0)
        {
        retval = gsd_index_buffer_free(&handle->buffer_index);
        if (retval != GSD_SUCCESS)
            {
            return retval;
            }
        }

    if (handle->write_buffer.reserved > 0)
        {
        retval = gsd_byte_buffer_free(&handle->write_buffer);
        if (retval != GSD_SUCCESS)
            {
            return retval;
            }
        }

    // keep a copy of the old header
    struct gsd_header old_header = handle->header;
    retval = gsd_initialize_file(handle->fd,
                                 old_header.application,
                                 old_header.schema,
                                 old_header.schema_version);

    if (retval != GSD_SUCCESS)
        {
        return retval;
        }

    return gsd_initialize_handle(handle);
    }

int gsd_close(struct gsd_handle* handle)
    {
    if (handle == NULL)
        {
        return GSD_ERROR_INVALID_ARGUMENT;
        }

    int retval;

    if (handle->open_flags != GSD_OPEN_READONLY)
        {
        retval = gsd_flush(handle);
        if (retval != GSD_SUCCESS)
            {
            return retval;
            }
        }

    // save the fd so we can use it after freeing the handle
    int fd = handle->fd;

    retval = gsd_index_buffer_free(&handle->file_index);
    if (retval != GSD_SUCCESS)
        {
        return retval;
        }

    if (handle->buffer_index.reserved > 0)
        {
        retval = gsd_index_buffer_free(&handle->buffer_index);
        if (retval != GSD_SUCCESS)
            {
            return retval;
            }
        }

    if (handle->write_buffer.reserved > 0)
        {
        retval = gsd_byte_buffer_free(&handle->write_buffer);
        if (retval != GSD_SUCCESS)
            {
            return retval;
            }
        }

    retval = gsd_name_id_map_free(&handle->name_map);
    if (retval != GSD_SUCCESS)
        {
        return retval;
        }

    if (handle->frame_names.data.reserved > 0)
        {
        handle->frame_names.n_names = 0;
        retval = gsd_byte_buffer_free(&handle->frame_names.data);
        if (retval != GSD_SUCCESS)
            {
            return retval;
            }
        }

    if (handle->file_names.data.reserved > 0)
        {
        handle->file_names.n_names = 0;
        retval = gsd_byte_buffer_free(&handle->file_names.data);
        if (retval != GSD_SUCCESS)
            {
            return retval;
            }
        }

    // close the file
    retval = close(fd);
    if (retval != 0)
        {
        return GSD_ERROR_IO;
        }

    return GSD_SUCCESS;
    }

int gsd_end_frame(struct gsd_handle* handle)
    {
    if (handle == NULL)
        {
        return GSD_ERROR_INVALID_ARGUMENT;
        }
    if (handle->open_flags == GSD_OPEN_READONLY)
        {
        return GSD_ERROR_FILE_MUST_BE_WRITABLE;
        }

    handle->buffer_frame++;
    handle->pending_index_entries = 0;

    return GSD_SUCCESS;
    }

int gsd_flush(struct gsd_handle* handle)
    {
    if (handle == NULL)
        {
        return GSD_ERROR_INVALID_ARGUMENT;
        }
    if (handle->open_flags == GSD_OPEN_READONLY)
        {
        return GSD_ERROR_FILE_MUST_BE_WRITABLE;
        }

    // flush the write buffer first so that all locations in the buffer_index are maintained.
    int retval = gsd_flush_write_buffer(handle);
    if (retval != GSD_SUCCESS)
        {
        return retval;
        }

    // flush the namelist buffer
    retval = gsd_flush_name_buffer(handle);
    if (retval != GSD_SUCCESS)
        {
        return retval;
        }

    // sync the data before writing the index
    retval = gsd_fsync(handle->fd);
    if (retval != 0)
        {
        return GSD_ERROR_IO;
        }

    // Write the buffer index to the file, excluding the index entries that are part of the current
    // frame.
    if (handle->pending_index_entries > handle->buffer_index.size)
        {
        return GSD_ERROR_INVALID_ARGUMENT;
        }
    uint64_t index_entries_to_write = handle->buffer_index.size - handle->pending_index_entries;

    if (index_entries_to_write > 0)
        {
        // ensure there is enough space in the index
        if ((handle->file_index.size + index_entries_to_write) > handle->file_index.reserved)
            {
            gsd_expand_file_index(handle, handle->file_index.size + index_entries_to_write);
            }

        // sort the index before writing
        retval = gsd_index_buffer_sort(&handle->buffer_index);
        if (retval != 0)
            {
            return retval;
            }

        // write the frame index entries to the file
        int64_t write_pos = handle->header.index_location
                            + sizeof(struct gsd_index_entry) * handle->file_index.size;

        size_t bytes_to_write = sizeof(struct gsd_index_entry) * index_entries_to_write;
        ssize_t bytes_written
            = gsd_io_pwrite_retry(handle->fd, handle->buffer_index.data, bytes_to_write, write_pos);

        if (bytes_written == -1 || bytes_written != bytes_to_write)
            {
            return GSD_ERROR_IO;
            }

        retval = gsd_fsync(handle->fd);
        if (retval != 0)
            {
            return GSD_ERROR_IO;
            }

#if !GSD_USE_MMAP
        // add the entries to the file index
        memcpy(handle->file_index.data + handle->file_index.size,
               handle->buffer_index.data,
               sizeof(struct gsd_index_entry) * index_entries_to_write);
#endif

        // update size of file index
        handle->file_index.size += index_entries_to_write;

        // Clear the frame index, keeping those in the current unfinished frame.
        if (handle->pending_index_entries > 0)
            {
            for (uint64_t i = 0; i < handle->pending_index_entries; i++)
                {
                handle->buffer_index.data[i]
                    = handle->buffer_index
                          .data[handle->buffer_index.size - handle->pending_index_entries + i];
                }
            }

        handle->buffer_index.size = handle->pending_index_entries;
        }

    handle->file_frame = handle->buffer_frame;

    return GSD_SUCCESS;
    }

int gsd_write_chunk(struct gsd_handle* handle,
                    const char* name,
                    enum gsd_type type,
                    uint64_t N,
                    uint32_t M,
                    uint8_t flags,
                    const void* data)
    {
    // validate input
    if (N > 0 && data == NULL)
        {
        return GSD_ERROR_INVALID_ARGUMENT;
        }
    if (M == 0)
        {
        return GSD_ERROR_INVALID_ARGUMENT;
        }
    if (handle->open_flags == GSD_OPEN_READONLY)
        {
        return GSD_ERROR_FILE_MUST_BE_WRITABLE;
        }
    if (flags != 0)
        {
        return GSD_ERROR_INVALID_ARGUMENT;
        }

    uint16_t id = gsd_name_id_map_find(&handle->name_map, name);
    if (id == UINT16_MAX)
        {
        // not found, append to the index
        int retval = gsd_append_name(&id, handle, name);
        if (retval != GSD_SUCCESS)
            {
            return retval;
            }

        if (id == UINT16_MAX)
            {
            // this should never happen
            return GSD_ERROR_NAMELIST_FULL;
            }
        }

    struct gsd_index_entry entry;
    // populate fields in the entry's data
    gsd_util_zero_memory(&entry, sizeof(struct gsd_index_entry));
    entry.frame = handle->buffer_frame;
    entry.id = id;
    entry.type = (uint8_t)type;
    entry.N = N;
    entry.M = M;
    entry.location = handle->file_size + handle->write_buffer.size;
    size_t size = N * M * gsd_sizeof_type(type);

    // flush the buffer if this entry won't fit
    if (size > (handle->maximum_write_buffer_size - handle->write_buffer.size))
        {
        int retval = gsd_flush_write_buffer(handle);
        if (retval != GSD_SUCCESS)
            {
            return retval;
            }
        }

    // add an entry to the buffer index
    struct gsd_index_entry* index_entry;
    int retval = gsd_index_buffer_add(&handle->buffer_index, &index_entry);
    if (retval != GSD_SUCCESS)
        {
        return retval;
        }
    *index_entry = entry;

    // decide whether to write this chunk to the buffer or straight to disk
    if (size < handle->maximum_write_buffer_size)
        {
        // add the data to the write buffer
        if (size > 0)
            {
            retval = gsd_byte_buffer_append(&handle->write_buffer, data, size);
            if (retval != GSD_SUCCESS)
                {
                return retval;
                }
            }
        }
    else
        {
        // write the data
        ssize_t bytes_written = gsd_io_pwrite_retry(handle->fd, data, size, entry.location);
        if (bytes_written == -1 || bytes_written != size)
            {
            return GSD_ERROR_IO;
            }

        // update the file_size in the handle
        handle->file_size += bytes_written;
        }

    handle->pending_index_entries++;
    return GSD_SUCCESS;
    }

uint64_t gsd_get_nframes(struct gsd_handle* handle)
    {
    if (handle == NULL)
        {
        return 0;
        }
    return handle->file_frame;
    }

const struct gsd_index_entry*
gsd_find_chunk(struct gsd_handle* handle, uint64_t frame, const char* name)
    {
    if (handle == NULL)
        {
        return NULL;
        }
    if (name == NULL)
        {
        return NULL;
        }
    if (frame >= gsd_get_nframes(handle))
        {
        return NULL;
        }

    // find the id for the given name
    uint16_t match_id = gsd_name_id_map_find(&handle->name_map, name);
    if (match_id == UINT16_MAX)
        {
        return NULL;
        }

    if (handle->header.gsd_version >= gsd_make_version(2, 0))
        {
        // gsd 2.0 files sort the entire index
        // binary search for the index entry
        ssize_t L = 0;
        ssize_t R = handle->file_index.size - 1;
        struct gsd_index_entry T;
        T.frame = frame;
        T.id = match_id;

        while (L <= R)
            {
            size_t m = (L + R) / 2;
            int cmp = gsd_cmp_index_entry(handle->file_index.data + m, &T);
            if (cmp == -1)
                {
                L = m + 1;
                }
            else if (cmp == 1)
                {
                R = m - 1;
                }
            else
                {
                return &(handle->file_index.data[m]);
                }
            }
        }
    else
        {
        // gsd 1.0 file: use binary search to find the frame and linear search to find the entry
        size_t L = 0;
        size_t R = handle->file_index.size;

        // progressively narrow the search window by halves
        do
            {
            size_t m = (L + R) / 2;

            if (frame < handle->file_index.data[m].frame)
                {
                R = m;
                }
            else
                {
                L = m;
                }
            } while ((R - L) > 1);

        // this finds L = the rightmost index with the desired frame
        int64_t cur_index;

        // search all index entries with the matching frame
        for (cur_index = L; (cur_index >= 0) && (handle->file_index.data[cur_index].frame == frame);
             cur_index--)
            {
            // if the frame matches, check the id
            if (match_id == handle->file_index.data[cur_index].id)
                {
                return &(handle->file_index.data[cur_index]);
                }
            }
        }

    // if we got here, we didn't find the specified chunk
    return NULL;
    }

int gsd_read_chunk(struct gsd_handle* handle, void* data, const struct gsd_index_entry* chunk)
    {
    if (handle == NULL)
        {
        return GSD_ERROR_INVALID_ARGUMENT;
        }
    if (data == NULL)
        {
        return GSD_ERROR_INVALID_ARGUMENT;
        }
    if (chunk == NULL)
        {
        return GSD_ERROR_INVALID_ARGUMENT;
        }

    size_t size = chunk->N * chunk->M * gsd_sizeof_type((enum gsd_type)chunk->type);
    if (size == 0)
        {
        return GSD_ERROR_FILE_CORRUPT;
        }
    if (chunk->location == 0)
        {
        return GSD_ERROR_FILE_CORRUPT;
        }

    // validate that we don't read past the end of the file
    if ((chunk->location + size) > (uint64_t)handle->file_size)
        {
        return GSD_ERROR_FILE_CORRUPT;
        }

    ssize_t bytes_read = gsd_io_pread_retry(handle->fd, data, size, chunk->location);
    if (bytes_read == -1 || bytes_read != size)
        {
        return GSD_ERROR_IO;
        }

    return GSD_SUCCESS;
    }

size_t gsd_sizeof_type(enum gsd_type type)
    {
    size_t val = 0;
    if (type == GSD_TYPE_UINT8)
        {
        val = sizeof(uint8_t);
        }
    else if (type == GSD_TYPE_UINT16)
        {
        val = sizeof(uint16_t);
        }
    else if (type == GSD_TYPE_UINT32)
        {
        val = sizeof(uint32_t);
        }
    else if (type == GSD_TYPE_UINT64)
        {
        val = sizeof(uint64_t);
        }
    else if (type == GSD_TYPE_INT8)
        {
        val = sizeof(int8_t);
        }
    else if (type == GSD_TYPE_INT16)
        {
        val = sizeof(int16_t);
        }
    else if (type == GSD_TYPE_INT32)
        {
        val = sizeof(int32_t);
        }
    else if (type == GSD_TYPE_INT64)
        {
        val = sizeof(int64_t);
        }
    else if (type == GSD_TYPE_FLOAT)
        {
        val = sizeof(float);
        }
    else if (type == GSD_TYPE_DOUBLE)
        {
        val = sizeof(double);
        }
    else if (type == GSD_TYPE_CHARACTER)
        {
        val = sizeof(char);
        }
    else
        {
        return 0;
        }
    return val;
    }

const char*
gsd_find_matching_chunk_name(struct gsd_handle* handle, const char* match, const char* prev)
    {
    if (handle == NULL)
        {
        return NULL;
        }
    if (match == NULL)
        {
        return NULL;
        }
    if (handle->file_names.n_names == 0)
        {
        return NULL;
        }

    // return nothing found if the name buffer is corrupt
    if (handle->file_names.data.data[handle->file_names.data.reserved - 1] != 0)
        {
        return NULL;
        }

    // determine search start index
    const char* search_str;
    if (prev == NULL)
        {
        search_str = handle->file_names.data.data;
        }
    else
        {
        // return not found if prev is not in range
        if (prev < handle->file_names.data.data)
            {
            return NULL;
            }
        if (prev >= (handle->file_names.data.data + handle->file_names.data.reserved))
            {
            return NULL;
            }

        if (handle->header.gsd_version < gsd_make_version(2, 0))
            {
            search_str = prev + GSD_NAME_SIZE;
            }
        else
            {
            search_str = prev + strlen(prev) + 1;
            }
        }

    size_t match_len = strlen(match);

    while (search_str < (handle->file_names.data.data + handle->file_names.data.reserved))
        {
        if (search_str[0] != 0 && 0 == strncmp(match, search_str, match_len))
            {
            return search_str;
            }

        if (handle->header.gsd_version < gsd_make_version(2, 0))
            {
            search_str += GSD_NAME_SIZE;
            }
        else
            {
            search_str += strlen(search_str) + 1;
            }
        }

    // searched past the end of the list, return NULL
    return NULL;
    }

int gsd_upgrade(struct gsd_handle* handle)
    {
    if (handle == NULL)
        {
        return GSD_ERROR_INVALID_ARGUMENT;
        }
    if (handle->open_flags == GSD_OPEN_READONLY)
        {
        return GSD_ERROR_INVALID_ARGUMENT;
        }
    if (handle->buffer_index.size > 0 || handle->frame_names.n_names > 0)
        {
        return GSD_ERROR_INVALID_ARGUMENT;
        }

    if (handle->header.gsd_version < gsd_make_version(2, 0))
        {
        if (handle->file_index.size > 0)
            {
            // make a copy of the file index
            struct gsd_index_buffer buf;
            gsd_util_zero_memory(&buf, sizeof(struct gsd_index_buffer));
            int retval = gsd_index_buffer_allocate(&buf, handle->file_index.size);
            if (retval != GSD_SUCCESS)
                {
                return retval;
                }
            memcpy(buf.data,
                   handle->file_index.data,
                   sizeof(struct gsd_index_entry) * handle->file_index.size);
            buf.size = handle->file_index.size;

            // sort the copy and write it back out to the file
            retval = gsd_index_buffer_sort(&buf);
            if (retval != GSD_SUCCESS)
                {
                gsd_index_buffer_free(&buf);
                return retval;
                }

            ssize_t bytes_written = gsd_io_pwrite_retry(handle->fd,
                                                        buf.data,
                                                        sizeof(struct gsd_index_entry) * buf.size,
                                                        handle->header.index_location);

            if (bytes_written == -1 || bytes_written != sizeof(struct gsd_index_entry) * buf.size)
                {
                gsd_index_buffer_free(&buf);
                return GSD_ERROR_IO;
                }

            retval = gsd_index_buffer_free(&buf);
            if (retval != GSD_SUCCESS)
                {
                return retval;
                }

            // sync the updated index
            retval = gsd_fsync(handle->fd);
            if (retval != 0)
                {
                return GSD_ERROR_IO;
                }
            }

        if (handle->file_names.n_names > 0)
            {
            // compact the name list without changing its size or position on the disk
            struct gsd_byte_buffer new_name_buf;
            gsd_util_zero_memory(&new_name_buf, sizeof(struct gsd_byte_buffer));
            int retval = gsd_byte_buffer_allocate(&new_name_buf, handle->file_names.data.reserved);
            if (retval != GSD_SUCCESS)
                {
                return retval;
                }

            const char* name = gsd_find_matching_chunk_name(handle, "", NULL);
            while (name != NULL)
                {
                retval = gsd_byte_buffer_append(&new_name_buf, name, strlen(name) + 1);
                if (retval != GSD_SUCCESS)
                    {
                    gsd_byte_buffer_free(&new_name_buf);
                    return retval;
                    }
                name = gsd_find_matching_chunk_name(handle, "", name);
                }

            if (new_name_buf.reserved != handle->file_names.data.reserved)
                {
                gsd_byte_buffer_free(&new_name_buf);
                return GSD_ERROR_FILE_CORRUPT;
                }

            // write the new names out to disk
            ssize_t bytes_written = gsd_io_pwrite_retry(handle->fd,
                                                        new_name_buf.data,
                                                        new_name_buf.reserved,
                                                        handle->header.namelist_location);

            if (bytes_written == -1 || bytes_written != new_name_buf.reserved)
                {
                gsd_byte_buffer_free(&new_name_buf);
                return GSD_ERROR_IO;
                }

            // swap in the re-organized name buffer
            retval = gsd_byte_buffer_free(&handle->file_names.data);
            if (retval != GSD_SUCCESS)
                {
                gsd_byte_buffer_free(&new_name_buf);
                return retval;
                }
            handle->file_names.data = new_name_buf;

            // sync the updated name list
            retval = gsd_fsync(handle->fd);
            if (retval != 0)
                {
                gsd_byte_buffer_free(&new_name_buf);
                return GSD_ERROR_IO;
                }
            }

        // GSD always writes files matching the current major and minor version.
        handle->header.gsd_version
            = gsd_make_version(GSD_CURRENT_FILE_VERSION_MAJOR, GSD_CURRENT_FILE_VERSION_MINOR);

        // write the new header out
        ssize_t bytes_written
            = gsd_io_pwrite_retry(handle->fd, &(handle->header), sizeof(struct gsd_header), 0);
        if (bytes_written != sizeof(struct gsd_header))
            {
            return GSD_ERROR_IO;
            }

        // sync the updated header
        int retval = gsd_fsync(handle->fd);
        if (retval != 0)
            {
            return GSD_ERROR_IO;
            }

        // remap the file index
        retval = gsd_index_buffer_free(&handle->file_index);
        if (retval != 0)
            {
            return retval;
            }

        retval = gsd_index_buffer_map(&handle->file_index, handle);
        if (retval != 0)
            {
            return retval;
            }
        }

    return GSD_SUCCESS;
    }

uint64_t gsd_get_maximum_write_buffer_size(struct gsd_handle* handle)
    {
    if (handle == NULL)
        {
        return 0;
        }
    return handle->maximum_write_buffer_size;
    }

int gsd_set_maximum_write_buffer_size(struct gsd_handle* handle, uint64_t size)
    {
    if (handle == NULL || size == 0)
        {
        return GSD_ERROR_INVALID_ARGUMENT;
        }

    handle->maximum_write_buffer_size = size;

    return GSD_SUCCESS;
    }

// undefine windows wrapper macros
#ifdef _WIN32
#undef lseek
#undef write
#undef read
#undef open
#undef ftruncate
#pragma warning(pop)

#endif
