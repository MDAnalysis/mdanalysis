/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2003 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/
/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: fastio.h,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.6 $       $Date: 2004/09/21 20:52:37 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *   This is a simple abstraction layer for system-dependent I/O calls
 * that allow plugins to do binary I/O using the fastest possible method.
 *
 * This code is intended for use by binary trajectory reader plugins that
 * work with multi-gigabyte data sets, reading only binary data.
 *
 ***************************************************************************/

#if defined(_MSC_VER)

/* Version for machines with plain old ANSI C  */

#include <stdio.h>

typedef FILE * fio_fd;
typedef size_t fio_size_t;  /* MSVC doesn't uinversally support ssize_t */
typedef void * fio_caddr_t; /* MSVC doesn't universally support caddr_t */

typedef struct {
  fio_caddr_t iov_base;
  int iov_len;
} fio_iovec;

#define FIO_READ  0x01
#define FIO_WRITE 0x02

#define FIO_SEEK_CUR SEEK_CUR
#define FIO_SEEK_SET SEEK_SET
#define FIO_SEEK_END SEEK_END

static int fio_open(const char *filename, int mode, fio_fd *fd) {
  char * modestr;
  FILE *fp;
 
  if (mode == FIO_READ) 
    modestr = "rb";

  if (mode == FIO_WRITE) 
    modestr = "wb";

  fp = fopen(filename, modestr);
  if (fp == NULL) {
    return -1;
  } else {
    *fd = fp;
    return 0;
  }
}

static int fio_fclose(fio_fd fd) {
  return fclose(fd);
}

static fio_size_t fio_fread(void *ptr, fio_size_t size, 
                            fio_size_t nitems, fio_fd fd) {
  return fread(ptr, size, nitems, fd);
}

static fio_size_t fio_readv(fio_fd fd, const fio_iovec * iov, int iovcnt) {
  int i;
  fio_size_t len = 0; 

  for (i=0; i<iovcnt; i++) {
    fio_size_t rc = fread(iov[i].iov_base, iov[i].iov_len, 1, fd);
    if (rc != 1)
      break;
    len += iov[i].iov_len;
  }

  return len;
}

static fio_size_t fio_fwrite(void *ptr, fio_size_t size, 
                             fio_size_t nitems, fio_fd fd) {
  return fwrite(ptr, size, nitems, fd);
}

static fio_size_t fio_fseek(fio_fd fd, fio_size_t offset, int whence) {
  return fseek(fd, offset, whence);
}

static fio_size_t fio_ftell(fio_fd fd) {
  return ftell(fd);
}

#else 

/* Version for UNIX machines */
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

typedef int fio_fd;
typedef ssize_t fio_size_t;

/* enable use of kernel readv() if available */
#if defined(__sun)
#define USE_KERNEL_READV 1
#endif
#define USE_KERNEL_READV 1

typedef void * fio_caddr_t;

#if defined(USE_KERNEL_READV)
#include <sys/uio.h>
typedef struct iovec fio_iovec;
#else

typedef struct {
  fio_caddr_t iov_base;
  int iov_len;
} fio_iovec;
#endif


#define FIO_READ  0x01
#define FIO_WRITE 0x02

#define FIO_SEEK_CUR SEEK_CUR
#define FIO_SEEK_SET SEEK_SET
#define FIO_SEEK_END SEEK_END

static int fio_open(const char *filename, int mode, fio_fd *fd) {
  int nfd;
  int oflag = 0;

  if (mode == FIO_READ) 
    oflag = O_RDONLY;

  if (mode == FIO_WRITE) 
    oflag = O_WRONLY | O_CREAT;

  nfd = open(filename, oflag, 0700);
  if (nfd < 0) {
    return -1;
  } else {
    *fd = nfd;
    return 0;
  }
}

static int fio_fclose(fio_fd fd) {
  return close(fd);
}

static fio_size_t fio_fread(void *ptr, fio_size_t size, 
                            fio_size_t nitems, fio_fd fd) {
  int i;
  fio_size_t len = 0; 
  int cnt = 0;

  for (i=0; i<nitems; i++) {
    fio_size_t rc = read(fd, ptr, size);
    if (rc != size)
      break;
    len += rc;
    cnt++;
  }

  return cnt;
}

static fio_size_t fio_readv(fio_fd fd, const fio_iovec * iov, int iovcnt) {
#if defined(USE_KERNEL_READV)
  return readv(fd, iov, iovcnt);
#else
  int i;
  fio_size_t len = 0; 

  for (i=0; i<iovcnt; i++) {
    fio_size_t rc = read(fd, iov[i].iov_base, iov[i].iov_len);
    if (rc != iov[i].iov_len)
      break;
    len += iov[i].iov_len;
  }

  return len;
#endif
}

static fio_size_t fio_fwrite(void *ptr, fio_size_t size, 
                             fio_size_t nitems, fio_fd fd) {
  int i;
  fio_size_t len = 0; 
  int cnt = 0;

  for (i=0; i<nitems; i++) {
    fio_size_t rc = write(fd, ptr, size);
    if (rc != size)
      break;
    len += rc;
    cnt++;
  }

  return cnt;
}

static fio_size_t fio_fseek(fio_fd fd, fio_size_t offset, int whence) {
  return lseek(fd, offset, whence);
}

static fio_size_t fio_ftell(fio_fd fd) {
  return lseek(fd, 0, SEEK_CUR);
}

#endif


/* higher level routines that are OS independent */

static int fio_write_int32(fio_fd fd, int i) {
  return (fio_fwrite(&i, 4, 1, fd) != 1);
}

static int fio_read_int32(fio_fd fd, int *i) {
  return (fio_fread(i, 4, 1, fd) != 1);
}

static int fio_write_str(fio_fd fd, const char *str) {
  int len = strlen(str);
  return (fio_fwrite((void *) str, len, 1, fd) != 1);
}

