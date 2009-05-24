/*  debyer -- program for calculation of diffration patterns
 *  Copyright (C) 2006-2007 Marcin Wojdyr
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  $Id: fileio.h 82 2008-10-15 19:30:45Z wojdyr $
 *
 *  class LineInput -- see description below
 */

#ifndef DEBYER_LINEIO_H_
#define DEBYER_LINEIO_H_

#include <cassert>
#include <string>
#include <cstdio>

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#ifdef HAVE_ZLIB
#  include <zlib.h>
#endif

#ifdef HAVE_BZLIB
#  include <bzlib.h>
#endif


// Line-oriented (fgets()-based) API for reading files.
// Provides random access to the first 1kb of the file (which is buffered),
// to allow guess file format.
// Can handle normal files, stdin, gzipped and bzip2-ed files.
class LineInput
{
public:
    static const int buffer_size = 1024;

    LineInput();
    ~LineInput();

    // if filename == "-", stdin is read
    // on error returns false and error message is stored in buffer
    bool init(const char* filename_);

    // after the first get_line() call, get_buffer() can't be called
    const char* get_buffer() { assert(line_number == 0); return buffer; }

    const char* get_error() { return buffer; }

    // returns pointer to 0-terminated array (without a new line character).
    // The string can be changed.
    const char* get_line();

    std::string const& get_filename() const { return filename; }
    std::string const& get_orig_filename() const { return orig_filename; }
    int get_line_number() const { return line_number; }

private:
    std::string filename;
    std::string orig_filename;
    int line_number;
    char *buffer;
    char* next_line;
    FILE *stream;
#ifdef HAVE_ZLIB
    gzFile gz_stream;
#endif
#ifdef HAVE_BZLIB
    BZFILE* bz_stream;
#endif

    int fill_buffer(size_t offset);

    // failure in init, store error message in buffer
    void failure(const char *msg, const char* fn=NULL)
    {
        snprintf(buffer, buffer_size, "%s%s\n", msg, (fn ? fn : ""));
    }
};

#endif // DEBYER_LINEIO_H_
