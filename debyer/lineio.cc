//  debyer -- program for calculation of diffration patterns
//  Copyright (C) 2006-2007 Marcin Wojdyr
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  $Id$

#include "lineio.h"
#include <cstring>

using namespace std;

LineInput::LineInput()
    : line_number(0), stream(NULL)
#ifdef HAVE_ZLIB
      , gz_stream(NULL)
#endif
#ifdef HAVE_BZLIB
      , bz_stream(NULL)
#endif
{
    buffer = new char[buffer_size];
}

LineInput::~LineInput()
{
    if (stream && stream != stdin)
        fclose(stream);
#ifdef HAVE_ZLIB
    if (gz_stream)
        gzclose(gz_stream);
#endif
#ifdef HAVE_BZLIB
    if (bz_stream)
        BZ2_bzclose(bz_stream);
#endif
    delete [] buffer;
}

// open file and fill buffer
bool LineInput::init(const char* filename_)
{
    assert(!stream);
#ifdef HAVE_ZLIB
    assert(!gz_stream);
#endif
#ifdef HAVE_BZLIB
    assert(!bz_stream);
#endif
    assert(filename_);

    filename = filename_;

    size_t len = strlen(filename_);
    bool gzipped = (len > 3 && strcmp(filename_ + len - 3, ".gz") == 0);
    bool bz2ed = (len > 4 && strcmp(filename_ + len - 4, ".bz2") == 0);

    if (gzipped)
        orig_filename = filename.substr(0, len - 3);
    else if (bz2ed)
        orig_filename = filename.substr(0, len - 4);
    else
        orig_filename = filename;

    if (gzipped) {
#ifdef HAVE_ZLIB
        gz_stream = gzopen(filename_, "rb");
        if (!gz_stream) {
            failure("Can not gzopen file: ", filename_);
            return false;
        }
#else
        failure("Program is compiled with disabled zlib support.");
        return false;
#endif //HAVE_ZLIB
    }
    else if (bz2ed) {
#ifdef HAVE_BZLIB
        bz_stream = BZ2_bzopen(filename_, "rb");
        if (!bz_stream) {
            failure("Can not BZ2_bzopen file: ", filename_);
            return false;
        }
#else
        failure("Program is compiled with disabled bzlib support.");
        return false;
#endif //HAVE_BZLIB
    }
    else if (len == 0 || strcmp(filename_, "-") == 0) {
        stream = stdin;
    }
    else {
        stream = fopen(filename_, "rb");
        if (!stream) {
            failure("Can not open file: ", filename_);
            return false;
        }
    }

    // read the first buffer_size bytes into buffer
    int n = fill_buffer(0);
    if (n < 0) {
        failure("Reading file failed.");
        return false;
    }
    next_line = buffer;

    return true;
}

int LineInput::fill_buffer(size_t offset)
{
    int n = -1;
    if (stream)
        n = fread(buffer+offset, 1, buffer_size-offset-1, stream);
#ifdef HAVE_ZLIB
    else if (gz_stream)
        n = gzread(gz_stream, buffer+offset, buffer_size-offset-1);
#endif
#ifdef HAVE_BZLIB
    else if (bz_stream)
        n = BZ2_bzread(bz_stream, buffer+offset, buffer_size-offset-1);
#endif
    if (n >= 0)
        buffer[offset + n] = 0;
    return n;
}

char* LineInput::get_line()
{
    if (next_line == NULL) {
        buffer[0] = 0;
        return NULL;
    }

    char *p = strchr(next_line, '\n');

    if (p) {
    // full next line was found in buffer
        *p = 0;
        ++line_number;
        char* line = next_line;
        next_line = p+1;
        return line;
    }

    int len = strlen(next_line);

    if (next_line + len != buffer + buffer_size - 1) {
    // 0 found before the end of the buffer, it's EOF
        if (*next_line == 0) {
            buffer[0] = 0;
            return NULL;
        }
        ++line_number;
        char* line = next_line;
        next_line = NULL;
        return line;
    }

    // need to refill buffer
    memmove(buffer, next_line, len);
    int n = fill_buffer(len);
    if (n < 0) { // error
        sprintf(buffer, "File reading error.\n");
        return NULL;
    }

    p = strchr(buffer+len, '\n');
    if (p) {
        *p = 0;
        next_line = p+1;
    }
    else if (strlen(buffer) == buffer_size - 1) {
        sprintf(buffer, "Line in the file is too long (>%d)\n", buffer_size-1);
        return NULL;
    }
    else {
        // it's EOF now
        buffer[0] = 0;
        next_line = NULL;
    }
    ++line_number;
    return buffer;
}

