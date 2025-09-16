# BSD 3-Clause License
#
# Copyright (c) 2025, Spill-Tea
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""Common string handling utilities to shuttle between C and python."""

from libc.string cimport memcpy
from libc.stdlib cimport free, malloc

cdef extern from "Python.h":
    Py_ssize_t PyUnicode_GET_LENGTH(object)
    bytes PyUnicode_AsUTF8String(object)
    Py_ssize_t PyBytes_GET_SIZE(object)
    char* PyBytes_AS_STRING(object)
    str PyUnicode_DecodeUTF8Stateful(char*, Py_ssize_t, char*, Py_ssize_t*)
    bytes PyBytes_FromStringAndSize(char*, Py_ssize_t)


ctypedef struct StringView:
    unsigned char* ptr
    Py_ssize_t size
    bint origin


cdef inline StringView construct(bytes s, Py_ssize_t length, bint isbytes):
    """Construct the StringView from a python bytes object."""
    cdef:
        char* buffer = PyBytes_AS_STRING(s)
        StringView view

    view.ptr = <unsigned char *> malloc((length + 1) * sizeof(unsigned char))
    memcpy(view.ptr, buffer, length + 1)
    view.ptr[length] = "\0"  # c string terminator
    view.size = length
    view.origin = isbytes

    return view


cdef inline StringView bytes_to_view(bytes b):
    """Construct StringView from python bytes object."""
    cdef Py_ssize_t length = PyBytes_GET_SIZE(b)

    return construct(b, length, True)


cdef inline StringView str_to_view(str s):
    """Construct StringView from python string object."""
    cdef:
        Py_ssize_t length = PyUnicode_GET_LENGTH(s)
        bytes temp = PyUnicode_AsUTF8String(s)

    return construct(temp, length, False)


cdef inline str to_str(StringView view):
    """Convert StringView back into a python string object, safely releasing memory."""
    cdef str obj = PyUnicode_DecodeUTF8Stateful(<char*> view.ptr, view.size, NULL, NULL)
    free(view.ptr)

    return obj


cdef inline bytes to_bytes(StringView view):
    """Convert StringView back into a python bytes object, safely releasing memory."""
    cdef bytes obj = PyBytes_FromStringAndSize(<char*> view.ptr, view.size)
    free(view.ptr)

    return obj


cdef inline object to_object(StringView view):
    """Convert StringView back into a python object, safely releasing memory."""
    if view.origin:
        return to_bytes(view)

    return to_str(view)
