# BSD 3-Clause License

# Copyright (c) 2025, Spill-Tea

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.

# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.

# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.

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

"""setup script to assist with compilation of extensions.

References:
    * https://cython.readthedocs.io/en/latest/src/userguide/source_files_and_compilation.html#compiler-directives
    * https://github.com/cython/cython/issues/2995
    * https://stackoverflow.com/a/58116368/16771898
    * https://medium.com/@xpl/protecting-python-sources-using-cython-dcd940bb188e
    * https://cython.readthedocs.io/en/latest/src/tutorial/parallelization.html#compilation

"""

import sys  # noqa: I001

from setuptools import Extension, setup

# NOTE: Import cython only after setuptools
from Cython.Compiler import Options
from Cython.Distutils import build_ext


# Primitive determination if package is being installed in editable mode (via pip)
if any(map(lambda x: x in sys.argv, ("editable_wheel", "-e", "--editable"))):
    Options.annotate = True

MAJOR_VERSION: str = str(sys.version_info[0])
extensions: list[Extension] = []
openmp: str = "/openmp" if sys.platform.startswith("win") else "-fopenmp"


# Oligonucleotides
extensions.append(
    Extension(
        "designer_dna._oligos",
        ["src/designer_dna/_oligos.pyx"],
        include_dirs=["src/designer_dna/headers"],
        # extra_compile_args=[openmp],
        # extra_link_args=[openmp],
    )
)
extensions.append(
    Extension(
        "designer_dna._oligonucleotides",
        ["src/designer_dna/_oligonucleotides.pyx"],
        include_dirs=["src/designer_dna/headers"],
        language="c++",
    )
)


# Add cython directive to specify python version target
directives: dict = {"language_level": MAJOR_VERSION}
for ext in extensions:
    if hasattr(ext, "cython_directives") and isinstance(ext.cython_directives, dict):
        ext.cython_directives.update(directives)
    else:
        ext.cython_directives = directives

# NOTE: Project metadata is captured from pyproject.toml
setup(
    ext_modules=extensions,
    cmdclass={"build_ext": build_ext},
)
