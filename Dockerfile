FROM ubuntu:latest AS compile-image

# Install basic stuff
RUN apt-get update && apt-get install -y \
    nano       \
    wget       \
    git        \
    cmake      \
    pkg-config

# Install compilers
RUN apt-get update && apt-get install -y \
    gcc-8      \
    gfortran-8

# Install GLib and GObject related packages
RUN apt-get update && apt-get install -y \
    gtk-doc-tools          \
    gobject-introspection  \
    libgirepository1.0-dev \
    libglib2.0-dev

# Install package building tools
RUN apt-get update && apt-get install -y \
    autoconf      \
    automake      \
    autotools-dev \
    libtool

# Install dependencies
RUN apt-get update && apt-get install -y \
    libgsl-dev       \
    libgmp-dev       \
    libmpfr-dev      \
    libcfitsio-dev   \
    libfftw3-dev     \
    libnlopt-dev     \
    liblapack-dev    \
    libopenblas-dev  \
    libhdf5-dev      \
    libflint-arb-dev

# NumCosmo (Creates a NumCosmo dir and copy everything from context to it)
RUN cd && mkdir NumCosmo usr
COPY . /root/NumCosmo/

# Set environment variables 
ENV CUBACORES=1
ENV OMP_NUM_THREADS=1
ENV OMP_THREAD_LIMIT=1

WORKDIR /root/NumCosmo

RUN NOCONFIGURE=yes ./autogen.sh
RUN CC=gcc-8 FC=gfortran-8 F90=gfortran-8 F77=gfortran-8 CFLAGS="-O3 -g -Wall" FCFLAGS="-O3 -g -Wall" ./configure --prefix=/usr --enable-opt-cflags 
RUN make -j12
RUN make install DESTDIR=/root

FROM ubuntu:latest AS runtime-image

# Install python
RUN apt-get update && apt-get install -y \
    python3-gi         \
    python3-numpy      \
    python3-matplotlib \
    python3-scipy

# Install dependencies (runtime only)
RUN apt-get update && apt-get install -y \
    libgfortran5     \
    libgsl23         \
    libgmp10         \
    libmpfr6         \
    libcfitsio5      \
    libfftw3-double3 \
    libfftw3-single3 \
    libnlopt0        \
    liblapack3       \
    libopenblas-base \
    libhdf5-100      \
    libflint-arb2 

COPY --from=compile-image /root/usr /usr
