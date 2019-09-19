FROM ubuntu:latest

RUN apt-get update
# Install basic stuff (language-pack-en powerline dconf-cli emacs zsh)
RUN apt-get install -y nano wget git cmake pkg-config
# Install compilers
RUN apt-get install -y gcc-8 gfortran-8
# Install G things
RUN apt-get install -y gtk-doc-tools gobject-introspection libgirepository1.0-dev libglib2.0-dev
# Install python
RUN apt-get install -y python-numpy python-gobject python-gi-cairo

# Install package building tools
RUN apt-get install -y autoconf automake autotools-dev libtool

# Install dependencies
RUN apt-get install -y libgsl-dev libgmp-dev libmpfr-dev libcfitsio-dev libfftw3-dev libnlopt-dev liblapack-dev libopenblas-dev libhdf5-dev libflint-arb-dev

# NumCosmo (clone from dropbox)
RUN cd && mkdir NumCosmo
COPY . /root/NumCosmo/

# Set environment variables 
ENV CUBACORES=1
ENV OMP_NUM_THREADS=1
ENV OMP_THREAD_LIMIT=1

RUN cd && cd NumCosmo 
RUN NOCONFIGURE=yes ./autogen.sh
RUN CC=gcc-8 FC=gfortran-8 F90=gfortran-8 CFLAGS="-O3 -g -Wall" FCFLAGS="-O3 -g -Wall" ./configure --prefix=/usr --enable-opt-cflags 
RUN make -j12
RUN make install
