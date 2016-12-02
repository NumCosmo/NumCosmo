FROM ubuntu:latest

RUN apt-get update
# Install basic stuff (language-pack-en powerline dconf-cli emacs zsh)
RUN apt-get install -y nano wget git cmake
# Install compilers
RUN apt-get install -y gcc-5 gfortran gfortran-5
# Install python
RUN apt-get install -y python-healpy
# Install G things
RUN apt-get install -y gtk-doc-tools gobject-introspection 

# Install dependencies
RUN apt-get install -y libglib2.0-dev libgsl-dev libgmp-dev libmpfr-dev libcfitsio-dev libfftw3-dev libnlopt-dev
# sundials
RUN cd && wget http://computation.llnl.gov/projects/sundials/download/sundials-2.7.0.tar.gz && \
	tar -xvf sundials-2.7.0.tar.gz && \
	mkdir sundials-build && cd sundials-build && cmake ../sundials-2.7.0 && make && make install

# NumCosmo (clone from dropbox)
RUN cd && mkdir NumCosmo && cd NumCosmo
COPY * ./
RUN ./autogen.sh && make -j12 && make install # gtkdocize && 
