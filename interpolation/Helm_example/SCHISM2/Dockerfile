#  Single-Stage Dockerfile
FROM ubuntu:18.04 

# Set up working directory
WORKDIR /app

# copy files for integration
COPY . .

# Set up environment variables
ENV LC_ALL=C 

#  HOME="/gpfs/work/jacobb/data/SINGULARITY/SCHISM/"

# Update and install necessary dependencies
RUN apt-get update && apt-get install -y \
    netcat \
    gnupg curl \
    git \
    gfortran g++ \
    libnetcdf-dev libnetcdff-dev \
    libmpich-dev \
    make \
    python-minimal \
    python3-minimal \
    nano \
    vim \	
    wget \
    checkinstall

	
# Install cmake 3.20.0
RUN wget https://github.com/Kitware/CMake/releases/download/v3.20.0/cmake-3.20.0.tar.gz && \
    tar -zvxf cmake-3.20.0.tar.gz && \
    cd cmake-3.20.0 && \
    ./configure -- -DCMAKE_USE_OPENSSL=OFF && \
    ./bootstrap && \
    make -j8 && \
    checkinstall --pkgname=cmake --pkgversion="3.20-custom" --default && \
    hash -r

# Clone SCHISM repository
RUN mkdir -p /gpfs/work/jacobb/data/SINGULARITY/SCHISM && \
    git clone -b v5.11 --single-branch --depth=1 https://github.com/schism-dev/schism.git /SCHISM

# old 
#RUN mkdir /SCHISM && \
#git clone -b v5.11 --single-branch --depth=1 https://github.com/schism-dev/schism.git /SCHISM


#Fails
# Create the SCHISM directory inside the container
#RUN mkdir /SCHISM 

# Copy the local SCHISM repository to the container
#COPY ./schism_code/* /SCHISM


# Configure SCHISM build
#RUN cd /SCHISM && \
#    ln -s cmake/SCHISM.local.ubuntu SCHISM.local.build && \
#    mkdir GBconfig && cd GBconfig && \
#    cmake -C ../cmake/SCHISM.local.build ../src -DTVD_LIM=SB -DCMAKE_Fortran_FLAGS=" -DPREC_EVAP=ON" -DPREC_EVAP=ON && \
#    make all

# Create symlink and build directory once (reused for both builds)
RUN cd /SCHISM && \
    ln -s cmake/SCHISM.local.ubuntu SCHISM.local.build && \
    mkdir -p GBconfig

# Configure SCHISM build # with atmo (creates pschism_PREC_EVAP_TVD-SB)
RUN cd /SCHISM/GBconfig && \
    cmake -C ../cmake/SCHISM.local.build ../src -DTVD_LIM=SB -DCMAKE_Fortran_FLAGS=" -DPREC_EVAP=ON" -DPREC_EVAP=ON && \
    make all

# Configure SCHISM build # without atmo (creates pschism_TVD-SB, reuses same directory)
RUN cd /SCHISM/GBconfig && \
    cmake -C ../cmake/SCHISM.local.build ../src -DTVD_LIM=SB && \
    make all



# Add a custom entry script
#RUN echo -e "cd /rundir\nSCHISM/GBconfig/bin/pschism_PREC_EVAP_TVD-SB 6" > /SCHISM/exec.sh && \
#    chmod u+rwx /SCHISM/exec.sh

RUN echo -e "cd /rundir\n mpirun -n 6 SCHISM/GBconfig/bin/pschism_PREC_EVAP_TVD-SB 6" > /SCHISM/exec.sh && \
    chmod u+rwx /SCHISM/exec.sh


#### test conda

# Set environment variables to prevent prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies for downloading Miniforge and running Conda
RUN apt-get update && \
    apt-get install -y wget bash bzip2 && \
    rm -rf /var/lib/apt/lists/*

# Download and install Miniforge
RUN wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" -O Miniforge3.sh && \
    bash Miniforge3.sh -b -p /opt/miniforge && \
    rm Miniforge3.sh

# Add Miniforge to PATH
ENV PATH="/opt/miniforge/bin:$PATH"

# Pre-initialize Conda to ensure it works in non-interactive mode
RUN conda init bash && echo ". /opt/miniforge/etc/profile.d/conda.sh" >> ~/.bashrc

# Create Conda environment a with required channels
RUN /bin/bash -c ". /opt/miniforge/etc/profile.d/conda.sh && \
    conda create -n forcing -c conda-forge -y"

# Activate environment and install additional packages
RUN /bin/bash -c ". /opt/miniforge/etc/profile.d/conda.sh && \
    conda activate forcing && \
	conda install conda-forge::copernicusmarine -y"


RUN /bin/bash -c ". /opt/miniforge/etc/profile.d/conda.sh && \
    conda activate forcing && \
	conda update --name forcing copernicusmarine --yes"

RUN /bin/bash -c ". /opt/miniforge/etc/profile.d/conda.sh && \
    conda activate forcing && \
    conda install -c conda-forge numpy matplotlib pandas cartopy scipy shapely networkx owslib ipython requests boto3 s3fs netcdf4 xarray copernicusmarine PyYAML reportlab cdsapi -y"



#### test


# Add labels
LABEL Author="benjamin.jacob@hereon.de" \
      Version="v0.0.1"

# Set environment variables
ENV LC_ALL=C \
    HOME="/gpfs/work/jacobb/data/SINGULARITY/SCHISM/"

# Copy files from the build stage
#COPY --from=build /SCHISM /SCHISM

# Add files and directories
RUN mkdir -p /rundir/outputs

# Expose the listen port
#EXPOSE $LISTEN_PORT

# Add test command
RUN grep -q NAME=\"Ubuntu\" /etc/os-release || echo "Container base is not Ubuntu."

# Default command for running the container
CMD ["/app/entrypoint.sh"]