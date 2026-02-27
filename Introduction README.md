# Scientific Computing – Docker Environment (AlmaLinux 9)

This repository contains the solutions for the *Scientific Computing* tasks.
All codes have been tested inside a clean **AlmaLinux 9** Docker container to guarantee
full reproducibility across different systems.

The Docker image installs all required compilers, libraries, and Python packages
needed to compile and run the tasks.

---

## Requirements

- Docker (Docker Desktop on Windows/macOS, or Docker Engine on Linux)

No other dependencies are required on the host system.

---

## Docker Image Overview

The custom Docker image is built on top of:

- `almalinux:9`

and includes:

- GCC / G++
- OpenMP (via GCC)
- OpenMPI
- GSL
- FFTW
- HDF5
- Python 3
- NumPy, SciPy, Matplotlib

---

## Building the Docker Image

### Create the Dockerfile (File) as follow:
```
FROM almalinux:9

# ------------------------------------------------------------
# Base system + repositories
# ------------------------------------------------------------
RUN dnf update -y && \
    dnf install -y dnf-plugins-core && \
    dnf config-manager --set-enabled crb && \
    dnf install -y epel-release && \
    dnf makecache

# ------------------------------------------------------------
# Compilers, build tools, editors + pkg-config
# ------------------------------------------------------------
RUN dnf install -y \
        gcc \
        gcc-c++ \
        make \
        cmake \
        git \
        vim \
        which \
        pkgconf-pkg-config

# ------------------------------------------------------------
# Python (interpreted language tasks) + scientific stack
# ------------------------------------------------------------
RUN dnf install -y python3 python3-pip && \
    pip3 install --no-cache-dir numpy scipy matplotlib

# ------------------------------------------------------------
# MPI (Task 09)
# ------------------------------------------------------------
RUN dnf install -y openmpi openmpi-devel

# ------------------------------------------------------------
# Scientific libraries
# ------------------------------------------------------------
# GSL (Task 05)
RUN dnf install -y gsl gsl-devel

# FFTW (Task 06)
RUN dnf install -y fftw fftw-devel

# ------------------------------------------------------------
# HDF5 (Task 03) - install BOTH serial and openmpi flavors
# ------------------------------------------------------------
# This avoids non-standard include paths (e.g. /usr/include/openmpi-x86_64/hdf5.h)
# and makes compilation possible via pkg-config on RHEL-like distros.
RUN dnf install -y \
        hdf5 hdf5-devel && \
    (dnf install -y hdf5-openmpi hdf5-openmpi-devel || true)

# ------------------------------------------------------------
# Clean up
# ------------------------------------------------------------
RUN dnf clean all && rm -rf /var/cache/dnf

# ------------------------------------------------------------
# Environment variables (MPI)
# ------------------------------------------------------------
ENV PATH=/usr/lib64/openmpi/bin:$PATH
ENV LD_LIBRARY_PATH=/usr/lib64/openmpi/lib

# ------------------------------------------------------------
# Work directory
# ------------------------------------------------------------
WORKDIR /work

```
From the root directory of this repository (where the `Dockerfile` is located):

```bash
docker build -t scicomp-alma9 .
```

Then run an interactive container and mount the repository directory:
```
docker run -it \
  --name scicomp_container \
  -v $(pwd):/work \
  scicomp-alma9
```

If you are working on Windows, using Powershell:
```
docker run -it `
  --name scicomp_container `
  -v ${PWD}:/work `
  scicomp-alma9
```
