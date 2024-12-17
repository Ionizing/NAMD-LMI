FROM ghcr.io/cross-rs/x86_64-unknown-linux-gnu:main-centos

RUN sed -i.bak \
<<<<<<< HEAD
    -e 's|^mirrorlist=|#mirrorlist=|g' \
    -e 's|^#baseurl=http://mirror.centos.org/centos|baseurl=https://vault.centos.org/centos|g' \
    /etc/yum.repos.d/CentOS-Base.repo
RUN yum makecache
=======
            -e 's|^mirrorlist=|#mirrorlist=|g' \
            -e 's|^#baseurl=http://mirror.centos.org/centos|baseurl=https://vault.centos.org/centos|g' \
            /etc/yum.repos.d/CentOS-Base.repo
RUN yum makecache
RUN yum update -y && yum group install -y 'Development Tools'

# install rust
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | bash -s -- -y --target x86_64-unknown-linux-gnu
ENV PATH="${HOME}/.cargo/bin:${PATH}"

# install blas
RUN yum install -y epel-release
RUN yum update -y && yum install -y \
    lapack-devel \
    openblas-devel \
    cmake3

# rename cmake3 to cmake
RUN alternatives --install /usr/local/bin/cmake cmake /usr/bin/cmake 10 \
--slave /usr/local/bin/ctest ctest /usr/bin/ctest \
--slave /usr/local/bin/cpack cpack /usr/bin/cpack \
--slave /usr/local/bin/ccmake ccmake /usr/bin/ccmake \
--family cmake
RUN alternatives --install /usr/local/bin/cmake cmake /usr/bin/cmake3 20 \
--slave /usr/local/bin/ctest ctest /usr/bin/ctest3 \
--slave /usr/local/bin/cpack cpack /usr/bin/cpack3 \
--slave /usr/local/bin/ccmake ccmake /usr/bin/ccmake3 \
--family cmake
>>>>>>> 0e2e146b9b61ca3a1351d824b2db7d3a9fbb9d3d

# install mkl
RUN yum-config-manager --add-repo https://yum.repos.intel.com/mkl/setup/intel-mkl.repo
RUN rpm --import https://yum.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
RUN yum update -y && yum install -y \
    intel-mkl-64bit-2020.0-088

# setup environment variables
ENV MKLROOT="/opt/intel/mkl"
ENV LIBRARY_PATH="/opt/intel/tbb/lib/intel64/gcc4.8:/opt/intel/compiler/lib/intel64:/opt/intel/mkl/lib/intel64"
ENV LD_LIBRARY_PATH="/opt/intel/lib/intel64:/opt/intel/mkl/lib/intel64:/opt/intel/tbb/lib/intel64"
ENV CPATH="/opt/intel/mkl/include"
ENV NLSPATH="/opt/intel/compilers_and_libraries/linux/compiler/lib/intel64/locale/%l_%t/%N:/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64_lin/locale/%l_%t/%N"
ENV PKG_CONFIG_PATH="/opt/intel/mkl/bin/pkgconfig"
