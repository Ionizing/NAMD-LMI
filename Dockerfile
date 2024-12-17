FROM ghcr.io/cross-rs/x86_64-unknown-linux-gnu:main-centos

RUN sed -i.bak \
    -e 's|^mirrorlist=|#mirrorlist=|g' \
    -e 's|^#baseurl=http://mirror.centos.org/centos|baseurl=https://vault.centos.org/centos|g' \
    /etc/yum.repos.d/CentOS-Base.repo
RUN yum makecache

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
