# Installation

There are two ways to install `NAMD-LMI`, downloading the pre-built binaries or
building from scratch.


## Download the binary from Github Release

We provide the pre-built binaries at [Github Release](https://github.com/Ionizing/NAMD-LMI/releases).

There are three platforms we support with pre-built binaries:

- `namd_lmi-<VERSION>-linux-x86_64.tar.gz`: Dynamically linked against `glibc`
  and other fundamental system libraries. This version can run on most Linux
  distros (e.g. CentOS 7, Ubuntu 18.04, etc).
- `rsgrad-<VERSION>-macos-x86_64.tar.gz`: Binary for macOS (Intel Chip).
- `rsgrad-<VERSION>-macos-aarch64.tar.gz`: Binary for macOS (Apple Silicon).
- `rsgrad-<VERSION>-windows-x86_64.zip`: Binary for Windows system. Windows 7
  and upper 64-bit versions are supported.


## Build from scratch

The build requires Internet accessibility and Rust toolchain, Intel-MKL libraries.

- If you haven't installed Rust toolchain, run `curl --proto '=https' --tlsv1.2 -sSf
https://sh.rustup.rs | sh` and follow the instructions to configure it. When you finished
the toolchain installation, you should be able to run `cargo` in command-line. Detailed
instructions can be found [here](https://rustup.rs).
- If you have installed Rust toolchain and Intel-MKL already, just run `cargo
install --git https://github.com/Ionizing/NAMD-LMI`. Several minutes later, the
`namd_lmi` binary should be installed to `~/.cargo/bin`, no need to modify the
`$PATH` (or `%PATH%` on Windows).
