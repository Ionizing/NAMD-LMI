# Installation

There are two ways to install `NAMD-LMI`, downloading the pre-built binaries or
building from scratch.


## Download the binary from Github Release

We provide the pre-built binaries at [Github Release](https://github.com/Ionizing/NAMD-LMI/releases).

There are several platforms we support with pre-built binaries, in which we
can identify them by their file names parts:

- **namd_lmi-<VERSION>**: the version of NAMD-LMI;
- **linux/macos/windows**: which operating system it runs on;
- **x86_64/aarch64**: which CPU architecture it runs on;
- **mkl-system/mkl-static/openblas-system/openblas-static**: which BLAS
  implementation and which link scheme, where `xxx-system` **requires you
  install corresponding BLAS implementation on system**, while `xxx-static`
  statically link against to BLAS (which means you needn't install anything
  else) but come with larger binary sizes.
- **.tar.gz/zip**: which compressor format it uses. You need to decompress
  the binary first.

Example:
- `namd_lmi-1.0.0-linux-x86_64-mkl-system.tar.gz` should run on an x86
  machine with 64-bit Linux installed, where the Intel-MKL should be installed.
- `namd_lmi-1.0.0-macos-x86_64-mkl-system.tar.gz` should run on macOS with Intel
  chip.
- `namd_lmi-1.0.0-macos-aarch64-mkl-system.tar.gz` should run on macOS with Apple
  Silicon.
- `namd_lmi-1.0.0-windows-x86_64-mkl-static.zip` should run on Windows 7 or
  upper with Intel/AMD or other x86 chips.


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

There are four features to link BLAS implementations:
- `intel-mkl-system`: dynamically link against pre-installed Intel-MKL library. You need
  to install Intel-MKL first.
- `intel-mkl-static`: statically link against Intel-MKL library. If no pre-installed
  MKL is found, it will download Intel-MKL 2020 from web by itself.
- `openblas-system`: dynamically link against pre-installed OpenBLAS. You need to install
  OpenBLAS first.
- `openblas-static`: statically link against OpenBLAS library. If no pre-installed OpenBLAS
  is found, it will download OpenBLAS source code and build it from scratch.


## Test the installation

Once you obtain the `namd_lmi` binary, just put it in your `PATH`, and then you
should be able to run it with the output like
```shell
$ namd_lmi
+----------------------------------------------------------------------+
|                                                                      |
|    _   _            __  __  _____           _       __  __  _____    |
|   | \ | |    /\    |  \/  ||  __ \         | |     |  \/  ||_   _|   |
|   |  \| |   /  \   | \  / || |  | | ______ | |     | \  / |  | |     |
|   | . ` |  / /\ \  | |\/| || |  | ||______|| |     | |\/| |  | |     |
|   | |\  | / ____ \ | |  | || |__| |        | |____ | |  | | _| |_    |
|   |_| \_|/_/    \_\|_|  |_||_____/         |______||_|  |_||_____|   |
|                                                                      |
+----------------------------------------------------------------------+

Welcome to use namd!
    current version:    0.1.0
    git hash:           d659586
    author(s):          Ionizing
    host:               x86_64-unknown-linux-gnu
    built time:         2024-11-19 15:34:05 +08:00


Usage: namd_lmi <COMMAND>

Commands:
  nac      Calculate non-adiabatic coupling (NAC) including `<j| d/dt |k>` and momentum matrix `<i| p |j>`
  hamil    Generate the Hamiltonian from NAC according to config file
  surfhop  Perform the surface-hopping process with given Hamiltonian file and config file
  help     Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help (see more with '--help')
  -V, --version  Print version
```
