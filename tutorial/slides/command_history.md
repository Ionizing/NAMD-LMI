# # NAMD-LMI tutorial command history



# Load environment

```sh
module load intel/2020u4
source /ssd/app/CMsoft/anaconda/bin/activate
```

# Copy example file

```sh
cp -rv /data/app/Hefeitest/NAMD-LMI_tutorial .
cd NAMD-LMI_tutorial/MoSe2_input_with_result
```

# VASP PART

## Primitive cell opt

```sh
cd 1x1
sbatch sub.sh
```

## Supercell opt

```sh
cd 2x3/opt
cp ../../1x1/CONTCAR POSCAR.prim
./supercell.py  # generate POSCAR from POSCAR.prim
sbatch sub.sh
```

## Supercell NVT

```sh
cd 2x3/nvt
cp ../opt/CONTCAR POSCAR
sbatch sub.sh
```

## Supercell NVE

```sh
cd 2x3/nve
cp ../nvt/CONTCAR POSCAR
sbatch sub.sh
```

## Supercell SCF calculation on NVE trajectory

```sh
cd 2x3/nve
./genstatic.py
cd 2x3/static_ncl
mv ../nve/run .
sbatch sub_vasp_namd
./mixWave.py
```

# NAMD-LMI PART

## Calculate NAC

```sh
cd 2x3/namd
cp /data/app/Hefeitest/NAMD-LMI_tutorial/MoSe2_input_with_result/2x3/namd/namd_lmi .
./namd_lmi nac -c 01_nac_config.toml
```

## Generate HAMIL.h5

```sh
cd 2x3/namd
./namd_lmi hamil -c 02_hamil_config_sigmaplus.toml
./hamil_plot_pchiral.py HAMIL_sigmaplus.h5

./namd_lmi hamil -c 02_hamil_config_sigmaminus.toml
```

## Run surface hopping

```sh
cd 2x3/namd
./namd_lmi surfhop -c 03_surfhop_config_sigmaplus.toml
cd excitation_sigmaplus
../surfhop_plot.py


cd 2x3/namd
./namd_lmi surfhop -c 03_surfhop_config_sigmaminus.toml
cd excitation_sigmaminus
../surfhop_plot.py
```
