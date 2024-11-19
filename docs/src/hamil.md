# Hamiltonian

The Hamiltonian of `namd_lmi` is simply a cropped coupling of the Non-adiabatic
couplings plus the external field.

## Help message

```shell
$ namd_lmi hamil --help
Generate the Hamiltonian from NAC according to config file

Usage: namd_lmi hamil [OPTIONS]

Options:
  -c, --config <CONFIG>
          Config file name.

          Aliases: "cfg", "conf".

          [default: hamil_config.toml]

      --generate <GENERATE>
          Generate auxiliary files for the calculation and analysis.

          The generation of Hamiltonian will not run if this flag is set.

          Alias: "gen".

          Possible values:
          - config-template:      Generate config template for Hamiltonian generation. Aliases: "config", "cfg" and "conf"
          - efield-template:      Generate script template for external electric field. Aliases: "efield", "ef"
          - postprocess-template: Generate post-process scripts for Hamiltonian analysis. Aliases: "post-process", "postprocess", "pp"

  -h, --help
          Print help (see a summary with '-h')
```

## Procedures

1. Generate a configuration template

	```shell
    $ ./run.sh hamil --generate conf
    2024-11-19 21:49:48 [ INFO] Global logger initialized with targets being stderr and "./globalrun.log"
    2024-11-19 21:49:48 [ INFO] Writing `02_hamil_config_template.toml` ...
    2024-11-19 21:49:48 [ INFO] Writing config to file "02_hamil_config_template.toml"
    2024-11-19 21:49:48 [ INFO] Time used: 1.624654ms
	```
