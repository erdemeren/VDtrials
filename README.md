### VDtrials ###
This is a collection of programs using [SCOREC](https://github.com/SCOREC/core) and [VDLIB](https://www.google.com) to demonstrate some of the capabilities of both libraries. Not all meshes shown in VD_trials.sh are uploaded here due to space constraints, but can be provided upon request.

A typical simulation can be started by running the following commands:

```
make disptess
./vddisp <tess_file> <ext_opts.txt> <eqn.txt> <time_opts.txt> <col_flag> <ins_flag> <topological_th_ratio> verb.txt
```

<ext_opts.txt> is contains instructions to extract information to a csv file at every sub-iteration, <eqn.txt> contains options for the equations of motion used, and <time_opts.txt> contains time integration related parameters. <col_flag> and <ins_flag> specify whether collapses or insertions are allowed and <topological_th_ratio> is the ratio of the topological length for considering transition to the reference length which is typically the median 1-stratum length.

# Examples: #

To run a 20 grain microstructure with local volume preserving boundary conditions and extracting the rate of change of volumes of grains: 
```
./vddisp ./mshfiles/gigi20.tess ./mshfiles/ext_opts_MSROCVOL.txt ./mshfiles/eqn_mason_shell_rk2.txt ./mshfiles/time_opts.txt 1 1 50. > verb.txt
```

VD_trials.sh contains options for compiling and running example programs to generate the results for articles [PhysRevMaterials.5.103802](https://doi.org/10.1103/PhysRevMaterials.5.103802), [PhysRevB.104.L140103](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.104.L140103), and [j.commatsci.2022.111632](https://doi.org/10.1016/j.commatsci.2022.111632).

To simulate the box in a box example case showing spurious insertions, call the following in the VDtrials directory:
```
./VD_trials.sh 21
```

## Microstructure generation ##
The library uses either tesselations (.tess files) generated using [Neper](https://neper.info/) or output meshes from SCOREC obtained during simulations as input files.

# Neper #
As an example, the following code can be used to generate a non-periodic tessellation composed of 1000 equi-axed grains.

```
neper --rcfile './.neperrc'  -T -n 1000 -id 1 -morpho  "gg" -o gigi1000 -format 'tess'
```

## Reference ##
More detailed information about the method is available on:
[Topological transitions during grain growth on a finite element mesh](https://doi.org/10.1103/PhysRevMaterials.5.103802).
The preprint is available on the arXiv:
[Topological transitions during grain growth on a finite element mesh](https://arxiv.org/abs/2101.12321).

## License ##
VDlib is licensed under either of

 * Apache License, Version 2.0, ([LICENSE-APACHE](LICENSE-APACHE) or https://www.apache.org/licenses/LICENSE-2.0)
 * MIT license ([LICENSE-MIT](LICENSE-MIT) or https://opensource.org/licenses/MIT)

at your option.

## Contact ##

* This code base is developed by Erdem Eren, ereren@ucdavis.edu 
* Other community or team contact: 
* Coordinator: Jeremy Mason jkmason@ucdavis.edu
