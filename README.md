# M<sup>4</sup>AGO sinking scheme 

Microstructure, Multiscale, Mechanistic, Marine Aggregates in the Global Ocean (M<sup>4</sup>AGO) sinking scheme


Both, the package and the description are in development, while the code provides the functionality as originally published in: 

Maerz et al. 2020: Microstructure and composition of marine aggregates as co-determinants for vertical particulate organic carbon transfer in the global ocean. Biogeosciences, 17, 1765-1803, https://doi.org/10.5194/bg-17-1765-2020

## Table of contents
[How To?](#howto)

  - [Use the package](#packageuse)
  
  - [Contribute](#contribute)

  - [Use the package as submodule](#aboutsubmodule)

[Licensing](#licensing)



<a name="howto"/>

## How To?

</a>

<a name="packageuse"/>

### Use the package

</a>

The package requires an interface, where the driving ocean biogeochemical (OBGC) model provides certain parameters and primary particle information generated from the model. 
Thus far, M<sup>4</sup>AGO is applied in variants of the HAMOCC (HAMburg Ocean Carbon Cycle) model and was originally developed for MPIOM as part of MPI-ESM. 
With the re-writing and separation into the interface functions and the M<sup>4</sup>AGO core module, it becomes easier to integrate and extend the M<sup>4</sup>AGO sinking scheme to other OBGC models.  

The module file `src/mo_ihamocc4m4ago.f90` holds the interface for iHAMOCC as part of NorESM and can serve as an example for including M<sup>4</sup>AGO into other OBGC models.

<a name="contribute"/>

### Contribute

</a>

Contributions are generally welcome!

In case you want to contribute, ideally get in touch and fork the repository, create your own branch in the forked repository and once new developments are ready, file a pull request.



<a name="aboutsubmodule"/>

### Use the package as submodule

</a>

Git submodules provides an easy way of integrating the package into existing code. The first few steps are:

Add the package as submodule: 
```
git submodule add https://github.com/your_account_name/M4AGO-sinking-scheme pckg_path/M4AGO-sinking-scheme
```
which will i) clone the repository into the packages path `pckg_path` and ii) add an entry to `.gitmodules` (in case other submodules are already present) or initialize `.gitmodules` (in case it is the first submodule).
Thereafter you can add and commit the new changes and push it to your OBGC repository. 

When installing the OBGC model including M<sup>4</sup>AGO, perform a:

```
git clone https://github.com/your_account_name/your_OBGC
cd your_OBGC
git submodule update --init --recursive
```
or 
```
git clone --recurse-submodules https://github.com/your_account_name/your_OBGC
```
which will recursively check out submodule packages while cloning your main OBGC.

<a name="licensing"/>

## Licensing

</a>

This software is provided under:

The 3-Clause BSD License
SPDX short identifier: BSD-3-Clause
See https://opensource.org/licenses/BSD-3-Clause




