project('M4AGO',['fortran','c'],
        license:'3-Clause BSD License',
	default_options:['warning_level=2'])

fcc=meson.get_compiler('fortran')

# Define global dependencies
netcdf = dependency('netcdf', language: 'fortran', version: '>=4.4.4')
omp    = dependency('openmp', required: get_option('openmp'))

# List of dependencies (created so that options can add dependencies)
deps = [netcdf, omp]


sources  = []
includes = []

# Add files to required sources to build the program
sources += files('M4AGO_driver.F90',
                 'mo_driver_routines.f90',
                 '../../src/mo_m4ago_HAMOCCinit.F90',
                 '../../src/mo_m4ago_HAMOCCPrimPart.F90',
		 '../../src/mo_m4ago_core.f90',
		 '../../src/mo_m4ago_kind.F90',
		 '../../src/mo_m4ago_types.f90',
		 '../../src/mo_m4ago_params.f90',
		 '../../src/mo_m4ago_control.f90',
		 '../../src/mo_m4ago_physics.f90')

subdir('tests')

executable('M4AGO',sources,include_directories:includes,dependencies:deps,link_language:'fortran')


