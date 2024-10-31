M4AGOrun = executable('M4AGO',
                      sources,
		      include_directories: includes,
		      dependencies:deps,
		      build_by_default:false,
		      link_language:'fortran')

test('run M4AGO',
     M4AGOrun,
     workdir:meson.current_build_dir())
