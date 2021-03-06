
import os

# Initialize the environment with path variables, CFLAGS, and so on
bioinfo_path = '#lib/bioinfo-libs'
commons_path = '#lib/common-libs'
#math_path = '#libs/math'

vars = Variables('buildvars.py')

compiler = ARGUMENTS.get('compiler', 'gcc')

env = Environment(tools = ['default', 'packaging'],
      		  CC = compiler,
                  variables = vars,
                  CFLAGS = '-std=c99 -D_XOPEN_SOURCE=600 -D_GNU_SOURCE -fopenmp ',
                  CPPPATH = ['#', '#src', '#include', bioinfo_path, commons_path, "%s/commons/argtable" % commons_path, "%s/commons/config" % commons_path ],
                  LIBPATH = ['#libs/common-libs/', commons_path],
                  LIBS = ['m', 'z'],
                  LINKFLAGS = ['-fopenmp'])
                  

if int(ARGUMENTS.get('debug', '0')) == 1:
    debug = 1
    env['CFLAGS'] += ' -O0 -g'
else:
    debug = 0
    env['CFLAGS'] += ' -O3 -g'

env['objects'] = []


# Targets

SConscript(['%s/SConscript' % bioinfo_path,
            '%s/SConscript' % commons_path
            ], exports = ['env', 'debug', 'compiler'])

env.Program('#bin/hpg-aligner',
             source = [Glob('src/*.c'),
	               Glob('src/build-index/*.c'),
	               Glob('src/dna/*.c'),
	               Glob('src/rna/*.c'),
                       "%s/libcommon.a" % commons_path,
                       "%s/libbioinfo.a" % bioinfo_path
                      ]
           )

# Create a Binary tarball
(distro, release) = os.popen('lsb_release -sir').read().strip().split('\n')
tb = env.Package(NAME          = 'hpg-aligner-'+distro+'_'+release+'-x86_64',
                VERSION        = '1.0.1',
                PACKAGEVERSION = 1,
                PACKAGETYPE    = 'targz',
                LICENSE         = 'gpl', 
                source         = ['#COPYING', '#README', '#bin/hpg-aligner' ] )

# Create a Source tarball
sourceFiles = [Glob("src/*.[ch]"), Glob("src/dna/*.[ch]"), Glob("src/rna/*.[ch]"), Glob("src/build-index/*.[ch]")]
sourceFiles = sourceFiles + [Glob("lib/common-libs/SConscript"), Glob("lib/common-libs/README"), Glob("lib/common-libs/COPYING"), Glob("lib/common-libs/*/*.[ch]"), Glob("lib/common-libs/*/*/*.[ch]")]
sourceFiles = sourceFiles + [Glob("lib/bioinfo-libs/SConscript"), Glob("lib/bioinfo-libs/README"), Glob("lib/bioinfo-libs/COPYING"), Glob("lib/bioinfo-libs/*/*/*.[ch]"), Glob("lib/bioinfo-libs/*/*/*/*.[ch]"), Glob("lib/bioinfo-libs/*/*/*/*/*.[ch]"), Glob("lib/bioinfo-libs/*/SConscript"), Glob("lib/bioinfo-libs/*/*/SConscript")]
tb = env.Package(NAME          = 'hpg-aligner-x86_64-src',
                VERSION        = '1.0.1',
                PACKAGEVERSION = 1,
                PACKAGETYPE    = 'src_targz',
                LICENSE         = 'gpl', 
                source         = ['#COPYING', '#README', '#SConstruct', sourceFiles ] )


'''
if 'debian' in COMMAND_LINE_TARGETS:
    SConscript("deb/SConscript", exports = ['env'] )
'''
