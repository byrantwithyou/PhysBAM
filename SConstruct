##############################################################################
# Copyright 2004-2008, Geoffrey Irving, Frank Losasso, Andrew Selle.
# This file is party of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
##############################################################################
import sys
import os
import glob
import commands
import re

default_cxx='g++'

### variables and base environment
variables=Variables('SConstruct.options')
variables.AddVariables(
    ('CXX','C++ compiler',default_cxx),
    EnumVariable('ARCH','Architecture (e.g. pentium4, opteron, nocona, powerpc)','',allowed_values=('pentium3','pentium4','opteron','nocona','powerpc','test','wine','')),
    EnumVariable('DEFAULT_ARCH','Architecture that doesn\'t need a suffix','',allowed_values=('pentium3','pentium4','opteron','nocona','powerpc','test','wine','darwin','')),
    EnumVariable('TYPE','Type of build (release, debug, profile, optdebug)','release',allowed_values=('release','debug','profile','optdebug')),
    ('cache','Cache directory to use',''),
    BoolVariable('shared','Build shared libraries',1),
    BoolVariable('wrapper','Build wrapper executable file',1),
    BoolVariable('USE_RAY_TRACING','Use ray tracing',0),
    BoolVariable('USE_COMPRESSIBLE','Use compressible flow',0),
    BoolVariable('USE_DEFORMABLES','Use deformable objects',0),
    BoolVariable('USE_DYNAMICS','Use dynamics',0),
    BoolVariable('USE_FLUIDS','Use fluids',0),
    BoolVariable('USE_GEOMETRY','Use geometry',0),
    BoolVariable('USE_INCOMPRESSIBLE','Use incompressible',0),
    BoolVariable('USE_RIGIDS','Use rigid bodies',0),
    BoolVariable('USE_SOLIDS','Use solids',0),
    BoolVariable('USE_HYBRID','Use hybrid methods',0),
    BoolVariable('USE_GRID_PDE','Use grid PDE',0),
    BoolVariable('USE_GRID_TOOLS','Use grid tools',0),
    BoolVariable('USE_TOOLS','Use tools',0),
    BoolVariable('USE_CORE','Use core',0),
    BoolVariable('USE_LEX_YACC','Use lex and yacc to allow parsing for symbolics',0),
    BoolVariable('compile_id_types_as_int','Treat ID types as int to avoid possible performance consequences',0),
    BoolVariable('fast_math','compile with -ffast-math',0),
    BoolVariable('install_programs','install programs into source directories',1),
    BoolVariable('compile_headers','compile headers without .cpp files',1),
    BoolVariable('install_headers','install Public_Library headers into $INSTALL_DIR/include/physbam',1),
    BoolVariable('use_rpath','use rpaths for dynamic libraries',1),
    ('CXXFLAGS_EXTRA','',[]),
    ('LINKFLAGS_EXTRA','',[]),
    ('CPPPATH_EXTRA','',[]),
    ('LIBPATH_EXTRA','',[]),
    ('RPATH_EXTRA','',[]),
    ('LIBS_EXTRA','',[]),
    ('INSTALL_PATH','Path to install libraries, binaries, and scripts',''))

### external libraries
arch=os.popen("uname -m").read()[:-1]

### Choose libraries
external_libraries={
    'boost':      {'default':1,'libs':['']},
    'boostregex': {'libs':['boost_regex']},
    'zlib':       {'default':1,'libs':['z']},
    'ffmpeg':     {'default':0,'flags':['USE_FFMPEG'],'libs':['libavformat','libavcodec','libavutil']},
    'boost_geometry':     {'default':0,'flags':['USE_BOOST_GEOMETRY'],'libs':[]},
    'geompack':     {'default':0,'flags':['USE_GEOMPACK'],'libs':['geompack']},
    'boost_serialization':     {'default':0,'flags':['USE_BOOST_SERIALIZATION'],'libs':['libboost_serialization']},
    'libjpeg':    {'default':1,'flags':['USE_LIBJPEG'],'libs':['jpeg']},
    'libpng':     {'default':1,'flags':['USE_LIBPNG'],'libs':['png']},
    'fftw':       {'default':0,'flags':['USE_FFTW'],'libs':['fftw3f','fftw3']},
    'ipopt':      {'default':0,'flags':['USE_IPOPT'],'libs':['ipopt']},
    'lapack':     {'default':0,'flags':['USE_LAPACK'],'libs':['lapack','lapacke','blas']},
    'mkl':        {'default':0,'flags':['USE_MKL'],'libs':[]},
    'gl2ps':      {'libs':['gl2ps']},
    'OpenGL':     {'libs':['GL','GLU','glut'],'libpath':['/usr/X11R6/lib','/usr/X11R6/lib64'],'cpppath':['/opt/X11/include']},
    'boostpython':{'flags':['USE_BOOSTPYTHON'],'cpppath':['/usr/include/python2.4'],'libs':['boost_python','python2.4']},
    'numpy':      {'flags':['USE_NUMPY'],'cpppath':[],'libs':[]},
    'lam':        {'flags':['USE_MPI'],'libs':['lammpio','lammpi++','mpi','lam','util','dl'],'linkflags':' -pthread'},
    'openmpi':    {'flags':['USE_MPI'],'libs':['mpi_cxx','open-rte','mpi','open-pal','util','dl','nsl'],'linkflags':' -pthread'},
    'mpich':      {'flags':['USE_MPI'],'libs':['lammpio','lammpi++','mpi','lam','util','dl'],'linkflags':' -pthread'},
    'pthreads':    {'default':0,'flags':['USE_PTHREADS'],'libs':[],'linkflags':' -pthread'},
    'openmp':    {'default':0,'flags':['USE_OPENMP'],'libs':[],'linkflags':' -fopenmp'},
    'renderman':  {'default':0,'flags':['USE_RENDERMAN'],'cpppath':[],'libs':[]}}
for name,lib in external_libraries.items():
    defaults={'default':0,'flags':'','linkflags':'','cpppath':[],'libpath':[]}
    for f in defaults.keys(): lib.setdefault(f,defaults[f])
    variables.AddVariables(BoolVariable('USE_'+name.upper(),'Use '+name,lib['default']),
                           (name+'_include','Include directory for '+name,0),
                           (name+'_libpath','Library directory for '+name,0),
                           (name+'_rpath','Extra rpath directory for '+name,0),
                           (name+'_libs','Libraries for'+name,0),
                           (name+'_linkflags','Linker flags for '+name,0))

### parse variables
if os.environ.has_key('PLATFORM') and os.environ['PLATFORM']=="wine":
    env=Environment(variables=variables,tools=["msvc","mslink","mslib"],platform="win32cross",CC="clwrap",CXX="clwrap",LINK="linkwrap",AR="libwrap")
    env["CXX"]="clwrap" # TODO: this /is/ necessary but why, maybe the default above?
else:
    env=Environment(variables=variables)
env.Replace(ENV=os.environ) # do this here to allow SConstruct.options to change environment
Help(variables.GenerateHelpText(env))
if env["USE_OPENMP"]:
    env.Append(CXXFLAGS_EXTRA='-fopenmp');

### improve performance
if 'Decider' in dir(env): # if Decider exists, use it to avoid deprecation warnings
    env.Decider('MD5-timestamp')
else:
    env.TargetSignatures('build')
env.SetOption('max_drift',100)
env.SetDefault(CPPPATH_HIDDEN=[]) # directories in CPPPATH_HIDDEN won't be searched for dependencies
env.SetDefault(CPPDEFINES=[])
env.Replace(_CPPINCFLAGS=env['_CPPINCFLAGS']+re.sub(r'\$\( (.*)\bCPPPATH\b(.*)\$\)',r'\1CPPPATH_HIDDEN\2',env['_CPPINCFLAGS']))

### avoid annoying bug in previous versions of scons
env.EnsureSConsVersion(0,96,92)

### avoid deprecation warnings about env.Copy
if 'AddMethod' in dir(env) and 'Clone' in dir(env):
    AddMethod(Environment,Environment.Clone,'Copy')

## override platform specific library names
if env['PLATFORM'].startswith('win32'):
    external_libraries['OpenGL']['libs']=['opengl32','glu32','glaux']
    external_libraries['ffmpeg']['libs']=['avcodec','avformat']
    env.Replace(compile_headers=0)
elif env['PLATFORM']=='darwin':
    opengl=external_libraries['OpenGL']
    opengl['linkflags']='-framework OpenGL -framework GLUT'
    opengl['libpath']=[]
    opengl['libs']=[]
    external_libraries['boostpython']['linkflags']='-framework Python'

### platform
if env['ARCH']=='':
    if os.environ.has_key('PLATFORM'): env['ARCH']=os.environ['PLATFORM']
    else: env['ARCH']='pentium4'

### build cache
if env['cache']!='': CacheDir(env['cache'])

### variant build setup
variant_build=os.path.join('build',env['ARCH'],env['TYPE'])
BuildDir(variant_build,'.',duplicate=0)
program_suffix=''
if env['ARCH']!=env['DEFAULT_ARCH']: program_suffix+='_'+env['ARCH']
if env['TYPE']!='release': program_suffix+='_'+env['TYPE']
if env['PLATFORM'].startswith('win32'): executable_suffix='.exe'
else: executable_suffix=''
library_suffix=''
if not env['shared']: library_suffix='_static'

### process external library related options
for key,lib in external_libraries.items():
    name=key
    lib['rpath']=[]
    if env[name+'_include']!=0: lib['cpppath']=env[name+'_include']
    if env[name+'_libpath']!=0: lib['libpath']=env[name+'_libpath']
    if env[name+'_rpath']!=0: lib['rpath']=env[name+'_rpath']
    if env[name+'_libs']!=0: lib['libs']=env[name+'_libs']
    if env[name+'_linkflags']!=0: lib['linkflags']=env[name+'_linkflags']

### extra flag options
env.Append(CXXFLAGS=env['CXXFLAGS_EXTRA'])
env.Append(LINKFLAGS=env['LINKFLAGS_EXTRA'])
env.Append(CPPPATH=env['CPPPATH_EXTRA'])
env.Append(LIBPATH=env['LIBPATH_EXTRA'])
env.Append(RPATH=env['RPATH_EXTRA'])
env.Append(LIBS=env['LIBS_EXTRA'])

### compiler flags
if env['CXX'].endswith('icc'):
    env.Append(CXXFLAGS=' -g')
    if env['TYPE']=='release' or env['TYPE']=='optdebug' or env['TYPE']=='profile':
        if env['ARCH']=='pentium4': env.Append(CXXFLAGS=' -O3 -xN')
        elif env['ARCH']=='nocona': env.Append(CXXFLAGS=' -O3 -xP')
        else: env.Append(CXXFLAGS=' -O3')
    env.Append(CXXFLAGS=' -w -vec-report0',LINKFLAGS=' -w -Kc++ -cxxlib-icc')
elif env['PLATFORM'].startswith('win32'):
    if env['PLATFORM']=="win32": env.Append(CXXFLAGS=' /Zi')
    if env['TYPE']=='release' or env['TYPE']=='optdebug' or env['TYPE']=='profile':
        env.Append(CXXFLAGS=' /O2') 
    env.Append(CXXFLAGS=' /GR /W3 /Wp64 /wd4996 /wd4355 /wd4150 /WL /EHsc',CPPDEFINES=["WIN32","_CRT_SECURE_NO_DEPRECATE","NOMINMAX"])
    if env['TYPE']=='debug': env.Append(CXXFLAGS=' /RTC1 /MDd',CCFLAGS=' /MDd',LINKFLAGS=' /DEBUG')
    else: env.Append(CXXFLAGS=' /MD',CCFLAGS=' /MD')
    env.Append(CXXFLAGS=' /WX',LINKFLAGS=' /wd422')
else: # assume g++...
    # machine flags
    if env['ARCH']=='pentium4': machine_flags=' -march=pentium4 -msse2 -m32 -Wa,--32'
    elif env['ARCH']=='pentium3': machine_flags=' -march=pentium3 -msse'
    elif env['ARCH']=='athlon': machine_flags=' -march=athlon-xp -msse'
    elif env['ARCH']=='nocona': machine_flags=' -march=nocona -msse3'
    elif env['ARCH']=='opteron': machine_flags=' -march=opteron -msse3'
    elif env['ARCH']=='powerpc': machine_flags=' -msse3'
    else: machine_flags=''
    env.Append(CXXFLAGS=machine_flags)
    # type specific flags
    if env['TYPE']=='release' or env['TYPE']=='optdebug' or env['TYPE']=='profile':
        optimization_flags=''
        if env['ARCH']=='pentium4': optimization_flags+=' -O2 -fexpensive-optimizations -falign-functions=4 -funroll-loops -fprefetch-loop-arrays'
        elif env['ARCH']=='pentium3': optimization_flags+=' -O2 -fexpensive-optimizations -falign-functions=4 -funroll-loops -fprefetch-loop-arrays'
        elif env['ARCH']=='opteron': optimization_flags+=' -O2'
        elif env['ARCH']=='nocona': optimization_flags+=' -O3 -funroll-loops'
        elif env['ARCH']=='powerpc': optimization_flags+=' -O2'
#        if env['fast_math']:optimization_flags+=' -ffast-math'
        optimization_flags+=' -fno-math-errno -fno-signed-zeros'
        env.Append(CXXFLAGS=optimization_flags)
        if env['TYPE']=='profile': env.Append(CXXFLAGS=' -pg',LINKFLAGS=' -pg')
    env.Append(CXXFLAGS=' -g3',LINKFLAGS=' -g')
    env.Append(CXXFLAGS=' -Wall -Werror -Winit-self -Woverloaded-virtual -Wstrict-aliasing=2 -std=gnu++14 -Wno-unknown-pragmas -Wno-strict-overflow -Wno-sign-compare')

if env['TYPE']=='release' or env['TYPE']=='profile' or env['TYPE']=='optdebug': env.Append(CPPDEFINES=['NDEBUG'])
else: env['USE_BOOSTIO']=1
if env['compile_id_types_as_int']: env.Append(CPPDEFINES=['COMPILE_ID_TYPES_AS_INT'])
if env['USE_LEX_YACC']: env.Append(CPPDEFINES=['USE_LEX_YACC'])

### teach scons about icecream environment variables
if env['CXX'].find('ice')!=-1:
    env.Replace(CXX='%s -DICECC_CXX=%s'%(env['CXX'],os.environ['ICECC_CXX']))

### darwin specific linking options
if env['PLATFORM']=='darwin':
    env.Replace(LDMODULESUFFIX='.so')

### library configuration
env.Append(CPPPATH=['#/Public_Library'])

### linker flags
extra_program_depends=[] # for win32 because dll's are actually needed for running while linking can be done with just .libs
def Link_Flags(env):
    if env['INSTALL_PATH']:
        public_library=os.path.join(env['INSTALL_PATH'],'lib')
    else:
        public_library=os.path.join('#'+variant_build,'Public_Library')
    env.Append(LIBPATH=public_library)
    if env['shared']:
        env.Append(RPATH=[Dir(public_library).abspath])
    if env['USE_OPENGL']:
        env.Append(LIBS=['PhysBAM_OpenGL'+library_suffix])
        env['USE_HYBRID']=1
    if env['USE_RAY_TRACING']:
        env.Append(LIBS=['PhysBAM_Ray_Tracing'+library_suffix])
        env['USE_DYNAMICS']=1
    if env['USE_COMPRESSIBLE'] or env['USE_INCOMPRESSIBLE'] or env['USE_DYNAMICS']:
        env.Append(LIBS=['PhysBAM_Dynamics'+library_suffix])
        env['USE_SOLIDS']=1
    if env['USE_HYBRID']:
        env.Append(LIBS=['PhysBAM_Hybrid_Methods'+library_suffix])
        env['USE_DEFORMABLES']=1
        env['USE_SOLIDS']=1
    if env['USE_SOLIDS']:
        env.Append(LIBS=['PhysBAM_Solids'+library_suffix])
        env['USE_DEFORMABLES']=1
    if env['USE_DEFORMABLES']:
        env.Append(LIBS=['PhysBAM_Deformables'+library_suffix])
        env['USE_RIGIDS']=1
    if env['USE_RIGIDS']:
        env.Append(LIBS=['PhysBAM_Rigids'+library_suffix])
        env['USE_GEOMETRY']=1
    if env['USE_GEOMETRY']:
        env.Append(LIBS=['PhysBAM_Geometry'+library_suffix])
        env['USE_GRID_PDE']=1
    if env['USE_GRID_PDE']:
        env.Append(LIBS=['PhysBAM_Grid_PDE'+library_suffix])
        env['USE_GRID_TOOLS']=1
    if env['USE_GRID_TOOLS']:
        env.Append(LIBS=['PhysBAM_Grid_Tools'+library_suffix])
        env['USE_TOOLS']=1
    if env['USE_TOOLS']:
        env.Append(LIBS=['PhysBAM_Tools'+library_suffix])
        env['USE_CORE']=1
    if env['USE_CORE']: env.Append(LIBS=['PhysBAM_Core'+library_suffix])
    for name,lib in external_libraries.items():
        if env['USE_'+name.upper()]:
            env.Append(LINKFLAGS=lib['linkflags'],LIBS=lib['libs'])
            env.PrependUnique(LIBPATH=lib['libpath'])
            if lib.has_key('rpath'): env.PrependUnique(RPATH=lib['rpath'])

### find SConscript files two levels down (for Projects and Tools)
def Find_SConscripts(env,dir):	
    for c in glob.glob(os.path.join(dir,"*","SConscript"))+(glob.glob(os.path.join(dir,"SConscript")))+glob.glob(os.path.join(dir,"*/*","SConscript")):
        env.SConscript(os.path.join(variant_build,c))

### find SConscript files two levels down (for Projects and Tools)
def Find_SConscripts_In_Subdirectories(env,dir):	
    for c in glob.glob(os.path.join(dir,"*","SConscript"))+glob.glob(os.path.join(dir,"*/*","SConscript")):
        env.SConscript(os.path.join(variant_build,c))

### find all .cpp files below current directory
def Find_Sources(dirs,ignore=[],sources=[],exclude=[]):
    local_sources=sources[:]
    def source_filter(test_file,ignore):
        if test_file.startswith('#') or (test_file.startswith('.') and not test_file.startswith('./')): return False
        test_results=map(lambda x: not test_file.endswith(x),ignore)
        test_results_exclude=map(lambda x: test_file.find(x)<0,exclude)
        type_ok=test_file.endswith(".h") or test_file.endswith(".cpp")
        if(env['USE_LEX_YACC'] and (test_file.endswith(".ll") or test_file.endswith(".yy"))): type_ok=True
        return type_ok and reduce(lambda x,y: x and y,test_results,True) and reduce(lambda x,y: x and y,test_results_exclude,True)
    build_directory=Dir('.').srcnode().abspath;build_directory_length=len(build_directory)+1
    for d in dirs:
        for root,dirs,files in os.walk(os.path.join(build_directory,d)):
            modified_root=root[build_directory_length:]
            local_sources.extend(filter(lambda filename: source_filter(filename,ignore),map(lambda file: os.path.join(modified_root,file),files)))
    cpps=filter(lambda s:s.endswith('.cpp') or s.endswith('.ll') or s.endswith('.yy'),local_sources)
    cpp_set=set(cpps)
    headers=filter(lambda s:s.endswith('.h') and not (s[:-2]+'.cpp' in cpp_set),local_sources)
    return cpps,headers

### convert sources into objects
def Automatic_Object_Helper(env,source,libraries):
    if type(source)!=str: return source # assume it's already an object
    if source.endswith('.so'): return source
    cppdefines_reversed=env['CPPDEFINES'][::-1]
    cpppath_reversed=env['CPPPATH_HIDDEN'][::-1]
    for lib in libraries:
        cppdefines_reversed.extend(lib['flags'][::-1])
        cpppath_reversed.extend(lib['cpppath'][::-1])
    if source.endswith('.h'):
        cpp=source+'.cpp'
        source_path=File(source).srcnode().path
        source_path=source_path.replace("Public_Library/","")
        env.Command(cpp,source,'echo \#include \\<'+source_path+'\\> > $TARGET')
    else: cpp=source
    if env['shared']: builder=env.SharedObject
    else: builder=env.StaticObject
    return builder(cpp,CPPDEFINES=cppdefines_reversed[::-1],CPPPATH_HIDDEN=cpppath_reversed[::-1])

def Automatic_Objects(env,sources):
    libraries=[external_libraries[name] for name in external_libraries.keys() if env['USE_'+name.upper()]]
    if type(sources)==list: return [Automatic_Object_Helper(env,source,libraries) for source in sources]
    else: return Automatic_Object_Helper(env,sources,libraries)

### automatic generation of library targets
def Automatic_Library(env,name,dirs=['.'],ignore=[],exclude=[],sources=[],link=False,prefix_install_with_lib=False,install=False,loadable_module=False):
    env_local=env.Copy()
    if env['PLATFORM']=='darwin' and not loadable_module:
        env_local.Append(LINKFLAGS='-undefined dynamic_lookup')
    cpps,headers=Find_Sources(dirs,ignore,sources,exclude)
    if env['compile_headers']:
        Automatic_Objects(env_local,headers)
    if env['PLATFORM'].startswith('win32') and env_local['TYPE']=='debug':
        env_local.Append(CXXFLAGS="  /Fd\"%s\""%os.path.join(Dir('.').path,name))
    objects=Automatic_Objects(env_local,cpps)
    if not objects: return None
    name=name+library_suffix
    if link: Link_Flags(env_local)
    if env_local['shared']:
        if env_local['use_rpath']==0 or env_local['PLATFORM']=='darwin': 
            env_local.Replace(RPATH=[])
        if loadable_module:
            ret=env_local.LoadableModule(name,source=objects)
            suffix=env_local['LDMODULESUFFIX']
        else:
            ret=env_local.SharedLibrary(name,source=objects)
            suffix=env_local['SHLIBSUFFIX']
        if env['INSTALL_PATH']:
            dest=os.path.join(env['INSTALL_PATH'],'lib','lib'+name+suffix)
            ret=env_local.InstallAs(dest,ret)
        elif install:
            prefix=""
            if prefix_install_with_lib==True: prefix="lib"
            shared_target_path=os.path.join(Dir('.').srcnode().abspath,prefix+name+suffix)
            env_local.InstallAs(shared_target_path,ret)
    else:
        ret=env_local.StaticLibrary(name,source=objects)
        if env['INSTALL_PATH']:
            dest=os.path.join(env['INSTALL_PATH'],'lib','lib'+name+env['LIBSUFFIX'])
            ret=env_local.InstallAs(dest,ret)
    return ret

def Automatic_Circular_DLLs(env,dirs):
    #print "IN AUTO %s"%repr(dirs)
    built_libraries={}
    for library in dirs:
        env_local=env.Copy()
        #print library
        library_name="PhysBAM_"+library
        if env_local['TYPE']=='debug':
            env_local.Append(CXXFLAGS="  /Fd\"%s\""%os.path.join(Dir('.').path,library_name+".pdb"))
        cpps,headers=Find_Sources([library],ignore=[],sources=[])
        if len(cpps)==0: continue
        objects=Automatic_Objects(env_local,cpps)
        def_file=env_local.DllDef(library_name+".def",objects)
        import_library=env_local.DllImportLibrary([library_name+".lib",library_name+".exp"],[def_file,objects])
        built_libraries[library]=(objects,def_file,import_library)
    
    for library,guys in built_libraries.items():
        objects,def_file,import_library=guys
        
        library_name="PhysBAM_"+library
        env_link=env.Copy()
        if env_link['TYPE']=='debug':
            env_link.Append(LINKFLAGS=" /PDB:\"%s\""%os.path.join(Dir('.').path,library_name+".pdb"))
        Link_Flags(env_link)
        env_link["no_import_lib"]=1
        other_libraries=map(lambda x: built_libraries[x][2],filter(lambda x: x!=library,built_libraries.keys()))
        other_libraries=filter(lambda x: not x.path.endswith(".exp"),Flatten(other_libraries))
        #print map(lambda x: x.path,other_libraries)
        env_link.Append(LIBS=other_libraries)
        dlls=env_link.SharedLibrary(library_name,source=[def_file,objects])
        extra_program_depends.extend(Flatten(dlls))

# Generate a wrapper on Posix 
wrapper_template_file=File('#Scripts/scons/wrapper_template')
wrapper_template=None
def Generate_Wrapper(source,target,env):
    global wrapper_template,wrapper_template_file
    if not wrapper_template:
        wrapper_template=open(wrapper_template_file.abspath).read()
    def relative_path(base,directory):
        b,d=os.path.abspath(base).split(os.sep),os.path.abspath(directory).split(os.sep)
        i=0
        while i<len(b) and i<len(d) and b[i]==d[i]: i+=1
        path=map(lambda x:'..',b[i:])+d[i:]
        if len(path)==0: return "."
        return os.path.join(*path)
    template_file,binary,paths=source
    wrapper_dir=os.path.dirname(target[0].abspath)
    paths=map(lambda x: relative_path(wrapper_dir,x),eval(paths.get_contents()))
    binary_path,binary_name=os.path.split(binary.path)
    open(target[0].abspath,"w").write(wrapper_template%(":".join(paths),relative_path(wrapper_dir,binary_path),binary_name))
    os.chmod(target[0].abspath,0755)

env.Append(BUILDERS={'Wrapper':Builder(action=Generate_Wrapper)})

# Generate DEF/DLL files for Windows
def DllDef_File(source,target,env):
    objs=source
    library_name=target[0].path.replace(".def","")
    print "Generating DEF file for %s"%library_name
    # use dumpbin to figure out what to export
    dumpbin="dumpbin"
    if env['PLATFORM']=="win32cross":
        dumpbin="dumpbinwrap"

    cmd="%s /symbols %s"%(dumpbin," ".join(map(lambda x: x.path,objs)))
    print cmd
    fp=os.popen(cmd)
    defs={}
    while 1:
        line=fp.readline()
        if line=="": break
        line=line.replace("\r","").replace("\n","")
        if line.find("External")!=-1:
            left,right=line.split(" | ")
            if left.find("External")==-1: continue # one more check!
            if left.find("UNDEF")!=-1: continue
            mangled_symbol=right.split(" ")[0]
            if mangled_symbol.find("@PhysBAM")!=-1: # only export PhysBAM namespace symbols
                defs[mangled_symbol]=True
    print "%s has %d symbols exported"%(library_name,len(defs.keys()))
    deffp=open(target[0].abspath,"w")
    deffp.write("LIBRARY %s\n"%os.path.basename(library_name))
    deffp.write("EXPORTS\n")
    for item in defs.keys():
        deffp.write(item+"\n")
    deffp.close()
    return None

env.Append(BUILDERS={"DllDef":Builder(action=DllDef_File)})
env.Append(BUILDERS={"DllImportLibrary":Builder(action="${AR} /out:${TARGETS[0]} /def:${SOURCES[0]} ${SOURCES[1:]} ")})
#env.Append(BUILDERS={"Dll":Builder(action="${AR} /out:${TARGETS[0]} /def:${SOURCES[0]} ${SOURCES[1:]} ")})

# Builder(action=DllImportLibrary)})


### automatic generation of program targets
def Automatic_Program(env,name='',sources=None):
    if name=='':
        dir=Dir('.').path
        name=os.path.basename(dir)
        if sources==None: sources=['main.cpp']
    elif sources==None: sources=[name+'.cpp']
    env_compile=env
    if env_compile['PLATFORM'].startswith('win32') and env_compile['TYPE']=='debug':
        env_compile=env_compile.Copy()
        env_compile.Append(CXXFLAGS="  /Fd\"%s\""%os.path.join(Dir('.').path,"vc80.pdb"))
    objects=Automatic_Objects(env_compile,sources)
    env_link=env_compile.Copy()
    Link_Flags(env_link)
    if env['PLATFORM']=='darwin':
        env_link.Append(LINKFLAGS='-bind_at_load')
    install_name=os.path.basename(name)+program_suffix+executable_suffix
    if env['INSTALL_PATH']:
        executable_target_path=os.path.join(env['INSTALL_PATH'],'bin',install_name)
    else:
        executable_target_path=os.path.join(Dir('.').srcnode().abspath,install_name)
    ### Mac OS X does not support -rdynamic    
    if sys.platform!="darwin" and not env['shared']: env_link.Append(LINKFLAGS='-rdynamic')
    if env['install_programs']:
        if env['shared'] and env['wrapper']:
            rpath_save=env_link['RPATH']
            env_link.Replace(RPATH=[])
            program=env_link.Program(name,objects)
            for i in extra_program_depends: env_link.Depends(program,i)
            if env_link['INSTALL_PATH']:
                program,=env_link.InstallAs(os.path.join(env['INSTALL_PATH'],'bin','physbam-raw',install_name),program)
            env_link.Wrapper(executable_target_path,source=[wrapper_template_file,program,Value(rpath_save)])
        else:
            program=env_link.Program(name,objects)
            for i in extra_program_depends: env_link.Depends(program,i)
            env_link.InstallAs(executable_target_path,program)
        env.Depends('.',executable_target_path)
    else:
        env_link.Program(name,objects)

### automatic generation of a global library that links to a bunch of sub libraries
def Automatic_Global_Library(env,name,source_libs):
    env_global=env.Copy()
    Link_Flags(env_global)
    if env_global["use_rpath"]==0 or env_global['PLATFORM']=='darwin': env_global.Replace(RPATH=[])
    return env_global.SharedLibrary(name,source_libs)

### build everything
Export('env Automatic_Library Automatic_Program Automatic_Circular_DLLs Automatic_Global_Library Find_Sources variant_build')
env.SConscript(variant_build+'/Public_Library/SConscript')
Find_SConscripts_In_Subdirectories(env,'Projects')
Find_SConscripts(env,'Tests')
Find_SConscripts(env,'Tools')
if os.path.exists('Scripts/SConscript'):
    env.SConscript(variant_build+'/Scripts/SConscript')
