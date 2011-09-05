#!/usr/bin/env python
######################################################################
# Copyright 2007, Geoffrey Irving, Andrew Selle.
# This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
######################################################################
# TODO: make force classes use custodian_and_ward when they objects/meshes

import os
import re
import sys
import optparse

######################################################################
# Arguments
######################################################################

parser=optparse.OptionParser('usage: %prog [options] <output dir>')
parser.add_option('--no-double',action='store_true',help='generate float-only bindings')
(options,args)=parser.parse_args()
if len(args)!=1: parser.error('expected exactly one argument')
output_directory,=args
use_double=not options.no_double

######################################################################

class_to_file_map={}
def register_class_include(file):
    match=re.match('(PhysBAM_\w+/)?(\w+)/(\w+)\.h$',file)
    if match: class_to_file_map[match.group(3)]=file

def parse_template(instance):
    start=instance.find('<')
    template=instance[:start].strip()
    arguments=[]
    nest=0
    for i in range(start,len(instance)):
        c=instance[i]
        if c=='<':
            nest=nest+1
        elif c=='>':
            nest=nest-1;
            assert(nest>=0)
            if nest==0:
                arguments.append(instance[start+1:i].strip())
        elif c==',':
            if nest==1:
                arguments.append(instance[start+1:i].strip())
                start=i
    assert(nest==0)
    for i,a in enumerate(arguments):
        try:
            arguments[i]=int(a)
        except ValueError:
            pass
    return template,arguments

######################################################################
# Class CODE
######################################################################
class CODE:
    def __init__(self,function,main_module=False):
        self.cpp_filename=os.path.join(output_directory,function+'.cpp')
        self.includes=set()
        self.routines=[]
        self.defs=[]
        self.reference_count=0
        if main_module==True:
            self.include('PhysBAM_Tools/Utilities/INTERRUPTS.h')
            self.definition="BOOST_PYTHON_MODULE(libphysbam_python_internal)"
            self.defs.append('scope().attr("__name__") = "physbam";\n')
            self.routines.append('void Check_For_Keyboard_Interrupts(){if(PyErr_CheckSignals()) boost::python::throw_error_already_set();}\n')
            self.defs.append('PhysBAM::Add_Interrupt_Checker(Check_For_Keyboard_Interrupts);\n')
            self.f=None
        else:
            self.f=function
            self.definition="void %s()"%self.f

    def prototype(self):
        return "void %s();\n"%self.f

    def call(self):
        return "%s();\n"%self.f

    ######################################################################
    # Function write
    ######################################################################
    def write(self):
        fp=open(self.cpp_filename,"w")
        fp.write("#include <boost/python.hpp>\n")
        fp.writelines(self.includes)
        if self.f:
            fp.write("""
namespace PhysBAM{}
using namespace PhysBAM;
//######################################################
""")

        fp.writelines(self.routines)
        fp.write("//######################################################\n")
        fp.write("%s{\nusing namespace boost::python;\n"%self.definition)
        fp.writelines(self.defs)
        fp.write("}\n")

    ######################################################################
    # Function typedef
    ######################################################################
    def typedef(self,**defs):
        for name,T in defs.items():
            self.defs.append('typedef %s %s;\n'%(T,name))

    ######################################################################
    # Function include
    ######################################################################
    def raw_include(self,file):
        self.includes.add('#include <%s>\n'%file)

    def include(self,file):
        register_class_include(file)
        self.raw_include(file)

    def include_class(self,cppclass,file=None):
        if file:
            class_to_file_map[cppclass]=file
            self.raw_include(file)
        elif isinstance(cppclass,str):
            if class_to_file_map.has_key(cppclass): self.raw_include(class_to_file_map[cppclass])
            elif cppclass.find('<')!=-1:
                template,arguments=parse_template(cppclass)
                if not class_to_file_map.has_key(template):
                    raise RuntimeError("Must specify file for unknown class %s"%template)
                self.raw_include(class_to_file_map[template])
                for a in arguments:
                    self.include_class(a)

    def local_include(self,file):
        self.includes.add('#include "%s"\n'%file)

    ######################################################################
    # Function member_function
    ######################################################################
    def member_function(self,cppclass,spec,overloads=None):
        # TODO: make an interface for specifying the member template using parse template
        if type(spec)!=tuple: spec=spec,
        name=spec[0]
        function="&%s::%s"%(cppclass,name)
        policy=""
        if len(spec)>=3:
            assert(len(spec)<=5)
            return_T,argument_T=spec[1:3]
            if len(spec)>=4: modifier=" "+spec[3]
            else: modifier=""
            if len(spec)==5: policy=","+spec[4]
            function="static_cast<%s (%s::*)(%s)%s>(%s)"%(return_T,cppclass,(",".join(argument_T)),modifier,function)
        elif len(spec)>=2: policy=","+spec[-1]
        if overloads:
            overload_min,overload_max=overloads
            overload_name='%s_overloads_%d'%(name,self.reference_count)
            self.reference_count+=1
            self.routines.append('BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(%s,%s,%d,%d)\n'%(overload_name,name,overload_min,overload_max))
            overload=',%s()'%overload_name
        else: overload=''
        self.defs.append("    .def(\"%s\",%s%s%s)\n"%(name,function,policy,overload))

    ######################################################################
    # Function static_function
    ######################################################################
    def static_function(self,cppclass,spec,decorate=True,overloads=None):
        if type(spec)!=tuple: spec=spec,
        name=spec[0]
        function="&%s::%s"%(cppclass,name)
        policy=""
        if len(spec)>=3:
            assert(len(spec)<=5)
            return_T,argument_T=spec[1:3]
            if len(spec)>=4: modifier=" "+spec[3]
            else: modifier=""
            if len(spec)==5: policy=","+spec[4]
            function="static_cast<%s (*)(%s)%s>(%s)"%(return_T,(",".join(argument_T)),modifier,function)
        elif len(spec)>=2: policy=","+spec[-1]
        if overloads:
            overload_min,overload_max=overloads
            overload_name='%s_overloads_%d'%(name,self.reference_count)
            self.reference_count+=1
            self.routines.append('BOOST_PYTHON_FUNCTION_OVERLOADS(%s,%s::%s,%d,%d)\n'%(overload_name,cppclass,name,overload_min,overload_max))
            overload=',%s()'%overload_name
        else: overload=''
        self.defs.append("    .def(\"%s\",%s%s%s)\n"%(name,function,policy,overload))
        if decorate:
            decoration='      .staticmethod("%s")\n'%name
            for i in range(len(self.defs)-1,-1,-1): # remove matching decorations
                if 'class_<' in self.defs[i]:
                    break
                elif self.defs[i]==decoration:
                    del self.defs[i]
                    break
            self.defs.append(decoration)

    ######################################################################
    # Function static_functions
    ######################################################################
    def static_functions(self,cppclass,specs):
        names=set()
        for spec in specs:
            if type(spec)!=tuple: spec=spec,
            names.add(spec[0])
            self.static_function(cppclass,spec,decorate=False)
        for n in names:
            self.defs.append("      .staticmethod(\"%s\")\n"%n)

    ######################################################################
    # Function member_variable
    ######################################################################
    def member_variable(self,cppclass,name):
        access="readwrite"
        if type(name)==tuple:
            name,access=name
        self.defs.append("    .def_%s(\"%s\",&%s::%s)\n"%(access,name,cppclass,name))

    ######################################################################
    # operators
    ######################################################################
    def equality_operators(self):
        self.defs.extend(["    .def(self==self)\n",
                          "    .def(self!=self)\n"])

    def vector_operators(self,scalar_type):
        self.equality_operators()
        self.defs.extend(["    .def(self+self)\n",
                          "    .def(self-self)\n",
                          "    .def(self*%s())\n"%scalar_type,
                          "    .def(%s()*self)\n"%scalar_type,
                          "    .def(self/%s())\n"%scalar_type,
                          "    .def(self*=%s())\n"%scalar_type,
                          "    .def(self/=%s())\n"%scalar_type,
                          "    .def(self+=self)\n",
                          "    .def(self-=self)\n",
                          "    .def(-self)\n"])

    def hack_vector_operators(self,scalar_type):
        self.vector_operators(scalar_type)
        self.defs.extend(["    .def(self*=self)\n",
                          "    .def(self/=self)\n"])

    ######################################################################
    # Function reference
    ######################################################################
    def reference(self,cppclass,member):
        self.local_include("UTILITIES.h")
        self.includes.add("#include <boost/type_traits/remove_pointer.hpp>\n")
        self.routines.append("boost::remove_pointer<typeof(((%s*)0)->%s)>::type& Get_Field_%d(%s& x){return Dereference_If_Pointer(x.%s);}\n"
            %(cppclass,member,self.reference_count,cppclass,member))
        self.defs.append("    .add_property(\"%s\",make_function(&::Get_Field_%d,return_internal_reference<>()))\n"%(member,self.reference_count))
        self.reference_count+=1

    ######################################################################
    # Function read_write
    ######################################################################
    def read_write(self,cppclass,T):
        self.defs.append("def(\"Read\",Read_Write<%s,%s >::Read);\n"%(cppclass,T))
        self.defs.append("def(\"Write\",Read_Write<%s,%s >::Write);\n"%(cppclass,T))

    ######################################################################
    # Function str
    ######################################################################
    def str(self,cppclass):
        self.defs.append("    .def(\"__str__\",&boost::lexical_cast<std::string,%s >)\n"%cppclass)

    ######################################################################
    # Function repr
    ######################################################################
    def repr(self,cppclass):
        self.local_include("UTILITIES.h")
        self.defs.append("    .def(\"__repr__\",&Repr<%s >)\n"%cppclass)

    ######################################################################
    # Function numpy
    ######################################################################
    def numpy(self,cppclass):
        self.local_include("NUMPY.h")
        self.defs.append('    .def("__array__",&As_Numpy<%s >)\n'%cppclass)
        self.defs.append('    .def("__array__",&As_Numpy_With_Context<%s >)\n'%cppclass)

    ######################################################################
    # Function exceptions
    ######################################################################
    def exceptions(self,exception_map):
        self.include('PhysBAM_Tools/Utilities/EXCEPTIONS.h')
        self.routines.append('\n'.join([
            'void Translate(const PHYSBAM_ERROR& error)',
            '{',
            '    PyObject* exc=PyExc_RuntimeError;',
            '         '+'\n    else '.join('if(typeid(error)==typeid(%s)) exc=PyExc_%s;'%(cpp,py) for cpp,py in exception_map.items()),
            '    PyErr_SetString(exc,error.what());',
            '}\n']))
        self.defs.append('register_exception_translator<PHYSBAM_ERROR>(&Translate);\n')
    
    ######################################################################
    # Function cloneable_base
    ######################################################################
    def cloneable_base(self):
        self.include('PhysBAM_Tools/Clone/CLONEABLE.h')
        self.defs.append('class_<CLONEABLE_BASE,boost::noncopyable>("CLONEABLE_BASE",no_init)\n')
        self.defs.append('    .def("__copy__",&CLONEABLE_BASE::Clone,return_value_policy<manage_new_object>())\n')
        self.defs.append(";\n")

    ######################################################################
    # Function numeric_limits_instance
    ######################################################################
    def numeric_limits_instance(self,name,T):
        self.includes.add('#include <limits>\n')

        cppclass="std::numeric_limits<%s>"%T
        self.defs.append("class_<%s >(\"%s\")\n"%(cppclass,name)) # class

        # make list of functions
        functions=['min','max','epsilon','round_error','infinity','quiet_NaN','signaling_NaN','denorm_min']
        for f in functions:
            self.static_function(cppclass,f)
        self.defs.append(";\n")

    ######################################################################
    # Function vector_instance
    ######################################################################
    def vector_instance(self,name,T,d):
        self.include("PhysBAM_Tools/Vectors/VECTOR.h")
        self.local_include("INDEXING.h")
        self.local_include("CONVERSIONS.h")

        dot_product_wrapper="template<class T_VECTOR> typename T_VECTOR::SCALAR Dot_Product_Wrapper(const T_VECTOR& u,const T_VECTOR& v){return T_VECTOR::Dot_Product(u,v);}\n"
        if dot_product_wrapper not in self.routines:
            self.routines.append(dot_product_wrapper)

        cppclass="VECTOR<%s,%d>"%(T,d)
        self.defs.append("class_<%s >(\"%s\")\n"%(cppclass,name)) # class
        self.defs.append("    .def(init<const %s&>())\n"%cppclass) # copy constructor
        self.defs.append("    .def(init<%s >())\n"%(",".join(d*[T]))) # constructor
        # members
        if d<=3:
            for i in ["x","y","z"][:d]:
                self.defs.append("    .def_readwrite(\"%s\",&%s::%s)\n"%(i,cppclass,i))
        # make list of functions
        functions=["Find","Contains","Remove_Index"]
        self.str(cppclass)
        self.repr(cppclass)
        self.numpy(cppclass)
        # add things if not bool
        if T!="bool" and (T=="float" or T=="double"):
            # add operators
            self.hack_vector_operators(T)
            non_bool_functions=["Min","Max","Max_Abs","Elements_Equal","Dominant_Axis",
                                "Componentwise_Min","Componentwise_Max","Sum","Average","Product",
                                "Horizontal_Vector","Sorted"]
            functions.extend(non_bool_functions)
            self.defs.append("    .def(\"Dot_Product\",&Dot_Product_Wrapper<%s >)\n      .staticmethod(\"Dot_Product\")\n"%cppclass)
            for f in ['Axis_Vector','All_Ones_Vector']:
                self.static_function(cppclass,f)
        # add things if float or double
        if T=="float" or T=="double":
            scalar_only_functions=["Magnitude_Squared","Magnitude","Lp_Norm","L1_Norm","Normalize",
                                   "Normalized","Orthogonal_Vector","Unit_Orthogonal_Vector",
                                   "Projected_On_Unit_Direction","Projected","Project_On_Unit_Direction","Project",
                                   "Projected_Orthogonal_To_Unit_Direction","Project_Orthogonal_To_Unit_Direction","Cross_Product",
                                   "Triple_Product","Angle_Between"]
            functions.extend(scalar_only_functions)
        # exclude things that 1d and 2d don't have
        dimensional_excludes={1:["Horizontal_Vector","Dominant_Axis","Remove_Index","Sorted","Max_Abs","Lp_Norm","Orthogonal_Vector",
                              "Unit_Orthogonal_Vector","Project_Orthogonal_To_Unit_Direction","Project_Orthogonal_To_Unit_Direction","Triple_Product","Cross_Product",
                              "Angle_Between","Projected_Orthogonal_To_Unit_Direction"],
                              2:["Orthogonal_Vector","Project_Orthogonal_To_Unit_Direction","Orthogonal_Vector","Unit_Orthogonal_Vector",
                                 "Project_Orthogonal_To_Unit_Direction","Triple_Product","Cross_Product"],
                              4:["Horizontal_Vector","Angle_Between","Orthogonal_Vector","Project_Orthogonal_To_Unit_Direction","Orthogonal_Vector","Unit_Orthogonal_Vector",
                                 "Project_Orthogonal_To_Unit_Direction","Triple_Product","Cross_Product"]}
        if dimensional_excludes.has_key(d): functions=filter(lambda x: not x in dimensional_excludes[d],functions)
        # actually add them
        for member in functions:
            self.member_function(cppclass,member)

        if not T.startswith("VECTOR"):
            self.defs.append("    .def(ARRAY_INDEXING_SUITE<%s >())\n"%cppclass)
            # terminate
            self.defs.append(";\n")
            self.defs.append("Register_Array_Conversion<%s >();\n"%cppclass)
        else:
            self.defs.append(";\n")

    ######################################################################
    # Function matrix_instances
    ######################################################################
    def matrix_instances(self,T,templates):
        self.include("PhysBAM_Tools/Vectors/ARITHMETIC_POLICY.h")
        self.include("PhysBAM_Tools/Matrices/VECTOR_POLICY.h")
        self.include("PhysBAM_Tools/Vectors/VECTOR.h")
        self.include("PhysBAM_Tools/Matrices/MATRIX.h")
        self.include("PhysBAM_Tools/Matrices/MATRIX_MXN.h")
        self.include("PhysBAM_Tools/Matrices/UPPER_TRIANGULAR_MATRIX_2X2.h")
        self.include("PhysBAM_Tools/Matrices/UPPER_TRIANGULAR_MATRIX_3X3.h")
        self.include("PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h")
        self.include("PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h")
        self.include("PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_2X2.h")
        self.include("PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_3X3.h")
        self.includes.add('#include <boost/python/tuple.hpp>\n')

        routine="""
template<class T_MATRIX> typename T_MATRIX::SCALAR Get_Item(const T_MATRIX& matrix,const boost::python::tuple& index)
{int i=boost::python::extract<int>(index[0]),j=boost::python::extract<int>(index[1]);
if(!(1<=i && i<=matrix.Rows() && 1<=j && j<=matrix.Columns())){
    PyErr_SetString(PyExc_IndexError,"Index out of range");
    boost::python::throw_error_already_set();}
return matrix.Valid_Index(i,j)?matrix(i,j):0;}

template<class T_MATRIX> typename MATRIX_INFO<T_MATRIX>::LEFT_VECTOR Get_Column(const T_MATRIX& matrix,const int i)
{typename MATRIX_INFO<T_MATRIX>::LEFT_VECTOR r(INITIAL_SIZE(matrix.Rows()));
for(int j=1;j<=matrix.Rows();j++) r(j)=matrix(j,i);
return r;}

template<class T_MATRIX> void Set_Item(T_MATRIX& matrix,const boost::python::tuple& index,const typename T_MATRIX::SCALAR& value)
{int i=boost::python::extract<int>(index[0]),j=boost::python::extract<int>(index[1]);
if(!(1<=i && i<=matrix.Rows() && 1<=j && j<=matrix.Columns())){
    PyErr_SetString(PyExc_IndexError,"Index out of range");
    boost::python::throw_error_already_set();}
else if(!matrix.Valid_Index(i,j)){
    PyErr_SetString(PyExc_IndexError,"Index is immutable (value is always zero)");
    boost::python::throw_error_already_set();}
matrix(i,j)=value;}

template<class T,int d> boost::python::tuple Solve_Eigenproblem_Wrapper(const SYMMETRIC_MATRIX<T,d>& matrix)
{typename VECTOR_POLICY<VECTOR<T,d> >::DIAGONAL_MATRIX D;MATRIX<T,d> V;
matrix.Solve_Eigenproblem(D,V);return boost::python::make_tuple(D,V);}

template<class T,int d> boost::python::tuple Fast_Solve_Eigenproblem_Wrapper(const SYMMETRIC_MATRIX<T,d>& matrix)
{typename VECTOR_POLICY<VECTOR<T,d> >::DIAGONAL_MATRIX D;MATRIX<T,d> V;
matrix.Fast_Solve_Eigenproblem(D,V);return boost::python::make_tuple(D,V);}

template<class T,int m,int n> boost::python::tuple Fast_Singular_Value_Decomposition_Wrapper(const MATRIX<T,m,n>& matrix)
{static const int k=m<n?m:n;MATRIX<T,m,k> U;typename VECTOR_POLICY<VECTOR<T,k> >::DIAGONAL_MATRIX D;MATRIX<T,n,k> V;
matrix.Fast_Singular_Value_Decomposition(U,D,V);return boost::python::make_tuple(U,D,V);}

"""
        if routine not in self.routines: self.routines.append(routine)

        def template_to_class(C,m,n,T2=None):
            if not T2: T2=T
            if C=='MATRIX_MXN': return '%s<%s>'%(C,T2)
            elif m==n: return '%s<%s,%d>'%(C,T2,m)
            else: return '%s<%s,%d,%d>'%(C,T2,m,n)

        mn_pattern=re.compile(r'\w+<\w+(?:,(\d+))?(?:,(\d+))?>$')
        def class_to_mn(c):
            match=mn_pattern.match(c)
            m,n=match.group(1),match.group(2)
            if not m: m=-1
            if not n: n=m
            return int(m),int(n)

        classes=[template_to_class(*t) for t in templates]

        for C,m,n in templates:
            dynamic=C=='MATRIX_MXN'
            regular=C.startswith('MATRIX')
            diagonal=C=='DIAGONAL_MATRIX' or m==n==1
            symmetric=C=='SYMMETRIC_MATRIX' or diagonal
            upper_triangular=C=='UPPER_TRIANGULAR_MATRIX' or diagonal

            cppclass=template_to_class(C,m,n)
            if dynamic: name='%s_%s'%(C,T[0])
            elif m==n: name='%s_%s%d'%(C,T[0],m)
            else: name='%s_%s%d%d'%(C,T[0],m,n)

            if dynamic: dofs=-1
            elif regular: dofs=m*n
            elif diagonal: dofs=m
            elif symmetric: dofs=m*(m+1)/2
            elif upper_triangular: dofs=m*(m+1)/2
            left_vector='MATRIX_INFO<%s >::LEFT_VECTOR'%cppclass
            right_vector='MATRIX_INFO<%s >::RIGHT_VECTOR'%cppclass

            functions=[]
            static_functions=[]

            # header
            self.defs.append("class_<%s >(\"%s\")\n"%(cppclass,name)) # class
            # constructors
            if dofs>0:
                self.defs.append("    .def(init<%s >())\n"%(",".join(dofs*[T])))
            if dynamic:
                self.defs.append("    .def(init<int,int>())\n")
                for other in classes:
                    if other!=cppclass:
                        self.defs.append("    .def(init<const %s&>())\n"%other)
            elif regular:
                self.defs.append("    .def(init<%s >())\n"%(",".join(n*['const %s&'%left_vector])))
                self.defs.append("    .def(init<const MATRIX_MXN<%s >&>())\n"%T)
            elif diagonal:
                self.defs.append("    .def(init<%s >())\n"%(",".join(['const VECTOR<%s,%d>&'%(T,m)])))
            self.defs.append("    .def(init<INITIAL_SIZE,INITIAL_SIZE>())\n")
            if T=='float': T_other='double'
            else: T_other='float'
            self.defs.append("    .def(init<const %s&>())\n"%template_to_class(C,m,n,T_other))
            self.defs.append("    .def(init<const %s&>())\n"%cppclass) # copy constructor
            # basics
            functions.extend(['Rows','Columns'])
            self.defs.append("    .def(\"__getitem__\",&Get_Item<%s >)\n"%cppclass)
            self.defs.append("    .def(\"__setitem__\",&Set_Item<%s >)\n"%cppclass)
            self.defs.append("    .def(\"Column\",&Get_Column<%s >)\n"%cppclass)
            self.str(cppclass)
            if regular:
                self.numpy(cppclass)
            if not dynamic: self.repr(cppclass)
            # arithmetic
            self.vector_operators(T)
            self.defs.append("    .def(self*%s())\n"%(right_vector))
            if m==n:
                self.defs.extend(["    .def(self*self)\n"])
                if not dynamic:
                    self.defs.extend(["    .def(self+%s())\n"%T,
                                      "    .def(self-%s())\n"%T,
                                      "    .def(%s()+self)\n"%T,
                                      "    .def(%s()-self)\n"%T,
                                      "    .def(self+=%s())\n"%T,
                                      "    .def(self-=%s())\n"%T])
                if (regular or upper_triangular) and not dynamic:
                    self.defs.append("    .def(self*=self)\n")
            for other in classes:
                m2,n2=class_to_mn(other)
                other_dynamic=other.startswith('MATRIX_MXN')
                if other!=cppclass and m==n==m2==n2:
                    self.defs.extend(["    .def(self+%s())\n"%other,
                                      "    .def(self-%s())\n"%other])
                if other!=cppclass and (n==m2 or dynamic or other_dynamic):
                    self.defs.append("    .def(self*%s())\n"%other)
                if symmetric and (not diagonal or m==1) and m==m2 and other.startswith('MATRIX'):
                    for center in classes:
                        if (center.startswith('SYMMETRIC') or center.startswith('DIAGONAL') or class_to_mn(center)==(1,1)) and class_to_mn(center)[0]==n2:
                            static_functions.append(('Conjugate',cppclass,['const %s&'%other,'const %s&'%center]))
                #if m==m2 or dynamic or other.startswith('MATRIX_MXN'):
                #    self.defs.append("    .def(\"Transpose_Times\",&Transpose_Times_Wrapper<%s,%s >)\n"%(cppclass,other))
                #if n==n2 or dynamic or other.startswith('MATRIX_MXN'):
                #    self.defs.append("    .def(\"Times_Transpose\",&Times_Transpose_Wrapper<%s,%s >)\n"%(cppclass,other))
            # other
            functions.extend(['Max_Abs','Frobenius_Norm','Frobenius_Norm_Squared'])
            if m==n and 0<=m<=3 or m==3 and n==2:
                functions.append('Cofactor_Matrix')
            if regular or symmetric:
                if m==2 and n==3:
                    functions.append(('Transposed','const TRANSPOSE<%s >::TYPE&'%cppclass,[],'const','return_value_policy<copy_const_reference>()'))
                else:
                    functions.append('Transposed')
                if m==n:
                    functions.append(('Transpose','void',[]))
            if m==n:
                functions.extend(['Trace'])
                if not dynamic:
                    functions.append(('Solve_Linear_System',right_vector,['const %s&'%(left_vector)],"const"))
                    functions.extend(['Inverse','Determinant'])
                    if diagonal or not upper_triangular:
                        functions.append('Robust_Solve_Linear_System')
                static_functions.extend(['Identity_Matrix'])
            if symmetric:
                functions.extend(['Positive_Definite','Positive_Semidefinite','Positive_Definite_Part'])
                for f in ['log','exp']:
                    self.defs.append("    .def(\"%s\",static_cast<%s (*)(const %s&)>(%s))\n"%(f,cppclass,cppclass,f))
                if not diagonal:
                    functions.extend(['Largest_Column'])
                    static_functions.extend([('Outer_Product',cppclass,['const VECTOR<%s,%d>&'%(T,m)]),
                        'Unit_Matrix','First_Eigenvector_From_Ordered_Eigenvalues','Last_Eigenvector_From_Ordered_Eigenvalues','Fast_Eigenvalues'])
                    self.defs.extend("    .def(\"Solve_Eigenproblem\",&Solve_Eigenproblem_Wrapper<%s,%d>)\n"%(T,m))
                    self.defs.append("    .def(\"Fast_Solve_Eigenproblem\",&Fast_Solve_Eigenproblem_Wrapper<%s,%d>)\n"%(T,m))
            if diagonal:
                functions.extend(['Sqrt','Min','Max','Abs','Sign','To_Vector','Clamp_Min','Clamp_Max'])
            if m==n and not dynamic and not upper_triangular:
                functions.append('Diagonal_Part')
            if regular and not dynamic:
                if m>=n: functions.extend(['Normal_Equations_Matrix'])
                static_functions.extend([('Outer_Product',cppclass,['const VECTOR<%s,%d>&'%(T,m),'const VECTOR<%s,%d>&'%(T,n)])])
                self.defs.append("    .def(\"Fast_Singular_Value_Decomposition\",Fast_Singular_Value_Decomposition_Wrapper<%s,%d,%d>)\n"%(T,m,n))
            if m==n and m>=2 and not dynamic and not diagonal and (regular or upper_triangular):
                functions.append('Simplex_Minimum_Altitude')
            # finish
            for member in functions:
                self.member_function(cppclass,member)
            self.static_functions(cppclass,static_functions)
            self.defs.append(";\n")

    ######################################################################
    # Function shape_helper
    ######################################################################
    def shape_helper(self,cppclass):
        self.include("PhysBAM_Geometry/Geometry/TRIANGULATED_SURFACE.h")
        c=self.reference_count
        self.reference_count+=1
        self.routines.append('bool Lazy_Inside_%d(const %s& shape,const %s::VECTOR_T& X){return shape.Lazy_Inside(X);}\n'%(c,cppclass,cppclass))
        for f in [('Inside','bool',['const TV&','T'],'const'),
                  ('Outside','bool',['const TV&','T'],'const'),
                  ('Lazy_Inside','bool',['const TV&'],'const'),
                  ('Lazy_Outside','bool',['const TV&'],'const'),
                  ('Normal','TV',['const TV&'],'const'),
                  'Boundary','Surface','Signed_Distance','Principal_Curvatures']:
            self.member_function(cppclass,f)
        if not cppclass.startswith('BOX<'):
            self.member_function(cppclass,'Bounding_Box')
        self.defs.append("    .def(\"__contains__\",Lazy_Inside_%d)\n"%c)
        self.repr(cppclass)

    ######################################################################
    # Function range
    ######################################################################
    def range(self,name,TV):
        self.defs.append('{\n')
        (T,d)=parse_template(TV)[1]
        d=int(d)
        self.include("PhysBAM_Tools/Math_Tools/RANGE.h")
        self.include("PhysBAM_Tools/Grids_Uniform/GRID_%dD.h"%d)
        self.typedef(T=T,TV=TV)

        cppclass="RANGE<%s >"%(TV)
        self.defs.append('class_<%s >("%s")\n'%(cppclass,name))
        self.defs.append("    .def(init<%s >())\n"%(",".join(2*d*[T])))
        self.defs.append("    .def(init<%s,%s >())\n"%(TV,TV))
        self.defs.append("    .def(init<%s >())\n"%(TV))
        # members
        for m in ["min_corner","max_corner"]:
            self.defs.append("    .def_readwrite(\"%s\",&%s::%s)\n"%(m,cppclass,m))
        # static functions
        statics=['Unit_Box','Zero_Box','Combine','Intersect',('Bounding_Box',cppclass,['const ARRAY<%s >&'%TV])]
        if T in ['float','double']:
            statics.extend(['Empty_Box','Full_Box'])
        self.static_functions(cppclass,statics)
        # member functions
        functions=['Empty','Edge_Lengths','Center','Minimum_Corner','Maximum_Corner','Size','Reset_Bounds','Enlarge_To_Include_Point','Enlarge_Nonempty_Box_To_Include_Point',
            'Enlarge_To_Include_Box','Thickened']
        self.str(cppclass)
        self.vector_operators(T)
        # actually add them
        for member in functions:
            self.member_function(cppclass,member)
        # terminate
        self.defs.append(";}\n")

    ######################################################################
    # Function box
    ######################################################################
    def box(self,name,TV):
        self.defs.append('{\n')
        (T,d)=parse_template(TV)[1]
        d=int(d)
        self.include("PhysBAM_Geometry/Geometry/BOX.h")
        self.include("PhysBAM_Tools/Grids_Uniform/GRID_%dD.h"%d)
        self.typedef(T=T,TV=TV)

        cppclass="BOX<%s >"%(TV)
        self.defs.append('class_<%s,bases<RANGE<TV> > >("%s")\n'%(cppclass,name))
        self.defs.append("    .def(init<%s >())\n"%(",".join(2*d*[T])))
        self.defs.append("    .def(init<%s,%s >())\n"%(TV,TV))
        self.defs.append("    .def(init<%s >())\n"%(TV))
        # member functions
        self.shape_helper(cppclass)
        functions=[]
        if d==3: functions.append('Surface_Area')
        self.str(cppclass)
        self.vector_operators(T)
        # actually add them
        for member in functions:
            self.member_function(cppclass,member)
        # terminate
        self.defs.append(";}\n")

    ######################################################################
    # Function plane
    ######################################################################
    def plane(self,name,TV):
        self.defs.append('{\n')
        (T,d)=parse_template(TV)[1]
        d=int(d)
        self.include("PhysBAM_Geometry/Geometry/PLANE.h")
        self.typedef(T=T)
        self.typedef(TV='VECTOR<T,3>')

        cppclass="PLANE<%s >"%T
        self.defs.append('class_<%s >("%s")\n'%(cppclass,name))
        self.defs.append("    .def(init<%s,%s >())\n"%(TV,TV))
        self.defs.append("    .def(init<%s,%s,%s >())\n"%(TV,TV,TV))
        # members
        for m in ["normal","x1"]:
            self.defs.append("    .def_readwrite(\"%s\",&%s::%s)\n"%(m,cppclass,m))
        # member functions
        self.shape_helper(cppclass)
        self.member_function(cppclass,('Normal',TV,[],'const'))
        # terminate
        self.defs.append(";}\n")
        self.read_write(cppclass,T)

    ######################################################################
    # Function sphere
    ######################################################################
    def sphere(self,name,TV):
        self.defs.append('{\n')
        (T,d)=parse_template(TV)[1]
        d=int(d)
        self.include("PhysBAM_Geometry/Geometry/RAY.h")
        self.include("PhysBAM_Geometry/Geometry/SPHERE.h")
        self.typedef(T=T,TV=TV)

        routine='template<class TV> SPHERE<TV> Bounding_Sphere_Wrapper(const ARRAY<TV>& X){return SPHERE<TV>::Bounding_Sphere(X);}\n'
        if routine not in self.routines: self.routines.append(routine)

        cppclass="SPHERE<%s >"%(TV)
        self.defs.append('class_<%s >("%s")\n'%(cppclass,name))
        self.defs.append("    .def(init<%s,%s >())\n"%(TV,T))
        # members
        for m in ["center","radius"]:
            self.defs.append("    .def_readwrite(\"%s\",&%s::%s)\n"%(m,cppclass,m))
        # member functions
        self.shape_helper(cppclass)
        functions=['Size']
        if d==2: functions.extend(['Circular_Segment_Area'])
        for member in functions:
            self.member_function(cppclass,member)
        # terminate
        self.defs.append(";}\n")
        self.read_write(cppclass,T)

    ######################################################################
    # Function ring
    ######################################################################
    def ring(self,name,T):
        self.defs.append('{\n')
        self.include("PhysBAM_Geometry/Geometry/RING.h")
        self.typedef(T=T)
        self.typedef(TV='VECTOR<T,3>')

        cppclass="RING<%s >"%(T)
        self.defs.append('class_<%s >("%s")\n'%(cppclass,name))
        self.defs.append("    .def(init<TV,TV,T,T>())\n")
        # members
        for m in ['plane1','plane2','height','outer_radius','inner_radius']:
            self.defs.append("    .def_readonly(\"%s\",&%s::%s)\n"%(m,cppclass,m))
        # member functions
        self.shape_helper(cppclass)
        # terminate
        self.defs.append(";}\n")
        self.read_write(cppclass,T)

    ######################################################################
    # Function cylinder
    ######################################################################
    def cylinder(self,name,T):
        self.defs.append('{\n')
        self.include("PhysBAM_Geometry/Geometry/CYLINDER.h")
        self.typedef(T=T)
        self.typedef(TV='VECTOR<T,3>')

        cppclass="CYLINDER<%s >"%(T)
        self.defs.append('class_<%s >("%s")\n'%(cppclass,name))
        self.defs.append("    .def(init<TV,TV,T>())\n")
        # members
        for m in ['plane1','plane2','height','radius']:
            self.defs.append("    .def_readonly(\"%s\",&%s::%s)\n"%(m,cppclass,m))
        # member functions
        self.shape_helper(cppclass)
        # terminate
        self.defs.append(";}\n")
        self.read_write(cppclass,T)

    ######################################################################
    # Function raw_array_instance
    ######################################################################
    def raw_array_instance(self,name,T):
        cppclass="RAW_ARRAY<%s >"%(T)
        self.include("PhysBAM_Tools/Arrays/RAW_ARRAY.h")
        self.include_class(T)
        self.local_include("INDEXING.h")
        self.local_include("CONVERSIONS.h")

        self.defs.append('class_<%s,boost::noncopyable>("%s",no_init)\n'%(cppclass,name)) # class
        # member functions
        self.str(cppclass)
        self.numpy(cppclass)
        self.equality_operators()
        self.defs.append("    .def(ARRAY_INDEXING_SUITE<%s >())\n"%cppclass)
        # trailer
        self.defs.append(";\n")

    ######################################################################
    # Function array
    ######################################################################
    def array_instance(self,name,T):
        self.raw_array_instance('R'+name,T)
        cppclass="ARRAY<%s >"%(T)
        self.include("PhysBAM_Tools/Arrays/ARRAY.h")
        self.include_class(T)
        self.local_include("INDEXING.h")
        self.local_include("CONVERSIONS.h")

        self.defs.append('class_<%s,bases<RAW_ARRAY<%s > > >("%s")\n'%(cppclass,T,name)) # class
        self.defs.append("    .def(init<const %s&>())\n"%cppclass) # copy constructor
        # member functions
        for member in ["Clean_Memory","Exchange_Arrays"]:
            self.member_function(cppclass,member)
        self.member_function(cppclass,("Resize","void",["const int","const bool","const bool"]),overloads=(1,3))
        self.str(cppclass)
        self.numpy(cppclass)
        # trailer
        self.defs.append(";\n")
        self.read_write(cppclass,T)
        self.defs.append("Register_Array_Conversion<%s >();\n"%cppclass)
        self.defs.append("implicitly_convertible<%s,RAW_ARRAY<const %s > >();\n"%(cppclass,T))

    ######################################################################
    # Function list_array
    ######################################################################
    def list_array_instance(self,name,T,sort=False):
        cppclass="LIST_ARRAY<%s >"%(T)
        self.include("PhysBAM_Tools/Arrays/SORT.h")
        self.include("PhysBAM_Tools/Arrays/LIST_ARRAY.h")
        self.include_class(T)
        self.local_include("INDEXING.h")
        self.local_include("CONVERSIONS.h")

        self.defs.append("class_<%s >(\"%s\")\n"%(cppclass,name)) # class
        self.defs.append("    .def(init<const %s&>())\n"%cppclass) # copy constructor
        # member functions
        for member in ["Preallocate","Append_Unique","Remove_End","Remove_Index_Lazy","Append",
                       "Remove_All","Clean_Memory","Insert","Pop","Exchange_Arrays",('Append_Elements','void',['const ARRAY<%s >&'%T])]:
            self.member_function(cppclass,member)
        self.member_function(cppclass,'Resize',overloads=(1,3))
        self.member_function(cppclass,'Exact_Resize',overloads=(1,2))
        self.str(cppclass)
        self.numpy(cppclass)
        self.equality_operators()
        self.defs.append("    .def(ARRAY_INDEXING_SUITE<%s >())\n"%cppclass)
        # trailer
        self.defs.append(";\n")
        if sort: self.defs.append("def(\"Sort\", static_cast<void (*)(%s&)>(Sort<%s >));\n"%(cppclass,cppclass))
        self.read_write(cppclass,T)
        self.defs.append("Register_Array_Conversion<%s >();\n"%cppclass)

    ######################################################################
    # Function arrays_nd_instance
    ######################################################################
    def arrays_nd_instance(self,name,T,T2,d,length=1):
        cppclass="ARRAYS_%dD<%s,%d>"%(d,T,length)
        self.include("PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_1D.h")
        self.include("PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_2D.h")
        self.include("PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_3D.h")
        routine="""\
template<class T_ARRAYS,class T_INDEX> const typename T_ARRAYS::ELEMENT& Get_Item(const T_ARRAYS& array,const T_INDEX& index)
{if(!array.Valid_Index(index)){
    PyErr_SetString(PyExc_IndexError,"Index out of range");
    boost::python::throw_error_already_set();}
return array(index);}

template<class T_ARRAYS,class T_INDEX> void Set_Item(T_ARRAYS& array,const T_INDEX& index,const typename T_ARRAYS::ELEMENT& value)
{if(!array.Valid_Index(index)){
    PyErr_SetString(PyExc_IndexError,"Index out of range");
    boost::python::throw_error_already_set();}
array(index)=value;}

"""
        if routine not in self.routines: self.routines.append(routine)
        self.defs.append('{\n')
        self.typedef(T=T)
        self.defs.append('static const int d=%d;\n'%d)
        self.defs.append('class_<%s >("%s")\n'%(cppclass,name))
        self.defs.append('    .def(init<const %s&>())\n'%cppclass)
        self.defs.append('    .def(init<%s>())\n'%','.join(['int']*(2*d)))
        self.defs.append('    .def("__getitem__",&Get_Item<%s,VECTOR<int,d> >,return_value_policy<copy_const_reference>())\n'%cppclass)
        self.defs.append('    .def("__setitem__",&Set_Item<%s,VECTOR<int,d> >)\n'%cppclass)
        self.member_function(cppclass,('Fill','void',['const T&']))
        if length>1:
            self.member_function(cppclass,('Fill','void',['const VECTOR<T,%d>&'%length]))
        self.member_function(cppclass,'Domain_Indices')
        if T in ['float','double']:
            self.static_function(cppclass,'Dot_Product')
        if T=='bool':
            self.member_function(cppclass,'Number_True')
        if T in ['int','float','double']:
            for f in ['Max','Maxabs','Maxmag','Min','Minmag','Sum','Sumabs']:
                self.member_function(cppclass,f)
        self.member_function(cppclass,('Resize','void',['int']*(2*d)+['bool','bool','const T&']),overloads=(2*d,2*d))
        self.member_function(cppclass,('Resize','void',['const RANGE<VECTOR<int,%d> >&'%d,'bool','bool','const T&']),overloads=(1,1))
        for f in ['m','n','mn'][:d]:
            self.member_variable(cppclass,(f,'readonly'))
        self.numpy(cppclass)
        self.defs.append(";\n")
        self.read_write(cppclass,T2)
        self.defs.append('}\n')

    ######################################################################
    # Function hashtable_instance
    ######################################################################
    def hashtable_instance(self,name,TKEY,T):
        cppclass="HASHTABLE<%s,%s >"%(TKEY,T)
        self.include("PhysBAM_Tools/Data_Structures/HASHTABLE.h")
        self.include("PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_HASHTABLE.h")
        self.local_include("INDEXING.h")

        self.defs.append("class_<%s >(\"%s\")\n"%(cppclass,name))
        functions=["Max_Size","Clean_Memory","Size","Next_Resize"]
                   #("Insert","%s&"%T,["const %s&"%TKEY,"const %s&"%T]),
                   #("Contains","bool",["const %s&"%TKEY],"return_internal_reference<>()")]
        for member in functions:
            self.member_function(cppclass,member)
        self.defs.append("    .def(HASHTABLE_INDEXING_SUITE<%s >())\n"%cppclass)
        self.defs.append(";\n")
        # wrap HASHTABLE_ENTRY_TEMPLATE element type
        entry_cppclass="HASHTABLE_ENTRY_TEMPLATE<%s,%s >"%(TKEY,T)
        self.defs.append("class_<%s >(\"%s\")\n"%(entry_cppclass,name+"_ENTRY"))
        for j in ["state","key","data"]:
            self.defs.append("    .def_readwrite(\"%s\",&%s::%s)\n"%(j,entry_cppclass,j))
        self.defs.append(";\n")
        # wrap ENTRY_STATE
        enum_cppclass="HASHTABLE_ENTRY_STATE"
        self.defs.append("enum_<%s >(\"%s\")\n"%(enum_cppclass,enum_cppclass))
        for j in ["ENTRY_FREE","ENTRY_ACTIVE","ENTRY_DELETED"]:
            self.defs.append("    .value(\"%s\",%s)\n"%(j,j))
        self.defs.append("    .export_values()")
        self.defs.append(";\n")
        self.read_write(cppclass,T)

    ######################################################################
    # Function pair_instance
    ######################################################################
    def pair_instance(self,name,T1,T2):
        cppclass="PAIR<%s,%s >"%(T1,T2)
        self.include("PhysBAM_Tools/Data_Structures/PAIR.h")

        self.defs.append("class_<%s >(\"%s\")\n"%(cppclass,name))
        for guy in ["x","y"]:
            self.defs.append("    .def_readwrite(\"%s\",&%s::%s)\n"%(guy,cppclass,guy))
        self.defs.append(";\n")

    ######################################################################
    # Function structure
    ######################################################################
    def structure(self,name,TV):
        self.include("PhysBAM_Geometry/Geometry/STRUCTURE.h")
        cppclass="STRUCTURE<%s >"%TV
        self.defs.append("class_<%s,boost::noncopyable>(\"%s\",no_init)\n"%(cppclass,name))
        for member in [("Append_Particles_And_Create_Copy","STRUCTURE<%s >* "%TV,["DEFORMABLE_GEOMETRY_PARTICLES<%s >&"%TV,"LIST_ARRAY<int>*"],"const","return_value_policy<manage_new_object>()")]:
            self.member_function(cppclass,member)
        self.defs.append(";\n")

    ######################################################################
    # Function solids_forces
    ######################################################################
    def solids_forces(self,name,TV):
        self.include("Forces_And_Torques/SOLIDS_FORCES.h")
        self.include("PhysBAM_Tools/Vectors/VECTOR_3D.h")
        self.include("Parallel_Computation/MPI_SOLIDS.h")
        cppclass="SOLIDS_FORCES<%s >"%TV
        self.defs.append("class_<%s,boost::noncopyable>(\"%s\",no_init)\n"%(cppclass,name))
        self.defs.append(";\n")

    ######################################################################
    # Function segmented_curve_instance
    ######################################################################
    def segmented_curve_instance(self,name,TV):
        (T,d)=parse_template(TV)[1]
        d=int(d)
        self.include("PhysBAM_Geometry/Geometry/SEGMENTED_CURVE_2D.h")
        if d==2: cppclass="SEGMENTED_CURVE_2D<%s >"%T
        else: cppclass="SEGMENTED_CURVE<%s >"%TV
        base_cppclass="STRUCTURE<%s >"%TV
        self.defs.append("class_<%s,bases<%s >,boost::noncopyable>(\"%s\",no_init)\n"%(cppclass,base_cppclass,name))
        # references to internal data
        self.reference(cppclass,"mesh")
        self.reference(cppclass,"particles")
        # member functions
        functions=["Update_Number_Nodes"]
        for member in functions:
            self.member_function(cppclass,member)
        self.defs.append(";\n")
        self.read_write(cppclass,T)

    ######################################################################
    # Function tetrahedron
    ######################################################################
    def tetrahedron(self,name,T):
        self.include("PhysBAM_Geometry/Geometry/TETRAHEDRON.h")
        TV="VECTOR<%s,3>"%T
        cppclass="TETRAHEDRON<%s >"%T
        self.defs.append("class_<%s,boost::noncopyable>(\"%s\",init<const %s&,const %s&,const %s&,const %s&>())\n"%(cppclass,name,TV,TV,TV,TV))
        for member in ["Create_Triangles",
                       ("Volume",T,[],"const"),
                       ("Signed_Volume",T,[],"const"),
                       ("Barycentric_Coordinates",TV,["const %s&"%TV],"const"),
                       ("Point_From_Barycentric_Coordinates",TV,["const %s&"%TV],"const"),
                       "Barycentric_Inside",
                       ("Center",TV,[],"const"),
                       ("Minimum_Edge_Length",T,[],"const"),
                       ("Size",T,[],"const"),
                       ("Signed_Size",T,[],"const"),
                       "Inside","Normal","Outside","Boundary","Thickened","Bounding_Box",
                       "Surface","Closest_Point","Minimum_Angle","Maximum_Angle","Minimum_Altitude",
                       ("Aspect_Ratio",T,[],"const"),
                       "Minimum_Dihedral_Angle",
                       "Maximum_Dihedral_Angle"]: # todo all collision and interaction and ray functions
            self.member_function(cppclass,member)
        for i in range(1,5):
            self.reference(cppclass,"triangle%d"%i)

        self.defs.append(";\n")

    ######################################################################
    # Function segment_3d
    ######################################################################
    def segment_3d(self,name,T):
        self.include("PhysBAM_Geometry/Geometry/SEGMENT_3D.h")
        TV="VECTOR<%s,3>"%T
        TV2="VECTOR<%s,2>"%T
        cppclass="SEGMENT_3D<%s >"%T
        self.defs.append("class_<%s,boost::noncopyable>(\"%s\",init<const %s&,const %s&>())\n"%(cppclass,name,TV,TV))
        for member in [("Length",T,[],"const"),
                       ("Size",T,[],"const"),
                       ("Barycentric_Coordinates",TV2,["const %s&"%TV],"const"),
                       ("Interpolation_Fraction",T,["const %s&"%TV],"const"),
                       ("Closest_Point_On_Segment"),
                       ("Closest_Point_On_Line"),
                       ("Distance_From_Point_To_Line"),
                       ("Shortest_Vector_Between_Lines"),
                       ("Shortest_Vector_Between_Segments")]:
            self.member_function(cppclass,member)
        for member in ["x1","x2"]:
            self.member_variable(cppclass,member)
        self.defs.append(";\n")

    ######################################################################
    # Function triangle_3d
    ######################################################################
    def triangle_3d(self,name,T):
        self.include("PhysBAM_Geometry/Geometry/TRIANGLE_3D.h")
        TV="VECTOR<%s,3>"%T
        cppclass="TRIANGLE_3D<%s >"%T
        self.defs.append("class_<%s,boost::noncopyable>(\"%s\",init<const %s&,const %s&,const %s&>())\n"%(cppclass,name,TV,TV,TV))
        for member in [("Area",T,[],"const"),
                       ("Size",T,[],"const"),
                       ("Aspect_Ratio",T,[],"const"),
                       ("Barycentric_Coordinates",TV,["const %s&"%TV],"const"),
                       ("Point_From_Barycentric_Coordinates",TV,["const %s&"%TV],"const"),
                       ("Center",TV,[],"const"),
                       "Signed_Distance",
                       "Incenter","Change_Size",
                       "Point_Inside_Triangle","Planar_Point_Inside_Triangle","Lazy_Planar_Point_Inside_Triangle",
                       ("Minimum_Edge_Length",T,[],"const"),
                       ("Maximum_Edge_Length",T,[],"const"),
                       ("Normal",TV,[],"const"),
                       "Bounding_Box",
                       "Region","Closest_Point","Distance_To_Triangle","Minimum_Angle","Maximum_Angle","Signed_Solid_Angle"]: # todo all collision and interaction and ray functions
            self.member_function(cppclass,member)
        for member in ["x1","x2","x3"]:
            self.member_variable(cppclass,member)
        self.defs.append(";\n")

    ######################################################################
    # Function triangulated_surface_instance
    ######################################################################
    def triangulated_surface_instance(self,name,T):
        self.include("PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h")
        self.include("PhysBAM_Tools/Data_Structures/HASHTABLE.h")
        self.include("PhysBAM_Geometry/Geometry/TRIANGULATED_SURFACE.h")
        self.include("PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.h")
        self.include("PhysBAM_Geometry/Particles/PARTICLE_SUBSET.h")
        self.include("PhysBAM_Geometry/Particles/DEFORMABLE_GEOMETRY_PARTICLES.h")
        cppclass="TRIANGULATED_SURFACE<%s >"%T
        base_cppclass="STRUCTURE<VECTOR<%s,3> >"%T
        self.defs.append("class_<%s,bases<%s >,boost::noncopyable>(\"%s\",no_init)\n"%(cppclass,base_cppclass,name))
        # references to internal data
        self.reference(cppclass,"mesh")
        self.reference(cppclass,"particles")
        # member functions
        for member in ["Update_Number_Nodes","Get_Element","Normal","Use_Vertex_Normals","Update_Vertex_Normals","Refresh_Auxiliary_Structures",
            "Initialize_Torus_Mesh_And_Particles","Initialize_Cylinder_Mesh_And_Particles","Update_Bounding_Box",("Update_Triangle_List",'void',[]),
            'Linearly_Subdivide','Loop_Subdivide','Root_Three_Subdivide','Total_Area','Make_Orientations_Consistent_With_Implicit_Surface',
            ("Discard_Valence_Zero_Particles_And_Renumber",'void',['ARRAY<int>&']),('Rescale','void',['const '+T]),('Rescale','void',['const '+T]*3)]:
            self.member_function(cppclass,member)
        for member in ["Inside","Outside","Calculate_Signed_Distance"]:
            self.member_function(cppclass,member,overloads=(1,2))
        self.member_function(cppclass,'Initialize_Hierarchy',overloads=(0,2))
        self.member_function(cppclass,'Initialize_Particle_Hierarchy',overloads=(1,3))
        self.member_function(cppclass,'Close_Surface',overloads=(3,4))
        self.member_function(cppclass,'Remove_Degenerate_Triangles',overloads=(0,1))
        # surface
        routine="""
template<class T> boost::python::object Surface_Wrapper(const TRIANGULATED_SURFACE<T>& surface,const VECTOR<T,3>& location,const T max_depth=0,const T thickness_over_2=0)
{int closest_triangle;T distance;VECTOR<T,3> result=surface.Surface(location,max_depth,thickness_over_2,&closest_triangle,&distance);
return boost::python::make_tuple(result,closest_triangle,distance);}

BOOST_PYTHON_FUNCTION_OVERLOADS(Surface_Wrapper_Overloads,Surface_Wrapper,2,4);

"""
        if routine not in self.routines: self.routines.append(routine)
        self.defs.append('    .def("Surface",&Surface_Wrapper<%s >,Surface_Wrapper_Overloads())\n'%T)
        # finish
        self.defs.append(";\n")
        self.read_write(cppclass,T)

    ######################################################################
    # Function triangulated_area_intsance
    ######################################################################
    def triangulated_area_instance(self,name,T):
        self.include("PhysBAM_Geometry/Geometry/TRIANGULATED_AREA.h")
        cppclass="TRIANGULATED_AREA<%s >"%T
        base_cppclass="STRUCTURE<VECTOR<%s,2> >"%T
        self.defs.append("class_<%s,bases<%s >,boost::noncopyable>(\"%s\",no_init)\n"%(cppclass,base_cppclass,name))
        # references to internal data
        self.reference(cppclass,"mesh")
        self.reference(cppclass,"particles")
        # member functions
        functions=["Update_Number_Nodes"]
        for member in functions:
            self.member_function(cppclass,member)
        self.defs.append(";\n")
        self.read_write(cppclass,T)

    ######################################################################
    # Function tetrahedralized_volume_instance
    ######################################################################
    def tetrahedralized_volume_instance(self,name,T):
        self.include("PhysBAM_Geometry/Geometry/TETRAHEDRALIZED_VOLUME.h")
        self.include("PhysBAM_Geometry/Spatial_Acceleration/TETRAHEDRON_HIERARCHY.h")
        cppclass="TETRAHEDRALIZED_VOLUME<%s >"%T
        base_cppclass="STRUCTURE<VECTOR<%s,3> >"%T
        self.defs.append("class_<%s,bases<%s >,boost::noncopyable>(\"%s\",no_init)\n"%(cppclass,base_cppclass,name))
        # references to internal data
        self.reference(cppclass,"mesh")
        self.reference(cppclass,"particles")
        self.reference(cppclass,"hierarchy")
        # member functions
        functions=["Update_Number_Nodes",("Discard_Valence_Zero_Particles_And_Renumber",'void',['ARRAY<int>&']),
        ('Rescale','void',['const '+T]),('Rescale','void',['const '+T]*3),("Initialize_Hierarchy","void",["const bool"])]
        for member in functions:
            self.member_function(cppclass,member)
        self.defs.append(";\n")
        self.read_write(cppclass,T)

    ######################################################################
    # Function segment_mesh
    ######################################################################
    def segment_mesh(self,name):
        self.include("PhysBAM_Geometry/Topology/SEGMENT_MESH.h")
        self.include("PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_MESH_OBJECT.h")
        cppclass="SEGMENT_MESH"
        self.defs.append("class_<%s >(\"%s\")\n"%(cppclass,name))
        self.member_variable(cppclass,'elements')
        for member in ["Initialize_Neighbor_Nodes","Initialize_Incident_Elements",'Initialize_Adjacent_Elements',"Initialize_Ordered_Loop_Nodes"]:
            self.member_function(cppclass,member)
        self.reference(cppclass,"ordered_loop_nodes")
        self.reference(cppclass,"neighbor_nodes")
        self.reference(cppclass,"incident_elements")
        self.reference(cppclass,"adjacent_elements")
        self.defs.append(";\n")

    ######################################################################
    # Function triangle_mesh
    ######################################################################
    def triangle_mesh(self,name):
        self.include("PhysBAM_Geometry/Topology/TRIANGLE_MESH.h")
        self.include("PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_MESH_OBJECT.h")
        cppclass="TRIANGLE_MESH"
        self.defs.append("class_<%s >(\"%s\")\n"%(cppclass,name))
        self.member_variable(cppclass,'elements')
        for member in ["Initialize_Neighbor_Nodes","Initialize_Incident_Elements","Initialize_Adjacent_Elements","Initialize_Herring_Bone_Mesh",
                       "Initialize_Square_Mesh","Initialize_Segment_Mesh",'Orientations_Consistent','Make_Orientations_Consistent',
                       ("Delete_Sorted_Elements","void",["const LIST_ARRAY<int>&","HASHTABLE<int,int>&"])]:
            self.member_function(cppclass,member)
        self.reference(cppclass,"segment_mesh")
        self.reference(cppclass,"neighbor_nodes")
        self.reference(cppclass,"incident_elements")
        self.reference(cppclass,"adjacent_elements")
        self.routines.append('LIST_ARRAY<int> Non_Manifold_Nodes(TRIANGLE_MESH& mesh){LIST_ARRAY<int> nodes;mesh.Non_Manifold_Nodes(nodes);return nodes;}\n')
        self.defs.append('    .def("Non_Manifold_Nodes",Non_Manifold_Nodes)\n')
        self.defs.append(";\n")

    ######################################################################
    # Function tetrahedron_mesh
    ######################################################################
    def tetrahedron_mesh(self,name):
        self.include("PhysBAM_Geometry/Topology/TETRAHEDRON_MESH.h")
        self.include("PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_MESH_OBJECT.h")
        cppclass="TETRAHEDRON_MESH"
        self.defs.append("class_<%s >(\"%s\")\n"%(cppclass,name))
        self.member_variable(cppclass,'elements')
        for member in ['Initialize_Neighbor_Nodes','Initialize_Incident_Elements','Initialize_Adjacent_Elements',"Initialize_Bending_Tetrahedrons"]:
            self.member_function(cppclass,member)
        self.reference(cppclass,"neighbor_nodes")
        self.reference(cppclass,"incident_elements")
        self.reference(cppclass,"adjacent_elements")
        self.defs.append(";\n")

    ######################################################################
    # Function single_particles
    ######################################################################
    def single_particles(self):
        self.include("PhysBAM_Tools/Particles/PARTICLES.h")
        name='SINGLE_PARTICLES'
        cppclass='SINGLE_PARTICLES<PARTICLES_BASE>'
        self.defs.append('class_<%s,boost::noncopyable>("%s",no_init)\n'%(cppclass,name))
        self.str(cppclass)
        self.defs.append(";\n")
        self.defs.append('register_ptr_to_python<shared_ptr<%s > >();\n'%cppclass)

    ######################################################################
    # Function particle_base
    ######################################################################
    def particle_base(self):
        self.include("PhysBAM_Tools/Particles/PARTICLES.h")
        self.include("boost/shared_ptr.hpp>\nusing boost::shared_ptr;//This file is convoluted.")
        name=cppclass='PARTICLES_BASE'
        self.defs.append('class_<%s,bases<CLONEABLE_BASE>,boost::noncopyable>("%s",no_init)\n'%(cppclass,name))
        self.routines.append("""\
shared_ptr<SINGLE_PARTICLES<PARTICLES_BASE> > Particle_Get_Item(PARTICLES_BASE& particles,const int index)
{if(!particles.Valid_Index(index)) throw INDEX_ERROR("Particle index out of range");
return shared_ptr<SINGLE_PARTICLES<PARTICLES_BASE> >(new SINGLE_PARTICLES<PARTICLES_BASE>(particles,index));}

void Particle_Set_Item(PARTICLES_BASE& particles,const int index,const SINGLE_PARTICLES<PARTICLES_BASE>& source)
{if(!particles.Valid_Index(index)) throw INDEX_ERROR("Destination index out of range");
if(!source.particles.Valid_Index(source.index)) throw INDEX_ERROR("Source index out of range");
particles(index)=source;}

void Particle_Del_Item(PARTICLES_BASE& particles,const int index)
{if(!particles.Valid_Index(index)) throw INDEX_ERROR("Particle deletion index out of range");
particles.Delete_Particle(index);}

PARTICLE_ATTRIBUTE_BASE& Particle_Get_Attr(PARTICLES_BASE& particles,const std::string& name)
{PARTICLE_ATTRIBUTE_BASE* attribute=particles.Get_Attribute<PARTICLE_ATTRIBUTE_BASE>(name);
if(!attribute) throw KEY_ERROR(str(boost::format("Invalid particle attribute name %s")%name));
return *attribute;} // TODO: this is not safe for dynamic attributes, since they can disappear without warning

void Particle_Set_Attr(PARTICLES_BASE& particles,const std::string& name,const boost::python::object& value)
{boost::python::object(Particle_Get_Attr(particles,name))[boost::python::_]=value;}

void Add_Particles(PARTICLES_BASE& particles,const int count)
{particles.Add_Particles(count);}
""")
        self.defs.append('    .def("__len__",&PARTICLES_BASE::Size)\n')
        self.defs.append('    .def("__getitem__",&Particle_Get_Item,with_custodian_and_ward_postcall<0,1>())\n')
        self.defs.append('    .def("__setitem__",&Particle_Set_Item)\n')
        self.defs.append('    .def("__delitem__",&Particle_Del_Item)\n')
        self.defs.append('    .def("__getattr__",&Particle_Get_Attr,return_internal_reference<>())\n')
        self.defs.append('    .def("__setattr__",&Particle_Set_Attr)\n')
        self.equality_operators()
        for member in ["Preallocate","Clean_Memory","Add_Particle","Add_Attributes","Number_Of_Attributes"]:
            self.member_function(cppclass,member)
        self.defs.append('    .def("Add_Particles",Add_Particles)\n')
        self.member_function(cppclass,('Add_Attribute','void',['const PARTICLE_ATTRIBUTE_BASE&']))
        self.member_function(cppclass,('Add_Attribute','void',['const SYMBOL','const PARTICLE_ATTRIBUTE_BASE&']))
        self.defs.append(";\n")
        self.read_write(cppclass,T)

    ######################################################################
    # Function particles
    ######################################################################
    def particles(self,name,TV):
        self.include("PhysBAM_Tools/Particles/PARTICLES.h")
        cppclass="PARTICLES<%s >"%(TV)
        self.defs.append("class_<%s,boost::noncopyable,bases<PARTICLES_BASE> >(\"%s\")\n"%(cppclass,name))
        self.defs.append(";\n")

    ######################################################################
    # Function simple_particles
    ######################################################################
    def simple_particles(self,name,cppclass):
        self.include_class(cppclass)
        self.defs.append('{\n')
        self.typedef(TV='%s::VECTOR_T'%cppclass)
        self.defs.append('class_<%s,boost::noncopyable,bases<PARTICLES<TV> > >("%s")\n'%(cppclass,name))
        if cppclass.startswith('DEFORMABLE_BODY_PARTICLES'):
            for member in ['Store_Velocity','Store_Mass','Delete_Particle']:
                self.member_function(cppclass,member)
        self.defs.append(";}\n")

    ######################################################################
    # Function box_hierarchy
    ######################################################################
    def box_hierarchy(self,name,TV):
        T=parse_template(TV)[1][0]
        self.include("PhysBAM_Geometry/Geometry/BOX.h")
        self.include("PhysBAM_Geometry/Spatial_Acceleration/BOX_HIERARCHY.h")
        cppclass="BOX_HIERARCHY<%s >"%TV

        routine="""
template<class TV> LIST_ARRAY<int> Intersection_List_Wrapper(const BOX_HIERARCHY<TV>& hierarchy,const TV& point,const typename TV::SCALAR thickness_over_two=0)
{LIST_ARRAY<int> list;hierarchy.Intersection_List(point,list,thickness_over_two);return list;}

BOOST_PYTHON_FUNCTION_OVERLOADS(Intersection_List_Overloads,Intersection_List_Wrapper,2,3);

"""
        if routine not in self.routines: self.routines.append(routine)

        self.defs.append("class_<%s,boost::noncopyable>(\"%s\")\n"%(cppclass,name))
        for m in ['leaves','root','parents','children','box_hierarchy','box_radius']:
            self.member_variable(cppclass,(m,'readonly'))
        self.defs.append("    .def(\"Intersection_List\",&Intersection_List_Wrapper<%s >,Intersection_List_Overloads())\n"%TV)
        self.member_function(cppclass,("Intersection_List","void",["const %s&"%TV,"LIST_ARRAY<int>&","const %s"%T],"const"))
        self.defs.append(";\n")

    ######################################################################
    # Function particle_hierarchy
    ######################################################################
    def particle_hierarchy(self,name,TV):
        T=parse_template(TV)[1][0]
        self.include("PhysBAM_Geometry/Spatial_Acceleration/PARTICLE_HIERARCHY.h")
        self.include("Particles/DEFORMABLE_BODY_PARTICLES.h")
        cppclass="PARTICLE_HIERARCHY<%s >"%TV
        particles='DEFORMABLE_BODY_PARTICLES<%s >'%TV
        base="BOX_HIERARCHY<%s >"%TV
        self.defs.append("class_<%s,bases<%s >,boost::noncopyable>(\"%s\",no_init)\n"%(cppclass,base,name))
        self.defs.append('    .def(init<const %s&>()[with_custodian_and_ward<1,2>()])\n'%particles)
        self.defs.append('    .def(init<const %s&,bool>()[with_custodian_and_ward<1,2>()])\n'%particles)
        self.defs.append('    .def(init<const %s&,bool,int>()[with_custodian_and_ward<1,2>()])\n'%particles)
        for m in ['particles_in_group','particles_per_group']:
            self.member_variable(cppclass,(m,'readonly'))
        self.defs.append(";\n")

    ######################################################################
    # Function tetrahedron_hierarchy
    ######################################################################
    def tetrahedron_hierarchy(self,name,T):
        TV="VECTOR<%s,3>"%T
        self.include("PhysBAM_Geometry/Spatial_Acceleration/TETRAHEDRON_HIERARCHY.h")
        cppclass="TETRAHEDRON_HIERARCHY<%s >"%T
        base_cppclass="BOX_HIERARCHY<%s >"%TV
        self.defs.append("class_<%s,bases<%s >,boost::noncopyable>(\"%s\",no_init)\n"%(cppclass,base_cppclass,name))
        self.defs.append(";\n")

    ######################################################################
    # Function particle_attribute_base
    ######################################################################
    def particle_attribute_base(self):
        self.defs.append('{\n')
        name=cppclass='PARTICLE_ATTRIBUTE_BASE'
        self.defs.append('class_<%s,boost::noncopyable>("%s",no_init)\n'%(cppclass,name))
        self.defs.append(";}\n")
        self.defs.append('register_ptr_to_python<shared_ptr<PARTICLE_ATTRIBUTE_BASE> >();\n')

    ######################################################################
    # Function particle_attribute
    ######################################################################
    def particle_attribute(self,name,T):
        self.local_include("INDEXING.h")
        self.defs.append('{\n')
        self.typedef(T=T)
        cppclass='PARTICLE_ATTRIBUTE<T>'
        self.defs.append("class_<%s,boost::noncopyable,bases<RAW_ARRAY<T>,PARTICLE_ATTRIBUTE_BASE> >(\"%s\")\n"%(cppclass,name))
        self.defs.append(";}\n")

    ######################################################################
    # Function particle_simple_attribute
    ######################################################################
    def particle_simple_attribute(self,name,cppclass):
        self.local_include("INDEXING.h")
        self.defs.append('{\n')
        self.typedef(T='%s::ELEMENT'%cppclass)
        self.defs.append("class_<%s,boost::noncopyable,bases<PARTICLE_ATTRIBUTE<T> > >(\"%s\")\n"%(cppclass,name))
        self.defs.append(";}\n")

    ######################################################################
    # Function particle_mass_attribute
    ######################################################################
    def particle_mass_attribute(self,name,T):
        self.include("Particles/PARTICLE_MASS_ATTRIBUTE.h")
        self.include("Deformable_Objects/SOFT_BINDINGS.h")
        self.defs.append('{\n')
        self.typedef(T=T)
        cppclass='PARTICLE_MASS_ATTRIBUTE<T>'
        self.defs.append("class_<%s,boost::noncopyable,bases<PARTICLE_ATTRIBUTE<T> > >(\"%s\",no_init)\n"%(cppclass,name))
        for d in [2,3]:
            self.member_function(cppclass,("Compute_Auxiliary_Attributes","void",["const SOFT_BINDINGS<VECTOR<%s,%d> >&"%(T,d)]))
        self.defs.append(";}\n")

    ######################################################################
    # Function levelset_instance
    ######################################################################
    def levelset_instance(self,name,TV):
        self.include("PhysBAM_Geometry/Level_Sets/LEVELSET_1D.h")
        self.include("PhysBAM_Geometry/Level_Sets/LEVELSET_2D.h")
        self.include("PhysBAM_Geometry/Level_Sets/LEVELSET_3D.h")
        self.include("PhysBAM_Geometry/Geometry/TRIANGULATED_SURFACE.h")
        self.include("PhysBAM_Tools/Matrices/MATRIX_1X1.h")
        self.include("PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h")
        self.include("PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h")
        (T,d)=parse_template(TV)[1]
        T_GRID="GRID_%sD<%s >"%(d,T)
        T_ARRAYS="GRID_ARRAYS_POLICY<%s >::ARRAYS_SCALAR"%T_GRID
        cppclass="LEVELSET_POLICY<%s >::LEVELSET"%T_GRID
        self.defs.append('{\n')
        self.typedef(T=T,TV=TV,T_GRID=T_GRID,T_ARRAYS=T_ARRAYS)
        self.defs.append('class_<%s,boost::noncopyable>("%s",init<T_GRID&,T_ARRAYS&>()[with_custodian_and_ward<1,2,with_custodian_and_ward<1,3> >()])\n'%(cppclass,name))
        # members
        for f in ['Phi','Phi_Secondary','Extended_Phi','Curvature','Normal','Extended_Normal','Hessian','Principal_Curvatures',
                  'Compute_Cell_Minimum_And_Maximum','Curvature']:
            self.member_function(cppclass,f)
        for f,o in [('Compute_Normals',(0,1)),
                    ('Fast_Marching_Method_Outside_Band',(1,3)),
                    ('Lazy_Inside',(1,2)),
                    ('Lazy_Inside_Extended_Levelset',(1,2)),
                    ('Lazy_Outside',(1,2)),
                    ('Lazy_Outside_Extended_Levelset',(1,2))]:
            self.member_function(cppclass,f,overloads=o)
        
        if d>1:
            self.member_function(cppclass,('Compute_Curvature','void',[T]),overloads=(0,1))
            self.member_function(cppclass,('Fast_Marching_Method','void',['const %s'%T,'const %s'%T,'const LIST_ARRAY<VECTOR<int,%d> >*'%d,'const bool']),overloads=(0,2))
        if d==3:
            self.member_function(cppclass,'Approximate_Surface_Area')
            self.member_function(cppclass,('Calculate_Triangulated_Surface_From_Marching_Tetrahedra','void',['const T_GRID&','TRIANGULATED_SURFACE<T>&'],'const'))
            self.member_function(cppclass,('Calculate_Triangulated_Surface_From_Marching_Tetrahedra','void',['TRIANGULATED_SURFACE<T>&','bool'],'const'),overloads=(1,2))

        # references to internal data
        self.reference(cppclass,"grid")
        self.reference(cppclass,"phi")
        self.defs.append(";}\n")
        self.read_write(cppclass,T)

    ######################################################################
    # Function levelset_implicit_object
    ######################################################################
    def levelset_implicit_object(self,name,TV):
        self.include("PhysBAM_Geometry/Implicit_Objects/LEVELSET_IMPLICIT_OBJECT.h")
        (T,d)=parse_template(TV)[1]
        T_GRID="GRID_%sD<%s >"%(d,T)
        T_ARRAYS="GRID_ARRAYS_POLICY<%s >::ARRAYS_SCALAR"%T_GRID
        cppclass="LEVELSET_IMPLICIT_OBJECT<%s >"%TV
        self.defs.append('class_<%s,boost::noncopyable>("%s",init<%s&,%s&>()[with_custodian_and_ward<1,2,with_custodian_and_ward<1,3> >()])\n'%(cppclass,name,T_GRID,T_ARRAYS))
        # references to internal data
        self.reference(cppclass,"levelset")
        # member functions
        for member in ['Extended_Phi','Normal','Extended_Normal']:
            self.member_function(cppclass,member)
        for member in [("Lazy_Inside","bool",["const %s&"%TV,"const %s"%T],"const"),("Compute_Cell_Minimum_And_Maximum","void",["bool"])]:
            self.member_function(cppclass,member)
        self.defs.append('    .def("__call__",&%s::operator())\n'%cppclass)
        self.defs.append(";\n")
        self.read_write(cppclass,T)

    ######################################################################
    # Function heat_uniform
    ######################################################################
    def heat_uniform(self,name,TV):
        self.include("Heat_Flows/HEAT_UNIFORM.h")
        (T,d)=parse_template(TV)[1]
        T_GRID="GRID_%sD<%s >"%(d,T)
        self.defs.append('{\n')
        self.typedef(T=T,TV=TV,T_GRID=T_GRID)
        T_ARRAYS="GRID_ARRAYS_POLICY<%s >::ARRAYS_SCALAR"%T_GRID
        cppclass="HEAT_UNIFORM<%s >"%T_GRID
        self.defs.append('class_<%s,boost::noncopyable>("%s",no_init)\n'%(cppclass,name))
        # members
        self.static_function(cppclass,('Smooth','void',['ARRAYS_%dD<T>&'%d,'int','const ARRAYS_%dD<T>*'%d]),overloads=(2,2))
        self.defs.append(";}\n")

    ######################################################################
    # Function dualcontour_instance
    ######################################################################
    def dualcontour_instance(self,name,T_GRID):
        self.include("PhysBAM_Geometry/Geometry/SEGMENTED_CURVE_2D.h")
        self.include("PhysBAM_Geometry/Geometry/TRIANGULATED_AREA.h")
        self.include("PhysBAM_Geometry/Geometry/TRIANGULATED_SURFACE.h")
        self.include("PhysBAM_Geometry/Dual_Contour/DUALCONTOUR_2D.h")
        self.include("PhysBAM_Geometry/Dual_Contour/DUALCONTOUR_3D.h")
        d=int(T_GRID[5])
        T=parse_template(T_GRID)[1][0]
        cppclass="DUALCONTOUR_%dD<%s>"%(d,T)
        self.defs.append('class_<%s,boost::noncopyable>("%s",no_init)\n'%(cppclass,name))
        if d==2: statics=['Create_Segmented_Curve_From_Levelset','Create_Triangulated_Area_From_Levelset']
        else: statics=['Create_Triangulated_Surface_From_Levelset']
        for f in statics:
            if d==2:
                self.routines.append("BOOST_PYTHON_FUNCTION_OVERLOADS(%s_overloads,%s::%s,1,2)\n"%(f,cppclass,f))
                self.defs.append('    .def("%s",%s::%s,%s_overloads()[return_value_policy<manage_new_object>()])\n'%(f,cppclass,f,f))
            else: self.defs.append('    .def("%s",%s::%s,return_value_policy<manage_new_object>())\n'%(f,cppclass,f))
            self.defs.append('    .staticmethod("%s")\n'%f)
        self.defs.append(";\n")

    ######################################################################
    # Function solids_fluid_example_uniform
    ######################################################################
    def solids_fluids_example_uniform(self,name,T,T_GRID):
        self.include("PhysBAM_Tools/Grids_Uniform/GRID_3D.h")
        self.local_include("SOLIDS_FLUIDS_EXAMPLE_UNIFORM_PYTHON.h")
        self.include("Solids_And_Fluids/SOLIDS_PARAMETERS.h")
        self.include("Rigid_Bodies/RIGID_BODY_3D.h")
        self.include("PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h")
        r=re.compile("GRID_([0-9])D")
        m=r.match(T_GRID)
        d=int(m.group(1))

        fp_cppclass="FLUIDS_PARAMETERS<%s >"%T_GRID

        self.defs.append("enum_<%s::TYPE >(\"FLUIDS_PARAMETERS_f%d\")\n"%(fp_cppclass,d))
        for j in ["NONE","SMOKE","FIRE","WATER","SPH"]:
            self.defs.append("    .value(\"%s\",%s::%s)\n"%(j,fp_cppclass,j))
        self.defs.append(";\n")

        #base_cppclass="SOLIDS_FLUIDS_EXAMPLE_UNIFORM<%s >"%(T_GRID)
        cppclass="SOLIDS_FLUIDS_EXAMPLE_UNIFORM_PYTHON<%s >"%(T_GRID)
        self.defs.append("class_<%s,boost::noncopyable>(\"%s\",init<const STREAM_TYPE,const int,const FLUIDS_PARAMETERS<%s >::TYPE>())\n"%(cppclass,name,T_GRID))
        #self.defs.append('    .def("Initialize_Bodies",&%s::Initialize_Bodies,&%s::default_Initialize_Bodies)'%(base_cppclass,cppclass))
        for i in ["frame_rate","last_frame","initial_time","first_frame","last_frame","frame_rate","restart","restart_frame","write_output_files",
                  "write_first_frame","write_last_frame","write_time","write_frame_title","frame_title","output_directory","data_directory","abort_when_dt_below"]:
            self.member_variable(cppclass,i)
        for member in ["Unite_All_Fragments","Collide_All_Bodies","Collide_Segmented_Curve","Collide_Rigid_Bodies","Set_Write_Substeps_Level","Zero_Nodes",
                       "Set_Rigid_Velocities","Set_Rigid_Positions"]:
            self.member_function(cppclass,member)
        self.reference(cppclass,"solids_parameters")
        self.defs.append(";\n")

    ######################################################################
    # Function binding
    ######################################################################
    def binding(self,name,T,TV,d):
        self.include("Deformable_Objects/BINDING.h")
        cppclass="BINDING<%s >"%TV
        self.defs.append("class_<%s,std::auto_ptr<%s >,boost::noncopyable>(\"%s\",no_init)\n"%(cppclass,cppclass,name))
        self.defs.append(";\n")

    ######################################################################
    # Function linear_binding
    ######################################################################
    def linear_binding(self,name,T,TV,d):
        self.include("Deformable_Objects/LINEAR_BINDING.h")
        self.include("Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h")

        base_cppclass="BINDING<%s >"%TV
        cppclass="LINEAR_BINDING<%s,%d>"%(TV,d)
        self.defs.append("implicitly_convertible<std::auto_ptr<%s > ,std::auto_ptr<%s > >();\n"%(cppclass,base_cppclass))
        self.defs.append("class_<%s,bases<%s >,std::auto_ptr<%s >,boost::noncopyable>(\"%s\","%(cppclass,base_cppclass,cppclass,name)
                         +"init<DEFORMABLE_BODY_PARTICLES<%s >&,const int,const VECTOR<int,%d>&,const VECTOR<%s,%d>&>())\n"%(TV,d,T,d))
        self.defs.append("    .def(init<DEFORMABLE_BODY_PARTICLES<%s >&,const int,const VECTOR<int,%d>&,const VECTOR<%s,%d>&>())\n"%(TV,d,T,d-1))
        self.defs.append(";\n")

    ######################################################################
    # Function super_fragment
    ######################################################################
    def super_fragment(self,name,TV):
        self.include("Deformable_Objects/FRAGMENT.h")
        cppclass="SUPER_FRAGMENT<%s >"%(TV)
        (T,d)=parse_template(TV)[1]
        d=int(d)

        self.defs.append("class_<%s,boost::noncopyable>(\"%s\",no_init)\n"%(cppclass,name))
        for data in [("id","readonly")]:
            self.member_variable(cppclass,data)

        self.defs.append(";\n")

    ######################################################################
    # Function binding_list
    ######################################################################
    def binding_list(self,name,T,TV):
        self.include("Deformable_Objects/BINDING_LIST.h")
        cppclass="BINDING_LIST<%s >"%TV
        self.routines.append('void Add_Binding_%s(BINDING_LIST<%s >& bl,std::auto_ptr<BINDING<%s > > b){bl.Add_Binding(b.get());b.release();}\n'%(name,TV,TV))

        self.defs.append("class_<%s,boost::noncopyable>(\"%s\",no_init)\n"%(cppclass,name))
        self.defs.append('    .def("Add_Binding",make_function(&::Add_Binding_%s))\n'%name)
        for member in ["Update_Binding_Index_From_Particle_Index","Binding_Index_From_Particle_Index",
                       #("Clamp_Particles_To_Embedded_Positions","void",["int"],"const"),
                       #("Clamp_Particles_To_Embedded_Positions","void",["int","ARRAY<%s >&"%TV],"const"),
                       #("Clamp_Particles_To_Embedded_Velocities","void",["int"],"const"),
                       # TODO: Other clamp Velocities requires TWIST
                       ("Embedded_Position",TV,["const int"],"const"),
                       ("Embedded_Position",TV,["const int","RAW_ARRAY<const %s >"%TV],"const"),
                       ("Embedded_Velocity",TV,["const int"],"const"),
                       ("Binding","BINDING<%s >*"%TV,["const int"],"const","return_internal_reference<>()"),
                       "V","Apply_Impulse",
                       ("One_Over_Effective_Mass",T,["const int"],"const"),
                       ("One_Over_Effective_Mass",T,["const int","const %s&"%TV],"const"),
                       ('Distribute_Force_To_Parents','void',['const SUPER_FRAGMENT<%s >&'%TV,'RAW_ARRAY<%s >'%TV],'const'),
                       ('Distribute_Force_To_Parents','void',['const SUPER_FRAGMENT<%s >&'%TV,'RAW_ARRAY<%s >'%TV,'RAW_ARRAY<TWIST<%s > >'%TV],'const'),
                       ("Clear_Hard_Bound_Particles","void",["RAW_ARRAY<%s >&"%T],"const"),
                       "Distribute_Mass_To_Parents",
                       "Parents","Dynamic_Parents"]:
            self.member_function(cppclass,member)
        self.defs.append(";\n")

    ######################################################################
    # Function solids_parameters
    ######################################################################
    def solids_parameters(self,name,T,d):
        self.include("Solids_And_Fluids/SOLIDS_PARAMETERS.h")
        self.include("Solids_And_Fluids/SOLID_BODY_COLLECTION.h")

        cppclass="SOLIDS_PARAMETERS<VECTOR<%s,%d> >"%(T,d)
        self.defs.append("class_<%s,boost::noncopyable>(\"%s\",no_init)\n"%(cppclass,name))

        self.reference(cppclass,"solid_body_collection")
        self.reference(cppclass,"triangle_repulsions")
        # TODO: collision_body_List, solids_evolution, fracture_evolution
        # Set_Evolution, Set_Fracture_Evolution, Initialize_Triangle_Collisions, Read_Output_Files, Write_Output_Files
        for data in ["write_deformable_body","verbose","verbose_dt","cfl","min_dt","gravity","gravity_direction","cg_iterations","cg_tolerance",
                     "spectral_analysis","lanczos_iterations","perform_collision_body_collisions","collision_tolerance","collide_with_interior","fracture",
                     "write_static_variables_every_frame","newton_tolerance","newton_iterations","use_partially_converged_result",
                     "use_spatial_partition_for_levelset_collision_objects","disable_multiple_levelset_collisions","maximum_levelset_collision_projection_velocity",
                     #"spatial_partition_voxel_size_heuristic"
                     "spatial_partition_number_of_cells","spatial_partition_voxel_size_scale_factor", # spatial_partition_union, mpi_solids
                     # TODO: all collision parameters
                     "perform_self_collision","repulsion_pair_update_count","repulsion_pair_update_frequency",
                     "total_collision_loops","topological_hierarchy_build_count","topological_hierarchy_build_frequency",
                     "check_initial_mesh_for_self_intersection","check_mesh_for_self_intersection","turn_off_all_collisions","allow_intersections",
                     "allow_intersections_tolerance","collisions_small_number","collisions_repulsion_thickness","clamp_repulsion_thickness",
                     "collisions_repulsion_clamp_fraction","collisions_collision_thickness",
                     "repulsions_youngs_modulus","collisions_final_repulsion_youngs_modulus","repulsions_limiter_fraction","collisions_final_repulsion_limiter_fraction",
                     "collisions_disable_repulsions_based_on_proximity_factor",
                     "collisions_output_repulsion_results","collisions_output_collision_results","collisions_output_number_checked","self_collision_friction_coefficient",
                     "collisions_nonrigid_collision_attempts","perform_per_time_step_repulsions","perform_per_collision_step_repulsions","output_interaction_pairs",
                     "perform_per_collision_step_repulsions","repulsion_pair_attractions_threshold",
                     "use_rigid_deformable_contact","use_push_out","use_epsilon_scaling","use_epsilon_scaling_for_level","use_shock_propagation",
                     "collision_iterations","contact_iterations","contact_project_iterations","cg_projection_iterations","use_trapezoidal_rule_for_velocities",
                     "throw_exception_on_backward_euler_failure","rigid_collisions_use_triangle_hierarchy","rigid_collisions_use_triangle_hierarchy_center_phi_test",
                     "rigid_collisions_use_edge_intersection","rigid_collisions_print_interpenetration_statistics",
#                     "rigid_collisions_spatial_partition_voxel_size_heuristic",
                     "rigid_collisions_spatial_partition_number_of_cells","rigid_collisions_spatial_partition_voxel_size_scale_factor","simulate_rigid_bodies","write_rigid_bodies",
                     "rigid_body_ether_viscosity","max_rigid_body_rotation_per_time_step","max_rigid_body_linear_movement_fraction_per_time_step","minimum_rigid_body_time_step_fraction",
                     "maximum_rigid_body_time_step_fraction","clamp_rigid_body_velocities","max_rigid_body_linear_velocity","max_rigid_body_angular_velocity","rigid_cfl",
                     "rigid_minimum_dt","rigid_maximum_dt"]:
            self.member_variable(cppclass,data)

        self.defs.append(";\n")

    ######################################################################
    # Function triangle_repulsions
    ######################################################################
    def triangle_repulsions(self,name,TV):
        self.include("Collisions_And_Interactions/TRIANGLE_REPULSIONS.h")
        T,d=parse_template(TV)[1]
        d=int(d)

        cppclass="TRIANGLE_REPULSIONS<VECTOR<%s,%d> >"%(T,d)
        self.defs.append("class_<%s,boost::noncopyable>(\"%s\",no_init)\n"%(cppclass,name))
        for i in ["hierarchy_repulsion_thickness_multiplier","compute_point_face_friction","compute_point_face_inelastic_collision_repulsion","compute_point_face_repulsion",
                  "compute_edge_edge_friction","compute_edge_edge_inelastic_collision_repulsion","compute_edge_edge_repulsion",
                  "output_repulsion_results","attractions_threshold"]:
            self.member_variable(cppclass,i)
        self.defs.append(";\n")

    ######################################################################
    # Function solid_Standard_tests
    ######################################################################
    def solids_standard_tests(self,name,TV):
        T,d=parse_template(TV)[1]
        d=int(d)
        self.include("Standard_Tests/SOLIDS_STANDARD_TESTS.h")
        self.include("Forces_And_Torques/GRAVITY.h")
        self.include("PhysBAM_Geometry/Geometry/SEGMENTED_CURVE_2D.h")
        self.include("PhysBAM_Geometry/Geometry/TETRAHEDRALIZED_VOLUME.h")
        self.include("PhysBAM_Geometry/Geometry/TRIANGULATED_SURFACE.h")
        self.include("Solids_And_Fluids/SOLIDS_PARAMETERS.h")
        self.include("Rigid_Bodies/RIGID_BODY_2D.h")
        self.include("Rigid_Bodies/RIGID_BODY_3D.h")
        self.local_include("SOLIDS_FLUIDS_EXAMPLE_UNIFORM_PYTHON.h")
        tri="GEOMETRY_POLICY<%s >::TRIANGULATED_OBJECT"%TV
        cppclass="SOLIDS_STANDARD_TESTS<%s >"%TV
        self.defs.append("class_<%s,boost::noncopyable>(\"%s\",init<SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<%s > >&>()[with_custodian_and_ward<1,2>()])\n"%(cppclass,name,TV))
        members=["Set_Initial_Particle_Configuration",
                 ("Add_Rigid_Body","RIGID_BODY<%s >&"%TV,["const std::string&","const %s"%T,"const %s"%T,"const bool","const bool"],"","return_internal_reference<>()"),
                 ("Add_Ground","RIGID_BODY<%s >&"%TV,["const %s friction"%T,"const %s height"%T,"const %s coefficient_of_restitution"%T,"const %s scale"%T],"",
                  "return_internal_reference<>()"),
                 "Add_Gravity",
                 ("Create_Segmented_Curve","SIMPLEX_POLICY<%s,1>::OBJECT&"%TV,["const int","const RIGID_BODY_STATE<%s >&"%TV],"","return_internal_reference<>()"),
                 ("Create_Segmented_Curve","SIMPLEX_POLICY<%s,1>::OBJECT&"%TV,["const std::string&","const RIGID_BODY_STATE<%s >&"%TV,"const bool","const bool"],
                  "","return_internal_reference<>()")]
        if d==3:
            members.extend([
                ("Create_Tetrahedralized_Volume","TETRAHEDRALIZED_VOLUME<%s >&"%T,["const std::string&","const RIGID_BODY_STATE<%s >&"%TV,"const bool","const bool","const %s"%T,"const %s"%T],
                 "","return_internal_reference<>()"),
                ("Create_Triangulated_Object","return_internal_reference<>()"),
                ("Create_Cloth_Panel","TRIANGULATED_SURFACE<%s >&"%T,["const int","const %s"%T,"const %s"%T,"const RIGID_BODY_STATE<%s >&"%TV,"LIST_ARRAY<int>*"],
                 "","return_internal_reference<>()"),"Make_Lathe_Chain"])

        for member in members:
            self.member_function(cppclass,member)
        self.defs.append(";\n")

    ######################################################################
    # Function solid_body_collection
    ######################################################################
    def solid_body_collection(self,name,TV):
        self.include("Solids_And_Fluids/SOLID_BODY_COLLECTION.h")

        cppclass="SOLID_BODY_COLLECTION<%s >"%TV
        self.routines.append('void Update_Fragments_%s(SOLID_BODY_COLLECTION<%s >& d){d.Update_Fragments();}\n'%(name,TV))

        self.defs.append("class_<%s,boost::noncopyable>(\"%s\",no_init)\n"%(cppclass,name))
        self.member_function(cppclass,("Add_Force","int",["SOLIDS_FORCES<%s >*"%TV]))
        self.defs.append('    .def("Update_Fragments",make_function(&::Update_Fragments_%s))\n'%name)
        for data in [("deformable_body_collection","readonly"),("rigid_body_collection","readonly"),("soft_bindings","readonly"),("binding_list","readonly")]:
            self.member_variable(cppclass,data)

        self.defs.append(";\n")
        
    ######################################################################
    # Function deformable_body_collection
    ######################################################################
    def deformable_body_collection(self,name,TV):
        self.include("Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h")

        cppclass="DEFORMABLE_BODY_COLLECTION<%s >"%TV

        self.defs.append("class_<%s,boost::noncopyable>(\"%s\",no_init)\n"%(cppclass,name))
        for data in [("particles","readonly"),("deformable_geometry","readonly")]:
            self.member_variable(cppclass,data)

        self.defs.append(";\n")

    ######################################################################
    # Function deformable_geometry_collection
    ######################################################################
    def deformable_geometry_collection(self,name,TV):
        self.include("PhysBAM_Geometry/Solids/DEFORMABLE_GEOMETRY_COLLECTION.h")
        self.include("PhysBAM_Geometry/Geometry/STRUCTURE.h")
        
        cppclass="DEFORMABLE_GEOMETRY_COLLECTION<%s >"%TV
        self.routines.append('void Add_Structure_%s(DEFORMABLE_GEOMETRY_COLLECTION<%s >& d,STRUCTURE<%s >* s){d.Add_Structure(s);}\n'%(name,TV,TV))

        self.defs.append("class_<%s,boost::noncopyable>(\"%s\",no_init)\n"%(cppclass,name))
        self.defs.append('    .def("Add_Structure",make_function(&::Add_Structure_%s))\n'%name)
        self.defs.append(";\n")

    ######################################################################
    # Function soft_bindings
    ######################################################################
    def soft_bindings(self,name,TV):
        self.include("Deformable_Objects/SOFT_BINDINGS.h")

        cppclass="SOFT_BINDINGS<%s >"%TV
        self.defs.append("class_<%s,boost::noncopyable>(\"%s\",no_init)\n"%(cppclass,name))
        for member in ["Add_Binding","Initialize_Binding_Mesh","Set_Mass_From_Effective_Mass"]:
            self.member_function(cppclass,member)
        self.reference(cppclass,"binding_mesh")
        self.defs.append(";\n")

    ######################################################################
    # Function gravity
    ######################################################################
    def gravity(self,name,T,TV):
        self.include("Forces_And_Torques/GRAVITY.h")
        self.include("Rigid_Bodies/RIGID_BODY_2D.h")
        self.include("Rigid_Bodies/RIGID_BODY_3D.h")
        base_cppclass="SOLIDS_FORCES<%s >"%TV
        cppclass="GRAVITY<%s >"%TV
        self.defs.append("class_<%s,bases<%s >,boost::noncopyable>(\"%s\",init<DEFORMABLE_BODY_PARTICLES<%s >&,RIGID_BODY_COLLECTION<%s >&,bool,bool,const %s,const %s&>())\n"%(cppclass,base_cppclass,name,TV,TV,T,TV))
        self.defs.append("    .def(init<DEFORMABLE_BODY_PARTICLES<%s >&,RIGID_BODY_COLLECTION<%s >&,LIST_ARRAY<int>*,LIST_ARRAY<int>*,const %s,const %s&>())\n"%(TV,TV,T,TV))
        for member in ["Set_Gravity"]:
            self.member_function(cppclass,member)

        self.defs.append(";\n")

    ######################################################################
    # Function linear_springs
    ######################################################################
    def linear_springs(self,name,T,TV):
        T,d=parse_template(TV)[1]
        d=int(d)

        self.include("Forces_And_Torques/LINEAR_SPRINGS.h")
        self.include("PhysBAM_Geometry/Geometry/TETRAHEDRALIZED_VOLUME.h")
        self.include("PhysBAM_Geometry/Geometry/TRIANGULATED_SURFACE.h")
        self.include("PhysBAM_Geometry/Geometry/SEGMENTED_CURVE_2D.h")
        base_cppclass="SOLIDS_FORCES<%s >"%TV
        cppclass="LINEAR_SPRINGS<%s >"%TV
        self.defs.append("class_<%s,bases<%s >,boost::noncopyable>(\"%s\",init<DEFORMABLE_BODY_PARTICLES<%s >&,SEGMENT_MESH&,int>())\n"%(cppclass,base_cppclass,name,TV))
        for member in ["Set_Restlength_From_Particles","Print_Restlength_Statistics","Clamp_Restlength",
                       ("Set_Overdamping_Fraction","void",["const %s"%T])]:
            self.member_function(cppclass,member)
        self.defs.append(";\n")
        count=0
        geoms=[]
        if d==3: geoms=["TRIANGULATED_SURFACE<%s >"%T,"TETRAHEDRALIZED_VOLUME<%s >"%T,"SIMPLEX_POLICY<%s,1>::OBJECT"%TV]
        elif d==2: geoms=["SEGMENTED_CURVE_2D<%s >"%T]
        for geom in geoms:
            count+=1
            self.defs.append("%s* (*Create_Edge_Springs_%s_%d)(%s&,%s,%s,bool,%s,bool,%s,bool,bool)=Create_Edge_Springs;\n"%(cppclass,name,count,geom,T,T,T,T))
            self.defs.append('def("Create_Edge_Springs",Create_Edge_Springs_%s_%d,return_value_policy<manage_new_object>());\n'%(name,count))

    ######################################################################
    # Function segment_adhesion
    ######################################################################
    def segment_adhesion(self,name,T,TV):
        T,d=parse_template(TV)[1]
        d=int(d)

        self.include("Forces_And_Torques/SEGMENT_ADHESION.h")
        self.include("PhysBAM_Geometry/Geometry/SEGMENTED_CURVE.h")
        base_cppclass="SOLIDS_FORCES<%s >"%TV
        cppclass="SEGMENT_ADHESION<%s >"%TV

        self.defs.append("class_<%s,bases<%s >,boost::noncopyable>(\"%s\",init<DEFORMABLE_BODY_PARTICLES<%s >&,SEGMENT_MESH&,LIST_ARRAY<HAIR_ID>&,HASHTABLE<VECTOR<int,4> >&>())\n"%(cppclass,base_cppclass,name,TV))
        self.defs.append("    .def(init<DEFORMABLE_BODY_PARTICLES<%s >&,SEGMENT_MESH&,LIST_ARRAY<HAIR_ID>&>())"%TV)

        for member in ["Set_Parameters","Set_Restlength","Write_State","Update_Springs","Update_Partitions"]:
            self.member_function(cppclass,member)
        self.defs.append(";\n")

    ######################################################################
    # Function linear_tet_springs
    ######################################################################
    def linear_tet_springs(self,name,T):
        TV="VECTOR<%s,3>"%T

        self.include("Particles/DEFORMABLE_BODY_PARTICLES.h")
        self.include("Forces_And_Torques/LINEAR_TET_SPRINGS.h")
        self.include("PhysBAM_Geometry/Geometry/TETRAHEDRALIZED_VOLUME.h")
        base_cppclass="SOLIDS_FORCES<%s >"%TV
        cppclass="LINEAR_TET_SPRINGS<%s >"%T
        self.defs.append("class_<%s,bases<%s >,boost::noncopyable>(\"%s\",init<DEFORMABLE_BODY_PARTICLES<%s >&,TETRAHEDRON_MESH&,bool>())\n"%(cppclass,base_cppclass,name,TV))
        for member in ["Print_Restlength_Statistics","Clamp_Restlength",
                       ("Set_Overdamping_Fraction","void",["const %s"%T])]:
            self.member_function(cppclass,member)
        self.defs.append(";\n")
        count=0
        geoms=[]
        geoms=["TETRAHEDRALIZED_VOLUME<%s >"%T]
        for geom in geoms:
            count+=1
            self.defs.append("%s* (*Create_Tet_Springs_%s_%d)(%s&,%s,%s,bool,%s,bool,%s,bool,%s,bool,bool)=Create_Tet_Springs;\n"%(cppclass,name,count,geom,T,T,T,T,T))
            self.defs.append('def("Create_Tet_Springs",Create_Tet_Springs_%s_%d,return_value_policy<manage_new_object>());\n'%(name,count))

    ######################################################################
    # Function linear_tet_springs
    ######################################################################
    def linear_altitude_springs_3d(self,name,T):
        TV="VECTOR<%s,3>"%T

        self.include("Forces_And_Torques/LINEAR_ALTITUDE_SPRINGS_3D.h")
        self.include("PhysBAM_Geometry/Geometry/TETRAHEDRALIZED_VOLUME.h")
        base_cppclass="SOLIDS_FORCES<%s >"%TV
        cppclass="LINEAR_ALTITUDE_SPRINGS_3D<%s >"%T
        self.defs.append("class_<%s,bases<%s >,boost::noncopyable>(\"%s\",init<DEFORMABLE_BODY_PARTICLES<%s >&,TETRAHEDRON_MESH&>())\n"%(cppclass,base_cppclass,name,TV))
        for member in ["Print_Restlength_Statistics","Clamp_Restlength",
                       ("Set_Overdamping_Fraction","void",["const %s"%T])]:
            self.member_function(cppclass,member)
        self.defs.append(";\n")
        count=0
        geoms=[]
        geoms=["TETRAHEDRALIZED_VOLUME<%s >"%T]
        for geom in geoms:
            count+=1
            self.defs.append("%s* (*Create_Altitude_Springs_%s_%d)(%s&,%s,%s,bool,%s,bool,%s,bool,%s,bool)=Create_Altitude_Springs;\n"%(cppclass,name,count,geom,T,T,T,T,T))
            self.defs.append('def("Create_Altitude_Springs",Create_Altitude_Springs_%s_%d,return_value_policy<manage_new_object>());\n'%(name,count))

    ######################################################################
    # Function triangle_bending_springs
    ######################################################################
    def triangle_bending_springs(self,name,T):
        self.include("Forces_And_Torques/TRIANGLE_BENDING_SPRINGS.h")
        self.include("PhysBAM_Geometry/Geometry/TRIANGULATED_SURFACE.h")
        self.include("PhysBAM_Geometry/Geometry/SEGMENTED_CURVE_2D.h")
        TV="VECTOR<%s,3>"%T
        base_cppclass="SOLIDS_FORCES<%s >"%TV
        cppclass="TRIANGLE_BENDING_SPRINGS<%s >"%T
        self.defs.append("class_<%s,bases<%s >,boost::noncopyable>(\"%s\",init<DEFORMABLE_BODY_PARTICLES<%s >&,TRIANGLE_MESH&>())\n"%(cppclass,base_cppclass,name,TV))
        self.defs.append(";\n")
        count=0
        for geom in ["TRIANGULATED_SURFACE<%s >"%T]:
            count+=1
            self.defs.append("%s* (*Create_Bending_Springs_%s_%d)(%s&,const %s,const %s,const bool,const %s,const bool,const %s,const bool,const bool)=Create_Bending_Springs;\n"%(cppclass,name,count,geom,T,T,T,T))
            self.defs.append('def("Create_Bending_Springs",Create_Bending_Springs_%s_%d,return_value_policy<manage_new_object>());\n'%(name,count))

    ######################################################################
    # Function triangle_bending_springs
    ######################################################################
    def triangle_bending_elements(self,name,T):
        self.include("Forces_And_Torques/TRIANGLE_BENDING_ELEMENTS.h")
        self.include("PhysBAM_Geometry/Geometry/TRIANGULATED_SURFACE.h")
        self.include("PhysBAM_Geometry/Geometry/SEGMENTED_CURVE_2D.h")
        TV="VECTOR<%s,3>"%T
        base_cppclass="SOLIDS_FORCES<%s >"%TV
        cppclass="TRIANGLE_BENDING_ELEMENTS<%s >"%T
        self.defs.append("class_<%s,bases<%s >,boost::noncopyable>(\"%s\",init<DEFORMABLE_BODY_PARTICLES<%s >&>())\n"%(cppclass,base_cppclass,name,TV))
        self.defs.append(";\n")
        count=0
        for geom in ["TRIANGULATED_SURFACE<%s >"%T]:
            count+=1
            self.defs.append("%s* (*Create_Bending_Elements_%s_%d)(%s&,const %s,const %s,const bool,const %s,const bool,const %s,const %s,const %s,const %s,const bool)=Create_Bending_Elements;\n"
                             %(cppclass,name,count,geom,T,T,T,T,T,T,T))
            self.defs.append('def("Create_Bending_Elements",Create_Bending_Elements_%s_%d,return_value_policy<manage_new_object>());\n'%(name,count))

    ######################################################################
    # Function segment_bending_springs
    ######################################################################
    def segment_bending_springs(self,name,TV):
        T,d=parse_template(TV)[1]
        d=int(d)

        self.include("Forces_And_Torques/SEGMENT_BENDING_SPRINGS.h")
        self.include("PhysBAM_Geometry/Geometry/TRIANGULATED_SURFACE.h")
        self.include("PhysBAM_Geometry/Geometry/SEGMENTED_CURVE.h")
        base_cppclass="SOLIDS_FORCES<%s >"%TV
        cppclass="SEGMENT_BENDING_SPRINGS<%s >"%TV
        self.defs.append("class_<%s,bases<%s >,boost::noncopyable>(\"%s\",init<DEFORMABLE_BODY_PARTICLES<%s >&,SEGMENT_MESH&>())\n"%(cppclass,base_cppclass,name,TV))
        self.defs.append(";\n")
        count=0
        for geom in ["SEGMENTED_CURVE<%s >"%TV]:
            count+=1
            self.defs.append("%s* (*Create_Segment_Bending_Springs_%s_%d)(%s&,const %s,const %s,const bool,const %s,const bool,const %s,const bool,const bool)=Create_Segment_Bending_Springs;\n"%(cppclass,name,count,geom,T,T,T,T))
            self.defs.append('def("Create_Segment_Bending_Springs",Create_Segment_Bending_Springs_%s_%d,return_value_policy<manage_new_object>());\n'%(name,count))

    ######################################################################
    # Function binding_springs
    ######################################################################
    def binding_springs(self,name,T,TV):
        self.include("Forces_And_Torques/BINDING_SPRINGS.h")
        base_cppclass="SOLIDS_FORCES<%s >"%TV
        cppclass="BINDING_SPRINGS<%s >"%(TV)
        self.defs.append("class_<%s,bases<%s >,boost::noncopyable>(\"%s\",no_init)\n"%(cppclass,base_cppclass,name))
        self.defs.append(";\n")
        count=1
        self.defs.append('    %s* (*Create_Edge_Binding_Springs_%s_%d) (DEFORMABLE_BODY_PARTICLES<%s >&,SEGMENT_MESH&,const %s,const %s,const bool)=Create_Edge_Binding_Springs;\n'%(cppclass,name,count,TV,T,T))
        self.defs.append('    def("Create_Edge_Binding_Springs",Create_Edge_Binding_Springs_%s_%d,return_value_policy<manage_new_object>());\n'%(name,count))

    ######################################################################
    # Function finite_volume
    ######################################################################
    def finite_volume(self,name,T,TV,d):
        self.include("Forces_And_Torques/FINITE_VOLUME.h")
        self.include("PhysBAM_Geometry/Geometry/TETRAHEDRALIZED_VOLUME.h")
        self.include("PhysBAM_Geometry/Geometry/TRIANGULATED_SURFACE.h")
        self.include("PhysBAM_Geometry/Geometry/SEGMENTED_CURVE_2D.h")
        base_cppclass="SOLIDS_FORCES<%s >"%TV
        cppclass="FINITE_VOLUME<%s,%d>"%(TV,d)
        self.defs.append("class_<%s,bases<%s >,boost::noncopyable>(\"%s\",no_init)\n"%(cppclass,base_cppclass,name))
        self.defs.append(";\n")
        count=0

        if d==2: geom ="TRIANGULATED_SURFACE<%s >"%T
        elif d==3: geom="TETRAHEDRALIZED_VOLUME<%s >"%T
        count+=1
        self.routines.append("%s* Create_Finite_Volume_%s_%d_Wrapper(%s& geom,std::auto_ptr<CONSTITUTIVE_MODEL<%s,%d> > constitutive_model)\n"%(cppclass,name,count,geom,T,d))
        self.routines.append("{%s* ret=Create_Finite_Volume(geom,constitutive_model.get());constitutive_model.release();return ret;}\n"%cppclass)
        self.defs.append('def("Create_Finite_Volume",&::Create_Finite_Volume_%s_%d_Wrapper,boost::python::return_value_policy<boost::python::manage_new_object>());\n'%(name,count))

    ######################################################################
    # Function constitutive_model
    ######################################################################
    def constitutive_model(self,name,T,d):
        self.include("Constitutive_Models/CONSTITUTIVE_MODEL.h")
        cppclass="CONSTITUTIVE_MODEL<%s,%d>"%(T,d)
        self.defs.append("class_<%s,std::auto_ptr<%s >,boost::noncopyable>(\"%s\",no_init)\n"%(cppclass,cppclass,name))
        self.defs.append(";\n")

    ######################################################################
    # Function neo_hookean
    ######################################################################
    def neo_hookean(self,name,T,d):
        self.include("Constitutive_Models/NEO_HOOKEAN.h")
        base_cppclass="CONSTITUTIVE_MODEL<%s,%d>"%(T,d)
        cppclass="NEO_HOOKEAN<%s,%d>"%(T,d)
        self.defs.append("implicitly_convertible<std::auto_ptr<%s > ,std::auto_ptr<%s > >();\n"%(cppclass,base_cppclass))
        self.defs.append("class_<%s,bases<%s >,std::auto_ptr<%s >,boost::noncopyable>(\"%s\",init<const %s,const %s,const %s,const %s >())\n"%(cppclass,base_cppclass,cppclass,name,T,T,T,T))
        self.defs.append(";\n")

    ######################################################################
    # Function solids_fluids_driver_uniform
    ######################################################################
    def solids_fluids_driver_uniform(self,name,T_GRID):
        self.include("Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_UNIFORM.h")
        self.include("PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS_3D.h")

        cppclass="SOLIDS_FLUIDS_DRIVER_UNIFORM<%s >"%(T_GRID)
        self.defs.append("class_<%s,boost::noncopyable>(\"%s\",init<SOLIDS_FLUIDS_EXAMPLE_UNIFORM<%s >&>())\n"%(cppclass,name,T_GRID))
        for data in ["Execute_Main_Program"]:
            self.member_function(cppclass,data)
        #self.defs.append(".staticmethod(\"Create\")\n")

        self.defs.append(";\n")

    ######################################################################
    # Function rotation
    ######################################################################
    def rotation(self,name,TV):
        self.include("PhysBAM_Tools/Matrices/ROTATION.h")
        T=parse_template(TV)[1][0]
        d=int(parse_template(TV)[1][1])
        T_SPIN="VECTOR<%s,%d>"%(T,d*(d-1)/2)

        cppclass="ROTATION<%s >"%TV
        self.defs.append('class_<%s >("%s",init<>())\n'%(cppclass,name))
        if d==3: self.defs.append("    .def(init<const %s,const %s&>())"%(T,TV))
        for i in ["Rotate","Inverse_Rotate"]:
            if d==1:
                self.member_function(cppclass,(i,'const %s&'%TV,['const %s&'%TV],'const','return_value_policy<copy_const_reference>()'))
            else:
                self.member_function(cppclass,(i,TV,['const %s&'%TV],'const'))
        for i in ["Inverse","Normalize","Normalized","Is_Normalized","Angle",
                  ("Euler_Angles",T_SPIN,[],"const"),
                  "Rotation_Vector","Rotation_Matrix","Scale_Angle"]:
            self.member_function(cppclass,i)
        for i in [("From_Euler_Angles","ROTATION<%s >"%TV,["const %s&"%T_SPIN]),
                  "Average_Rotation","Spherical_Linear_Interpolation","From_Rotation_Vector","From_Rotated_Vector"]:
            self.static_function(cppclass,i)
        self.str(cppclass)
        self.repr(cppclass)
        self.defs.extend(["    .def(self*=self)\n"])
        self.defs.extend(["    .def(self*self)\n"])
        self.defs.append(";\n")

    ######################################################################
    # Function twist
    ######################################################################
    def twist(self,name,TV):
        self.include("PhysBAM_Tools/Vectors/TWIST.h")
        cppclass="TWIST<%s >"%TV
        self.defs.append('class_<%s >("%s",init<const %s,const %s::SPIN&>())\n'%(cppclass,name,TV,TV))
        for member in ["linear","angular"]:
            self.member_variable(cppclass,member)
        self.str(cppclass)
        self.repr(cppclass)
        self.defs.append(";\n")

    ######################################################################
    # Function rigid_body
    ######################################################################
    def rigid_body(self,name,TV):
        self.include("Rigid_Bodies/RIGID_BODY.h")
        self.include("PhysBAM_Geometry/Topology/TRIANGLE_MESH.h")
        self.include("PhysBAM_Geometry/Geometry/MESH_OBJECT.h")
        self.include("PhysBAM_Geometry/Geometry/SEGMENT_2D.h")
        self.include("PhysBAM_Geometry/Geometry/TRIANGLE_3D.h")

        cppclass="RIGID_BODY<%s >"%(TV)
        self.defs.append("class_<%s,boost::noncopyable>(\"%s\",no_init)\n"%(cppclass,name))
        variables=["rigid_body_particle_index","surface_roughness",
                   # TODO: structures, body_forces, fluid_forces, example_forces_and_velocities, example_forces_and_velocities_default, force_accumulator, implicit_object, simplicial_object,
                   # TODO: moving_simplex_hierarchy, oriented_box, axis_aligned_bounding_box, bounding_box_up_to_date
                   "coefficient_of_restitution","coefficient_of_rolling_friction","name","is_static","is_kinematic","is_temporarily_static", # TODO: sAved_states
                   "thin_shell","CFL_initialized","bounding_box_radius"]
        for v in variables:
            self.member_variable(cppclass,v)
        functions=[('Add_Structure','void',['STRUCTURE<%s >&'%TV]),
                   "Set_Name","Set_Mass","Set_Inertia_Tensor",
                   "Rescale","Set_Surface_Roughness","Set_Coefficient_Of_Restitution", # TODO: Find_Structure, REmove_Structure Coefficient_Of_Restitution, COefficient_Of_Friction
                   "Set_Coefficient_Of_Rolling_Friction","Has_Infinite_Inertia","Is_Simulated", # TODO: Coefficient_Of_Rolling_Friction, Object_Space_Ray
                   "Object_Space_Point","Object_Space_Vector","World_Space_Point","World_Space_Vector",
                   ("Frame","FRAME<%s >&"%TV,[],"","return_internal_reference<>()"),
                   ("Twist","TWIST<%s >&"%TV,[],"","return_internal_reference<>()"),
                   ("Angular_Momentum","%s::SPIN&"%TV,[],"","return_internal_reference<>()")] # TODO: more wrapping
        for member in functions:
            self.member_function(cppclass,member)
        self.defs.append(";\n")

    ######################################################################
    # Function rigid_body_2d
    ######################################################################
    def rigid_body_2d(self,name,T):
        self.include("PhysBAM_Geometry/Geometry/TRIANGULATED_SURFACE.h")
        self.include("Rigid_Bodies/RIGID_BODY_2D.h")
        TV="VECTOR<%s,2>"%T
        base_cppclass="RIGID_BODY_BASE<%s >"%TV
        cppclass="RIGID_BODY<%s >"%(TV)
        #self.defs.append("implicitly_convertible<%s,%s >();\n"%(cppclass,base_cppclass))
        self.defs.append("class_<%s,bases<%s >,boost::noncopyable>(\"%s\",no_init)\n"%(cppclass,base_cppclass,name))
        functions=["Simplex_Inside","Simplex_Outside",
                   ("World_Space_Simplex_Bounding_Box","RANGE<%s >"%TV,["const int"],"const"),# TODO: wrap this guy("World_Space_Simplex_Bounding_Box","BOX_3D<T>",["const int","RIGID_BODY_STATE<TV>&"])
                   ]
        # TODO: "Implicit_Geometry_Intersection", TODO: more wrapping
        for member in functions:
            self.member_function(cppclass,member)
        self.defs.append(";\n")


    ######################################################################
    # Function rigid_body_3d
    ######################################################################
    def rigid_body_3d(self,name,T):
        self.include("PhysBAM_Geometry/Geometry/TRIANGULATED_SURFACE.h")
        self.include("Rigid_Bodies/RIGID_BODY_3D.h")
        TV="VECTOR<%s,3>"%T
        base_cppclass="RIGID_BODY_BASE<%s >"%TV
        cppclass="RIGID_BODY<%s >"%(TV)
        #self.defs.append("implicitly_convertible<%s,%s >();\n"%(cppclass,base_cppclass))
        self.defs.append("class_<%s,bases<%s >,boost::noncopyable>(\"%s\",no_init)\n"%(cppclass,base_cppclass,name))
        functions=["Simplex_Normal","Simplex_Inside","Simplex_Outside","Simplex_Boundary",
                   ("World_Space_Simplex_Bounding_Box","RANGE<%s >"%TV,["const int"],"const"),# TODO: wrap this guy("World_Space_Simplex_Bounding_Box","BOX_3D<T>",["const int","RIGID_BODY_STATE<TV>&"])
                   ("Simplex_World_Space_Point_From_Barycentric_Coordinates",TV,["const int","const %s&"%TV],"const")]
        # TODO: "Implicit_Geometry_Intersection", TODO: more wrapping
        for member in functions:
            self.member_function(cppclass,member)
        self.defs.append(";\n")


    ######################################################################
    # Function frame
    ######################################################################
    def frame(self,name,TV):
        self.include("PhysBAM_Tools/Matrices/FRAME.h")
        cppclass="FRAME<%s >"%TV
        T,d=parse_template(TV)[1]

        self.defs.append("class_<%s >(\"%s\")\n"%(cppclass,name))
        self.defs.append("    .def(init<const %s&>())\n"%TV)
        if int(d)==3:
            self.defs.append("    .def(init<const ROTATION<%s >&>())\n"%(TV))
            self.defs.append("    .def(init<const %s&,const ROTATION<%s >&>())\n"%(TV,TV))
        self.defs.extend(["    .def(self*self)\n",
                          "    .def(self*=self)\n",
                          "    .def(self*%s())\n"%TV])
        for member in ["t","r"]:
            self.member_variable(cppclass,member)
        members=["Inverse","Invert","Interpolation"]
        for member in members:
            self.member_function(cppclass,member)
        self.str(cppclass)
        self.repr(cppclass)
        self.defs.append(";\n")

    ######################################################################
    # Function rigid_body_state
    ######################################################################
    def rigid_body_state(self,name,TV):
        self.include("Rigid_Bodies/RIGID_BODY_STATE.h")
        cppclass="RIGID_BODY_STATE<%s >"%TV

        self.defs.append("class_<%s >(\"%s\")\n"%(cppclass,name))
        self.defs.append("    .def(init<const FRAME<%s >&>())\n"%TV)
        self.defs.append("    .def(init<const FRAME<%s >&,const TWIST<%s >&>())\n"%(TV,TV))
        members=["Compute_Velocity_Between_States","Object_Space_Point","Object_Space_Vector",
                 "World_Space_Point","World_Space_Vector","World_Space_Vector",
                 # TODO: "World_Space_Intertia_Tensor","World_Space_Rigid_Mass"
                 # TODO: "World_Space_Intertia_Tensor_Inverse","World_Space_Rigid_Mass_Inverse"
                 # TODO: "Update_Angular_Velocity",
                 # TODO: "Update_Angular_Momentum","Pointwise_Object_Velocity"
                 ]
        for member in members:
            self.member_function(cppclass,member)
        self.defs.append(";\n")


    ######################################################################
    # Function robust_simplex_interactions
    ######################################################################
    def robust_simplex_interactions(self,name,TV):
        self.include("PhysBAM_Tools/Vectors/VECTOR_3D.h")
        self.include("PhysBAM_Tools/Data_Structures/PAIR.h")
        self.include("Collisions_And_Interactions/ROBUST_SIMPLEX_INTERACTIONS.h")
        cppclass="ROBUST_SIMPLEX_INTERACTIONS<%s >"%TV

        self.defs.append("class_<%s >(\"%s\")\n"%(cppclass,name))
        #members=[("Intersection","bool",["const VECTOR<%s,3>&,const"%TV,"VECTOR<%s,2>&"%TV,"bool*"])]
        members=[("Intersection_Test","VECTOR<bool,2>",["const VECTOR<%s,3>&"%TV,"const VECTOR<%s,2>&"%TV]),
                 ("Intersection_Test","VECTOR<bool,2>",["const VECTOR<%s,3>&"%TV,"const VECTOR<%s,1>&"%TV])]
        #members=[]
        for member in members:
            self.member_function(cppclass,member)
        self.defs.append(";\n")

    ######################################################################
    # Function simplex_interactions
    ######################################################################
    def simplex_interactions(self,name,T):
        self.include("PhysBAM_Tools/Vectors/VECTOR_3D.h")
        self.include("PhysBAM_Tools/Data_Structures/PAIR.h")
        self.include("Collisions_And_Interactions/SIMPLEX_INTERACTIONS.h")
        cppclass="SIMPLEX_INTERACTIONS<%s >"%T

        self.defs.append("class_<%s >(\"%s\")\n"%(cppclass,name))
        #members=[("Intersection","bool",["const VECTOR<%s,3>&,const"%TV,"VECTOR<%s,2>&"%TV,"bool*"])]
        members=["Two_Segment_Intersection_Barycentric_Coordinates"]
        #members=[]
        for member in members:
            self.member_function(cppclass,member)
            self.defs.append(".staticmethod(\"%s\")"%member)
        self.defs.append(";\n")

    ######################################################################
    # Function random_numbers
    ######################################################################
    def random_numbers(self,name,T):
        self.include("PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h")
        self.include("PhysBAM_Tools/Math_Tools/RANGE.h")
        self.include("PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_2X2.h")
        self.include("PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_3X3.h")
        self.include("PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h")
        self.include("PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h")
        self.include("PhysBAM_Tools/Matrices/UPPER_TRIANGULAR_MATRIX_2X2.h")
        self.include("PhysBAM_Tools/Matrices/UPPER_TRIANGULAR_MATRIX_3X3.h")
        self.include("PhysBAM_Tools/Matrices/MATRIX.h")
        self.include("PhysBAM_Tools/Matrices/MATRIX_MXN.h")
        self.include("PhysBAM_Tools/Matrices/ROTATION.h")

        routine="""
template<class T,class T_VECTOR> void Fill_Uniform_Vector_Wrapper(RANDOM_NUMBERS<T>& random,T_VECTOR& v,const T a,const T b)
{random.Fill_Uniform_Vector(v,a,b);}

template<class T,class T_MATRIX> void Fill_Uniform_Matrix_Wrapper(RANDOM_NUMBERS<T>& random,T_MATRIX& m,const T a,const T b)
{random.Fill_Uniform_Matrix(m,a,b);}

template<class T>
boost::python::object Get_Direction_Wrapper(RANDOM_NUMBERS<T>& random,const int d)
{if(d==1) return boost::python::object(random.template Get_Direction<VECTOR<T,1> >());
else if(d==2) return boost::python::object(random.template Get_Direction<VECTOR<T,2> >());
else if(d==3) return boost::python::object(random.template Get_Direction<VECTOR<T,3> >());
else{PyErr_SetString(PyExc_ValueError,"Dimension must be 1, 2, or 3");boost::python::throw_error_already_set();}
PHYSBAM_FATAL_ERROR();} // silence noreturn warnings

template<class T>
boost::python::object Get_Rotation_Wrapper(RANDOM_NUMBERS<T>& random,const int d)
{if(d==1) return boost::python::object(random.template Get_Rotation<VECTOR<T,1> >());
else if(d==2) return boost::python::object(random.template Get_Rotation<VECTOR<T,2> >());
else if(d==3) return boost::python::object(random.template Get_Rotation<VECTOR<T,3> >());
else{PyErr_SetString(PyExc_ValueError,"Dimension must be 1, 2, or 3");boost::python::throw_error_already_set();}
PHYSBAM_FATAL_ERROR();} // silence noreturn warnings

"""
        if routine not in self.routines: self.routines.append(routine)

        cppclass="RANDOM_NUMBERS<%s>"%T
        self.defs.append("class_<%s,boost::noncopyable>(\"%s\")\n"%(cppclass,name))
        members=['Set_Seed','Get_Uniform_Integer','Get_Gaussian','Get_Number']
        members.append(('Get_Uniform_Number',T,['const '+T]*2))
        for d in range(4):
            members.append(('Get_Uniform_Vector','VECTOR<%s,%d>'%(T,d),['const VECTOR<%s,%d>&'%(T,d)]*2))
            self.defs.append("    .def(\"Fill_Uniform_Vector\",&Fill_Uniform_Vector_Wrapper<%s,VECTOR<%s,%d> >)\n"%(T,T,d))
        for d in [1,2,3]:
            members.append(('Get_Uniform_Vector','VECTOR<%s,%d>'%(T,d),['const RANGE<VECTOR<%s,%d> >&'%(T,d)]))
            self.defs.append("    .def(\"Get_Direction\",&Get_Direction_Wrapper<%s >)\n"%T)
            self.defs.append("    .def(\"Get_Rotation\",&Get_Rotation_Wrapper<%s >)\n"%T)
        matrices=(['MATRIX_MXN<%s>'%T]
            +['%s<%s,%d>'%(M,T,d) for M in ['DIAGONAL_MATRIX','SYMMETRIC_MATRIX','UPPER_TRIANGULAR_MATRIX'] for d in [2,3]]
            +['MATRIX<%s,%d,%d>'%(T,m,n) for m in [1,2,3] for n in [1,2,3]])
        for m in matrices:
            self.defs.append("    .def(\"Fill_Uniform_Matrix\",&Fill_Uniform_Matrix_Wrapper<%s,%s >)\n"%(T,m))
        for member in members:
            self.member_function(cppclass,member)
        self.defs.append(";\n")
        self.read_write(cppclass,T)

    ######################################################################
    # Function sobol
    ######################################################################
    def sobol(self,name,TV):
        self.include("PhysBAM_Tools/Random_Numbers/SOBOL.h")
        self.defs.append('{\n')
        self.typedef(TV=TV)
        cppclass='SOBOL<TV>'
        self.defs.append("class_<%s,boost::noncopyable>(\"%s\",init<const RANGE<TV>&>())\n"%(cppclass,name))
        self.member_function(cppclass,'Get_Vector')
        self.defs.append(";}\n")

    ######################################################################
    # Function piecewise_constant_pdf
    ######################################################################
    def piecewise_constant_pdf(self,name,T):
        self.include("PhysBAM_Tools/Random_Numbers/PIECEWISE_CONSTANT_PDF.h")
        cppclass="PIECEWISE_CONSTANT_PDF<%s>"%T
        self.defs.append("class_<%s >(\"%s\")\n"%(cppclass,name))
        self.defs.append("    .def(init<int>())\n")
        for m in [('pdf','readonly'),('cdf','readonly')]:
            self.member_variable(cppclass,m)
        members=["Initialize","Compute_Cumulative_Distribution_Function","Sample"]
        for member in members:
            self.member_function(cppclass,member)
        self.defs.append(";\n")

    ######################################################################
    # Function cutting
    ######################################################################
    def cutting_2d(self,name,TV,d):
        self.include("PhysBAM_Geometry/Geometry/TRIANGULATED_AREA.h")
        self.include("PhysBAM_Geometry/Geometry/SEGMENTED_CURVE_2D.h")
        self.include("Fracture/CUTTING_GEOMETRY_2D.h")
        cppclass="CUTTING_GEOMETRY_2D<%s,%d>"%(TV,d)
        self.defs.append("class_<%s,boost::noncopyable>(\"%s\",init<bool>())\n"%(cppclass,name))
        for member in ["Initialize_Original_Embedding","Cut_Material"]:
            self.member_function(cppclass,member)
        self.reference(cppclass,"cutting")
        self.defs.append(";\n")

    def cutting_3d(self,name,TV,d):
        self.include("PhysBAM_Geometry/Geometry/TETRAHEDRALIZED_VOLUME.h")
        self.include("PhysBAM_Geometry/Geometry/TRIANGULATED_SURFACE.h")
        self.include("Fracture/CUTTING_GEOMETRY_3D.h")
        cppclass="CUTTING_GEOMETRY_3D<%s,%d>"%(TV,d)
        self.defs.append("class_<%s,boost::noncopyable>(\"%s\",init<bool>())\n"%(cppclass,name))
        for member in ["Initialize_Original_Embedding","Cut_Material"]:
            self.member_function(cppclass,member)
        self.reference(cppclass,"cutting")
        self.defs.append(";\n")

    ######################################################################
    # Function point_repulsion
    ######################################################################
    def point_repulsion(self,name,T):
        cppclass="POINT_REPULSION<%s>"%T
        self.include("Meshing/POINT_REPULSION.h")

        self.defs.append("class_<%s,boost::noncopyable>(\"%s\",init<TRIANGULATED_SURFACE<%s>&,int>())\n"%(cppclass,name,T))
        for member in ["Seed_Points","Move_Points","Get_Points","Get_Triangles"]:
            self.member_function(cppclass,member)
        for m in [('points','readonly')]:
            self.member_variable(cppclass,m)

        self.defs.append(";\n")

    ######################################################################
    # Function point_repulsion_data
    ######################################################################
    def point_repulsion_data(self,name,T):
        cppclass="POINT_REPULSION_DATA<%s>"%T
        class_to_file_map["POINT_REPULSION_DATA"]="Meshing/POINT_REPULSION.h"
        self.defs.append("class_<%s >(\"%s\")\n"%(cppclass,name))

        for member in ["position","triangle","neighbors"]:
            self.member_variable(cppclass,(member,'readonly'))

        self.defs.append(";\n")

    ######################################################################
    # Function body_motion_sequence
    ######################################################################
    def body_motion_sequence(self,name,T):
        cppclass="BODY_MOTION_SEQUENCE<%s>"%T
        self.include("Motion/BODY_MOTION_SEQUENCE.h")

        self.defs.append("class_<%s >(\"%s\")\n"%(cppclass,name))
        #for member in [("Initialize","void",["const int","const GRID_1D<%s>&"%T]),("Set_Frame_Rate","void",["const %s"%T,"const %s"%T])]:
        for member in ["Initialize","Set_Frame_Rate","Update_Name_Lookup","Set_Targeted","Set_Transform","Resize","Rescale"]:
            self.member_function(cppclass,member)
        for m in ['time_grid','trajectories','valid','names','base_position','bone_hierarchy','name_to_track_index','saved_frame']:
            self.member_variable(cppclass,m)

        self.defs.append(";\n")
        self.read_write(cppclass,T)

    ######################################################################
    # Function bone
    ######################################################################
    def bone(self,name,T):
        cppclass="BONE<%s>"%T
        self.include("Motion/BONE.h")
        self.include("Read_Write/Motion/READ_WRITE_BONE.h")

        self.defs.append("class_<%s >(\"%s\")\n"%(cppclass,name))
        for m in ['transform','length','targeted_transform']:
            self.member_variable(cppclass,m)

        self.defs.append(";\n")
        self.read_write(cppclass,T)

    ######################################################################
    # Function uniform_grid
    ######################################################################
    def uniform_grid(self,name,T,d):
        self.include("PhysBAM_Tools/Grids_Uniform/GRID_%dD.h"%d)
        cppclass="GRID_%dD<%s >"%(d,T)
        TV="VECTOR<%s,%d>"%(T,d)
        self.defs.append('class_<%s >("%s")\n'%(cppclass,name))
        for MAC in ['',',const bool']:
            self.defs.append('    .def(init<%s%s>())\n'%(','.join(["const int"]*d+["const %s"%T]*2*d),MAC))
            self.defs.append('    .def(init<%s,const RANGE<%s >&%s>())\n'%(','.join(["const int"]*d),TV,MAC))
            self.defs.append('    .def(init<const VECTOR<int,%d>&,const RANGE<%s >&%s>())\n'%(d,TV,MAC))
        for i in ["x","y","z"][:d]:
            for m in [i+'min',i+'max']:
                self.defs.append("    .def_readwrite(\"%s\",&%s::%s)\n"%(m,cppclass,m))
        for member in ["DX","One_Over_DX","Numbers_Of_Cells","Numbers_Of_Nodes","Is_MAC_Grid","Domain",'Domain_Indices',
                       'Set_To_Double_Resolution_Grid','Cell_Diagonal','Cell_Size','Face','Block_Index','Counts','Domain_Indices',
                       'Get_MAC_Grid','Get_Regular_Grid','Get_Face_Grid',
                       'Get_Regular_Grid_At_MAC_Positions','Get_MAC_Grid_At_Regular_Positions','Remove_Dimension',
                       ('Outside','bool',["const %s&"%TV],'const'),
                       ("Clamp_To_Cell","VECTOR<int,%d>"%d,["const %s&"%TV],"const"),
                       ("Clamp_To_Cell","VECTOR<int,%d>"%d,["const %s&"%TV,'int'],"const"),
                       ('Closest_Node',"VECTOR<int,%d>"%d,["const %s&"%TV],"const")]+['Get_%s_Face_Grid'%a for a in 'XYZ'[:d]]:
            self.member_function(cppclass,member)
        for s in ['X','Center','Node']+['%s_Face'%a for a in 'XYZ'[:d]]:
            self.member_function(cppclass,(s,"VECTOR<%s,%d>"%(T,d),["const int"]*d,"const"))
            self.member_function(cppclass,(s,"VECTOR<%s,%d>"%(T,d),["const VECTOR<int,%d>&"%d],"const"))
        if d>1:
            for f in ['Get_Horizontal_Grid','Remove_Dimension']:
                self.member_function(cppclass,f)
        self.static_function(cppclass,'Create_Grid_Given_Cell_Size',overloads=(3,4))
        self.static_function(cppclass,'Create_Even_Sized_Grid_Given_Cell_Size',overloads=(3,4))
        self.str(cppclass)
        self.repr(cppclass)
        for member in ['Node_Indices','Cell_Indices','Block_Indices','Face_Indices']:
            self.member_function(cppclass,member,overloads=(0,1))
        self.defs.append(";\n")
        self.read_write(cppclass,T)

    ######################################################################
    # Function element_id
    ######################################################################
    def element_id(self,class_definition,name):
        #self.include("PhysBAM_Tools/Data_Structures/ELEMENT_ID.h")
        self.include(class_definition)
        self.defs.append('class_<%s >("%s")\n'%(name,name))
        self.defs.append('    .def(init<int>())\n')
        self.member_function(name,"Value")
        self.defs.append(";\n")

    ######################################################################
    # Function registration
    ######################################################################
    def registration(self,name,T):
        TV='VECTOR<%s,3>'%T
        cppclass="REGISTRATION<%s >"%(T)
        self.include("PhysBAM_Geometry/Geometry/REGISTRATION.h")
        self.include("PhysBAM_Tools/Vectors/VECTOR.h")
        self.include("PhysBAM_Tools/Matrices/MATRIX.h")
        self.include("PhysBAM_Tools/Data_Structures/PAIR.h")
        self.routines.append("""
template<class T> boost::python::object Affine_Registration_Wrapper(const ARRAY<VECTOR<T,3> >& X1,const ARRAY<VECTOR<T,3> >& X2)
{PAIR<MATRIX<T,3>,VECTOR<T,3> > result=REGISTRATION<T>::Affine_Registration(X1,X2);return boost::python::make_tuple(result.x,result.y);}
""")
        self.defs.append('class_<%s,boost::noncopyable>("%s",no_init)\n'%(cppclass,name))
        self.static_function(cppclass,('Rigid_Registration','FRAME<%s >'%TV,['const ARRAY<%s >&'%TV]*2+['%s*'%T]),overloads=(2,2))
        self.static_function(cppclass,('Rigid_Registration','FRAME<%s >'%TV,['const ARRAY<%s >&'%T]+['const ARRAY<%s >&'%TV]*2+['%s*'%T]),overloads=(3,3))
        self.static_function(cppclass,('Rigid_Momenta','TWIST<%s >'%TV,['const ARRAY<%s >&'%T]+['const ARRAY<%s >&'%TV]*2))
        self.defs.append('    .def("Affine_Registration",&Affine_Registration_Wrapper<%s >).staticmethod("Affine_Registration")\n'%T)
        self.defs.append(";\n")

    ######################################################################
    # Function union_find
    ######################################################################
    def union_find(self,name,ID):
        cppclass="UNION_FIND<%s >"%(ID)
        self.include("PhysBAM_Tools/Data_Structures/UNION_FIND.h")
        self.include("Deformable_Objects/FRAGMENT_ID.h")
        self.defs.append('class_<%s,boost::noncopyable>("%s")\n'%(cppclass,name))
        self.defs.append('    .def(init<%s >())\n'%ID)
        for guy in ['Initialize','Size','Clear_Connectivity','Is_Root','Find','Merge',('Union',ID,[ID,ID])]:
            self.member_function(cppclass,guy)
        self.defs.append(";\n")

    ######################################################################
    # Function union_find
    ######################################################################
    def sparse_union_find(self,name,ID):
        cppclass="SPARSE_UNION_FIND<%s >"%(ID)
        self.include("PhysBAM_Tools/Data_Structures/SPARSE_UNION_FIND.h")
        self.include("PhysBAM_Tools/Data_Structures/DATA_STRUCTURES_FORWARD.h")
        self.defs.append('class_<%s,boost::noncopyable>("%s")\n'%(cppclass,name))
        self.defs.append('    .def(init<%s >())\n'%ID)
        for guy in ['Initialize','Size','Clear_Connectivity','Is_Root','Find','Merge',('Union',ID,[ID,ID])]:
            self.member_function(cppclass,guy)
        self.defs.append(";\n")

    ######################################################################
    # Function particle_connectivity
    ######################################################################
    def particle_connectivity(self,name,TV):
        cppclass="PARTICLE_CONNECTIVITY<%s >"%(TV)
        self.include("Deformable_Objects/PARTICLE_CONNECTIVITY.h")
        self.defs.append("class_<%s,boost::noncopyable>(\"%s\",no_init)\n"%(cppclass,name))
        for guy in [("Union","int",["const int","const int"])]:
            self.member_function(cppclass,guy)
        self.defs.append(";\n")

    ######################################################################
    # Function file_utilities
    ######################################################################
    def file_utilities(self):
        self.include("PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h")
        self.routines.append("""\
static STREAM_TYPE stream_type_float=STREAM_TYPE(float());
static STREAM_TYPE stream_type_double=STREAM_TYPE(double());
""")
        self.defs.append("""\
class_<STREAM_TYPE>("STREAM_TYPE",no_init);
scope().attr("stream_type_float")=stream_type_float;
scope().attr("stream_type_double")=stream_type_double;
def("Safe_Open_Output",FILE_UTILITIES::Safe_Open_Output);
def("Safe_Open_Input",FILE_UTILITIES::Safe_Open_Input);
class_<std::istream,boost::noncopyable>("istream",no_init);
class_<std::ostream,boost::noncopyable>("ostream",no_init);
register_ptr_to_python<std::auto_ptr<std::istream> >();
register_ptr_to_python<std::auto_ptr<std::ostream> >();
class_<TYPED_ISTREAM,boost::noncopyable>("TYPED_ISTREAM",init<std::istream&,const STREAM_TYPE>()[with_custodian_and_ward<1,2>()],return_value_policy<reference_existing_object>());
class_<TYPED_OSTREAM,boost::noncopyable>("TYPED_OSTREAM",init<std::ostream&,const STREAM_TYPE>()[with_custodian_and_ward<1,2>()],return_value_policy<reference_existing_object>());
""")

    ######################################################################
    # Function permutations
    ######################################################################
    def permutations(self,T):
        self.include("PhysBAM_Tools/Math_Tools/permutation.h")
        for i in ["permute_two","permute_two_inverse","permute_three","permute_three_inverse","permute_four","permute_four_inverse"]:
            self.defs.append("def(\"%s\",%s<%s >);\n"%(i,i,T))

    ######################################################################
    # Function log
    ######################################################################
    def log(self):
        self.include("PhysBAM_Tools/Log/LOG.h")
        self.defs.append("""
    def("Initialize_Logging",LOG::Initialize_Logging);
    def("Finish_Logging",LOG::Finish_Logging);
    def("Stop_Time",LOG::Stop_Time);
    def("Time",static_cast<void(*)(const std::string&)>(LOG::Time));
    def("Stat",LOG::Stat<float>);
    def("Stat",LOG::Stat<int>);
    def("Stat",LOG::Stat<double>);

class_<LOG::SCOPE,boost::noncopyable>("LOG_SCOPE",init<const std::string>());
""")

        self.defs.append(';\n')

    ######################################################################
    # Function progress_indicator
    ######################################################################
    def progress_indicator(self):
        self.include("PhysBAM_Tools/Log/PROGRESS_INDICATOR.h")
        self.routines.append("BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Progress_overloads,Progress,0,1)\n")
        self.defs.append("""class_<PROGRESS_INDICATOR>("PROGRESS_INDICATOR")
    .def(init<int>())
    .def("Initialize",&PROGRESS_INDICATOR::Initialize)
    .def("Progress",&PROGRESS_INDICATOR::Progress,Progress_overloads())
;
""")

    ######################################################################
    # Function image
    ######################################################################
    def image(self,name,T):
        self.include("PhysBAM_Tools/Images/IMAGE.h")
        self.defs.append('{\n')
        self.typedef(T=T)
        cppclass='IMAGE<T>'
        self.routines.append("BOOST_PYTHON_FUNCTION_OVERLOADS(Write_overloads,IMAGE<%s>::Write<3>,2,4)\n"%T)
        self.defs.append("class_<%s >(\"%s\")\n"%(cppclass,name))
        for f in ['Is_Supported']:
            self.static_function(cppclass,f)
        self.defs.append('    .def("Read",&IMAGE<T>::Write<3>).staticmethod("Read")\n')
        self.defs.append('    .def("Write",&IMAGE<T>::Write<3>,Write_overloads())\n      .staticmethod("Write")\n')
        self.defs.append(';}\n')

    ######################################################################
    # Function misc
    ######################################################################
    def misc(self):
        self.local_include("TEST.h")
        self.include("Particles/DEFORMABLE_BODY_PARTICLES.h")
        self.include("PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h")
        self.include("PhysBAM_Tools/Log/DEBUG_UTILITIES.h")
        self.defs.append('def("Initialize_Incident_Old",Initialize_Incident_Old);\n')
        self.defs.append('def("Initialize_Incident_New",Initialize_Incident_New);\n')
        self.defs.append('def("Initialize_Incident_Hash",Initialize_Incident_Hash);\n')
        self.defs.append('def("Scale_Restlength",Scale_Restlength<VECTOR<float,3> >);\n')
        self.defs.append('def("Project",Project<VECTOR<float,3> >);\n')
        self.defs.append("""\
def("Set_Floating_Point_Exception",PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling);
def("Backtrace",PROCESS_UTILITIES::Backtrace);
def("Set_Backtrace",PROCESS_UTILITIES::Set_Backtrace);
""")



codes=[]

######################################################################
# INSTANTIATIONS
######################################################################

register_class_include("PhysBAM_Geometry/Geometry/BOX.h")
register_class_include("PhysBAM_Tools/Math_Tools/RANGE.h")
register_class_include("PhysBAM_Tools/Matrices/FRAME.h")
register_class_include("PhysBAM_Tools/Vectors/TWIST.h")
register_class_include("Rigid_Bodies/RIGID_BODY_MASS.h")
register_class_include("Particles/DEFORMABLE_BODY_PARTICLES.h")
register_class_include("Particles/PARTICLE_LEVELSET_PARTICLES.h")
register_class_include("Particles/PARTICLE_LEVELSET_REMOVED_PARTICLES.h")
register_class_include("Particles/RIGID_BODY_PARTICLES.h")

misc_code=CODE("misc")
misc_code.misc()
misc_code.progress_indicator()
misc_code.image('IMAGE_f','float')
for T in ['float','double']:
    if not use_double and T=='double': continue
    misc_code.numeric_limits_instance('numeric_limits_'+T[0],T)
misc_code.registration('REGISTRATION_f','float')
misc_code.cloneable_base()
codes.append(misc_code)

# vectors
for T in ["float","double","bool","int"]:
    code=CODE('vector_'+T)
    codes.append(code)
    if not use_double and T=='double': continue
    for d in [1,2,3,4]:
        code.vector_instance("V%s%d"%(T[0],d),T,d)
for T in ["float"]:
    code=CODE('vector_vector')
    codes.append(code)
    for d1 in [2,3,4]:
        for d2 in [1,2,3,4]:
            code.vector_instance("V_V%s%d_%d"%(T[0],d1,d2),"VECTOR<%s,%d>"%(T,d1),d2)

# matrices
matrix_templates=([('MATRIX_MXN',-1,-1)]
    +[(M,d,d) for M in ['SYMMETRIC_MATRIX','DIAGONAL_MATRIX','UPPER_TRIANGULAR_MATRIX'] for d in [2,3]]
    +[('MATRIX',m,n) for m in [2,3] for n in [2,3] if m!=n]
    +[('MATRIX',m,m) for m in [0,1,2,3]])
for T in ['float','double']:
    matrix_code=CODE('matrices_'+T)
    codes.append(matrix_code)
    if not use_double and T=='double': continue
    matrix_code.matrix_instances(T,matrix_templates)

# arrays
array_code=CODE("arrays")
array_code.array_instance("A_Vf1","VECTOR<float,1>")
array_code.array_instance("A_Vf2","VECTOR<float,2>")
array_code.array_instance("A_Vf3","VECTOR<float,3>")
array_code.array_instance("A_BOX_Vf2","BOX<VECTOR<float,2> >")
array_code.array_instance("A_BOX_Vf3","BOX<VECTOR<float,3> >")
array_code.array_instance("A_FRAME_Vf2","FRAME<VECTOR<float,2> >")
array_code.array_instance("A_FRAME_Vf3","FRAME<VECTOR<float,3> >")
array_code.array_instance("A_TWIST_Vf2","TWIST<VECTOR<float,2> >")
array_code.array_instance("A_TWIST_Vf3","TWIST<VECTOR<float,3> >")
array_code.array_instance("A_RIGID_BODY_MASS_Vf2","RIGID_BODY_MASS<VECTOR<float,2> >")
array_code.array_instance("A_RIGID_BODY_MASS_Vf3","RIGID_BODY_MASS<VECTOR<float,3> >")
array_code.array_instance("A_f","float")
array_code.array_instance("A_i","int")
array_code.array_instance("A_b","bool")
array_code.array_instance("A_ushort","unsigned short")
array_code.array_instance("A_string","std::string")
array_code.array_instance("A_A_i","ARRAY<int>")
codes.append(array_code)
# list arrays
list_array_code=CODE("list_arrays")
list_array_code.list_array_instance("LA_bool","bool")
list_array_code.list_array_instance("LA_i","int",True)
list_array_code.list_array_instance("LA_string","std::string")
list_array_code.list_array_instance("LA_f","float")
list_array_code.list_array_instance("LA_Vf3","VECTOR<float,3>")
list_array_code.list_array_instance("LA_Vi2","VECTOR<int,2>")
list_array_code.list_array_instance("LA_Vi3","VECTOR<int,3>")
list_array_code.list_array_instance("LA_Vi4","VECTOR<int,4>")
list_array_code.list_array_instance("LA_LA_i","LIST_ARRAY<int>")
list_array_code.list_array_instance("LA_LA_f","LIST_ARRAY<float>")
list_array_code.array_instance("A_LA_string","LIST_ARRAY<std::string>")
codes.append(list_array_code)
# multidimensional arrays
for T in ['float','double']:
    arrays_nd_code=CODE('arrays_nd_'+T)
    codes.append(arrays_nd_code)
    if not use_double and T=='double': continue
    arrays_nd_code.arrays_nd_instance("ARRAYS_1D_%s"%T[0],T,T,1)
    arrays_nd_code.arrays_nd_instance("ARRAYS_2D_%s"%T[0],T,T,2)
    arrays_nd_code.arrays_nd_instance("ARRAYS_2D_V%s3"%T[0],"VECTOR<%s,3>"%T,T,2)
    arrays_nd_code.arrays_nd_instance("ARRAYS_3D_%s"%T[0],T,T,3)
arrays_nd_code=CODE('arrays_nd_bool')
arrays_nd_code.arrays_nd_instance("ARRAYS_1D_b","bool","bool",1)
arrays_nd_code.array_instance("A_ARRAYS_1D_b","ARRAYS_1D<bool> ")
codes.append(arrays_nd_code)

element_id_code=CODE("element_id")
element_id_code.element_id("Articulated_Rigid_Bodies/JOINT_ID.h","JOINT_ID")
element_id_code.element_id("PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_ID.h","COLLISION_GEOMETRY_ID")
element_id_code.element_id("Deformable_Objects/FRAGMENT_ID.h","FRAGMENT_ID")
element_id_code.element_id("PhysBAM_Tools/Data_Structures/ELEMENT_ID.h","INITIAL_SIZE")
element_id_code.element_id("Deformable_Objects/FRAGMENT_ID.h","SUPER_FRAGMENT_ID")
element_id_code.element_id("Parallel_Computation/PARTITION_ID.h","PARTITION_ID")
element_id_code.element_id("Deformable_Objects/HAIR_ID.h","HAIR_ID")
codes.append(element_id_code)

solid_body_collection_code=CODE("solid_body_collection")
solid_body_collection_code.solid_body_collection("SOLID_BODY_COLLECTION_Vf2","VECTOR<float,2>")
solid_body_collection_code.solid_body_collection("SOLID_BODY_COLLECTION_Vf3","VECTOR<float,3>")
solid_body_collection_code.soft_bindings("SOFT_BINDING_Vf2","VECTOR<float,2>")
solid_body_collection_code.soft_bindings("SOFT_BINDING_Vf3","VECTOR<float,3>")
solid_body_collection_code.binding("BINDING_Vf2_3","float","VECTOR<float,2>",3)
solid_body_collection_code.binding("BINDING_Vf3_3","float","VECTOR<float,3>",3)
solid_body_collection_code.linear_binding("LINEAR_BINDING_Vf3_3","float","VECTOR<float,3>",3)
solid_body_collection_code.linear_binding("LINEAR_BINDING_Vf3_4","float","VECTOR<float,3>",4)
solid_body_collection_code.binding_list("BINDING_LIST_Vf2","float","VECTOR<float,2>")
solid_body_collection_code.binding_list("BINDING_LIST_Vf3","float","VECTOR<float,3>")
solid_body_collection_code.super_fragment("SUPER_FRAGMENT_Vf2","VECTOR<float,2>")
solid_body_collection_code.super_fragment("SUPER_FRAGMENT_Vf3","VECTOR<float,3>")
solid_body_collection_code.pair_instance("PAIR_SUPER_FRAGMENT_ID_SUPER_FRAGMENT_ID","SUPER_FRAGMENT_ID","SUPER_FRAGMENT_ID")
solid_body_collection_code.list_array_instance("LA_SUPER_FRAGMENT_ID","SUPER_FRAGMENT_ID")
solid_body_collection_code.list_array_instance("LA_PAIR_SUPER_FRAGMENT_ID_SUPER_FRAGMENT_ID","PAIR<SUPER_FRAGMENT_ID,SUPER_FRAGMENT_ID>")
codes.append(solid_body_collection_code)

deformable_body_collection_code=CODE("deformable_body_collection")
deformable_body_collection_code.deformable_body_collection("DEFORMABLE_BODY_COLLECTION_Vf2","VECTOR<float,2>")
deformable_body_collection_code.deformable_body_collection("DEFORMABLE_BODY_COLLECTION_Vf3","VECTOR<float,3>")
codes.append(deformable_body_collection_code)

deformable_geometry_collection_code=CODE("deformable_geometry_collection")
deformable_geometry_collection_code.deformable_geometry_collection("DEFORMABLE_GEOMETRY_COLLECTION_Vf2","VECTOR<float,2>")
deformable_geometry_collection_code.deformable_geometry_collection("DEFORMABLE_GEOMETRY_COLLECTION_Vf3","VECTOR<float,3>")
codes.append(deformable_geometry_collection_code)

forces_code=CODE("forces")
forces_code.solids_forces("SOLIDS_FORCES_Vf2","VECTOR<float,2>")
forces_code.solids_forces("SOLIDS_FORCES_Vf3","VECTOR<float,3>")
forces_code.gravity("GRAVITY_Vf2","float","VECTOR<float,2>")
forces_code.gravity("GRAVITY_Vf3","float","VECTOR<float,3>")
forces_code.linear_springs("LINEAR_SPRINGS_Vf2","float","VECTOR<float,2>")
forces_code.linear_springs("LINEAR_SPRINGS_Vf3","float","VECTOR<float,3>")
forces_code.segment_adhesion("SEGMENT_ADHESION_Vf3","float","VECTOR<float,3>")
forces_code.linear_tet_springs("LINEAR_TET_SPRINGS_f","float")
forces_code.linear_altitude_springs_3d("LINEAR_ALTITUDE_SPRINGS_3D_f","float")
forces_code.triangle_bending_springs("TRIANGLE_BENDING_SPRINGS_f","float")
forces_code.triangle_bending_elements("TRIANGLE_BENDING_ELEMENTS_f","float")
forces_code.segment_bending_springs("SEGMENT_BENDING_SPRINGS_f","VECTOR<float,3>")
forces_code.finite_volume("FINITE_VOLUME","float","VECTOR<float,3>",3)
forces_code.binding_springs("BINDING_SPRINGS_Vf3","float","VECTOR<float,3>")
forces_code.constitutive_model("CONSTITUTIVE_MODEL_f3","float",3)
if use_double: forces_code.constitutive_model("CONSTITUTIVE_MODEL_f3","double",3)
forces_code.neo_hookean("NEO_HOOKEAN_f3","float",3)
codes.append(forces_code)

union_find=CODE("union_find")
union_find.union_find("UNION_FIND_int","int")
union_find.union_find("UNION_FIND_FRAGMENT_ID","FRAGMENT_ID")
union_find.sparse_union_find("SPARSE_UNION_FIND_int","int")
union_find.sparse_union_find("SPARSE_UNION_FIND_FRAGMENT_ID","FRAGMENT_ID")
union_find.particle_connectivity("PARTICLE_CONNECTIVITY_Vf2","VECTOR<float,2>")
union_find.particle_connectivity("PARTICLE_CONNECTIVITY_Vf3","VECTOR<float,3>")
codes.append(union_find)

geometry_code=CODE("geometry")
for T in ['float','int']:
    for d in [1,2,3]:
        geometry_code.range('RANGE_V%s%d'%(T[0],d),'VECTOR<%s,%d>'%(T,d))
for T in ['float']:
    for d in [1,2,3]:
        geometry_code.box('BOX_V%s%d'%(T[0],d),'VECTOR<%s,%d>'%(T,d))
geometry_code.tetrahedron("TETRAHEDRON_f","float")
geometry_code.segment_3d("SEGMENT_3D_f","float")
geometry_code.triangle_3d("TRIANGLE_3D_f","float")
geometry_code.sphere("SPHERE_Vf2","VECTOR<float,2>")
geometry_code.sphere("SPHERE_Vf3","VECTOR<float,3>")
geometry_code.plane("PLANE_f","VECTOR<float,3>")
geometry_code.ring("RING_f","float")
geometry_code.cylinder("CYLINDER_f","float")
codes.append(geometry_code)

triangulated_surface_code=CODE("triangulated_surface")
triangulated_surface_code.structure("STRUCTURE_Vf2","VECTOR<float,2>")
triangulated_surface_code.structure("STRUCTURE_Vf3","VECTOR<float,3>")
triangulated_surface_code.segmented_curve_instance("SEGMENTED_CURVE_Vf2","VECTOR<float,2>")
triangulated_surface_code.segmented_curve_instance("SEGMENTED_CURVE_Vf3","VECTOR<float,3>")
triangulated_surface_code.triangulated_surface_instance("TRIANGULATED_SURFACE_f","float")
triangulated_surface_code.triangulated_area_instance("TRIANGULATED_AREA_f","float")
triangulated_surface_code.tetrahedralized_volume_instance("TETRAHEDRALIZED_VOLUME_f","float")
triangulated_surface_code.segment_mesh("SEGMENT_MESH")
triangulated_surface_code.triangle_mesh("TRIANGLE_MESH")
triangulated_surface_code.tetrahedron_mesh("TETRAHEDRON_MESH")
codes.append(triangulated_surface_code)

triangulated_surface_code=CODE("hierarchies")
triangulated_surface_code.box_hierarchy("BOX_HIERARCHY_Vf3","VECTOR<float,3>")
triangulated_surface_code.particle_hierarchy("PARTICLE_HIERARCHY_Vf3","VECTOR<float,3>")
triangulated_surface_code.tetrahedron_hierarchy("TETRAHEDRON_HIERARCHY_f","float")
codes.append(triangulated_surface_code)

particles_code=CODE("particles")
# base attributes
particles_code.particle_attribute_base()
particles_code.particle_attribute('PARTICLE_ATTRIBUTE_f','float')
particles_code.particle_attribute('PARTICLE_ATTRIBUTE_int','int')
particles_code.particle_attribute('PARTICLE_ATTRIBUTE_ushort','unsigned short')
particles_code.particle_attribute('PARTICLE_ATTRIBUTE_Vf1','VECTOR<float,1>')
particles_code.particle_attribute('PARTICLE_ATTRIBUTE_Vf2','VECTOR<float,2>')
particles_code.particle_attribute('PARTICLE_ATTRIBUTE_Vf3','VECTOR<float,3>')
particles_code.particle_attribute('PARTICLE_ATTRIBUTE_FRAME_Vf2','FRAME<VECTOR<float,2> >')
particles_code.particle_attribute('PARTICLE_ATTRIBUTE_FRAME_Vf3','FRAME<VECTOR<float,3> >')
particles_code.particle_attribute('PARTICLE_ATTRIBUTE_TWIST_Vf2','TWIST<VECTOR<float,2> >')
particles_code.particle_attribute('PARTICLE_ATTRIBUTE_TWIST_Vf3','TWIST<VECTOR<float,3> >')
particles_code.particle_attribute('PARTICLE_ATTRIBUTE_RIGID_BODY_MASS_Vf2','RIGID_BODY_MASS<VECTOR<float,2> >')
particles_code.particle_attribute('PARTICLE_ATTRIBUTE_RIGID_BODY_MASS_Vf3','RIGID_BODY_MASS<VECTOR<float,3> >')
# named attributes
particles_code.particle_mass_attribute("PARTICLE_MASS_ATTRIBUTE_f","float")
particles_code.particle_simple_attribute('PARTICLE_POSITION_ATTRIBUTE_Vf2','PARTICLE_POSITION_ATTRIBUTE<VECTOR<float,2> >')
particles_code.particle_simple_attribute('PARTICLE_POSITION_ATTRIBUTE_Vf3','PARTICLE_POSITION_ATTRIBUTE<VECTOR<float,3> >')
particles_code.particle_simple_attribute('PARTICLE_VELOCITY_ATTRIBUTE_Vf2','PARTICLE_VELOCITY_ATTRIBUTE<VECTOR<float,2> >')
particles_code.particle_simple_attribute('PARTICLE_VELOCITY_ATTRIBUTE_Vf3','PARTICLE_VELOCITY_ATTRIBUTE<VECTOR<float,3> >')
particles_code.particle_simple_attribute('PARTICLE_ENERGY_ATTRIBUTE_f','PARTICLE_ENERGY_ATTRIBUTE<float>')
particles_code.particle_simple_attribute('PARTICLE_RADIUS_ATTRIBUTE_f','PARTICLE_RADIUS_ATTRIBUTE<float>')
particles_code.particle_simple_attribute('PARTICLE_DENSITY_ATTRIBUTE_f','PARTICLE_DENSITY_ATTRIBUTE<float>')
particles_code.particle_simple_attribute('PARTICLE_TEMPERATURE_ATTRIBUTE_f','PARTICLE_TEMPERATURE_ATTRIBUTE<float>')
particles_code.particle_simple_attribute('PARTICLE_ID_ATTRIBUTE_int','PARTICLE_ID_ATTRIBUTE<int>')
particles_code.particle_simple_attribute('PARTICLE_VORTICITY_ATTRIBUTE_Vf2','PARTICLE_VORTICITY_ATTRIBUTE<VECTOR<float,2> >')
particles_code.particle_simple_attribute('PARTICLE_VORTICITY_ATTRIBUTE_Vf3','PARTICLE_VORTICITY_ATTRIBUTE<VECTOR<float,3> >')
particles_code.particle_simple_attribute('PARTICLE_ANGULAR_MOMENTUM_ATTRIBUTE_Vf2','PARTICLE_ANGULAR_MOMENTUM_ATTRIBUTE<VECTOR<float,2> >')
particles_code.particle_simple_attribute('PARTICLE_ANGULAR_MOMENTUM_ATTRIBUTE_Vf3','PARTICLE_ANGULAR_MOMENTUM_ATTRIBUTE<VECTOR<float,3> >')
particles_code.particle_simple_attribute('PARTICLE_QUANTIZED_COLLISION_DISTANCE_ATTRIBUTE','PARTICLE_QUANTIZED_COLLISION_DISTANCE_ATTRIBUTE')
particles_code.particle_simple_attribute('PARTICLE_STORED_PHI_ATTRIBUTE_f','PARTICLE_STORED_PHI_ATTRIBUTE<float>')
particles_code.particle_simple_attribute('PARTICLE_MATERIAL_VOLUME_ATTRIBUTE_f','PARTICLE_MATERIAL_VOLUME_ATTRIBUTE<float>')
particles_code.particle_simple_attribute('PARTICLE_FRAME_ATTRIBUTE_Vf2','PARTICLE_FRAME_ATTRIBUTE<VECTOR<float,2> >')
particles_code.particle_simple_attribute('PARTICLE_FRAME_ATTRIBUTE_Vf3','PARTICLE_FRAME_ATTRIBUTE<VECTOR<float,3> >')
particles_code.particle_simple_attribute('PARTICLE_TWIST_ATTRIBUTE_Vf2','PARTICLE_TWIST_ATTRIBUTE<VECTOR<float,2> >')
particles_code.particle_simple_attribute('PARTICLE_TWIST_ATTRIBUTE_Vf3','PARTICLE_TWIST_ATTRIBUTE<VECTOR<float,3> >')
particles_code.particle_simple_attribute('PARTICLE_RIGID_BODY_MASS_ATTRIBUTE_Vf2','PARTICLE_RIGID_BODY_MASS_ATTRIBUTE<VECTOR<float,2> >')
particles_code.particle_simple_attribute('PARTICLE_RIGID_BODY_MASS_ATTRIBUTE_Vf3','PARTICLE_RIGID_BODY_MASS_ATTRIBUTE<VECTOR<float,3> >')
# particles
particles_code.particle_base()
particles_code.single_particles()
particles_code.particles("PARTICLE_Vf2","VECTOR<float,2>")
particles_code.particles("PARTICLE_Vf3","VECTOR<float,3>")
particles_code.simple_particles("DEFORMABLE_BODY_PARTICLES_Vf2","DEFORMABLE_BODY_PARTICLES<VECTOR<float,2> >")
particles_code.simple_particles("DEFORMABLE_BODY_PARTICLES_Vf3","DEFORMABLE_BODY_PARTICLES<VECTOR<float,3> >")
particles_code.simple_particles("PARTICLE_LEVELSET_PARTICLES_Vf2","PARTICLE_LEVELSET_PARTICLES<VECTOR<float,2> >")
particles_code.simple_particles("PARTICLE_LEVELSET_PARTICLES_Vf3","PARTICLE_LEVELSET_PARTICLES<VECTOR<float,3> >")
particles_code.simple_particles("PARTICLE_LEVELSET_REMOVED_PARTICLES_Vf2","PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<float,2> >")
particles_code.simple_particles("PARTICLE_LEVELSET_REMOVED_PARTICLES_Vf3","PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<float,3> >")
particles_code.simple_particles("RIGID_BODY_PARTICLES_Vf2","RIGID_BODY_PARTICLES<VECTOR<float,2> >")
particles_code.simple_particles("RIGID_BODY_PARTICLES_Vf3","RIGID_BODY_PARTICLES<VECTOR<float,3> >")
codes.append(particles_code)

levelset_code=CODE("levelset")
levelset_code.levelset_instance("LEVELSET_UNIFORM_Vf1","VECTOR<float,1>")
levelset_code.levelset_instance("LEVELSET_UNIFORM_Vf2","VECTOR<float,2>")
levelset_code.levelset_instance("LEVELSET_UNIFORM_Vf3","VECTOR<float,3>")
levelset_code.levelset_implicit_object("LEVELSET_IMPLICIT_OBJECT_Vf3","VECTOR<float,3>")
codes.append(levelset_code)

heat_code=CODE("heat")
heat_code.heat_uniform("HEAT_UNIFORM_Vf1","VECTOR<float,1>")
heat_code.heat_uniform("HEAT_UNIFORM_Vf2","VECTOR<float,2>")
heat_code.heat_uniform("HEAT_UNIFORM_Vf3","VECTOR<float,3>")
codes.append(heat_code)

dualcontour_code=CODE("dualcontour")
dualcontour_code.dualcontour_instance("DUALCONTOUR_UNIFORM_Vf2","GRID_2D<float>")
dualcontour_code.dualcontour_instance("DUALCONTOUR_UNIFORM_Vf3","GRID_3D<float>")
codes.append(dualcontour_code)

solids_fluid_example_code=CODE("solids_fluid_example")
solids_fluid_example_code.solids_fluids_example_uniform("SOLIDS_FLUIDS_EXAMPLE_UNIFORM_f2","float","GRID_2D<float>")
solids_fluid_example_code.solids_fluids_example_uniform("SOLIDS_FLUIDS_EXAMPLE_UNIFORM_f3","float","GRID_3D<float>")
solids_fluid_example_code.solids_fluids_driver_uniform("SOLIDS_FLUIDS_DRIVER_UNIFORM_f2","GRID_2D<float>")
solids_fluid_example_code.solids_fluids_driver_uniform("SOLIDS_FLUIDS_DRIVER_UNIFORM_f3","GRID_3D<float>")
solids_fluid_example_code.solids_parameters("SOLIDS_PARAMETERS_Vf2","float",2)
solids_fluid_example_code.solids_parameters("SOLIDS_PARAMETERS_Vf3","float",3)
solids_fluid_example_code.triangle_repulsions("TRIANGLE_REPULSIONS_Vf2","VECTOR<float,2>")
solids_fluid_example_code.triangle_repulsions("TRIANGLE_REPULSIONS_Vf3","VECTOR<float,3>")
codes.append(solids_fluid_example_code)

solids_standard_tests_code=CODE("solids_standard_tests")
solids_standard_tests_code.solids_standard_tests("SOLIDS_STANDARD_TESTS_Vf2","VECTOR<float,2> ")
solids_standard_tests_code.solids_standard_tests("SOLIDS_STANDARD_TESTS_Vf3","VECTOR<float,3> ")
codes.append(solids_standard_tests_code)

data_structures_code=CODE("data_structures")
data_structures_code.hashtable_instance("HASHTABLE_i_i","int","int")
#data_structures_code.hashtable_instance("HASHTABLE_Vf4","VECTOR<int,4>","void")
data_structures_code.pair_instance("PAIR_i_f","int","float")
codes.append(data_structures_code)

rigid_bodies_code=CODE("rigid_bodies")
for T in ['float','double']:
    if not use_double and T=='double': continue
    for d in [1,2,3]:
        rigid_bodies_code.rotation('ROTATION_V%s%d'%(T[0],d),'VECTOR<%s,%d>'%(T,d))
rigid_bodies_code.twist("TWIST_Vf2","VECTOR<float,2>")
rigid_bodies_code.twist("TWIST_Vf3","VECTOR<float,3>")
rigid_bodies_code.rigid_body("RIGID_BODY_Vf3_BASE","VECTOR<float,3>")
rigid_bodies_code.rigid_body("RIGID_BODY_Vf2_BASE","VECTOR<float,2>")
rigid_bodies_code.rigid_body_2d("RIGID_BODY_Vf2","float")
rigid_bodies_code.rigid_body_3d("RIGID_BODY_Vf3","float")
rigid_bodies_code.frame("FRAME_Vf2","VECTOR<float,2>")
rigid_bodies_code.frame("FRAME_Vf3","VECTOR<float,3>")
rigid_bodies_code.rigid_body_state("RIGID_BODY_STATE_Vf2","VECTOR<float,2>")
rigid_bodies_code.rigid_body_state("RIGID_BODY_STATE_Vf3","VECTOR<float,3>")
rigid_bodies_code.list_array_instance("LA_FRAME_Vf3","FRAME<VECTOR<float,3> >")
rigid_bodies_code.list_array_instance("LA_LA_FRAME_Vf3","LIST_ARRAY<FRAME<VECTOR<float,3> > >")
rigid_bodies_code.list_array_instance("LA_HAIR_ID","HAIR_ID")

codes.append(rigid_bodies_code)

simplex_interactions_code=CODE("simplex_interactions")
simplex_interactions_code.simplex_interactions("SIMPLEX_INTERACTIONS_f","float")
simplex_interactions_code.robust_simplex_interactions("ROBUST_SIMPLEX_INTERACTIONS_Vf2","VECTOR<float,2>")
codes.append(simplex_interactions_code)

for T in ['float','double']:
    grids_code=CODE('grids_'+T)
    codes.append(grids_code)
    if not use_double and T=='double': continue
    grids_code.uniform_grid("GRID_V%s1"%T[0],T,1)
    grids_code.uniform_grid("GRID_V%s2"%T[0],T,2)
    grids_code.uniform_grid("GRID_V%s3"%T[0],T,3)

random_numbers_code=CODE("random_numbers")
random_numbers_code.random_numbers("RANDOM_NUMBERS_f","float")
if use_double: random_numbers_code.random_numbers("RANDOM_NUMBERS_d","double")
random_numbers_code.piecewise_constant_pdf("PIECEWISE_CONSTANT_PDF_f","float")
for d in [1,2,3]:
    random_numbers_code.sobol('SOBOL_Vf%d'%d,'VECTOR<float,%d>'%d)
codes.append(random_numbers_code)

permutation_code=CODE("permutation")
permutation_code.permutations("int")
codes.append(permutation_code)

meshing_code=CODE("meshing")
meshing_code.point_repulsion("POINT_REPULSION_f","float")
meshing_code.point_repulsion_data("POINT_REPULSION_DATA_float","float")
#meshing_code.list_array_instance("LA_POINT_REPULSION_DATA","POINT_REPULSION_DATA<float>")
codes.append(meshing_code)

cutting_code=CODE("cutting")
cutting_code.cutting_2d("CUTTING_GEOMETRY_2D_Vf2_2","VECTOR<float,2>",2)
cutting_code.cutting_2d("CUTTING_GEOMETRY_2D_Vf3_2","VECTOR<float,3>",2)
cutting_code.cutting_3d("CUTTING_GEOMETRY_3D_Vf3_3","VECTOR<float,3>",3)
codes.append(cutting_code)

motion_code=CODE("motion")
motion_code.body_motion_sequence("BODY_MOTION_SEQUENCE_f","float")
motion_code.bone("BONE_f","float");
motion_code.array_instance("A_BONE_f","BONE<float>")
motion_code.arrays_nd_instance("ARRAYS_1D_BONE_f","BONE<float> ","float",1)
motion_code.array_instance("A_ARRAYS_1D_BONE_f","ARRAYS_1D<BONE<float> >")
codes.append(motion_code)

file_utilities_code=CODE("file_utilities")
file_utilities_code.file_utilities()
codes.append(file_utilities_code)

log_code=CODE("log")
log_code.log()
codes.append(log_code)

exceptions_code=CODE("exceptions")
exception_map={
    'READ_ERROR':'IOError',
    'FILESYSTEM_ERROR':'OSError',
    'LOOKUP_ERROR':'LookupError',
    'INDEX_ERROR':'IndexError',
    'KEY_ERROR':'KeyError',
    'TYPE_ERROR':'TypeError',
    'VALUE_ERROR':'ValueError',
    'NOT_IMPLEMENTED_ERROR':'NotImplementedError',
    'ASSERTION_ERROR':'AssertionError'}
exceptions_code.exceptions(exception_map)
codes.append(exceptions_code)

main_code=CODE("physbam",True)
main_code.includes.add('#define PHYSBAM_IMPORT_NUMPY\n#include "NUMPY.h"\n')
main_code.defs.append('PhysBAM::Import_Numpy();\n')
for i in codes:
    i.write()
    main_code.defs.append(i.call())
    main_code.routines.append(i.prototype())
main_code.write()
