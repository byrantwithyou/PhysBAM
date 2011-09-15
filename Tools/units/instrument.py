#!/usr/bin/python

import os
import re
import sys
import fileinput

files=sys.argv[1:]
invert_mapping=0
if '-v' in files:
    invert_mapping=1
    files.remove('-v')

ignore_file=re.compile(r'\b(CVS|build|Tools/units|External_Libraries|(READ_WRITE|READ_WRITE_FUNCTIONS|TYPED_STREAM|TYPE_UTILITIES|Hash|STATIC_ASSERT|BLOCK_DYADIC|TIMER|cbrt|min|max|clamp|LOG_ENTRY|ARRAY_DIFFERENCE|DEBUG_UTILITIES)\.(h|cpp))')
forward_file=re.compile(r'SCALAR_POLICY|FORWARD')

header_line=re.compile(r'#include|namespace PhysBAM{|#ifdef (?!__)')
units_include=re.compile(r'#include <Tools/units/\w+.h>')
include_quantity='#include <Tools/units/QUANTITY.h>\n'
forward_quantity='#include <Tools/units/QUANTITY_FORWARD.h>\n'
opengl_file=re.compile('OpenGL')
opengl_wrappers='#include <Tools/units/OPENGL_WRAPPERS.h>\n'

quantify=re.compile(r'\b(float\b|double\b|GLfloat\b|GLdouble\b|QUANTITY<(?:float|double|GLfloat|GLdouble)>)(.|$)')
unquantify=re.compile(r'\bQUANTITY<(float|double|GLfloat|GLdouble)>( >)?')
fix_casts=re.compile(r'(?<!\w)\((int|unsigned char|unsigned short|[A-Z]\w*)\)(\(|\w|-\d)')
unfix_casts=re.compile(r'\((?:typename )?CASTER<([\w ]+)>::TYPE\)')
template_context_guess=re.compile(r'T|[A-Z]*_MODE')

def instrument_number(m):
    s=m.group()
    if s[0]=='Q': return s
    q=m.expand(r'QUANTITY<\1>')
    next=m.group(2)
    if next=='>': q+=' '
    return q+next

def uninstrument_number(m):
    if m.group(2)==' >': return m.group(1)+'>'
    else: return m.group(1)

def fix_cast(m):
    prefix=''
    if template_context_guess.match(m.group(1)): prefix='typename '
    return '('+prefix+m.expand(r'CASTER<\1>::TYPE)\2')

def instrument_file(file):
    if not invert_mapping:
        header=0
        for line in fileinput.input([file],inplace=1):
            if not header and header_line.match(line):
                if not units_include.match(line):
                    if forward_file.search(file): print forward_quantity,
                    else: print include_quantity,
                    if opengl_file.search(file): print opengl_wrappers,
                header=1
            line=quantify.sub(instrument_number,line)
            line=fix_casts.sub(fix_cast,line)
            print line,
    else:
        for line in fileinput.input([file],inplace=1):
            if units_include.match(line): pass
            else:
                line=unquantify.sub(uninstrument_number,line)
                line=unfix_casts.sub(r'(\1)',line)
                print line,

def instrument_top(file):
    if ignore_file.search(os.path.abspath(file)): return
    if os.path.isdir(file):
        #print 'recursing into',file
        map(instrument_top,map(lambda x:os.path.join(file,x),os.listdir(file)))
    else:
        if file.endswith('.h') or file.endswith('.cpp'): instrument_file(file)

map(instrument_top,files)
