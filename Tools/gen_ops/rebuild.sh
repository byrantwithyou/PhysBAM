#!/bin/bash

VF=../../Public_Library/Core/Vectors/SMALL_VECTOR_OPS.h
MF=../../Public_Library/Core/Matrices/SMALL_MATRIX_OPS.h

scons
./gen-code < matrix_list.txt > SMALL_MATRIX_OPS.h
./gen-code < vector_list.txt > SMALL_VECTOR_OPS.h

cat <<EOF > $VF
#ifndef __SMALL_VECTOR_OPS__
#define __SMALL_VECTOR_OPS__

#include <Core/Vectors/VECTOR.h>

namespace PhysBAM{
EOF

cat SMALL_VECTOR_OPS.h >> $VF

cat <<EOF >>  $VF
}

#endif
EOF



cat <<EOF > $MF
#ifndef __SMALL_MATRIX_OPS__
#define __SMALL_MATRIX_OPS__

#include <Core/Matrices/MATRIX.h>

namespace PhysBAM{
EOF

cat SMALL_MATRIX_OPS.h >> $MF

cat <<EOF >>  $MF
}

#endif
EOF



