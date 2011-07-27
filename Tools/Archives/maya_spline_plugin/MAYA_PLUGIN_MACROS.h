//+
// Copyright (C) Alias Systems,  a division  of  Silicon Graphics Limited and/or
// its licensors ("Alias").  All rights reserved.  These coded instructions,
// statements, computer programs, and/or related material (collectively, the
// "Material") contain unpublished information proprietary to Alias, which is
// protected by Canadian and US federal copyright law and by international
// treaties.  This Material may not be disclosed to third parties, or be copied
// or duplicated, in whole or in part, without the prior written consent of
// Alias.  ALIAS HEREBY DISCLAIMS ALL WARRANTIES RELATING TO THE MATERIAL,
// INCLUDING, WITHOUT LIMITATION, ANY AND ALL EXPRESS OR IMPLIED WARRANTIES OF
// NON-INFRINGEMENT, MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.
// IN NO EVENT SHALL ALIAS BE LIABLE FOR ANY DAMAGES WHATSOEVER, WHETHER DIRECT,
// INDIRECT, SPECIAL, OR PUNITIVE, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
// OR OTHER TORTIOUS ACTION, OR IN EQUITY, ARISING OUT OF OR RELATED TO THE
// ACCESS TO, USE OF, OR RELIANCE UPON THE MATERIAL.
//-

////////////////////////////////////////////////////////////////////////////////
//
// api_macros.h
//
// Description:
//    Convenience macros for error checking and attribute creation,
//
////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
//
// Error checking
//
//    MCHECKERROR       - check the status and print the given error message
//    MCHECKERRORNORET  - same as above but does not return
//
//////////////////////////////////////////////////////////////////////

#define MCHECKERROR(STAT,MSG)       \
    if ( MS::kSuccess != STAT ) {   \
        cerr << MSG << endl;        \
            return MS::kFailure;    \
    }

#define MCHECKERRORNORET(STAT,MSG)  \
    if ( MS::kSuccess != STAT ) {   \
        cerr << MSG << endl;        \
    }

//////////////////////////////////////////////////////////////////////
//
// Attribute creation
//
//       MAKE_TYPED_ATTR   - creates and adds a typed attribute
//       MAKE_NUMERIC_ATTR - creates and adds a numeric attribute
//       ADD_ATTRIBUTE     - adds the given attribute
//       ATTRIBUTE_AFFECTS - calls attributeAffects
//
//////////////////////////////////////////////////////////////////////

#define MAKE_TYPED_ATTR( NAME, LONGNAME, SHORTNAME, TYPE, DEFAULT )         \
                                                                            \
    MStatus NAME##_stat;                                                    \
    MFnTypedAttribute NAME##_fn;                                            \
    NAME = NAME##_fn.create( LONGNAME, SHORTNAME, TYPE, DEFAULT );          \
    NAME##_fn.setHidden( true );                                            \
    NAME##_stat = addAttribute( NAME );                                     \
    MCHECKERROR(NAME##_stat, "addAttribute error");

#define MAKE_NUMERIC_ATTR( NAME, LONGNAME, SHORTNAME, TYPE, DEFAULT,        \
                            ARRAY, BUILDER, KEYABLE )                       \
                                                                            \
    MStatus NAME##_stat;                                                    \
    MFnNumericAttribute NAME##_fn;                                          \
    NAME = NAME##_fn.create( LONGNAME, SHORTNAME, TYPE, DEFAULT );          \
    MCHECKERROR(NAME##_stat, "numeric attr create error");                    \
    NAME##_fn.setArray( ARRAY );                                            \
    NAME##_fn.setUsesArrayDataBuilder( BUILDER );                           \
    NAME##_fn.setHidden( ARRAY );                                            \
    NAME##_fn.setKeyable( KEYABLE );                                        \
    NAME##_stat = addAttribute( NAME );                                     \
    MCHECKERROR(NAME##_stat, "addAttribute error");

#define ADD_ATTRIBUTE( ATTR )                                               \
    MStatus ATTR##_stat;                                                    \
    ATTR##_stat = addAttribute( ATTR );                                     \
    MCHECKERROR( ATTR##_stat, "addAttribute: ATTR" )

#define ATTRIBUTE_AFFECTS( IN, OUT )                                        \
    MStatus IN##OUT##_stat;                                                 \
    IN##OUT##_stat = attributeAffects( IN, OUT );                           \
    MCHECKERROR(IN##OUT##_stat,"attributeAffects:" #IN "->" #OUT);

