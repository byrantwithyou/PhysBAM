%{

#include <Tools/Symbolics/PROGRAM_DEFINITIONS.h>
#include <Tools/Symbolics/PROGRAM_YACC.hpp>

inline int yyerror(const char *msg) {
    fprintf(stderr,"Error:%s\n",msg); return 0;
}

%}

%option header-file="Public_Library/Tools/Symbolics/PROGRAM_LEX.hpp"
%option noyywrap

%%

"+=" { return TOKEN_ADDEQ; }
"-=" { return TOKEN_SUBEQ; }
"*=" { return TOKEN_MULEQ; }
"/=" { return TOKEN_DIVEQ; }
"^=" { return TOKEN_POWEQ; }
"==" { return TOKEN_EQ; }
"!=" { return TOKEN_NE; }
"<=" { return TOKEN_LE; }
">=" { return TOKEN_GE; }
"&&" { return TOKEN_AND; }
"||" { return TOKEN_OR; }
"<" { return '<'; }
">" { return '>'; }
"(" { return '('; }
")" { return ')'; }
"[" { return '['; }
"]" { return ']'; }
"{" { return '{'; }
"}" { return '}'; }
"=" { return '='; }
"+" { return '+'; }
"-" { return '-'; }
"*" { return '*'; }
"/" { return '/'; }
"%" { return '%'; }
"^" { return '^'; }
"!" { return '!'; }
"?" { return '?'; }
":" { return ':'; }
";" { return ';'; }
"," { return ','; }
[a-zA-Z][0-9a-zA-Z_]* { yylval.node=new PhysBAM::PROGRAM_PARSE_NODE(TOKEN_IDENT,::PhysBAM::parse_identifiers.Append(yytext),0,0); return TOKEN_IDENT; }
([0-9]+\.?|[0-9]*\.[0-9]+)([Ee][+-]?[0-9]+)? { yylval.node=new PhysBAM::PROGRAM_PARSE_NODE(TOKEN_NUMBER,::PhysBAM::parse_constants.Append(atof(yytext)),0,0); return TOKEN_NUMBER; }
[ \r\n\t]* { }
. { printf("Enexpected token: %s\n", yytext); return 0; }

%%

inline void Eliminate_Warning()
{
    (void) yyunput;
    (void) yyinput;
}
