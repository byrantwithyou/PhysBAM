%{

#include <Tools/Symbolics/PROGRAM_DEFINITIONS.h>
#include <Tools/Symbolics/PROGRAM_LEX.hpp>
#include <Tools/Symbolics/PROGRAM_YACC.hpp>

using namespace PhysBAM;

inline int yyerror(const char *msg) {
    fprintf(stderr,"Error:%s\n",msg); return 0;
}

%}

%defines "Public_Library/Tools/Symbolics/PROGRAM_YACC.hpp"

%union {
    PhysBAM::PROGRAM_PARSE_NODE *node;
}

%token TOKEN_ADDEQ
%token TOKEN_SUBEQ
%token TOKEN_MULEQ
%token TOKEN_DIVEQ
%token TOKEN_POWEQ
%token TOKEN_EQ
%token TOKEN_NE
%token TOKEN_LE
%token TOKEN_GE
%token TOKEN_AND
%token TOKEN_OR
%token TOKEN_LIST
%token TOKEN_FUNC
%token TOKEN_INDEX
%token TOKEN_IF
%token TOKEN_ELSE
%token <node> TOKEN_NUMBER
%token <node> TOKEN_IDENT
%nonassoc IF_NO_ELSE
%nonassoc TOKEN_ELSE

%type <node> and_expression or_expression program
%type <node> statement_list expression statement assignment_expression
%type <node> equals_expression conditional_expression inner_expression
%type <node> comparison_expression postfix_expression prefix_expression
%type <node> addition_expression multiplication_expression function_argument_list
%start program

%%

program
    : statement_list { parse_root = $1; }
    ;

statement_list
    : expression { $$ = $1; }
    | statement { $$ = $1; }  
    | statement statement_list { $$ = new PROGRAM_PARSE_NODE(TOKEN_LIST,0,$1,$2); }
    ;

statement
    : expression ';' { $$ = $1; }
    | '{' statement_list '}' { $$ = $2; }
    | TOKEN_IF '(' expression ')' statement %prec IF_NO_ELSE { $$ = new PROGRAM_PARSE_NODE('?',0,$3,new PROGRAM_PARSE_NODE(':',0,$5,0)); }
    | TOKEN_IF '(' expression ')' statement TOKEN_ELSE statement { $$ = new PROGRAM_PARSE_NODE('?',0,$3,new PROGRAM_PARSE_NODE(':',0,$5,$7)); }
    ;

expression
    : assignment_expression { $$ = $1; }
    | expression ',' assignment_expression { $$ = new PROGRAM_PARSE_NODE(',',0,$1,$3); }
    ;

assignment_expression
    : conditional_expression { $$ = $1; }
    | conditional_expression '=' assignment_expression { $$ = new PROGRAM_PARSE_NODE('=',0,$1,$3); }
    | conditional_expression TOKEN_ADDEQ assignment_expression { $$ = new PROGRAM_PARSE_NODE(TOKEN_ADDEQ,0,$1,$3); }
    | conditional_expression TOKEN_SUBEQ assignment_expression { $$ = new PROGRAM_PARSE_NODE(TOKEN_SUBEQ,0,$1,$3); }
    | conditional_expression TOKEN_MULEQ assignment_expression { $$ = new PROGRAM_PARSE_NODE(TOKEN_MULEQ,0,$1,$3); }
    | conditional_expression TOKEN_DIVEQ assignment_expression { $$ = new PROGRAM_PARSE_NODE(TOKEN_DIVEQ,0,$1,$3); }
    | conditional_expression TOKEN_POWEQ assignment_expression { $$ = new PROGRAM_PARSE_NODE(TOKEN_POWEQ,0,$1,$3); }
    ;

conditional_expression
    : or_expression { $$ = $1; }
    | or_expression '?' expression ':' conditional_expression { $$ = new PROGRAM_PARSE_NODE('?',0,$1,new PROGRAM_PARSE_NODE(':',0,$3,$5)); }
    ;

or_expression
    : and_expression { $$ = $1; }
    | or_expression TOKEN_OR and_expression { $$ = new PROGRAM_PARSE_NODE(TOKEN_OR,0,$1,$3); }
    ;

and_expression
    : equals_expression { $$ = $1; }
    | and_expression TOKEN_AND equals_expression { $$ = new PROGRAM_PARSE_NODE(TOKEN_AND,0,$1,$3); }
    ;

equals_expression
    : comparison_expression { $$ = $1; }
    | equals_expression TOKEN_EQ comparison_expression { $$ = new PROGRAM_PARSE_NODE(TOKEN_EQ,0,$1,$3); }
    | equals_expression TOKEN_NE comparison_expression { $$ = new PROGRAM_PARSE_NODE(TOKEN_NE,0,$1,$3); }
    ;

comparison_expression
    : addition_expression { $$ = $1; }
    | comparison_expression '<' addition_expression { $$ = new PROGRAM_PARSE_NODE('<',0,$1,$3); }
    | comparison_expression '>' addition_expression { $$ = new PROGRAM_PARSE_NODE('>',0,$1,$3); }
    | comparison_expression TOKEN_LE addition_expression { $$ = new PROGRAM_PARSE_NODE(TOKEN_LE,0,$1,$3); }
    | comparison_expression TOKEN_GE addition_expression { $$ = new PROGRAM_PARSE_NODE(TOKEN_GE,0,$1,$3); }
    ;

addition_expression
    : multiplication_expression { $$ = $1; }
    | addition_expression '+' multiplication_expression { $$ = new PROGRAM_PARSE_NODE('+',0,$1,$3); }
    | addition_expression '-' multiplication_expression { $$ = new PROGRAM_PARSE_NODE('-',0,$1,$3); }
    ;

multiplication_expression
    : prefix_expression { $$ = $1; }
    | multiplication_expression '*' prefix_expression { $$ = new PROGRAM_PARSE_NODE('*',0,$1,$3); }
    | multiplication_expression '/' prefix_expression { $$ = new PROGRAM_PARSE_NODE('/',0,$1,$3); }
    | multiplication_expression '%' prefix_expression { $$ = new PROGRAM_PARSE_NODE('%',0,$1,$3); }
    ;

prefix_expression
    : postfix_expression { $$ = $1; }
    | postfix_expression '^' prefix_expression { $$ = new PROGRAM_PARSE_NODE('^',0,$1,$3); }
    | '-' prefix_expression { $$ = new PROGRAM_PARSE_NODE('-',0,$2,0); }
    | '+' prefix_expression { $$ = new PROGRAM_PARSE_NODE('+',0,$2,0); }
    | '!' prefix_expression { $$ = new PROGRAM_PARSE_NODE('!',0,$2,0); }
    ;

postfix_expression
    : inner_expression { $$ = $1; }
    | TOKEN_IDENT '(' function_argument_list ')' { $$ = new PROGRAM_PARSE_NODE(TOKEN_FUNC,0,$1,$3); }
    | postfix_expression '[' expression ']' { $$ = new PROGRAM_PARSE_NODE(TOKEN_INDEX,0,$1,$3); }
    ;

function_argument_list
    :  { $$ = 0; }
    | assignment_expression { $$ = new PROGRAM_PARSE_NODE(TOKEN_LIST,0,$1,0); }
    | assignment_expression ',' function_argument_list { $$ = new PROGRAM_PARSE_NODE(TOKEN_LIST,0,$1,$3); }
    ;

inner_expression
    : '(' expression ')' { $$ = $2; }
    | TOKEN_IDENT { $$ = $1; }
    | TOKEN_NUMBER { $$ = $1; }
    ;

%%
