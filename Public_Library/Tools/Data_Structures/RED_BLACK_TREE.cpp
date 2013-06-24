#include <Tools/Data_Structures/RED_BLACK_TREE.h>
using namespace PhysBAM;
//#####################################################################
// Function Rotate_Right
//#####################################################################
void RED_BLACK_TREE_CORE::
Rotate_Right(NODE* node)
{
    assert(node->L);
    NODE* tmp=node->L;
    node->L=tmp->R;
    if(tmp->R) tmp->R->P=node;
    if(tmp->P) tmp->P=node->P;
    if(!tmp->P) root=tmp;
    else if(node==node->P->R) node->P->R=tmp;
    else node->P->L=tmp;
    tmp->R=node;
    node->P=tmp;
    if(func) func(tmp,0,rot_r);
}
//#####################################################################
// Function Rotate_Left
//#####################################################################
void RED_BLACK_TREE_CORE::
Rotate_Left(NODE* node)
{
    assert(node->R);
    NODE* tmp=node->R;
    node->R=tmp->L;
    if(tmp->L) tmp->L->P=node;
    if(tmp->P) tmp->P=node->P;
    if(!tmp->P) root=tmp;
    else if(node==node->P->L) node->P->L=tmp;
    else node->P->R=tmp;
    tmp->L=node;
    node->P=tmp;
    if(func) func(tmp,0,rot_l);
}
//#####################################################################
// Function Exchange_Nodes
//#####################################################################
void RED_BLACK_TREE_CORE::
Exchange_Nodes(NODE* a,NODE* b)
{
    assert(a!=b);
    exchange(a->L,b->L);
    exchange(a->R,b->R);
    exchange(a->P,b->P);
    exchange(a->black,b->black);
    if(a->P==a){a->P=b;if(b->L==b) b->L=a;else b->R=a;}
    if(b->P==b){b->P=a;if(a->L==a) a->L=b;else a->R=b;}
    if(!a->P) root=a;
    else if(a->P->L==b) a->P->L=a;
    else if(a->P->R==b) a->P->R=a;
    if(!b->P) root=b;
    else if(b->P->L==a) b->P->L=b;
    else if(b->P->R==a) b->P->R=b;
    if(a->L) a->L->P=a;
    if(a->R) a->R->P=a;
    if(b->L) b->L->P=b;
    if(b->R) b->R->P=b;
    if(func) func(a,b,exch);
}
//#####################################################################
// Function Assert_Valid_Helper
//#####################################################################
int RED_BLACK_TREE_CORE::
Assert_Valid_Helper(NODE* node) const
{
    if(!node) return 0;
    PHYSBAM_ASSERT(!node->L || node->L->P==node);
    PHYSBAM_ASSERT(!node->R || node->R->P==node);
    PHYSBAM_ASSERT(node->black || !node->L || node->L->black);
    PHYSBAM_ASSERT(node->black || !node->R || node->R->black);
    int L=Assert_Valid_Helper(node->L);
    int R=Assert_Valid_Helper(node->R);
    PHYSBAM_ASSERT(L==R);
    return L+node->black;
}
//#####################################################################
// Function Insert
//#####################################################################
void RED_BLACK_TREE_CORE::
Root_Insert(NODE* node)
{
    size++;
    node->Reset();
    root=node;
    root->black=true;
    if(func) func(node,0,insert);
}
//#####################################################################
// Function Insert
//#####################################################################
void RED_BLACK_TREE_CORE::
Insert(NODE* node,NODE* parent,bool left)
{
    size++;
    node->Reset();
    node->P=parent;
    if(left) parent->L=node;
    else parent->R=node;
    if(func) func(node,0,insert);
    Insert_Fixup(node);
}
//#####################################################################
// Function Insert_Fixup
//#####################################################################
void RED_BLACK_TREE_CORE::
Insert_Fixup(NODE* node)
{
    while(node->P && !node->P->black && node->P->P){
        NODE* N=node->P;
        if(N==N->P->L){
            NODE* uncle=N->P->R;
            if(uncle && !uncle->black){ // case 1
                N->black=true;
                uncle->black=true;
                N->P->black=false;
                node=N->P;}
            else{
                if(node==N->R){ // case 2
                    node=N;
                    Rotate_Left(node);
                    N=node->P;}
                N->black=true; // case 3
                N->P->black=false;
                Rotate_Right(N->P);}}
        else{
            NODE* uncle=N->P->L;
            if(uncle && !uncle->black){ // case 1
                N->black=true;
                uncle->black=true;
                N->P->black=false;
                node=N->P;}
            else{
                if(node==N->L){ // case 2
                    node=N;
                    Rotate_Right(node);
                    N=node->P;}
                N->black=true; // case 3
                N->P->black=false;
                Rotate_Left(N->P);}}}
    root->black=true;
}
//#####################################################################
// Function Remove
//#####################################################################
void RED_BLACK_TREE_CORE::
Remove(NODE* node)
{
    if(node->L && node->R) Exchange_Nodes(node,node->R->Far_Left());
    assert(!node->L || !node->R);
    NODE* child=node->L?node->L:node->R?node->R:0,*P=node->P;
    if(!child && node->black){ // Use node as a sentinal
        Remove_Fixup(node);
        if(func) func(node,0,remove);
        if(!P) root=0;
        else if(node==P->L) P->L=0;
        else P->R=0;}
    else{
        if(child) child->P=P;
        if(func) func(node,0,remove);
        if(!P) root=child;
        else if(node==P->L) P->L=child;
        else P->R=child;
        if(node->black) Remove_Fixup(child);}
    node->Reset();
    size--;
}
//#####################################################################
// Function Remove_Fixup
//#####################################################################
void RED_BLACK_TREE_CORE::
Remove_Fixup(NODE* node)
{
    while(node->P && node->black){
        NODE* N=node->P;
        if(node==N->L){
            NODE *sibbling=N->R;
            if(sibbling && !sibbling->black){ // case 1
                sibbling->black=true;
                N->black=false;
                Rotate_Left(N);
                sibbling=N->R;}
            if((!sibbling->L || sibbling->L->black) && (!sibbling->R || sibbling->R->black)){ // case 2
                sibbling->black=false;
                node=N;}
            else{
                if(!sibbling->R || sibbling->R->black){ // case 3
                    sibbling->L->black=true;
                    sibbling->black=false;
                    Rotate_Right(sibbling);
                    sibbling=N->R;}
                sibbling->black=N->black; // case 4
                N->black=true;
                sibbling->R->black=true;
                Rotate_Left(N);
                node=root;}}
        else{
            NODE *sibbling=N->L;
            if(sibbling && !sibbling->black){ // case 1
                sibbling->black=true;
                N->black=false;
                Rotate_Right(N);
                sibbling=N->L;}
            if((!sibbling->R || sibbling->R->black) && (!sibbling->L || sibbling->L->black)){ // case 2
                sibbling->black=false;
                node=N;}
            else{
                if(!sibbling->L || sibbling->L->black){ // case 3
                    sibbling->R->black=true;
                    sibbling->black=false;
                    Rotate_Left(sibbling);
                    sibbling=N->L;}
                sibbling->black=N->black; // case 4
                N->black=true;
                sibbling->L->black=true;
                Rotate_Right(N);
                node=root;}}}
    node->black=true;
}
