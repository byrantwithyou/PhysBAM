//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RED_BLACK_TREE
//#####################################################################
#ifndef __RED_BLACK_TREE__
#define __RED_BLACK_TREE__

#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
#include <cassert>
namespace PhysBAM{

template<class K,class T_DATA> struct RED_BLACK_TREE_NODE_PAYLOAD{RED_BLACK_TREE_NODE_PAYLOAD(const K& k,const T_DATA& t):key(k),data(t){} K key;T_DATA data;};
template<class K> struct RED_BLACK_TREE_NODE_PAYLOAD<K,void>{RED_BLACK_TREE_NODE_PAYLOAD(const K& k):key(k){} K key;};

template<class K,class T=void> // K=key,T=data
class RED_BLACK_TREE
{
    struct UNUSABLE{};
    typedef typename IF<IS_SAME<T,void>::value,UNUSABLE,T>::TYPE& T_REF;
public:
    struct NODE:public RED_BLACK_TREE_NODE_PAYLOAD<K,T>
    {
        bool black;
        NODE *left,*right,*parent;

        NODE(const K& k,const T_REF t)
            :RED_BLACK_TREE_NODE_PAYLOAD<K,T>(k,t),black(false),left(0),right(0),parent(0)
        {}

        NODE(const K& k)
            :RED_BLACK_TREE_NODE_PAYLOAD<K,T>(k),black(false),left(0),right(0),parent(0)
        {}

        ~NODE()
        {}

        void Reset()
        {left=right=parent=0;black=false;}

        NODE* Far_Left()
        {NODE* node=this;while(node->left) node=node->left;return node;}

        NODE* Far_Right()
        {NODE* node=this;while(node->right) node=node->right;return node;}

        NODE* Far_Left_Leaf()
        {NODE* node=this;while(node->left || node->right) node=node->left?node->left:node->right;return node;}

        NODE* Far_Right_Leaf()
        {NODE* node=this;while(node->right || node->left) node=node->right?node->right:node->left;return node;}

        NODE* Inorder_Next()
        {
            if(right) return right->Far_Left();
            NODE* save=this,*node=parent;
            while(node && node->right==save){save=node;node=node->parent;}
            return node;
        }

        NODE* Inorder_Previous()
        {
            if(left) return left->Far_Right();
            NODE* save=this,*node=parent;
            while(node && node->left==save){save=node;node=node->parent;}
            return node;
        }

        NODE* Preorder_Next()
        {
            if(left) return left;
            if(right) return right;
            NODE *save=this,*node=parent;
            while(node && (node->right==save || !node->right)){save=node;node=node->parent;}
            if(node) node=node->right;
            return node;
        }

        NODE* Preorder_Previous()
        {
            if(!parent || !parent->left || parent->left==this) return parent;
            NODE* node=parent->left;
            return node->Far_Right_Leaf();
        }

        NODE* Postorder_Next()
        {
            if(!parent || !parent->right || parent->right==this) return parent;
            NODE* node=parent->right;
            return node->Far_Left_Leaf();
        }

        NODE* Postorder_Previous()
        {
            if(right) return right;
            if(left) return left;
            NODE *save=this,*node=parent;
            while(node && (node->left==save || !node->left)){save=node;node=node->parent;}
            if(node) node=node->left;
            return node;
        }
    };

    struct PREORDER_ITERATOR
    {
        NODE* node;

        PREORDER_ITERATOR(NODE* n)
            :node(n)
        {}

        PREORDER_ITERATOR(RED_BLACK_TREE& tree)
            :node(tree.root)
        {}

        bool Valid() const
        {return node!=0;}

        void Next()
        {node=node->Preorder_Next();}

        K& Key()
        {return node->key;}

        const K& Key() const
        {return node->key;}

        K& Data()
        {return node->data;}

        const K& Data() const
        {return node->data;}
    };

    struct POSTORDER_ITERATOR
    {
        NODE* node;

        POSTORDER_ITERATOR(NODE* n)
            :node(n)
        {if(node) node=node->Far_Left_Leaf();}

        POSTORDER_ITERATOR(RED_BLACK_TREE& tree)
            :node(tree.root)
        {if(node) node=node->Far_Left_Leaf();}

        bool Valid() const
        {return node!=0;}

        void Next()
        {node=node->Postorder_Next();}

        K& Key()
        {return node->key;}

        const K& Key() const
        {return node->key;}

        K& Data()
        {return node->data;}

        const K& Data() const
        {return node->data;}
    };

    struct INORDER_ITERATOR
    {
        NODE* node;

        INORDER_ITERATOR(NODE* n)
            :node(n)
        {if(node) node=node->Far_Left();}

        INORDER_ITERATOR(RED_BLACK_TREE& tree)
            :node(tree.root)
        {if(node) node=node->Far_Left();}

        bool Valid() const
        {return node!=0;}

        void Next()
        {node=node->Inorder_Next();}

        K& Key()
        {return node->key;}

        const K& Key() const
        {return node->key;}

        T_REF Data()
        {return node->data;}

        const T_REF Data() const
        {return node->data;}
    };

    NODE *root;
    int size;
    
    RED_BLACK_TREE()
        :root(0),size(0)
    {}
    
    ~RED_BLACK_TREE()
    {POSTORDER_ITERATOR iterator(*this);while(iterator.Valid()){root=iterator.node;iterator.Next();delete root;}}

    int Size() const
    {return size;}

    void Assert_Valid() const
    {
        if(!root) return;
        PHYSBAM_ASSERT(!root->parent && root->black);
        Assert_Valid_Helper(root);
    }

    void Assert_Valid_Tree() const
    {
        if(!root) return;
        PHYSBAM_ASSERT(!root->parent);
        Assert_Valid_Tree_Helper(root);
    }

    void Print() const
    {Print(root);LOG::cout<<std::endl;}

    NODE* Find(const K& key) const
    {int result;NODE* node=Find(key,result);return result?0:node;}

    NODE* Insert(const K& key,const T_REF data)
    {NODE* node=new NODE(key,data);Insert(node);return node;}

    NODE* Insert(const K& key)
    {NODE* node=new NODE(key);Insert(node);return node;}

    NODE* Insert_Unique(const K& key,const T_REF data)
    {NODE* node=new NODE(key,data);if(Insert(node,true)) return node;delete node;return 0;}

    NODE* Insert_Unique(const K& key)
    {NODE* node=new NODE(key);if(Insert(node,true)) return node;delete node;return 0;}

    bool Remove(const K& key)
    {int result;NODE* found=Find(key,result);if(!result) Remove(found);return !result;}

    NODE* First() const
    {return root?root->Far_Left():0;}

    NODE* Last() const
    {return root?root->Far_Right():0;}

    void Exchange_Nodes(NODE* a,NODE* b); // Warning: only use if you know what you are doing
protected:
    void Print(NODE* node) const
    {
        if(!node){LOG::cout<<"null";return;}
        LOG::cout<<"("<<(node->black?"*":"")<<node->key<<" ";Print(node->left);LOG::cout<<" ";Print(node->right);LOG::cout<<")";
    }
    void Assert_Valid_Tree_Helper(NODE* node) const
    {
        if(!node) return;
        PHYSBAM_ASSERT(!node->left || node->left->parent==node);
        PHYSBAM_ASSERT(!node->right || node->right->parent==node);
    }
public:
    bool Insert(NODE* node,bool unique=false);
    void Insert_Fixup(NODE* node);
    void Remove(NODE* node);
    void Remove_Fixup(NODE* node);
    NODE* Find(const K& key,int& result) const; // result: -1 -> less than returned.  0 -> equal to returned.  1 -> greater than returned.
    template<class K2,class T_COMPARISON> NODE* Find(const K2& key,int& result,const T_COMPARISON& compare) const;
protected:
    void Rotate_Right(NODE* node);
    void Rotate_Left(NODE* node);
    int Assert_Valid_Helper(NODE* node) const;
};
//#####################################################################
// Function Rotate_Right
//#####################################################################
template<class K,class T> void RED_BLACK_TREE<K,T>::
Rotate_Right(NODE* node)
{
    assert(node->left);
    NODE* tmp=node->left;
    node->left=tmp->right;
    if(tmp->right) tmp->right->parent=node;
    if(tmp->parent) tmp->parent=node->parent;
    if(!tmp->parent) root=tmp;
    else if(node==node->parent->right) node->parent->right=tmp;
    else node->parent->left=tmp;
    tmp->right=node;
    node->parent=tmp;
}
//#####################################################################
// Function Rotate_Left
//#####################################################################
template<class K,class T> void RED_BLACK_TREE<K,T>::
Rotate_Left(NODE* node)
{
    assert(node->right);
    NODE* tmp=node->right;
    node->right=tmp->left;
    if(tmp->left) tmp->left->parent=node;
    if(tmp->parent) tmp->parent=node->parent;
    if(!tmp->parent) root=tmp;
    else if(node==node->parent->left) node->parent->left=tmp;
    else node->parent->right=tmp;
    tmp->left=node;
    node->parent=tmp;
}
//#####################################################################
// Function Find
//#####################################################################
template<class K,class T> typename RED_BLACK_TREE<K,T>::NODE* RED_BLACK_TREE<K,T>::
Find(const K& key,int& result) const
{
    if(!root){result=-1;return 0;}
    NODE* node=root;
    while(1){
        if(key<node->key){
            if(!node->left){result=-1;return node;}
            node=node->left;}
        else if(key>node->key){
            if(!node->right){result=1;return node;}
            node=node->right;}
        else{result=0;return node;}}
}
//#####################################################################
// Function Find
//#####################################################################
template<class K,class T> template<class K2,class T_COMPARISON> typename RED_BLACK_TREE<K,T>::NODE* RED_BLACK_TREE<K,T>::
Find(const K2& key,int& result,const T_COMPARISON& compare) const // int compare(K,K2) -> -1,0,1
{
    NODE* node=root;
    while(1){
        int cmp=compare(node->key,key);
        if(cmp>0){
            if(!node->left){result=-1;return node;}
            node=node->left;}
        else if(cmp<0){
            if(!node->right){result=1;return node;}
            node=node->right;}
        else{result=0;return node;}}
}
//#####################################################################
// Function Exchange_Nodes
//#####################################################################
template<class K,class T> void RED_BLACK_TREE<K,T>::
Exchange_Nodes(NODE* a,NODE* b)
{
    assert(a!=b);
    exchange(a->left,b->left);
    exchange(a->right,b->right);
    exchange(a->parent,b->parent);
    exchange(a->black,b->black);
    if(a->parent==a){a->parent=b;if(b->left==b) b->left=a;else b->right=a;}
    if(b->parent==b){b->parent=a;if(a->left==a) a->left=b;else a->right=b;}
    if(!a->parent) root=a;
    else if(a->parent->left==b) a->parent->left=a;
    else if(a->parent->right==b) a->parent->right=a;
    if(!b->parent) root=b;
    else if(b->parent->left==a) b->parent->left=b;
    else if(b->parent->right==a) b->parent->right=b;
    if(a->left) a->left->parent=a;
    if(a->right) a->right->parent=a;
    if(b->left) b->left->parent=b;
    if(b->right) b->right->parent=b;
}
//#####################################################################
// Function Assert_Valid_Helper
//#####################################################################
template<class K,class T> int RED_BLACK_TREE<K,T>::
Assert_Valid_Helper(NODE* node) const
{
    if(!node) return 0;
    PHYSBAM_ASSERT(!node->left || node->left->key<=node->key);
    PHYSBAM_ASSERT(!node->right || node->right->key>=node->key);
    PHYSBAM_ASSERT(!node->left || node->left->parent==node);
    PHYSBAM_ASSERT(!node->right || node->right->parent==node);
    PHYSBAM_ASSERT(node->black || !node->left || node->left->black);
    PHYSBAM_ASSERT(node->black || !node->right || node->right->black);
    int left=Assert_Valid_Helper(node->left);
    int right=Assert_Valid_Helper(node->right);
    PHYSBAM_ASSERT(left==right);
    return left+node->black;
}
//#####################################################################
// Function Insert
//#####################################################################
template<class K,class T> bool RED_BLACK_TREE<K,T>::
Insert(NODE* node,bool unique)
{
    size++;
    node->Reset();
    if(!root){root=node;root->black=true;return true;}
    NODE *parent=0,*tmp=root;
    while(tmp){
        parent=tmp;
        tmp=(node->key<tmp->key)?tmp->left:tmp->right;}
    node->parent=parent;
    if(node->key<parent->key) parent->left=node;
    else if(!unique || parent->key<node->key) parent->right=node;
    else return false;
    Insert_Fixup(node);
    return true;
}
//#####################################################################
// Function Insert_Fixup
//#####################################################################
template<class K,class T> void RED_BLACK_TREE<K,T>::
Insert_Fixup(NODE* node)
{
    while(node->parent && !node->parent->black && node->parent->parent){
        NODE* parent=node->parent;
        if(parent==parent->parent->left){
            NODE* uncle=parent->parent->right;
            if(uncle && !uncle->black){ // case 1
                parent->black=true;
                uncle->black=true;
                parent->parent->black=false;
                node=parent->parent;}
            else{
                if(node==parent->right){ // case 2
                    node=parent;
                    Rotate_Left(node);
                    parent=node->parent;}
                parent->black=true; // case 3
                parent->parent->black=false;
                Rotate_Right(parent->parent);}}
        else{
            NODE* uncle=parent->parent->left;
            if(uncle && !uncle->black){ // case 1
                parent->black=true;
                uncle->black=true;
                parent->parent->black=false;
                node=parent->parent;}
            else{
                if(node==parent->left){ // case 2
                    node=parent;
                    Rotate_Right(node);
                    parent=node->parent;}
                parent->black=true; // case 3
                parent->parent->black=false;
                Rotate_Left(parent->parent);}}}
    root->black=true;
}
//#####################################################################
// Function Remove
//#####################################################################
template<class K,class T> void RED_BLACK_TREE<K,T>::
Remove(NODE* node)
{
    if(node->left && node->right){
        NODE* temp=node->right;while(temp->left) temp=temp->left;
        Exchange_Nodes(node,temp);}
    assert(!node->left || !node->right);
    NODE* child=node->left?node->left:node->right?node->right:0,*parent=node->parent;
    if(!child && node->black){ // Use node as a sentinal
        Remove_Fixup(node);
        if(!parent) root=0;
        else if(node==parent->left) parent->left=0;
        else parent->right=0;}
    else{
        if(child) child->parent=parent;
        if(!parent) root=child;
        else if(node==parent->left) parent->left=child;
        else parent->right=child;
        if(node->black) Remove_Fixup(child);}
    node->Reset();
    size--;
    delete node;
}
//#####################################################################
// Function Remove_Fixup
//#####################################################################
template<class K,class T> void RED_BLACK_TREE<K,T>::
Remove_Fixup(NODE* node)
{
    while(node->parent && node->black){
        NODE* parent=node->parent;
        if(node==parent->left){
            NODE *sibbling=parent->right;
            if(sibbling && !sibbling->black){ // case 1
                sibbling->black=true;
                parent->black=false;
                Rotate_Left(parent);
                sibbling=parent->right;}
            if((!sibbling->left || sibbling->left->black) && (!sibbling->right || sibbling->right->black)){ // case 2
                sibbling->black=false;
                node=parent;}
            else{
                if(!sibbling->right || sibbling->right->black){ // case 3
                    sibbling->left->black=true;
                    sibbling->black=false;
                    Rotate_Right(sibbling);
                    sibbling=parent->right;}
                sibbling->black=parent->black; // case 4
                parent->black=true;
                sibbling->right->black=true;
                Rotate_Left(parent);
                node=root;}}
        else{
            NODE *sibbling=parent->left;
            if(sibbling && !sibbling->black){ // case 1
                sibbling->black=true;
                parent->black=false;
                Rotate_Right(parent);
                sibbling=parent->left;}
            if((!sibbling->right || sibbling->right->black) && (!sibbling->left || sibbling->left->black)){ // case 2
                sibbling->black=false;
                node=parent;}
            else{
                if(!sibbling->left || sibbling->left->black){ // case 3
                    sibbling->right->black=true;
                    sibbling->black=false;
                    Rotate_Left(sibbling);
                    sibbling=parent->left;}
                sibbling->black=parent->black; // case 4
                parent->black=true;
                sibbling->left->black=true;
                Rotate_Right(parent);
                node=root;}}}
    node->black=true;
}
}
#endif
