//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RED_BLACK_TREE
//#####################################################################
#ifndef __RED_BLACK_TREE__
#define __RED_BLACK_TREE__

#include <Tools/Log/DEBUG_UTILITIES.h>
#include <Tools/Log/LOG.h>
#include <Tools/Utilities/TYPE_UTILITIES.h>
#include <cassert>
namespace PhysBAM{

template<class K,class T=void> class RED_BLACK_TREE;

class RED_BLACK_TREE_CORE
{
    struct UNUSABLE{};
public:
    struct NODE
    {
        NODE *L,*R,*P;
        bool black;

        NODE() :L(0),R(0),P(0),black(false) {}
        ~NODE() {}

        void Reset()
        {L=R=P=0;black=false;}

        NODE* Far_Left()
        {NODE* node=this;while(node->L) node=node->L;return node;}

        NODE* Far_Right()
        {NODE* node=this;while(node->R) node=node->R;return node;}

        NODE* Far_Left_Leaf()
        {NODE* node=this;while(node->L || node->R) node=node->L?node->L:node->R;return node;}

        NODE* Far_Right_Leaf()
        {NODE* node=this;while(node->R || node->L) node=node->R?node->R:node->L;return node;}

        NODE* Inorder_Next()
        {
            if(R) return R->Far_Left();
            NODE* save=this,*node=P;
            while(node && node->R==save){save=node;node=node->P;}
            return node;
        }

        NODE* Inorder_Prev()
        {
            if(L) return L->Far_Right();
            NODE* save=this,*node=P;
            while(node && node->L==save){save=node;node=node->P;}
            return node;
        }

        NODE* Preorder_Next()
        {
            if(L) return L;
            if(R) return R;
            NODE *save=this,*node=P;
            while(node && (node->R==save || !node->R)){save=node;node=node->P;}
            if(node) node=node->R;
            return node;
        }

        NODE* Preorder_Prev()
        {
            if(!P || !P->L || P->L==this) return P;
            NODE* node=P->L;
            return node->Far_Right_Leaf();
        }

        NODE* Postorder_Next()
        {
            if(!P || !P->R || P->R==this) return P;
            NODE* node=P->R;
            return node->Far_Left_Leaf();
        }

        NODE* Postorder_Prev()
        {
            if(R) return R;
            if(L) return L;
            NODE *save=this,*node=P;
            while(node && (node->L==save || !node->L)){save=node;node=node->P;}
            if(node) node=node->L;
            return node;
        }
    };

    NODE *root;
    int size;

    enum FUNC_CODE {insert,remove,rot_l,rot_r,exch};
    void (*func)(NODE* A,NODE* B,int code);

    RED_BLACK_TREE_CORE()
        :root(0),size(0),func(0)
    {}
    
    ~RED_BLACK_TREE_CORE()
    {}

    void Assert_Valid() const
    {
        if(!root) return;
        PHYSBAM_ASSERT(!root->P && root->black);
        Assert_Valid_Helper(root);
    }

    void Assert_Valid_Tree() const
    {
        if(!root) return;
        PHYSBAM_ASSERT(!root->P);
        Assert_Valid_Tree_Helper(root);
    }

    NODE* First() const
    {return root?root->Far_Left():0;}

    NODE* Last() const
    {return root?root->Far_Right():0;}

    void Exchange_Nodes(NODE* a,NODE* b); // Warning: only use if you know what you are doing
    void Assert_Valid_Tree_Helper(NODE* node) const
    {
        if(!node) return;
        PHYSBAM_ASSERT(!node->L || node->L->P==node);
        PHYSBAM_ASSERT(!node->R || node->R->P==node);
    }
    void Root_Insert(NODE* node);
    void Insert(NODE* node,NODE* parent,bool left);
    void Insert_Fixup(NODE* node);
    void Remove(NODE* node);
    void Remove_Fixup(NODE* node);
    void Rotate_Right(NODE* node);
    void Rotate_Left(NODE* node);
    int Assert_Valid_Helper(NODE* node) const;
};

template<class K,class T_DATA> struct RED_BLACK_TREE_NODE_PAYLOAD{RED_BLACK_TREE_NODE_PAYLOAD(const K& k,const T_DATA& t):key(k),data(t){} K key;T_DATA data;};
template<class K> struct RED_BLACK_TREE_NODE_PAYLOAD<K,void>{RED_BLACK_TREE_NODE_PAYLOAD(const K& k):key(k){} K key;};

template<class K,class T> void (*Set_T_Callback_Helper(int(*)[sizeof(&T::Callback)]))(typename RED_BLACK_TREE<K,T>::NODE*,typename RED_BLACK_TREE<K,T>::NODE*,int) {return &T::Callback;}
void (*Set_T_Callback_Helper(...))(RED_BLACK_TREE_CORE::NODE*,RED_BLACK_TREE_CORE::NODE*,int) {return 0;}

template<class K,class T> // K=key,T=data (=void)
class RED_BLACK_TREE
{
    struct UNUSABLE{};
    typedef typename IF<IS_SAME<T,void>::value,UNUSABLE,T>::TYPE T_OR_UNUSABLE;
    typedef RED_BLACK_TREE_CORE CORE;
public:
    CORE core;

    struct NODE:public CORE::NODE,public RED_BLACK_TREE_NODE_PAYLOAD<K,T>
    {
        NODE(const K& k,const T_OR_UNUSABLE& t): RED_BLACK_TREE_NODE_PAYLOAD<K,T>(k,t) {}
        explicit NODE(const K& k): RED_BLACK_TREE_NODE_PAYLOAD<K,T>(k) {}
        ~NODE() {}
    };

    struct ITERATOR
    {
        NODE* node;
        explicit ITERATOR(NODE* n=0) :node(n) {}
        explicit ITERATOR(RED_BLACK_TREE& tree) :node(static_cast<NODE*>(tree.core.root)) {}
        bool Valid() const {return node!=0;}
        const K& Key() const {return node->key;}
        K& Data() {return node->data;}
        const K& Data() const {return node->data;}
        void Parent() {node=static_cast<NODE*>(node->P);}
        void Left() {node=static_cast<NODE*>(node->L);}
        void Right() {node=static_cast<NODE*>(node->R);}
    };

    struct PREORDER_ITERATOR:public ITERATOR
    {
        using ITERATOR::node;
        PREORDER_ITERATOR(ITERATOR it): ITERATOR(it) {}
        explicit PREORDER_ITERATOR(NODE* n=0): ITERATOR(n) {}
        explicit PREORDER_ITERATOR(RED_BLACK_TREE& tree): ITERATOR(tree) {}
        void Next() {node=static_cast<NODE*>(node->Preorder_Next());}
        void Prev() {node=static_cast<NODE*>(node->Preorder_Prev());}
    };

    struct POSTORDER_ITERATOR:public ITERATOR
    {
        using ITERATOR::node;
        POSTORDER_ITERATOR(ITERATOR it): ITERATOR(it) {}
        explicit POSTORDER_ITERATOR(NODE* n=0): ITERATOR(n) {if(node) node=static_cast<NODE*>(node->Far_Left_Leaf());}
        explicit POSTORDER_ITERATOR(RED_BLACK_TREE& tree): ITERATOR(tree) {if(node) node=static_cast<NODE*>(node->Far_Left_Leaf());}
        void Next() {node=static_cast<NODE*>(node->Postorder_Next());}
        void Prev() {node=static_cast<NODE*>(node->Postorder_Prev());}
    };

    struct INORDER_ITERATOR:public ITERATOR
    {
        using ITERATOR::node;
        INORDER_ITERATOR(ITERATOR it): ITERATOR(it) {}
        explicit INORDER_ITERATOR(NODE* n=0): ITERATOR(n) {if(node) node=static_cast<NODE*>(node->Far_Left());}
        explicit INORDER_ITERATOR(RED_BLACK_TREE& tree): ITERATOR(tree) {if(node) node=static_cast<NODE*>(node->Far_Left());}
        void Next() {node=static_cast<NODE*>(node->Inorder_Next());}
        void Prev() {node=static_cast<NODE*>(node->Inorder_Prev());}
    };

    RED_BLACK_TREE() {core.func=reinterpret_cast<void (*)(CORE::NODE*,CORE::NODE*,int)>(Set_T_Callback_Helper(0));}

    ~RED_BLACK_TREE() {Clean_Memory();}

    void Remove_All() {Clean_Memory();}

    void Clean_Memory()
    {
        POSTORDER_ITERATOR it(*this);
        while(it.Valid()){
            core.root=it.node;
            it.Next();
            delete static_cast<NODE*>(core.root);}
    }

    int Size() const {return core.size;}

    void Assert_Valid() const
    {
        if(!core.root) return;
        Assert_Valid_Helper(core.root);
        core.Assert_Valid_Helper(core.root);
    }

    ITERATOR Find(const K& key) const
    {int result;ITERATOR found=Find(key,result);return result?ITERATOR():found;}

    ITERATOR Insert(const K& key,const T_OR_UNUSABLE& data)
    {NODE* node=new NODE(key,data);Insert(node);return ITERATOR(node);}

    ITERATOR Insert(const K& key)
    {NODE* node=new NODE(key);Insert(node);return ITERATOR(node);}

    ITERATOR Insert_Unique(const K& key,const T_OR_UNUSABLE& data)
    {NODE* node=new NODE(key,data);if(Insert(node,true)) return ITERATOR(node);delete node;return ITERATOR(0);}

    ITERATOR Insert_Unique(const K& key)
    {NODE* node=new NODE(key);if(Insert(node,true)) return ITERATOR(node);delete node;return ITERATOR(0);}

    bool Remove(const K& key)
    {int result;ITERATOR found=Find(key,result);if(!result) Remove(found);return !result;}

    void Remove(ITERATOR it)
    {
        core.Remove(it.node);
        delete it.node;
    }

    ITERATOR Find(const K& key,int& result) const // result: -1 -> less than returned.  0 -> equal to returned.  1 -> greater than returned.
    {
        if(!core.root){result=-1;return ITERATOR();}
        NODE* node=static_cast<NODE*>(core.root);
        while(1){
            if(key<node->key){
                if(!node->L){result=-1;return ITERATOR(node);}
                node=static_cast<NODE*>(node->L);}
            else if(key>node->key){
                if(!node->R){result=1;return ITERATOR(node);}
                node=static_cast<NODE*>(node->R);}
            else{result=0;return ITERATOR(node);}}
    }

    ITERATOR Lower_Bound(const K& key) const
    {
        if(!core.root) return ITERATOR();
        NODE* node=static_cast<NODE*>(core.root);
        while(1){
            if(node->key<key){
                if(!node->R) return ITERATOR(node->Inorder_Next());
                node=static_cast<NODE*>(node->R);}
            else{
                if(!node->L) return ITERATOR(node);
                node=static_cast<NODE*>(node->L);}}
    }

    ITERATOR Upper_Bound(const K& key) const
    {
        if(!core.root) return ITERATOR();
        NODE* node=static_cast<NODE*>(core.root);
        while(1){
            if(key<node->key){
                if(!node->L) return ITERATOR(node);
                node=static_cast<NODE*>(node->L);}
            else{
                if(!node->R) return ITERATOR(node->Inorder_Next());
                node=static_cast<NODE*>(node->R);}}
    }

    void Print(std::ostream& o=LOG::cout) const
    {
        if(core.root) Print(o,core.root);
        else o<<"null";
    }

protected:
    bool Insert(NODE* node,bool unique=false)
    {
        if(!core.root){core.Root_Insert(node);return true;}
        NODE *P=0,*tmp=static_cast<NODE*>(core.root);
        while(tmp){
            P=tmp;
            tmp=static_cast<NODE*>((node->key<tmp->key)?tmp->L:tmp->R);}
        bool left=node->key<P->key;
        if(!left && unique && !(P->key<node->key)) return false;
        core.Insert(node,P,left);
        return true;
    }

    void Print_Data_Helper(std::ostream& o,CORE::NODE* node,const UNUSABLE*) const {}
    void Print_Data_Helper(std::ostream& o,CORE::NODE* node,const T*) const {o<<" : "<<static_cast<NODE*>(node)->data;}

    void Print(std::ostream& o,CORE::NODE* node) const
    {
        o<<'(';
        if(node->L){
            Print(o,node->L);
            o<<' ';}
        if(node->black) o<<'*';
        o<<static_cast<NODE*>(node)->key;
        Print_Data_Helper(o,node,(T_OR_UNUSABLE*)0);
        if(node->R){
            o<<' ';
            Print(o,node->R);}
        o<<')';
    }

    void Assert_Valid_Helper(CORE::NODE* node) const
    {
        if(node->L){
            PHYSBAM_ASSERT(static_cast<NODE*>(node->L)->key<=static_cast<NODE*>(node)->key);
            Assert_Valid_Helper(node->L);}
        if(node->R){
            PHYSBAM_ASSERT(static_cast<NODE*>(node->R)->key>=static_cast<NODE*>(node)->key);
            Assert_Valid_Helper(node->R);}
    }
};
}
#endif
