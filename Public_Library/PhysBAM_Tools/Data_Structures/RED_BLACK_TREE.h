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

        NODE* Inorder_Previous()
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

        NODE* Preorder_Previous()
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

        NODE* Postorder_Previous()
        {
            if(R) return R;
            if(L) return L;
            NODE *save=this,*node=P;
            while(node && (node->L==save || !node->L)){save=node;node=node->P;}
            if(node) node=node->L;
            return node;
        }
    };

    struct PREORDER_ITERATOR
    {
        NODE* node;

        PREORDER_ITERATOR(NODE* n)
            :node(n)
        {}

        PREORDER_ITERATOR(RED_BLACK_TREE_CORE& tree)
            :node(tree.root)
        {}

        bool Valid() const
        {return node!=0;}

        void Next()
        {node=node->Preorder_Next();}
    };

    struct POSTORDER_ITERATOR
    {
        NODE* node;

        POSTORDER_ITERATOR(NODE* n)
            :node(n)
        {if(node) node=node->Far_Left_Leaf();}

        POSTORDER_ITERATOR(RED_BLACK_TREE_CORE& tree)
            :node(tree.root)
        {if(node) node=node->Far_Left_Leaf();}

        bool Valid() const
        {return node!=0;}

        void Next()
        {node=node->Postorder_Next();}
    };

    struct INORDER_ITERATOR
    {
        NODE* node;

        INORDER_ITERATOR(NODE* n)
            :node(n)
        {if(node) node=node->Far_Left();}

        INORDER_ITERATOR(RED_BLACK_TREE_CORE& tree)
            :node(tree.root)
        {if(node) node=node->Far_Left();}

        bool Valid() const
        {return node!=0;}

        void Next()
        {node=node->Inorder_Next();}
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

template<class K,class T=void> // K=key,T=data
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
        NODE(const K& k): RED_BLACK_TREE_NODE_PAYLOAD<K,T>(k) {}
        ~NODE() {}
    };

    struct PREORDER_ITERATOR:public CORE::PREORDER_ITERATOR
    {
        PREORDER_ITERATOR(NODE* n): CORE::PREORDER_ITERATOR(n) {}
        PREORDER_ITERATOR(RED_BLACK_TREE& tree): CORE::PREORDER_ITERATOR(tree.core.root) {}
        const K& Key() const {return static_cast<NODE*>(node)->key;}
        K& Data() {return static_cast<NODE*>(node)->data;}
        const K& Data() const {return static_cast<NODE*>(node)->data;}
    };

    struct POSTORDER_ITERATOR:public CORE::POSTORDER_ITERATOR
    {
        POSTORDER_ITERATOR(NODE* n): CORE::POSTORDER_ITERATOR(n) {}
        POSTORDER_ITERATOR(RED_BLACK_TREE& tree): CORE::POSTORDER_ITERATOR(tree.core.root) {}
        const K& Key() const {return static_cast<NODE*>(node)->key;}
        K& Data() {return static_cast<NODE*>(node)->data;}
        const K& Data() const {return static_cast<NODE*>(node)->data;}
    };

    struct INORDER_ITERATOR:public CORE::INORDER_ITERATOR
    {
        INORDER_ITERATOR(NODE* n): CORE::INORDER_ITERATOR(n) {}
        INORDER_ITERATOR(RED_BLACK_TREE& tree): CORE::INORDER_ITERATOR(tree.core.root) {}
        const K& Key() const {return static_cast<NODE*>(node)->key;}
        T_OR_UNUSABLE& Data() {return static_cast<NODE*>(node)->data;}
        const T_OR_UNUSABLE& Data() const {return static_cast<NODE*>(node)->data;}
    };

    struct CONST_PREORDER_ITERATOR:public CORE::PREORDER_ITERATOR
    {
        CONST_PREORDER_ITERATOR(NODE* n): CORE::PREORDER_ITERATOR(n) {}
        CONST_PREORDER_ITERATOR(const RED_BLACK_TREE& tree): CORE::PREORDER_ITERATOR(const_cast<CORE::NODE*>(tree.core.root)) {}
        const K& Key() const {return static_cast<NODE*>(node)->key;}
        const K& Data() const {return static_cast<NODE*>(node)->data;}
    };

    struct CONST_POSTORDER_ITERATOR:public CORE::POSTORDER_ITERATOR
    {
        CONST_POSTORDER_ITERATOR(NODE* n): CORE::POSTORDER_ITERATOR(n) {}
        CONST_POSTORDER_ITERATOR(const RED_BLACK_TREE& tree): CORE::POSTORDER_ITERATOR(const_cast<CORE::NODE*>(tree.core.root)) {}
        const K& Key() const {return static_cast<NODE*>(node)->key;}
        const K& Data() const {return static_cast<NODE*>(node)->data;}
    };

    struct CONST_INORDER_ITERATOR:public CORE::INORDER_ITERATOR
    {
        CONST_INORDER_ITERATOR(NODE* n): CORE::INORDER_ITERATOR(n) {}
        CONST_INORDER_ITERATOR(const RED_BLACK_TREE& tree): CORE::INORDER_ITERATOR(const_cast<CORE::NODE*>(tree.core.root)) {}
        const K& Key() const {return static_cast<NODE*>(node)->key;}
        const T_OR_UNUSABLE& Data() const {return static_cast<NODE*>(node)->data;}
    };
    
    RED_BLACK_TREE() {}

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

    void Set_Callback(void (*func)(NODE* A,NODE* B,int code))
    {
        core.func=reinterpret_cast<void (*)(CORE::NODE*,CORE::NODE*,int)>(func);
    }
    
    void Assert_Valid() const
    {
        if(!core.root) return;
        Assert_Valid_Helper(core.root);
        core.Assert_Valid_Helper(core.root);
    }

    NODE* Find(const K& key) const
    {int result;NODE* node=Find(key,result);return result?0:node;}

    NODE* Insert(const K& key,const T_OR_UNUSABLE& data)
    {NODE* node=new NODE(key,data);Insert(node);return node;}

    NODE* Insert(const K& key)
    {NODE* node=new NODE(key);Insert(node);return node;}

    NODE* Insert_Unique(const K& key,const T_OR_UNUSABLE& data)
    {NODE* node=new NODE(key,data);if(Insert(node,true)) return node;delete node;return 0;}

    NODE* Insert_Unique(const K& key)
    {NODE* node=new NODE(key);if(Insert(node,true)) return node;delete node;return 0;}

    bool Remove(const K& key)
    {int result;NODE* found=Find(key,result);if(!result) Remove(found);return !result;}

    void Remove(NODE* node)
    {
        core.Remove(node);
        delete node;
    }

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

    NODE* Find(const K& key,int& result) const // result: -1 -> less than returned.  0 -> equal to returned.  1 -> greater than returned.
    {
        if(!core.root){result=-1;return 0;}
        NODE* node=static_cast<NODE*>(core.root);
        while(1){
            if(key<node->key){
                if(!node->L){result=-1;return node;}
                node=static_cast<NODE*>(node->L);}
            else if(key>node->key){
                if(!node->R){result=1;return node;}
                node=static_cast<NODE*>(node->R);}
            else{result=0;return node;}}
    }

    void Print_Data_Helper(std::ostream& o,CORE::NODE* node,const UNUSABLE*) const {}
    void Print_Data_Helper(std::ostream& o,CORE::NODE* node,const T*) const {o<<" : "<<static_cast<NODE*>(node)->data;}

    void Print(std::ostream& o=LOG::cout) const
    {
        if(core.root) Print(o,core.root);
        else o<<"null";
    }

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

protected:
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
