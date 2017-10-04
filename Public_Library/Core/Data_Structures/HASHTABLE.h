//#####################################################################
// Copyright 2002-2008, Robert Bridson, Kevin Der, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra, Neil Molino, Craig Schroeder, Andrew Selle, Tamar Shinar, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HASHTABLE
//#####################################################################
#ifndef __HASHTABLE__
#define __HASHTABLE__

#include <Core/Arrays/ARRAY.h>
#include <Core/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <Core/Data_Structures/HASH_REDUCE.h>
#include <Core/Data_Structures/PAIR.h>
#include <Core/Math_Tools/integer_log.h>
#include <Core/Utilities/EXCEPTIONS.h>
namespace PhysBAM{

enum HASHTABLE_ENTRY_STATE{ENTRY_FREE,ENTRY_ACTIVE,ENTRY_DELETED};
template<class TK,class T> struct HASHTABLE_ENTRY_TEMPLATE
{
    HASHTABLE_ENTRY_STATE state;
    TK key;
    T data;
};
template<class TK> struct HASHTABLE_ENTRY_TEMPLATE<TK,void>
{
    HASHTABLE_ENTRY_STATE state;
    TK key;
    TK& operator*(){return key;}
    const TK& operator*() const {return key;}
    TK* operator->(){return &key;}
    const TK* operator->() const {return &key;}
};

template<class TK,class T> // T = void
class HASHTABLE
{
private:
    struct UNUSABLE{};
    typedef HASHTABLE_ENTRY_TEMPLATE<TK,T> ENTRY; // don't store data if T is void
    typedef typename conditional<is_same<T,void>::value,UNUSABLE,T>::type T_UNLESS_VOID;
public:
    typedef TK KEY;
    typedef T ELEMENT;
    typedef ENTRY value_type; // for stl
    typedef int HAS_UNTYPED_READ_WRITE;
#ifdef NDEBUG
    static const bool is_debug=false;
#else
    static const bool is_debug=true;
#endif

    ARRAY<ENTRY> table;
    int number_of_entries;
    int next_resize;

    HASHTABLE(const int estimated_max_number_of_entries=5)
    {
        Initialize_New_Table(estimated_max_number_of_entries);
    }

    ~HASHTABLE()
    {}

    void Clean_Memory()
    {Initialize_New_Table(5);} // can't resize to zero since table.m must be a power of two

    int Size() const
    {return number_of_entries;}

    int Max_Size() const
    {return table.m;}

    int Next_Resize() const
    {return next_resize;}

    void Initialize_New_Table(const int estimated_max_number_of_entries_input)
    {next_resize=max(5,estimated_max_number_of_entries_input);
    int estimated_table_entries=(unsigned int)(next_resize*4/3+1); // choose so load is .75
    int number_of_lists=next_power_of_two(estimated_table_entries);
    table.Resize(number_of_lists,false,false); // TODO: only resize if significantly reducing the size
    Remove_All();}

    void Resize_Table(const int estimated_max_number_of_entries_input=0)
    {int estimated_max_number_of_entries=estimated_max_number_of_entries_input;if(!estimated_max_number_of_entries) estimated_max_number_of_entries=3*number_of_entries/2;
    ARRAY<ENTRY> old_table;old_table.Exchange(table);Initialize_New_Table(estimated_max_number_of_entries);
    for(int h=0;h<old_table.m;h++){ENTRY& entry=old_table(h);if(entry.state==ENTRY_ACTIVE) Insert(entry);}}

private:
    int Next_Index(const int h) const // linear probing
    {return ((h+1)&(table.m-1));} // power of two so mod is dropping high order bits

    int Hash_Index(const TK& v) const
    {return (Hash(v)&(table.m-1));} // power of two so mod is dropping high order bits

    // pair.x = location in hash where this entry is/should be be placed.
    // pair.y = true if the key is already present at pair.x
    // If check_exists=false, don't check for the key being present
    template<bool check_exists>
    PAIR<int,bool> Search(const TK& v) const
    {
        for(int h=Hash_Index(v);;h=Next_Index(h)){
            if(table(h).state==ENTRY_FREE) return PAIR<int,bool>(h,false);
            if(check_exists && table(h).state==ENTRY_ACTIVE && table(h).key==v) return PAIR<int,bool>(h,true);}
    }

    // pair.x = location in hash where this entry is
    // pair.y = true if the key was already present at pair.x
    // If check_exists=false, don't check for the key being present
    // In debug mode, checks anyway so it can assert.
    template<bool check_exists>
    PAIR<int,bool> Insert_Helper(const TK& v)
    {
        if(number_of_entries>next_resize) Resize_Table();
        auto pr=Search<check_exists||is_debug>(v);
        if(!check_exists && is_debug) assert(!pr.y);
        if(!check_exists || !pr.y){number_of_entries++;ENTRY& entry=table(pr.x);entry.key=v;entry.state=ENTRY_ACTIVE;}
        return pr;
    }

    void Insert(const HASHTABLE_ENTRY_TEMPLATE<TK,void>& entry)
    {Insert(entry.key);}

    template<class T2> typename enable_if<!is_same<T2,void>::value>::type
    Insert(const HASHTABLE_ENTRY_TEMPLATE<TK,T2>& entry)
    {Insert(entry.key,entry.data);}
public:

    void Insert(const TK& v) // assumes no entry with v exists
    {STATIC_ASSERT((is_same<T,void>::value));Insert_Helper<false>(v);}

    T_UNLESS_VOID& Insert(const TK& v,const T_UNLESS_VOID& value) // assumes no entry with v exists
    {STATIC_ASSERT((!is_same<T,void>::value));return table(Insert_Helper<false>(v).x).data=value;}

    // inserts the default if key not found
    T_UNLESS_VOID& Get_Or_Insert(const TK& v,const T_UNLESS_VOID& default_value=T_UNLESS_VOID())
    {auto pr=Insert_Helper<true>(v);T_UNLESS_VOID& data=table(pr.x).data;if(!pr.y) data=default_value;return data;}

    T* Get_Pointer(const TK& v) // returns NULL if key not found
    {auto pr=Search<true>(v);return pr.y?&table(pr.x).data:0;}

    const T* Get_Pointer(const TK& v) const // returns NULL if key not found
    {auto pr=Search<true>(v);return pr.y?&table(pr.x).data:0;}

    T_UNLESS_VOID& Get(const TK& v) // fails if key not found
    {auto pr=Search<true>(v);if(pr.y) return table(pr.x).data;throw KEY_ERROR("HASHTABLE::Get");}

    const T_UNLESS_VOID& Get(const TK& v) const // fails if key not found
    {auto pr=Search<true>(v);if(pr.y) return table(pr.x).data;throw KEY_ERROR("HASHTABLE::Get");}

    T Get_Default(const TK& v,const T_UNLESS_VOID& default_value=T_UNLESS_VOID()) const // returns default_value if key not found
    {auto pr=Search<true>(v);if(pr.y) return table(pr.x).data;return default_value;}

    bool Contains(const TK& v) const
    {return Search<true>(v).y;}

    bool Get(const TK& v,T_UNLESS_VOID& value) const
    {auto pr=Search<true>(v);if(pr.y) value=table(pr.x).data;return pr.y;}

    bool Set(const TK& v) // insert entry if doesn't already exists, returns whether it added a new entry
    {STATIC_ASSERT((is_same<T,void>::value));return !Insert_Helper<true>(v).y;}

    bool Set(const TK& v,const T_UNLESS_VOID& value) // if v doesn't exist insert value, else sets its value, returns whether it added a new entry
    {STATIC_ASSERT((!is_same<T,void>::value));auto pr=Insert_Helper<true>(v);if(!pr.y) table(pr.x).data=value;return !pr.y;}

    template<class T_ARRAY>
    void Set_All(const T_ARRAY& array)
    {STATIC_ASSERT(is_same<typename T_ARRAY::ELEMENT,TK>::value);
    for(typename T_ARRAY::INDEX i(0);i<array.Size();i++) Set(array(i));}

    template<class T_KEY,class T_ARRAY>
    void Set_All(const T_KEY& key,const T_ARRAY& array)
    {STATIC_ASSERT(is_same<typename T_ARRAY::ELEMENT,T>::value);
    STATIC_ASSERT(is_same<typename T_KEY::ELEMENT,TK>::value);
    for(typename T_ARRAY::INDEX i(0);i<array.Size();i++) Set(key(i),array(i));}

    bool Delete_If_Present(const TK& v)
    {auto pr=Search<true>(v);if(pr.y){table(pr.x).state=ENTRY_DELETED;number_of_entries--;next_resize--;}return pr.y;}

    void Delete(const TK& v)
    {if(!Delete_If_Present(v)) throw KEY_ERROR("HASHTABLE::Delete");}

    void Remove_All()
    {for(int i=0;i<table.m;i++) table(i).state=ENTRY_FREE;number_of_entries=0;}

    void Exchange(const TK& x,const TK& y) // Exchange values at entries x and y; valid if x or y (or both) are not present; efficient for array values
    {bool a=Contains(x),b=Contains(y);if(a || b){exchange(Get_Or_Insert(x),Get_Or_Insert(y));if(!a || !b){Delete(a?x:y);}}}

    void Exchange(HASHTABLE& hash)
    {table.Exchange(hash.table);exchange(number_of_entries,hash.number_of_entries);exchange(next_resize,hash.next_resize);}

    static void Exchange_Hashtables(HASHTABLE& hash1,HASHTABLE& hash2)
    {hash1.Exchange(hash2);}

    void Print_Table(std::ostream& output) const
    {output<<"Entry Count: "<<number_of_entries<<" Elements to resize at: "<<next_resize<<std::endl;
    for(int h=0;h<table.m;h++)
         output<<h<<":"<<(table(h).state==ENTRY_ACTIVE?"ACTIVE":(table(h).state==ENTRY_DELETED?"DELETED":"FREE"))<<" key="<<table(h).key<<" value="<<table(h).data<<std::endl;}

    template<class FUNC>
    void Map(FUNC function) // function(key,data)
    {for(int h=0;h<table.m;h++) if(table(h).state==ENTRY_ACTIVE) function(table(h).key,table(h).data);}

    void Append_Keys(ARRAY<TK>& keys) const
    {keys.Preallocate(Size());for(int h=0;h<table.m;h++) if(table(h).state==ENTRY_ACTIVE) keys.Append(table(h).key);}

    void Get_Keys(ARRAY<TK>& keys) const
    {keys.Remove_All();Append_Keys(keys);}

    void Append_Data(ARRAY<T_UNLESS_VOID>& data) const
    {data.Preallocate(data.m+Size());for(int h=0;h<table.m;h++) if(table(h).state==ENTRY_ACTIVE) data.Append(table(h).data);}

    void Get_Data(ARRAY<T_UNLESS_VOID>& data) const
    {data.Remove_All();Append_Data(data);}

    template<class RW> void Read(typename conditional<is_same<T,void>::value,std::istream&,UNUSABLE>::type input) // void version
    {int entries;Read_Binary<RW>(input,entries);Initialize_New_Table(entries);
    for(int i=0;i<entries;i++){TK key;Read_Binary<RW>(input,key);Insert(key);}}

    template<class RW> void Read(typename conditional<is_same<T,void>::value,UNUSABLE,std::istream&>::type input) // non-void version
    {int entries;Read_Binary<RW>(input,entries);Initialize_New_Table(entries);
    for(int i=0;i<entries;i++){TK key;T value;Read_Binary<RW>(input,key,value);Insert(key,value);}}

    template<class RW> void Write(typename conditional<is_same<T,void>::value,std::ostream&,UNUSABLE>::type output) const // void version
    {Write_Binary<RW>(output,number_of_entries);
    for(int h=0;h<table.m;h++) if(table(h).state==ENTRY_ACTIVE) Write_Binary<RW>(output,table(h).key);}

    template<class RW> void Write(typename conditional<is_same<T,void>::value,UNUSABLE,std::ostream&>::type output) const // non-void version
    {Write_Binary<RW>(output,number_of_entries);
    for(int h=0;h<table.m;h++) if(table(h).state==ENTRY_ACTIVE) Write_Binary<RW>(output,table(h).key,table(h).data);}

    template<bool is_const>
    class HASHTABLE_ITERATOR
    {
        typedef typename std::conditional<is_const,const HASHTABLE,HASHTABLE>::type IT_HASHTABLE;
        typedef typename std::conditional<is_const,const T_UNLESS_VOID,T_UNLESS_VOID>::type IT_T_UNLESS_VOID;
        typedef typename std::conditional<is_const,const ENTRY,ENTRY>::type IT_ENTRY;
    public:

        IT_HASHTABLE& hashtable;
        int index;
        enum INVALID_ENUM {END_ITERATOR};
        
        HASHTABLE_ITERATOR(IT_HASHTABLE& hashtable)
            :hashtable(hashtable),index(0)
        {
            Advance_Until_Valid();
        }

        HASHTABLE_ITERATOR(IT_HASHTABLE& hashtable,INVALID_ENUM)
            :hashtable(hashtable),index(hashtable.table.m)
        {}
        
        const TK& Key() const PHYSBAM_ALWAYS_INLINE
        {return hashtable.table(index).key;}

        IT_T_UNLESS_VOID& Data() PHYSBAM_ALWAYS_INLINE
        {return hashtable.table(index).data;}

        const T_UNLESS_VOID& Data() const PHYSBAM_ALWAYS_INLINE
        {return hashtable.table(index).data;}

        bool Valid() const PHYSBAM_ALWAYS_INLINE
        {return index<hashtable.table.m;}

        void Next() PHYSBAM_ALWAYS_INLINE
        {assert(Valid());index++;Advance_Until_Valid();}

        void Advance_Until_Valid()
        {for(;index<hashtable.table.m;index++) if(hashtable.table(index).state==ENTRY_ACTIVE) return;}

        void Prev() PHYSBAM_ALWAYS_INLINE
        {assert(Valid());index--;Reverse_Until_Valid();}

        void Reverse_Until_Valid()
        {for(;index>=0;index--) if(hashtable.table(index).state==ENTRY_ACTIVE) return;}

        // stl
        HASHTABLE_ITERATOR& operator++()
        {Next();return *this;}

        HASHTABLE_ITERATOR operator++(int)
        {HASHTABLE_ITERATOR it(*this);Next();return it;}

        // stl
        HASHTABLE_ITERATOR& operator--()
        {Prev();return *this;}

        HASHTABLE_ITERATOR operator--(int)
        {HASHTABLE_ITERATOR it(*this);Prev();return it;}

        auto* operator->()
        {return &Dereference_Operator((T*)0);}

        const auto* operator->() const
        {return &Dereference_Operator((T*)0);}

        auto& operator*()
        {return Dereference_Operator((T*)0);}

        const auto& operator*() const
        {return Dereference_Operator((T*)0);}

        const TK& Dereference_Operator(const void*) const
        {return Key();}

        const ENTRY& Dereference_Operator(const T_UNLESS_VOID*) const
        {return hashtable.table(index);}

        IT_ENTRY& Dereference_Operator(const T_UNLESS_VOID*)
        {return hashtable.table(index);}

        template<bool isc> bool operator==(const HASHTABLE_ITERATOR<isc>& it) const
        {return index==it.index;}

        template<bool isc> bool operator!=(const HASHTABLE_ITERATOR<isc>& it) const
        {return index!=it.index;}

        template<bool isc> bool operator<(const HASHTABLE_ITERATOR<isc>& it) const
        {return index<it.index;}

        template<bool isc> bool operator>(const HASHTABLE_ITERATOR<isc>& it) const
        {return index>it.index;}

        template<bool isc> bool operator<=(const HASHTABLE_ITERATOR<isc>& it) const
        {return index<=it.index;}

        template<bool isc> bool operator>=(const HASHTABLE_ITERATOR<isc>& it) const
        {return index>=it.index;}
    };

    typedef HASHTABLE_ITERATOR<false> ITERATOR;
    typedef HASHTABLE_ITERATOR<true> CONST_ITERATOR;
    typedef ITERATOR iterator; // for stl
    typedef CONST_ITERATOR const_iterator; // for stl

    // stl
    iterator begin()
    {return iterator(*this);}

    const_iterator begin() const
    {return const_iterator(*this);}

    iterator end()
    {return iterator(*this,iterator::END_ITERATOR);}

    const_iterator end() const
    {return const_iterator(*this,const_iterator::END_ITERATOR);}
};
template<class K,class T>
std::ostream& operator<<(std::ostream& output,const HASHTABLE<K,T>& hashtable)
{
    output<<"(";
    bool first=true;
    for(const auto& it:hashtable){
        if(!first) output<<" ";
        first=false;
        output<<it.key<<":"<<it.data;}
    output<<")";
    return output;
}

template<class K>
std::ostream& operator<<(std::ostream& output,const HASHTABLE<K,void>& hashtable)
{
    output<<"(";
    bool first=true;
    for(const auto& data:hashtable){
        if(!first) output<<" ";
        first=false;
        output<<data;}
    output<<")";
    return output;
}
template<class T,class T_ARRAY>
void Get_Unique(ARRAY<T>& uniq,const ARRAY_BASE<T,T_ARRAY>& array)
{
    const T_ARRAY& self=array.Derived();
    HASHTABLE<T> hash(Value(self.Size())*3/2);
    uniq.Remove_All();
    for(int i=0;i<self.Size();i++)
        if(hash.Set(self(i)))
            uniq.Append(self(i));
}
template<class T>
void Prune_Duplicates(ARRAY<T>& array)
{
    HASHTABLE<T> hash(Value(array.Size())*3/2);
    int j=0;
    for(int i=0;i<array.Size();i++)
        if(hash.Set(array(i)))
            array(j++)=array(i);
    array.Resize(j);
}
}
#endif
