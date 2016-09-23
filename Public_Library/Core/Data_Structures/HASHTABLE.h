//#####################################################################
// Copyright 2002-2008, Robert Bridson, Kevin Der, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra, Neil Molino, Craig Schroeder, Andrew Selle, Tamar Shinar, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HASHTABLE
//#####################################################################
#ifndef __HASHTABLE__
#define __HASHTABLE__

#include <Core/Arrays/ARRAY.h>
#include <Core/Arrays/PROJECTED_ARRAY.h>
#include <Core/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <Core/Data_Structures/HASH_REDUCE.h>
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
    typedef FIELD_PROJECTOR<HASHTABLE_ENTRY_TEMPLATE<TK,T>,HASHTABLE_ENTRY_STATE,&HASHTABLE_ENTRY_TEMPLATE<TK,T>::state> T_FIELD_PROJECTOR;
    typedef int HAS_UNTYPED_READ_WRITE;

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

    bool Contains(const TK& v,const int h) const
    {for(int i=h;;i=Next_Index(i))
        if(table(i).state==ENTRY_FREE) return false;
        else if(table(i).state==ENTRY_ACTIVE && table(i).key==v) return true;}

    void Insert(const HASHTABLE_ENTRY_TEMPLATE<TK,void>& entry)
    {Insert(entry.key);}

    template<class T2> typename enable_if<!is_same<T2,void>::value>::type
    Insert(const HASHTABLE_ENTRY_TEMPLATE<TK,T2>& entry)
    {Insert(entry.key,entry.data);}

public:

    void Insert(const TK& v) // assumes no entry with v exists
    {STATIC_ASSERT((is_same<T,void>::value));
    if(number_of_entries>next_resize) Resize_Table();
    number_of_entries++;
    int h=Hash_Index(v);assert(!Contains(v,h));
    for(;table(h).state!=ENTRY_FREE;h=Next_Index(h));
    ENTRY& entry=table(h);entry.key=v;entry.state=ENTRY_ACTIVE;}

    T_UNLESS_VOID& Insert(const TK& v,const T_UNLESS_VOID& value) // assumes no entry with v exists
    {STATIC_ASSERT((!is_same<T,void>::value));
    if(number_of_entries>next_resize) Resize_Table();
    number_of_entries++;
    int h=Hash_Index(v);assert(!Contains(v,h));
    for(;table(h).state!=ENTRY_FREE;h=Next_Index(h));
    ENTRY& entry=table(h);entry.key=v;entry.state=ENTRY_ACTIVE;entry.data=value;return entry.data;}

    T_UNLESS_VOID& Get_Or_Insert(const TK& v,const T_UNLESS_VOID& default_value=T_UNLESS_VOID()) // inserts the default if key not found
    {int h=Hash_Index(v);
    for(;table(h).state!=ENTRY_FREE;h=Next_Index(h))
        if(table(h).state==ENTRY_ACTIVE && table(h).key==v) return table(h).data;
    if(number_of_entries>next_resize){return Insert(v, default_value);}
    number_of_entries++;
    ENTRY& entry=table(h);entry.key=v;entry.state=ENTRY_ACTIVE;entry.data=default_value;return entry.data;}

    T* Get_Pointer(const TK& v) // returns NULL if key not found
    {for(int h=Hash_Index(v);table(h).state!=ENTRY_FREE;h=Next_Index(h))
        if(table(h).state==ENTRY_ACTIVE && table(h).key==v) return &table(h).data;
    return 0;}

    const T* Get_Pointer(const TK& v) const // returns NULL if key not found
    {for(int h=Hash_Index(v);table(h).state!=ENTRY_FREE;h=Next_Index(h))
        if(table(h).state==ENTRY_ACTIVE && table(h).key==v) return &table(h).data;
    return 0;}

    T_UNLESS_VOID& Get(const TK& v) // fails if key not found
    {if(T_UNLESS_VOID* data=Get_Pointer(v)) return *data;throw KEY_ERROR("HASHTABLE::Get");}

    const T_UNLESS_VOID& Get(const TK& v) const // fails if key not found
    {if(const T_UNLESS_VOID* data=Get_Pointer(v)) return *data;throw KEY_ERROR("HASHTABLE::Get");}

    T Get_Default(const TK& v,const T_UNLESS_VOID& default_value=T_UNLESS_VOID()) const // returns default_value if key not found
    {if(const T* data=Get_Pointer(v)) return *data;return default_value;}

    bool Contains(const TK& v) const
    {return Contains(v,Hash_Index(v));}

    bool Get(const TK& v,T_UNLESS_VOID& value) const
    {for(int h=Hash_Index(v);table(h).state!=ENTRY_FREE;h=Next_Index(h))
         if(table(h).state==ENTRY_ACTIVE && table(h).key==v){value=table(h).data;return true;}
    return false;}

    bool Set(const TK& v) // insert entry if doesn't already exists, returns whether it added a new entry
    {STATIC_ASSERT((is_same<T,void>::value));
    if(number_of_entries>next_resize) Resize_Table();
    int h=Hash_Index(v);
    for(;table(h).state!=ENTRY_FREE;h=Next_Index(h))
        if(table(h).state==ENTRY_ACTIVE && table(h).key==v) return false;
    number_of_entries++;
    ENTRY& entry=table(h);entry.key=v;entry.state=ENTRY_ACTIVE;return true;}

    bool Set(const TK& v,const T_UNLESS_VOID& value) // if v doesn't exist insert value, else sets its value, returns whether it added a new entry
    {STATIC_ASSERT((!is_same<T,void>::value));
    if(number_of_entries>next_resize) Resize_Table(); // if over load average, have to grow (must do this before computing hash index)
    int h=Hash_Index(v);
    for(;table(h).state!=ENTRY_FREE;h=Next_Index(h))
        if(table(h).state==ENTRY_ACTIVE && table(h).key==v){table(h).data=value;return false;}
    number_of_entries++;
    ENTRY& entry=table(h);entry.key=v;entry.state=ENTRY_ACTIVE;entry.data=value;return true;}

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
    {for(int h=Hash_Index(v);table(h).state!=ENTRY_FREE;h=Next_Index(h)) // reduce as still are using entries for deletions
        if(table(h).state==ENTRY_ACTIVE && table(h).key==v){table(h).state=ENTRY_DELETED;number_of_entries--;next_resize--;return true;}
    return false;}

    void Delete(const TK& v)
    {if(!Delete_If_Present(v)) throw KEY_ERROR("HASHTABLE::Delete");}

    void Remove_All()
    {PROJECTED_ARRAY<ARRAY<HASHTABLE_ENTRY_TEMPLATE<TK,T>,int>,T_FIELD_PROJECTOR> projected_array=table.template Project<HASHTABLE_ENTRY_STATE,&HASHTABLE_ENTRY_TEMPLATE<TK,T>::state>();
        projected_array.Fill(ENTRY_FREE);number_of_entries=0;}
    //{table.template Project<HASHTABLE_ENTRY_STATE,&HASHTABLE_ENTRY_TEMPLATE<TK,T>::state>().Fill(ENTRY_FREE);number_of_entries=0;}

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

    void Apply_Function_To_All_Entries(void (*function)(TK&,T_UNLESS_VOID&))
    {for(int h=0;h<table.m;h++) if(table(h).state==ENTRY_ACTIVE) function(table(h).key,table(h).data);}

    void Delete_Pointers_Stored_In_Table() // of course, only valid if pointers are stored in table
    {for(int h=0;h<table.m;h++) if(table(h).state==ENTRY_ACTIVE){delete table(h).data;table(h).data=0;}}

    void Reset_List_Arrays_Stored_In_Table() // of course, only works if pointers to ARRAY are stored in table
    {for(int h=0;h<table.m;h++) if(table(h).state==ENTRY_ACTIVE){table(h).data->Remove_All();}}

    void Append_Keys(ARRAY<TK>& keys) const
    {keys.Preallocate(Size());for(int h=0;h<table.m;h++) if(table(h).state==ENTRY_ACTIVE) keys.Append(table(h).key);}

    void Get_Keys(ARRAY<TK>& keys) const
    {keys.Remove_All();Append_Keys(keys);}

    template<class T_ARRAY0,class T_ARRAY1>
    void Get_Complementary_Keys(const T_ARRAY0& keys_universe,T_ARRAY1& keys_complementary) const
    {STATIC_ASSERT(is_same<typename T_ARRAY0::ELEMENT,TK>::value && is_same<typename T_ARRAY1::ELEMENT,TK>::value);
    keys_complementary.Remove_All();keys_complementary.Preallocate(keys_universe.Size()-Size());
    for(typename T_ARRAY0::INDEX i(0);i<keys_universe.Size();i++) if(!Contains(keys_universe(i))) keys_complementary.Append(keys_universe(i));}

    void Append_Data(ARRAY<T_UNLESS_VOID>& data) const
    {data.Preallocate(Size());for(int h=0;h<table.m;h++) if(table(h).state==ENTRY_ACTIVE) data.Append(table(h).data);}

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
        if(!first){output<<" ";first=false;}
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
        if(!first){output<<" ";first=false;}
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
