#ifndef __pack__
#define __pack__
#include <cstring>
#include <exception>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;

size_t pack_size(unsigned char c)
{
    return 1;
}

size_t pack_size(int i)
{
    return sizeof(int);
}

size_t pack_size(size_t i)
{
    return sizeof(size_t);
}

size_t pack_size(const string& s)
{
    return pack_size(s.size())+s.size();
}

template<class T>
size_t pack_size(const vector<T>& v)
{
    size_t k = pack_size(v.size());
    for(size_t i = 0; i < v.size(); i++)
        k += pack_size(v[i]);
    return k;
}

size_t pack(unsigned char * buff, unsigned char c)
{
    *buff=c;
    return 1;
}

size_t pack(unsigned char * buff, int i)
{
    *reinterpret_cast<int*>(buff) = i;
    return sizeof(int);
}

size_t pack(unsigned char * buff, size_t i)
{
    *reinterpret_cast<size_t*>(buff) = i;
    return sizeof(size_t);
}

size_t pack(unsigned char * buff, const string& s)
{
    size_t k = pack(buff,s.size());
    memcpy(buff+k,s.c_str(),s.size());
    return k+s.size();
}

template<class T>
size_t pack(unsigned char * buff, const vector<T>& v)
{
    size_t k = pack(buff,v.size());
    for(size_t i = 0; i < v.size(); i++)
        k += pack(buff + k, v[i]);
    return k;
}

size_t unpack(const unsigned char * buff, size_t size, unsigned char& c)
{
    if(size < 1) throw std::runtime_error("buffer too short");
    c = *buff;
    return 1;
}

size_t unpack(const unsigned char * buff, size_t size, int& i)
{
    if(size < sizeof(int)) throw std::runtime_error("buffer too short");
    i = *reinterpret_cast<const int*>(buff);
    return sizeof(int);
}

size_t unpack(const unsigned char * buff, size_t size, size_t& i)
{
    if(size < sizeof(size_t)) throw std::runtime_error("buffer too short");
    i = *reinterpret_cast<const size_t*>(buff);
    return sizeof(size_t);
}

size_t unpack(const unsigned char * buff, size_t size, string& s)
{
    size_t n, k=unpack(buff, size, n);
    if(size < k + n) throw std::runtime_error("buffer too short");
    s = string(buff + k, buff + k + n);
    return k + n;
}

template<class T>
size_t unpack(const unsigned char * buff, size_t size, vector<T>& v)
{
    size_t n, k = unpack(buff, size, n);
    v.resize(n);
    for(size_t i = 0; i < v.size(); i++)
        k += unpack(buff + k, size - k, v[i]);
    return k;
}

struct memory_layout
{
    size_t message_offset;
    size_t message_max_length;
    size_t message_length;
    int next_job_id;
    int priority;
    int filled;
};

#endif
