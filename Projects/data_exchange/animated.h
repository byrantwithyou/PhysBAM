#ifndef __DATA_EXCHANGE_ANIMATED__
#define __DATA_EXCHANGE_ANIMATED__

#include <vector>

namespace data_exchange{

template<class T>
struct animated
{
    float uniform_dt;
    float start_time;
    std::vector<float> times;
    std::vector<T> values;

    animated(): uniform_dt(0), start_time(0) {}
    animated(float start, float dt): uniform_dt(dt), start_time(start) {}

    void uniform_sampling(float start, float dt)
    {
        assert(times.empty());
        uniform_dt = dt;
        start_time = start;
    }

    void insert_uniform(const T& value)
    {
        assert(times.empty());
        values.push_back(value);
    }

    void insert_nonuniform(float time, const T& value)
    {
        assert(uniform_dt==0);
        times.push_back(time);
        values.push_back(value);
    }

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & uniform_dt & start_time & times & values;
    }
};

}

#endif
