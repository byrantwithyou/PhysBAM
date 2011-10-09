extern "C"
{
    void set_int(physbam_object * parent, const char * attribute, int x);
    int get_int(const physbam_object * parent, const char * attribute);

    void set_float(physbam_object * parent, const char * attribute, float x);
    float get_float(const physbam_object * parent, const char * attribute);

    int get_array_length(const physbam_object * parent, const char * attribute);

    void set_int_array(physbam_object * parent, const char * attribute, const int * x, int length, int start = 0);
    void get_int_array(const physbam_object * parent, const char * attribute, int * x, int length, int start = 0);

    void set_float_array(physbam_object * parent, const char * attribute, const float * x, int length, int start = 0);
    void get_float_array(const physbam_object * parent, const char * attribute, float * x, int length, int start = 0);
}
