#include <cstdio>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <regex>
#include <cassert>

using namespace std;

struct node;

enum object_type {invalid, scalar, vec, mat, sym_mat, diag_mat, upper_mat, id_mat};
object_type object_type_map[128];

enum prec {prec_add, prec_mul, prec_pre};

struct arg_data
{
    string type;
    string name;
    bool is_ref;
};
struct size_data
{
    string name;
    int a,b;
};
struct type_data
{
    string size0,size1;
    object_type type;
};

struct func_info;

struct node
{
    vector<node*> A;
    virtual ~node(){for(auto p:A) delete p;}

    virtual pair<string,prec> gen(func_info& func) const = 0;
};

struct func_info
{
    string ret_type;
    string op_name;
    vector<arg_data> prototype;
    node* ret_expr=0;
    string ret_ind0,ret_ind1;
    map<string,type_data> types;
    map<string,type_data> var_map;
    map<string,int> int_map;
    
    ~func_info()
    {
        delete ret_expr;
    }
};

string print_element(func_info& func, const string& var, const vector<int>& v)
{
    const auto& t=func.var_map[var];
    if(v.size()==0)
    {
        assert(t.type==scalar);
        return var;
    }
    char buff[1000];
    if(v.size()==1)
    {
        assert(t.type==vec);
        sprintf(buff,"%s.array[%i]",var.c_str(),v[0]);
        return buff;
    }
    if(v.size()==2)
    {
        if(var=="I")
        {
            if(v[0]!=v[1]) return "0";
            return "1";
        }
        int m=func.int_map[t.size0];
        assert(v[0]<m);
        if(t.type==mat) sprintf(buff,"%s.x[%i]",var.c_str(),v[1]*m+v[0]);
        else if(t.type==sym_mat)
        {
            int a=v[0],b=v[1];
            if(b>a) swap(a,b);
            sprintf(buff,"%s.x[%i]",var.c_str(),(2*m-b-1)*b/2+a);
        }
        else if(t.type==diag_mat)
        {
            if(v[0]!=v[1]) return "0";
            sprintf(buff,"%s.x[%i]",var.c_str(),v[0]);
        }
        else if(t.type==upper_mat)
        {
            if(v[0]>v[1]) return "0";
            sprintf(buff,"%s.x[%i]",var.c_str(),v[1]*(v[1]+1)/2+v[0]);
        }
        else
        {
            assert(false);
        }
        return buff;
    }
    if(v.size()==3)
    {
        if(var=="e")
        {
            if(v[0]==v[1] || v[0]==v[2] || v[1]==v[2] || v[0]>=3 || v[1]>=3|| v[2]>=3) return "0";
            if((v[0]==0 && v[1]==1) || (v[0]==1 && v[1]==2) || (v[0]==2 && v[1]==0)) return "1";
            return "-1";
        }
    }
    assert(false);
}

struct node_op : node
{
    string op;

    virtual pair<string,prec> gen(func_info& func) const override
    {
        if(op=="~")
        {
            assert(A.size()==1);
            auto p=A[0]->gen(func);
            if(p.second<prec_pre) return {"-("+p.first+")",prec_pre};
            return {"-"+p.first,prec_pre};
        }
        else if(op=="+" || op=="-" || op=="*" || op=="/")
        {
            assert(A.size()==2);
            auto p=A[0]->gen(func);
            auto q=A[1]->gen(func);
            prec r=prec_add;
            if(op=="/")
            {
                if(p.first=="0" || q.first=="1") return p;
                r=prec_mul;
                if(p.second<prec_mul) p.first="("+p.first+")";
                if(q.second<prec_pre) q.first="("+q.first+")";
            }
            else if(op=="*")
            {
                if(p.first=="0" || q.first=="1") return p;
                if(q.first=="0" || p.first=="1") return q;
                r=prec_mul;
                if(p.second<prec_mul) p.first="("+p.first+")";
                if(q.second<prec_mul) q.first="("+q.first+")";
            }
            else if(op=="-")
            {
                if(q.first=="0") return p;
                if(p.first=="0")
                {
                    if(p.second<prec_pre) return {"-("+p.first+")",prec_pre};
                    return {"-"+p.first,prec_pre};
                }
                if(q.second<prec_mul) q.first="("+q.first+")";
            }
            else
            {
                if(p.first=="0") return q;
                if(q.first=="0") return p;
            }
            return {p.first+op+q.first,r};
        }
        else
        {
            assert(false);
            return {"",prec_add};
        }
    }
};

struct node_sum : node
{
    string var,size;

    virtual pair<string,prec> gen(func_info& func) const override
    {
        int n=func.int_map[size];
        if(!n) return {"0",prec_pre};
        func.int_map[var]=0;
        auto p=A[0]->gen(func);
        if(n==1) return p;
        for(int i=1;i<n;i++)
        {
            func.int_map[var]=i;
            auto q=A[0]->gen(func);
            if(p.first=="0") p=q;
            else if(q.first!="0") p.first+="+"+q.first;
        }
        p.second=prec_add;
        return p;
    }
};

struct node_var : node
{
    string var;
    vector<string> index;

    virtual pair<string,prec> gen(func_info& func) const override
    {
        vector<int> indices;
        for(const auto&s:index) indices.push_back(func.int_map[s]);
        return {print_element(func, var, indices), prec_pre};
    }
};

struct node_func : node
{
    string var;

    virtual pair<string,prec> gen(func_info& func) const override
    {
        auto p=A[0]->gen(func);
        pair<string,prec> r={var+"("+p.first,prec_pre};
        for(size_t i=1;i<A.size();i++)
            r.first+=","+A[i]->gen(func).first;
        r.first+=")";
        return r;
    }
};

// return-type function-name arg1 arg2 arg3 ...

regex re_arg("([A-Za-z0-9_]+)(&?):([A-Za-z0-9_]+)");
regex re_type("([A-Za-z0-9_]+)=([TVMSDU])(?:\\(([A-Za-z0-9_]+)(,([A-Za-z0-9_]+))?\\))?");
regex re_size("([A-Za-z0-9_]+)=\\(([0-9]),([0-9])\\)");
regex re_func("([A-Za-z0-9_]+)\\(([A-Za-z0-9_]+)(?:,([A-Za-z0-9_]+))?(?:,([A-Za-z0-9_]+))?\\)");
regex re_ret("return(?:\\(([A-Za-z0-9_]+)(?:,([A-Za-z0-9_]+))?\\))?");

node* parse_node(const vector<string>& v)
{
    vector<node*> st;
    for(const auto& s:v)
    {
        smatch m;
        if(s=="~")
        {
            auto N=new node_op;
            N->op=s;
            N->A.push_back(st.back());
            st.pop_back();
            st.push_back(N);
        }
        else if(s=="+" || s=="-" || s=="*" || s=="/")
        {
            auto N=new node_op;
            N->op=s;
            N->A.push_back(st[st.size()-2]);
            N->A.push_back(st.back());
            st.pop_back();
            st.pop_back();
            st.push_back(N);
        }
        else if(regex_match(s,m,re_func))
        {
            if(m.str(1)=="sum")
            {
                auto N=new node_sum;
                N->A.push_back(st.back());
                st.pop_back();
                N->var=m.str(2);
                N->size=m.str(3);
                st.push_back(N);
            }
            else
            {
                auto N=new node_var;
                N->var=m.str(1);
                if(m.str(2)!="") N->index.push_back(m.str(2));
                if(m.str(3)!="") N->index.push_back(m.str(3));
                if(m.str(4)!="") N->index.push_back(m.str(4));
                st.push_back(N);
            }
        }
        else
        {
            if(s=="exp" || s=="log" || s=="abs" || s=="sqrt")
            {
                auto N=new node_func;
                N->var=s;
                N->A.push_back(st.back());
                st.pop_back();
                st.push_back(N);
            }
            else if(s=="min" || s=="max")
            {
                auto N=new node_func;
                N->var=s;
                N->A.push_back(st[st.size()-2]);
                N->A.push_back(st.back());
                st.pop_back();
                st.pop_back();
                st.push_back(N);
            }
            else
            {
                auto N=new node_var;
                N->var=s;
                st.push_back(N);
            }
        }
    }
    assert(st.size()==1);
    return st[0];
}

func_info func;

string type_to_str(const type_data& rt, map<string,int>& int_map)
{
    if(rt.type==scalar) return "T";
    int m=int_map[rt.size0];
    char buff[1000];
    if(rt.type==vec)
    {
        sprintf(buff,"VECTOR<T,%i>",m);
        return buff;
    }
    if(rt.type==sym_mat)
    {
        sprintf(buff,"SYMMETRIC_MATRIX<T,%i>",m);
        return buff;
    }
    if(rt.type==diag_mat)
    {
        sprintf(buff,"DIAGONAL_MATRIX<T,%i>",m);
        return buff;
    }
    if(rt.type==upper_mat)
    {
        sprintf(buff,"UPPER_TRIANGULAR_MATRIX<T,%i>",m);
        return buff;
    }
    assert(rt.type==mat);
    int n=int_map[rt.size1];
    if(m==n) sprintf(buff,"MATRIX<T,%i>",m);
    else sprintf(buff,"MATRIX<T,%i,%i>",m,n);
    return buff;
}

void gen_op()
{
    const type_data& ret_type=func.types[func.ret_type];
    string ret_name=type_to_str(ret_type,func.int_map);
    printf("template<class T> inline %s\n%s(",
        ret_name.c_str(),
        func.op_name.c_str());
    for(size_t i=0;i<func.prototype.size();i++)
    {
        if(i) printf(", ");
        const type_data& td=func.types[func.prototype[i].type];
        if(td.type!=scalar && !func.prototype[i].is_ref) printf("const ");
        printf("%s",type_to_str(td,func.int_map).c_str());
        if(td.type!=scalar || func.prototype[i].is_ref) printf("&");
        printf(" %s",func.prototype[i].name.c_str());
    }
    printf(")\n{\n");
    if(func.ret_expr)
    {
        if(ret_type.type==scalar)
            printf("    return %s;\n",func.ret_expr->gen(func).first.c_str());
        else
        {
            bool comma=false;
            printf("    return %s(",ret_name.c_str());
            int m=func.int_map[ret_type.size0];
            if(ret_type.type==vec)
            {
                for(int i=0;i<m;i++)
                {
                    func.int_map[func.ret_ind0]=i;
                    if(comma) printf(",");
                    else comma=true;
                    printf("%s",func.ret_expr->gen(func).first.c_str());
                }
            }
            else if(ret_type.type==diag_mat)
            {
                for(int i=0;i<m;i++)
                {
                    func.int_map[func.ret_ind0]=i;
                    func.int_map[func.ret_ind1]=i;
                    if(comma) printf(",");
                    else comma=true;
                    printf("%s",func.ret_expr->gen(func).first.c_str());
                }
            }
            else if(ret_type.type==sym_mat)
            {
                for(int j=0;j<m;j++)
                {
                    func.int_map[func.ret_ind1]=j;
                    for(int i=j;i<m;i++)
                    {
                        func.int_map[func.ret_ind0]=i;
                        if(comma) printf(",");
                        else comma=true;
                        printf("%s",func.ret_expr->gen(func).first.c_str());
                    }
                }
            }
            else if(ret_type.type==upper_mat)
            {
                for(int j=0;j<m;j++)
                {
                    func.int_map[func.ret_ind1]=j;
                    for(int i=0;i<=j;i++)
                    {
                        func.int_map[func.ret_ind0]=i;
                        if(comma) printf(",");
                        else comma=true;
                        printf("%s",func.ret_expr->gen(func).first.c_str());
                    }
                }
            }
            else
            {
                int n=func.int_map[ret_type.size1];
                for(int j=0;j<n;j++)
                {
                    func.int_map[func.ret_ind1]=j;
                    for(int i=0;i<m;i++)
                    {
                        func.int_map[func.ret_ind0]=i;
                        if(comma) printf(",");
                        else comma=true;
                        printf("%s",func.ret_expr->gen(func).first.c_str());
                    }
                }
            }
            printf(");\n");
        }
    }
    printf("}\n\n");
}

void gen_int_map(const vector<size_data>& sizes,int i)
{
    if(i>=(int)sizes.size())
    {
        return gen_op();
    }
    for(int j=sizes[i].a;j<=sizes[i].b;j++)
    {
        func.int_map[sizes[i].name]=j;
        gen_int_map(sizes,i+1);
    }
}

int main(int argc, char* argv[])
{
    object_type_map['T']=scalar;
    object_type_map['V']=vec;
    object_type_map['M']=mat;
    object_type_map['S']=sym_mat;
    object_type_map['D']=diag_mat;
    object_type_map['U']=upper_mat;
    object_type_map['I']=id_mat;

    string line;
    while(getline(cin,line))
    {
        istringstream ss(line);
        string token;
        ss>>token;
        smatch m;
        if(token=="op")
        {
            func=func_info();
            func.types["T"]={"","",scalar};
            for(int i=0;i<10;i++)
            {
                char buff[3]={'-',(char)(i+'0')};
                func.int_map[buff+1]=i;
                func.var_map[buff+1]={"","",scalar};
                func.var_map[buff]={"","",scalar};
            }
            ss>>func.ret_type;
            ss>>func.op_name;
            while(ss>>token)
            {
                bool found=regex_match(token,m,re_arg);
                assert(found);
                arg_data ad;
                ad.type=m.str(1);
                ad.name=m.str(3);
                ad.is_ref=m.str(2)=="&";
                func.prototype.push_back(ad);
            }
            continue;
        }
        else if(token=="gen")
        {
            vector<size_data> sizes;
            while(ss>>token)
            {
                if(regex_match(token,m,re_type))
                {
                    type_data t;
                    t.type=object_type_map[(int)m.str(2).c_str()[0]];
                    t.size0=m.str(3);
                    t.size1=m.str(5);
                    func.types[m.str(1)]=t;
                }
                else if(regex_match(token,m,re_size))
                {
                    size_data t;
                    t.name=m.str(1);
                    t.a=atoi(m.str(2).c_str());
                    t.b=atoi(m.str(3).c_str());
                    sizes.push_back(t);
                }
                else
                {
                    printf("bad token: %s\n",token.c_str());
                    assert(false);
                }
            }
            for(const auto&i:func.prototype)
                func.var_map[i.name]=func.types[i.type];
            gen_int_map(sizes,0);
        }
        else if(regex_match(token,m,re_ret))
        {
            func.ret_ind0=m.str(1);
            func.ret_ind1=m.str(2);
            vector<string> v;
            while(ss>>token) v.push_back(token);
            func.ret_expr=parse_node(v);
        }
    }


    
    return 0;
}
