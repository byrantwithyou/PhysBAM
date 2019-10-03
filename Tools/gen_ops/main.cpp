#include <cassert>
#include <cstdio>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <regex>

using namespace std;

struct node;

enum object_type {invalid, scalar, vec, mat, sym_mat, diag_mat, upper_mat, id_mat};
object_type object_type_map[128];

// sym, diag, upper
FILE * out_file[2][2][2];
FILE * vec_file;
FILE * cur_file;

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

struct node_result
{
    string ret, defs;
    prec pr;
};

struct node
{
    vector<node*> A;
    virtual ~node(){for(auto p:A) delete p;}

    virtual node_result gen() const = 0;
    virtual bool must_equal(string& e, const string& s) const = 0;
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
    map<string,string> name_remap;

    string get_str(const string& name)
    {
        auto it=int_map.find(name);
        if(it!=int_map.end() && it->second>=0) return std::to_string(it->second);
        auto is=name_remap.find(name);
        if(is!=name_remap.end()) return is->second;
        return name;
    }

    ~func_info()
    {
        delete ret_expr;
    }
} func;

string print_element(const string& var, const vector<int>& v, const vector<std::string>& vs)
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
        sprintf(buff,"%s.array[%s]",var.c_str(),vs[0].c_str());
        return buff;
    }
    if(v.size()==2)
    {
        if(var=="I")
        {
            if(v[1]>=0 && v[0]>=0)
            {
                if(v[0]!=v[1]) return "0";
                return "1";
            }
            else
            {
                if(vs[0]==vs[1]) return "1";
                sprintf(buff,"(%s==%s)",vs[0].c_str(),vs[1].c_str());
                return buff;
            }
        }
        int m=func.int_map[t.size0];
        std::string m_str=func.get_str(t.size0);
        assert(m<0 || v[0]<m);
        if(t.type==mat)
        {
            if(v[1]>=0 && v[0]>=0 && m>=0) sprintf(buff,"%s.x[%i]",var.c_str(),v[1]*m+v[0]);
            else sprintf(buff,"%s(%s,%s)",var.c_str(),vs[0].c_str(),vs[1].c_str());
        }
        else if(t.type==sym_mat)
        {
            if(v[1]>=0 && v[0]>=0)
            {
                int a=v[0],b=v[1];
                if(b>a) swap(a,b);
                sprintf(buff,"%s.x[%i]",var.c_str(),(2*m-b-1)*b/2+a);
            }
            else sprintf(buff,"%s(%s,%s)",var.c_str(),vs[0].c_str(),vs[1].c_str());
        }
        else if(t.type==diag_mat)
        {
            if(v[1]>=0 && v[0]>=0)
            {
                if(v[0]!=v[1]) return "0";
                sprintf(buff,"%s.x[%s]",var.c_str(),vs[0].c_str());
            }
            else
            {
                if(vs[0]==vs[1]) sprintf(buff,"%s.x[%s]",var.c_str(),vs[0].c_str());
                else sprintf(buff,"(%s==%s?%s.x[%s]:0)",
                    vs[0].c_str(),vs[1].c_str(),var.c_str(),vs[0].c_str());
            }
        }
        else if(t.type==upper_mat)
        {
            if(v[1]>=0 && v[0]>=0)
            {
                if(v[0]>v[1]) return "0";
                sprintf(buff,"%s.x[%i]",var.c_str(),v[1]*(v[1]+1)/2+v[0]);
            }
            else
            {
                sprintf(buff,"(%s<=%s?%s(%s,%s):0)",
                    vs[0].c_str(),vs[1].c_str(),var.c_str(),vs[0].c_str(),vs[1].c_str());
            }
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

    virtual bool must_equal(string& e, const string& s) const override
    {
        if(op=="/") assert(!A[1]->must_equal(e,s));
        if(op=="~" || op=="/") return A[0]->must_equal(e,s);
        if(op=="+" || op=="-")
        {
            string f;
            bool a=A[0]->must_equal(f,s);
            bool b=A[1]->must_equal(e,s);
            if(a && b && e==f) return true;
            return false;
        }
        return A[0]->must_equal(e,s) || A[1]->must_equal(e,s);
    }
    
    virtual node_result gen() const override
    {
        if(op=="~")
        {
            assert(A.size()==1);
            auto p=A[0]->gen();
            if(p.pr<prec_pre) return {"-("+p.ret+")",p.defs,prec_pre};
            return {"-"+p.ret,p.defs,prec_pre};
        }
        else if(op=="+" || op=="-" || op=="*" || op=="/")
        {
            assert(A.size()==2);
            auto p=A[0]->gen();
            auto q=A[1]->gen();
            p.defs=q.defs=p.defs+q.defs;
            prec r=prec_add;
            if(op=="/")
            {
                if(p.ret=="0" || q.ret=="1") return p;
                r=prec_mul;
                if(p.pr<prec_mul) p.ret="("+p.ret+")";
                if(q.pr<prec_pre) q.ret="("+q.ret+")";
            }
            else if(op=="*")
            {
                if(p.ret=="0" || q.ret=="1") return p;
                if(q.ret=="0" || p.ret=="1") return q;
                r=prec_mul;
                if(p.pr<prec_mul) p.ret="("+p.ret+")";
                if(q.pr<prec_mul) q.ret="("+q.ret+")";
            }
            else if(op=="-")
            {
                if(q.ret=="0") return p;
                if(p.ret=="0")
                {
                    if(p.pr<prec_pre) return {"-("+p.ret+")",p.defs,prec_pre};
                    return {"-"+p.ret,p.defs,prec_pre};
                }
                if(q.pr<prec_mul) q.ret="("+q.ret+")";
            }
            else
            {
                if(p.ret=="0") return q;
                if(q.ret=="0") return p;
            }
            return {p.ret+op+q.ret,p.defs,r};
        }
        else
        {
            assert(false);
            return {"","",prec_add};
        }
    }
};

struct node_sum : node
{
    string var,size;

    virtual bool must_equal(string& e, const string& s) const override
    {
        return A[0]->must_equal(e,s);
    }
    
    virtual node_result gen() const override
    {
        int n=func.int_map[size];
        if(n>=0)
        {
            if(!n) return {"0","",prec_pre};
            func.int_map[var]=0;
            auto p=A[0]->gen();
            if(n==1) return p;
            for(int i=1;i<n;i++)
            {
                func.int_map[var]=i;
                auto q=A[0]->gen();
                if(p.ret=="0") p=q;
                else if(q.ret!="0") p.ret+="+"+q.ret;
                p.defs+=q.defs;
            }
            p.pr=prec_add;
            func.int_map.erase(var);
            return p;
        }
        else
        {
            string eq_st;
            if(A[0]->must_equal(eq_st,var) && func.int_map.find(eq_st)!=func.int_map.end())
            {
                func.name_remap[var]=eq_st;
                func.int_map[var]=-1;
                auto p=A[0]->gen();
                func.int_map.erase(var);
                func.name_remap.erase(var);
                return p;
            }

            func.int_map[var]=-1;
            auto p=A[0]->gen();
            func.int_map.erase(var);
            string loop="T sum_"+var+"=0;";
            loop+="for(int "+var+"=0;"+var+"<"+func.get_str(size)+";"+var+"++){"+
                p.defs+
                "sum_"+var+"+="+p.ret+";}\n";
            return {"sum_"+var,loop,prec_pre};
        }
    }
};

struct node_var : node
{
    string var;
    vector<string> index;

    virtual bool must_equal(string& e, const string& s) const override
    {
        const auto& t=func.var_map[var];
        if(t.type!=id_mat && t.type!=diag_mat) return false;
        assert(index.size()==2);
        if(index[0]==s) e=index[1];
        else if(index[1]==s) e=index[0];
        else return false;
        return true;
    }
    
    virtual node_result gen() const override
    {
        vector<int> indices;
        vector<string> mapped;
        for(const auto&s:index)
        {
            indices.push_back(func.int_map[s]);
            mapped.push_back(func.get_str(s));
        }
        return {print_element(var, indices, mapped), "", prec_pre};
    }
};

struct node_func : node
{
    string var;

    virtual bool must_equal(string& e, const string& s) const override
    {
        return false;
    }
    
    virtual node_result gen() const override
    {
        auto p=A[0]->gen();
        node_result r={var+"("+p.ret,p.defs,prec_pre};
        for(size_t i=1;i<A.size();i++)
        {
            auto q=A[i]->gen();
            r.ret+=","+q.ret;
            r.defs+=q.defs;
        }
        r.ret+=")";
        return r;
    }
};

// return-type function-name arg1 arg2 arg3 ...

regex re_arg("([A-Za-z0-9_]+)(&?):([A-Za-z0-9_]+)");
regex re_type("([A-Za-z0-9_]+)=([TVMSDU])(?:\\(([A-Za-z0-9_]+)(,([A-Za-z0-9_]+))?\\))?");
regex re_size("([A-Za-z0-9_]+)=\\(([0-9]),([0-9])\\)");
regex re_temp("([A-Za-z0-9_]+)=temp");
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

string type_to_str(const type_data& rt, map<string,int>& int_map)
{
    if(rt.type==scalar) return "T";
    std::string m_str=func.get_str(rt.size0);
    char buff[1000];
    if(rt.type==vec)
    {
        sprintf(buff,"VECTOR<T,%s>",m_str.c_str());
        return buff;
    }
    if(rt.type==sym_mat)
    {
        sprintf(buff,"SYMMETRIC_MATRIX<T,%s>",m_str.c_str());
        return buff;
    }
    if(rt.type==diag_mat)
    {
        sprintf(buff,"DIAGONAL_MATRIX<T,%s>",m_str.c_str());
        return buff;
    }
    if(rt.type==upper_mat)
    {
        sprintf(buff,"UPPER_TRIANGULAR_MATRIX<T,%s>",m_str.c_str());
        return buff;
    }
    assert(rt.type==mat);
    std::string n_str=func.get_str(rt.size1);
    if(m_str==n_str) sprintf(buff,"MATRIX<T,%s>",m_str.c_str());
    else sprintf(buff,"MATRIX<T,%s,%s>",m_str.c_str(),n_str.c_str());
    return buff;
}

void gen_op()
{
    const type_data& ret_type=func.types[func.ret_type];
    string ret_name=type_to_str(ret_type,func.int_map);
    bool is_template=false;
    fprintf(cur_file,"template<class T");
    for(const auto& i:func.int_map)
        if(i.second==-1)
        {
            fprintf(cur_file,",int %s",i.first.c_str());
            is_template=true;
        }
    fprintf(cur_file,"> inline %s\n%s(",
        ret_name.c_str(),
        func.op_name.c_str());
    for(size_t i=0;i<func.prototype.size();i++)
    {
        if(i) fprintf(cur_file,", ");
        const type_data& td=func.types[func.prototype[i].type];
        if(td.type!=scalar && !func.prototype[i].is_ref) fprintf(cur_file,"const ");
        fprintf(cur_file,"%s",type_to_str(td,func.int_map).c_str());
        if(td.type!=scalar || func.prototype[i].is_ref) fprintf(cur_file,"&");
        fprintf(cur_file," %s",func.prototype[i].name.c_str());
    }
    fprintf(cur_file,")\n{\n");
    if(func.ret_expr)
    {
        if(ret_type.type==scalar)
        {
            auto p=func.ret_expr->gen();
            if(p.defs!="") fprintf(cur_file,"%s",p.defs.c_str());
            fprintf(cur_file,"    return %s;\n",p.ret.c_str());
        }
        else if(!is_template)
        {
            bool comma=false;
            string defs="";
            string body="";
            int m=func.int_map[ret_type.size0];
            if(ret_type.type==vec)
            {
                for(int i=0;i<m;i++)
                {
                    func.int_map[func.ret_ind0]=i;
                    if(comma) body+=",";
                    else comma=true;
                    auto p=func.ret_expr->gen();
                    body+=p.ret;
                    defs+=p.defs;
                }
                func.int_map.erase(func.ret_ind0);
            }
            else if(ret_type.type==diag_mat)
            {
                for(int i=0;i<m;i++)
                {
                    func.int_map[func.ret_ind0]=i;
                    func.int_map[func.ret_ind1]=i;
                    if(comma) body+=",";
                    else comma=true;
                    auto p=func.ret_expr->gen();
                    body+=p.ret;
                    defs+=p.defs;
                }
                func.int_map.erase(func.ret_ind0);
                func.int_map.erase(func.ret_ind1);
            }
            else if(ret_type.type==sym_mat)
            {
                for(int j=0;j<m;j++)
                {
                    func.int_map[func.ret_ind1]=j;
                    for(int i=j;i<m;i++)
                    {
                        func.int_map[func.ret_ind0]=i;
                        if(comma) body+=",";
                        else comma=true;
                        auto p=func.ret_expr->gen();
                        body+=p.ret;
                        defs+=p.defs;
                    }
                }
                func.int_map.erase(func.ret_ind0);
                func.int_map.erase(func.ret_ind1);
            }
            else if(ret_type.type==upper_mat)
            {
                for(int j=0;j<m;j++)
                {
                    func.int_map[func.ret_ind1]=j;
                    for(int i=0;i<=j;i++)
                    {
                        func.int_map[func.ret_ind0]=i;
                        if(comma) body+=",";
                        else comma=true;
                        auto p=func.ret_expr->gen();
                        body+=p.ret;
                        defs+=p.defs;
                    }
                }
                func.int_map.erase(func.ret_ind0);
                func.int_map.erase(func.ret_ind1);
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
                        if(comma) body+=",";
                        else comma=true;
                        auto p=func.ret_expr->gen();
                        body+=p.ret;
                        defs+=p.defs;
                    }
                }
                func.int_map.erase(func.ret_ind0);
                func.int_map.erase(func.ret_ind1);
            }
            if(defs!="") fprintf(cur_file,"%s",defs.c_str());
            fprintf(cur_file,"    return %s(%s);\n",ret_name.c_str(),body.c_str());
        }
        else
        {
            fprintf(cur_file,"    %s r;\n",ret_name.c_str());
            std::string m_str=func.get_str(ret_type.size0);
            const char* i0=func.ret_ind0.c_str(),*i1=func.ret_ind1.c_str();
            if(ret_type.type==vec)
            {
                func.int_map[func.ret_ind0]=-1;
                auto p=func.ret_expr->gen();
                func.int_map.erase(func.ret_ind0);
                fprintf(cur_file,"    for(int %s=0;%s<%s;%s++){\n",i0,i0,m_str.c_str(),i0);
                if(p.defs!="") fprintf(cur_file,"%s",p.defs.c_str());
                fprintf(cur_file,"        r.array[%s]=%s;}\n",i0,p.ret.c_str());
            }
            else if(ret_type.type==diag_mat)
            {
                func.int_map[func.ret_ind0]=-1;
                func.int_map[func.ret_ind1]=-1;
                func.name_remap[func.ret_ind1]=func.ret_ind0;
                auto p=func.ret_expr->gen();
                func.int_map.erase(func.ret_ind0);
                func.int_map.erase(func.ret_ind1);
                func.name_remap.erase(func.ret_ind1);
                fprintf(cur_file,"    for(int %s=0;%s<%s;%s++){\n",i0,i0,m_str.c_str(),i0);
                if(p.defs!="") fprintf(cur_file,"%s",p.defs.c_str());
                fprintf(cur_file,"        r.x(%s)=%s;}\n",i0,p.ret.c_str());
            }
            else if(ret_type.type==sym_mat)
            {
                func.int_map[func.ret_ind0]=-1;
                func.int_map[func.ret_ind1]=-1;
                auto p=func.ret_expr->gen();
                fprintf(cur_file,"    for(int %s=0;%s<%s;%s++)\n",i0,i0,m_str.c_str(),i0);
                fprintf(cur_file,"    for(int %s=0;%s<%s;%s++){\n",i1,i1,i0,i1);
                if(p.defs!="") fprintf(cur_file,"%s",p.defs.c_str());
                fprintf(cur_file,"        r.Element_Lower(%s,%s)=%s;}\n",i0,i1,p.ret.c_str());
                func.int_map.erase(func.ret_ind0);
                func.int_map.erase(func.ret_ind1);
            }
            else if(ret_type.type==upper_mat)
            {
                func.int_map[func.ret_ind0]=-1;
                func.int_map[func.ret_ind1]=-1;
                auto p=func.ret_expr->gen();
                func.int_map.erase(func.ret_ind0);
                func.int_map.erase(func.ret_ind1);
                fprintf(cur_file,"    for(int %s=0;%s<%s;%s++)\n",i0,i0,m_str.c_str(),i0);
                fprintf(cur_file,"    for(int %s=%s;%s<%s;%s++){\n",i1,i0,i1,m_str.c_str(),i1);
                if(p.defs!="") fprintf(cur_file,"%s",p.defs.c_str());
                fprintf(cur_file,"        r(%s,%s)=%s;}\n",i0,i1,p.ret.c_str());
            }
            else
            {
                std::string n_str=func.get_str(ret_type.size1);
                func.int_map[func.ret_ind0]=-1;
                func.int_map[func.ret_ind1]=-1;
                auto p=func.ret_expr->gen();
                fprintf(cur_file,"    for(int %s=0;%s<%s;%s++)\n",i0,i0,m_str.c_str(),i0);
                func.int_map.erase(func.ret_ind0);
                func.int_map.erase(func.ret_ind1);
                fprintf(cur_file,"    for(int %s=0;%s<%s;%s++){\n",i1,i1,n_str.c_str(),i1);
                if(p.defs!="") fprintf(cur_file,"%s",p.defs.c_str());
                fprintf(cur_file,"        r(%s,%s)=%s;}\n",i0,i1,p.ret.c_str());
            }
            fprintf(cur_file,"    return r;\n");
        }
    }
    fprintf(cur_file,"}\n\n");
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
    func.int_map.erase(sizes[i].name);
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

    for(int s=0;s<2;s++)
        for(int d=0;d<2;d++)
            for(int u=0;u<2;u++)
            {
                std::string f = "SMALL_MATRIX_OPS_M";
                if(s) f+="S";
                if(d) f+="D";
                if(u) f+="U";
                f+=".h";
                out_file[s][d][u] = fopen(f.c_str(), "w");
            }
    vec_file = fopen("SMALL_VECTOR_OPS.h", "w");

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
            bool want_generic=false;
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
                else if(regex_match(token,m,re_temp))
                {
                    size_data t;
                    t.name=m.str(1);
                    t.a=-1;
                    t.b=-1;
                    sizes.push_back(t);
                }
                else if(token=="full")
                {
                    want_generic=true;
                }
                else
                {
                    fprintf(cur_file,"bad token: %s\n",token.c_str());
                    assert(false);
                }
            }
            for(const auto&i:func.prototype)
                func.var_map[i.name]=func.types[i.type];

            bool is_m=false,is_s=false,is_d=false,is_u=false;
            for(const auto& t:func.types)
            {
                if(t.second.type==mat) is_m=true;
                if(t.second.type==sym_mat) is_s=true;
                if(t.second.type==diag_mat) is_d=true;
                if(t.second.type==upper_mat) is_u=true;
            }
            if(!is_m && !is_s && !is_d && !is_u) cur_file=vec_file;
            else cur_file = out_file[is_s][is_d][is_u];

            gen_int_map(sizes,0);
            if(want_generic)
            {
                for(auto& i:sizes)
                {
                    i.a=-1;
                    i.b=-1;
                }
                gen_int_map(sizes,0);
            }
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
