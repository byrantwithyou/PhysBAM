#include <vector>
#include <cstdio>
#include <cassert>
#include <map>
#include <string>
#include <cstring>
#include <cstdlib>

template<int b>
struct basis_vector
{
    int data[b];
    bool operator<(const basis_vector& v) const
    {
        for(int i=0;i<b;i++)
            if(data[i]!=v.data[i])
                return data[i]<v.data[i];
        return false;
    }

    basis_vector operator+(const basis_vector& v) const
    {
        basis_vector r;
        for(int i=0;i<b;i++) r.data[i]=data[i]+v.data[i];
        return r;
    }

    basis_vector operator-() const
    {
        basis_vector r;
        for(int i=0;i<b;i++) r.data[i]=-data[i];
        return r;
    }
};

template<int b>
void print(const basis_vector<b>& v)
{
    for(int i=0;i<b;i++)
        printf("%i%s",v.data[i],i==b-1?"":",");
}

template<int m, int n, int b>
struct main_table
{
    int data[m][n];
    std::vector<basis_vector<b> > entries;
    std::map<basis_vector<b>,int> basis_index;
    bool register_neg=false;
 
    int get_hash(const basis_vector<b>& v)
    {
        auto it=basis_index.find(v);
        if(it!=basis_index.end()) return it->second;
        else
        {
            int h=entries.size();
            if(register_neg) basis_index[-v]=h+1;
            basis_index[v]=h;
            entries.push_back(v);
            if(register_neg) entries.push_back(-v);
            return h;
        }
    }
    
    void set(int i,int j,const basis_vector<b>& v)
    {
        data[i][j]=get_hash(v);
    }

    int hash_sum(int g, int h)
    {
        return get_hash(entries[g]+entries[h]);
    }
};

template<int b>
basis_vector<b> parse_vector(const char* s, const std::map<std::string,int>& basis_map)
{
    basis_vector<b> v={{}};
    if(*s=='0' && s[1]==0) return v;
    while(*s)
    {
        char * e=0;
        int num=strtol(s,&e,10);
        if(e==s)
        {
            if(*s=='-')
            {
                num=-1;
                s++;
            }
            else if(*s=='+')
            {
                num=1;
                s++;
            }
            else num=1;
        }
        else
        {
            assert(*e=='*');
            s=e+1;
        }
        assert(num);
        int basis_len=strcspn(s,"-+");
        auto it=basis_map.find(std::string(s,basis_len));
        assert(it!=basis_map.end());
        v.data[it->second]=num;
        s+=basis_len;
    }
    return v;
}

template<int m, int n, int b>
void parse_main_table(const char* maple, main_table<m,n,b>& table, const std::map<std::string,int>& basis_map)
{
    FILE* F = fopen(maple,"r");
    int c,len=0;
    char buff[10000];
    while((c=fgetc(F))!=-1) if(c=='[') break;
    while((c=fgetc(F))!=-1)
    {
        if(c==')' || c==';') break;
        if(!isspace(c) && c!='\\' && c!=']' && c!='[')
            buff[len++]=c;
    }
    buff[len++]=',';
    buff[len]=0;
    fclose(F);

    char* s=buff;
    for(int i=0;i<m;i++)
    {
        for(int j=0;j<n;j++)
        {
            assert(*s);
            char* e=strchr(s,',');
            *e=0;
            basis_vector<b> v=parse_vector<b>(s, basis_map);
            table.set(i,j,v);
            s=e+1;
        }
    }
    assert(!*s);
}

void print_array_rec(int* data, int* sizes, int num_sizes, int depth, int total_size, int num_flat)
{
    if(num_sizes==0)
    {
        printf("%i",*data);
        return;
    }
    int n=sizes[0];
    int stride=total_size/n;
    if(num_sizes>=num_flat) printf("%*s",4*depth,"");
    putchar('{');
    if(num_sizes>num_flat) putchar('\n');
    for(int i=0;i<n;i++)
    {
        print_array_rec(data, sizes+1, num_sizes-1, depth+1, stride, num_flat);
        if(i!=n-1)
        {
            putchar(',');
            if(num_sizes>num_flat) putchar('\n');
            data+=stride;
        }
    }
    if(num_sizes>num_flat) printf("\n%*s",4*depth,"");
    putchar('}');
}

void print_array(int* data, int* sizes, int num_sizes, const char* name, int num_flat=1)
{
    printf("static int %s",name);
    int total_size=1;
    for(int i=0;i<num_sizes;i++)
    {
        printf("[%i]",sizes[i]);
        total_size*=sizes[i];
    }
    printf("=\n");
    print_array_rec(data,sizes,num_sizes,0,total_size, num_flat);
    printf(";\n\n");
}


template<int m, int n, int b>
void print(const main_table<m,n,b>& table, const char* name)
{
    int sizes[]={m,n},uniq[]={(int)table.entries.size(),b};
    char buff[100];
    sprintf(buff,"main_table_%s",name);
    print_array((int*)table.data,sizes,2,buff);
    sprintf(buff,"unique_entries_%s",name);
    print_array((int*)&table.entries[0],uniq,2,buff);
}

int main(int argc, char* argv[])
{
    assert(argc==3);
    const char* visc_maple=argv[1];
    const char* pres_maple=argv[2];

    std::map<std::string,int> basis_visc,basis_pres;

    std::vector<std::string> vars;
    char buff[100];
    for(int j=1;j<=2;j++)
        for(int i=2;i<=3;i++)
        {
            sprintf(buff,"x%i%i",i,j);
            basis_pres[buff]=vars.size();
            vars.push_back(buff);
        }
    for(int k=0,i=0;i<4;i++)
        for(int j=i;j<4;j++)
        {
            sprintf(buff,"%s*%s",vars[i].c_str(),vars[j].c_str());
            basis_visc[buff]=k;
            sprintf(buff,"%s*%s",vars[j].c_str(),vars[i].c_str());
            basis_visc[buff]=k;
            if(i==j)
            {
                sprintf(buff,"%s^2",vars[i].c_str());
                basis_visc[buff]=k;
            }
            k++;
        }

    main_table<12,12,10> table_v;
    table_v.get_hash({});
    parse_main_table(visc_maple, table_v, basis_visc);

    int vertex_table_visc[2][2][7];
    int vertex_table_visc_sizes[3]={2,2,7};
    memset(vertex_table_visc,-1,sizeof(vertex_table_visc));

    for(int i=0;i<2;i++)
        for(int j=0;j<2;j++)
        {
            int z=table_v.data[6*i+0][6*j+0];
            int x=table_v.data[6*i+1][6*j+1];
            int y=table_v.data[6*i+2][6*j+2];
            int xy=table_v.hash_sum(x,y);
            int xyz=table_v.hash_sum(xy,z);
            int xyz2=table_v.hash_sum(xyz,xyz);
            vertex_table_visc[i][j][1]=z;
            vertex_table_visc[i][j][2]=xy;
            vertex_table_visc[i][j][3]=xyz;
            vertex_table_visc[i][j][6]=xyz2;
        }

    // missing-vertex-tri-0
    // missing-vertex-tri-1
    // {vertex A, vertex B, middle}-row
    // dim-row
    // {vertex A, vertex B, middle}-col
    // dim-col
    int edge_table_visc[3][3][3][2][3][2];
    int edge_table_visc_sizes[]={3,3,3,2,3,2};
    memset(edge_table_visc,-1,sizeof(edge_table_visc));

    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
        {
            int di[3]={(i+1)%3,(i+2)%3,i+3};
            int dj[3]={(j+2)%3,(j+1)%3,j+3};
            for(int m=0;m<3;m++)
                for(int k=0;k<2;k++)
                    for(int n=0;n<3;n++)
                        for(int l=0;l<2;l++)
                        {
                            if(m<2 && m==n) continue;
                            int u=table_v.data[6*k+di[m]][6*l+di[n]];
                            int v=table_v.data[6*k+dj[m]][6*l+dj[n]];
                            edge_table_visc[i][j][m][k][n][l]=table_v.hash_sum(u,v);
                        }
        }
    
    print(table_v,"visc");

    print_array((int*)vertex_table_visc, vertex_table_visc_sizes, 3, "vertex_table_visc");

    print_array((int*)edge_table_visc, edge_table_visc_sizes, 6, "edge_table_visc", 3);
    
    main_table<3,12,4> table_p;
    table_p.register_neg=true;
    table_p.get_hash({});
    parse_main_table(pres_maple, table_p, basis_pres);

    // axis, {no x, pos x, neg x}, {no z, pos z, neg z}
    int vertex_table_pres[2][3][3];
    int vertex_table_pres_sizes[]={2,3,3};
    memset(vertex_table_pres,-1,sizeof(vertex_table_pres));

    for(int i=0;i<2;i++)
    {
        int z=table_p.data[0][6*i+0];
        int x=table_p.data[1][6*i+1];
        int y=table_p.data[2][6*i+2];
        int xny=table_p.hash_sum(x,y^1);
        int xnyz=table_p.hash_sum(xny,z);
        int xnynz=table_p.hash_sum(xny,z^1);
        vertex_table_pres[i][0][0]=-1;
        vertex_table_pres[i][0][1]=z;
        vertex_table_pres[i][0][2]=z^1;
        vertex_table_pres[i][1][0]=xny;
        vertex_table_pres[i][1][1]=xnyz;
        vertex_table_pres[i][1][2]=xnynz;
        vertex_table_pres[i][2][0]=xny^1;
        vertex_table_pres[i][2][1]=xnynz^1;
        vertex_table_pres[i][2][2]=xnyz^1;
    }

    // missing-vertex-tri-0
    // missing-vertex-tri-1
    // {vertex A, vertex B}-row
    // {vertex A, vertex B, middle}-col
    // dim-col
    int edge_table_pres[3][3][2][3][2];
    int edge_table_pres_sizes[]={3,3,2,3,2};
    memset(edge_table_pres,-1,sizeof(edge_table_pres));

    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
        {
            int di[3]={(i+1)%3,(i+2)%3,i+3};
            int dj[3]={(j+2)%3,(j+1)%3,j+3};
            for(int m=0;m<2;m++)
                for(int n=0;n<3;n++)
                    for(int l=0;l<2;l++)
                    {
                        int u=table_p.data[di[m]][6*l+di[n]];
                        int v=table_p.data[dj[m]][6*l+dj[n]]^1;
                        edge_table_pres[i][j][m][n][l]=table_p.hash_sum(u,v);
                    }
        }
    
    print(table_p,"pres");

    print_array((int*)vertex_table_pres, vertex_table_pres_sizes, 3, "vertex_table_pres");

    print_array((int*)edge_table_pres, edge_table_pres_sizes, 5, "edge_table_pres", 3);
    
    return 0;
}
/*

my @reg_vert=(-1) x (7*4);
for $i (0..1)
{
    for $j (0..1)
    {
        my $z = $table[$hashes[6*$i + 6*$j*12]];
        my $x = $table[$hashes[6*$i + 6*$j*12 + 13]];
        my $y = $table[$hashes[6*$i + 6*$j*12 + 13*2]];
        my $xy = &add($x,$y);
        my $xyz = &add($xy,$z);
        my $xyz2 = &add($xyz,$xyz);
        $reg_vert[7*$i + 14*$j + 1] = &get_hash($z);
        $reg_vert[7*$i + 14*$j + 2] = &get_hash($xy);
        $reg_vert[7*$i + 14*$j + 3] = &get_hash($xyz);
        $reg_vert[7*$i + 14*$j + 6] = &get_hash($xyz2);
    }
}

my @edge_table=(-1) x (9*4);
for $i (0..2)
{
    for $j (0..2)
    {
        my $ia=($i+1)%3;
        my $ib=($i+2)%3;
        my $im=$i+3;
        my $ja=($j+1)%3;
        my $jb=($j+2)%3;
        my $jm=$j+3;
        my $vv=&get_hash(&add($table[$hashes[12*$ia+$ib]],$table[$hashes[12*$jb+$ja]]));
        $edge_table[4*3*$i+4*$j+0]=$vv;
        
    }
}

print "static int fem_unique_entry_table[$next_index][10]=\n{\n";
for $i (0..$#table)
{
    print "    {$table[$i]}";
    if($i<$#table){print ",";}
    print "\n";
}
print "};\n\n";

print "static int fem_regular_vertex[2][2][7]={\n";
for $i (0..1)
{
    print "    {\n";
    for $j (0..1)
    {
        print "        {";
        for $k (0..6)
        {
            my $a=$reg_vert[7*$i + 14*$j + $k];
            print "$a";
            if($k!=6){print ",";}
        }
        print "}";
        if($j!=1){print ",";}
        print "\n";
    }
    print "    }";
    if($i!=1){print ",";}
    print "\n";
}
print "};\n";
*/
