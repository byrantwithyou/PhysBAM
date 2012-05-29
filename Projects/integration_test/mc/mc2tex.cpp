//#####################################################################
// Copyright 2012, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <fstream>
#include <iostream>

using namespace PhysBAM;

typedef double T;
typedef VECTOR<T,2> TV2;
typedef VECTOR<T,3> TV3;
typedef VECTOR<int,2> TV_INT2;
typedef VECTOR<int,3> TV_INT3;

TV2 Project(const TV3& X)
{return TV2(X.y-0.35*X.x,X.z-0.25*X.x);}

TV3 X[27]={

    /* CUBE POINTS */ 
    
    /* 0 */ TV3( 0, 0, 0 ),
    /* 0 */ TV3( 1, 0, 0 ),
    /* 0 */ TV3( 0, 1, 0 ),
    /* 0 */ TV3( 1, 1, 0 ),
    /* 0 */ TV3( 0, 0, 1 ),
    /* 0 */ TV3( 1, 0, 1 ),
    /* 0 */ TV3( 0, 1, 1 ),
    /* 0 */ TV3( 1, 1, 1 ),

    /* POINTS (+8) */ 

    /* 00 */ TV3( 0.5, 0.0, 0.0 ),
    /* 01 */ TV3( 0.5, 1.0, 0.0 ),
    /* 02 */ TV3( 0.5, 0.0, 1.0 ),
    /* 03 */ TV3( 0.5, 1.0, 1.0 ),

    /* 04 */ TV3( 0.0, 0.5, 0.0 ),
    /* 05 */ TV3( 1.0, 0.5, 0.0 ),
    /* 06 */ TV3( 0.0, 0.5, 1.0 ),
    /* 07 */ TV3( 1.0, 0.5, 1.0 ),

    /* 08 */ TV3( 0.0, 0.0, 0.475 ),
    /* 09 */ TV3( 1.0, 0.0, 0.525 ),
    /* 10 */ TV3( 0.0, 1.0, 0.525 ),
    /* 11 */ TV3( 1.0, 1.0, 0.475 ),
    
    /* 12 */ TV3( 0.0, 0.5, 0.5 ),
    /* 13 */ TV3( 1.0, 0.5, 0.5 ),

    /* 14 */ TV3( 0.5, 0.0, 0.5 ),
    /* 15 */ TV3( 0.5, 1.0, 0.5 ),

    /* 16 */ TV3( 0.5, 0.5, 0.0 ),
    /* 17 */ TV3( 0.5, 0.5, 1.0 ),
    
    /* 18 */ TV3( 0.5, 0.5, 0.5 )};

bool Visible(const TV3& x0,const TV3& x1)
{
    if(x0.x==1 && x1.x==1) return true;
    if(x0.y==1 && x1.y==1) return true;
    if(x0.z==1 && x1.z==1) return true;
    if(x0.x==0 && x1.x==0) return false;
    if(x0.y==0 && x1.y==0) return false;
    if(x0.z==0 && x1.z==0) return false;
    PHYSBAM_FATAL_ERROR();
}

struct OBJECT
{
    enum WORKAROUND{CENTER_SEGMENT,CUBE_SEGMENT,SEGMENT,POINT,CUBE_POINT};

    int v0,v1,c0,c1;
    int type;
    
    TV3 Centroid() const
    {
        switch(type){
            case POINT: return X[v0+8];
            case SEGMENT:case CENTER_SEGMENT: return (X[v0+8]+X[v1+8])/2;
            case CUBE_POINT: return X[v0];
            case CUBE_SEGMENT: return (X[v0]+X[v1])/2;
            default: PHYSBAM_FATAL_ERROR();}
    }
    
    bool operator<(const OBJECT& s) const
    {
        T c=Centroid().x;
        T sc=s.Centroid().x;
        if(c<sc) return true;
        if(c>sc) return false;
        return(type<s.type);
    }

    bool operator==(const OBJECT& s) const
    {return (v0==s.v0 && v1==s.v1 && c0==s.c0 && c1==s.c1 && type==s.type);}

    void Draw(std::ofstream& fout)
    {
        switch(type){
            case POINT:{
                TV2 xp=Project(X[v0+8]);
                fout<<"\\qdisk("<<xp.x<<","<<xp.y<<"){0.025}\n";
                break;}
            case SEGMENT:{
                TV3 x0=X[v0+8];
                TV3 x1=X[v1+8];
                TV2 xp0=Project(x0);
                TV2 xp1=Project(x1);
                TV2 c=(xp1+xp0)/2;
                bool visible=Visible(x0,x1);
                TV2 t=xp1-xp0;
                t.Normalize();
                TV2 n=t.Rotate_Clockwise_90();

                T shift=0.0115;
                T linewidth=0.0175;
                std::string dash=",dash=0.015 0.01";

                T space=-0.005;
                T length=0.04;
                T width=0.025;
                T back=0.005;

                T sign[2];
                if(visible) {sign[0]=-1;sign[1]=1;}
                else {sign[0]=1;sign[1]=-1;}

                int color[2];
                color[0]=c0;
                color[1]=c1;

                for(int s=0;s<2;s++){
                    int spikes=color[s]+1;
                    int spaces=color[s];
                    T start=(spikes*length+spaces*space)/2;

                for(int i=0;i<spikes;i++){
                    fout<<"\\pspolygon[fillstyle=solid,fillcolor=color"<<color[s]<<",linecolor=color"<<color[s]<<"]";
                    fout<<"("<<c.x+sign[s]*shift*n.x+start*t.x-i*(space+length)*t.x<<","
                             <<c.y+sign[s]*shift*n.y+start*t.y-i*(space+length)*t.y<<")";
                    fout<<"("<<c.x+sign[s]*shift*n.x+start*t.x-i*(space+length)*t.x-(length-back)*t.x<<","
                             <<c.y+sign[s]*shift*n.y+start*t.y-i*(space+length)*t.y-(length-back)*t.y<<")";
                    fout<<"("<<c.x+sign[s]*shift*n.x+start*t.x-i*(space+length)*t.x-length*t.x+sign[s]*width*n.x<<","
                             <<c.y+sign[s]*shift*n.y+start*t.y-i*(space+length)*t.y-length*t.y+sign[s]*width*n.y<<")";}

                fout<<"\\psline[linecolor=color"<<color[s]<<",";
                if(Visible(x0,x1)) fout<<"linewidth="<<linewidth;
                else fout<<"linewidth="<<linewidth<<",linestyle=dashed"<<dash;
                fout<<"]{c-c}("<<xp0.x+sign[s]*shift*n.x<<","<<xp0.y+sign[s]*shift*n.y<<")("<<xp1.x+sign[s]*shift*n.x<<","<<xp1.y+sign[s]*shift*n.y<<")\n";}

                break;}
            case CUBE_POINT:{
                TV2 xp=Project(X[v0]);
                fout<<"\\psset{linecolor=color"<<c0<<"}";
                fout<<"\\qdisk("<<xp.x<<","<<xp.y<<"){0.03}\n";
                fout<<"\\psset{linecolor=black}";

                T alpha=0.425;
                T start=(alpha*c0)/2;
                T r=0.05;
                int sectors=3;
                T offset=-M_PI/2;

                for(int i=0;i<sectors;i++){
                    T phi=2*M_PI*i/sectors-start+offset;
                    for(int j=0;j<c0+1;j++){
                        fout<<"\\psline[linecolor=color"<<c0<<",";
                        fout<<"linewidth=0.012";
                        fout<<"]{c-c}("<<xp.x<<","<<xp.y<<")("<<xp.x+cos(phi+j*alpha)*r<<","<<xp.y+sin(phi+j*alpha)*r<<")\n";}}

                break;}
            case CUBE_SEGMENT:{
                TV3 x0=X[v0];
                TV3 x1=X[v1];
                TV2 xp0=Project(x0);
                TV2 xp1=Project(x1);
                
                T linewidth=0.005;

                fout<<"\\psline[linecolor=black,";
                if(Visible(x0,x1)) fout<<"linewidth="<<linewidth;
                else fout<<"linewidth="<<linewidth<<",linestyle=dashed,dash=0.01 0.01";
                fout<<"]{c-c}("<<xp0.x<<","<<xp0.y<<")("<<xp1.x<<","<<xp1.y<<")\n";
                break;}
            case CENTER_SEGMENT:{
                TV3 x0=X[v0+8];
                TV3 x1=X[v1+8];
                TV2 xp0=Project(x0);
                TV2 xp1=Project(x1);
                
                T linewidth=0.003;

                fout<<"\\psline[linecolor=gray,";
                fout<<"linewidth="<<linewidth<<",linestyle=dashed,dash=0.005 0.005";
                fout<<"]{c-c}("<<xp0.x<<","<<xp0.y<<")("<<xp1.x<<","<<xp1.y<<")\n";
                
                T r=0.006;
                
                fout<<"\\qdisk("<<xp0.x<<","<<xp0.y<<"){"<<r<<"}\n";
                fout<<"\\qdisk("<<xp1.x<<","<<xp1.y<<"){"<<r<<"}\n";
                break;}
            default: PHYSBAM_FATAL_ERROR();}
    }
};

ARRAY<OBJECT> objects;

int main(int argc, char* argv[])
{
    PARSE_ARGS parse_args;
    parse_args.Add_String_Argument("-i","in.txt","input file");
    parse_args.Add_String_Argument("-o","out.tex","output file");
    parse_args.Parse(argc,argv);
    std::string input_file=parse_args.Get_String_Value("-i");
    std::string output_file=parse_args.Get_String_Value("-o");

    std::ifstream fin(input_file.c_str());
    std::ofstream fout(output_file.c_str());

    int v0,v1,c0,c1;
    for(int i=0;i<8 && fin>>c0;i++){
        OBJECT cp={i,-1,c0,-1,OBJECT::CUBE_POINT};objects.Append(cp);}
    while(fin>>v0>>v1>>c0>>c1){
        OBJECT s={v0,v1,c0,c1,OBJECT::SEGMENT};objects.Append(s);
        OBJECT p0={v0,-1,c0,c1,OBJECT::POINT};objects.Append(p0);
        OBJECT p1={v1,-1,c0,c1,OBJECT::POINT};objects.Append(p1);}
    
    OBJECT cs00={0,1,-1,-1,OBJECT::CUBE_SEGMENT};objects.Append(cs00);
    OBJECT cs01={2,3,-1,-1,OBJECT::CUBE_SEGMENT};objects.Append(cs01);
    OBJECT cs02={4,5,-1,-1,OBJECT::CUBE_SEGMENT};objects.Append(cs02);
    OBJECT cs03={6,7,-1,-1,OBJECT::CUBE_SEGMENT};objects.Append(cs03);

    OBJECT cs10={0,2,-1,-1,OBJECT::CUBE_SEGMENT};objects.Append(cs10);
    OBJECT cs11={1,3,-1,-1,OBJECT::CUBE_SEGMENT};objects.Append(cs11);
    OBJECT cs12={4,6,-1,-1,OBJECT::CUBE_SEGMENT};objects.Append(cs12);
    OBJECT cs13={5,7,-1,-1,OBJECT::CUBE_SEGMENT};objects.Append(cs13);

    OBJECT cs20={0,4,-1,-1,OBJECT::CUBE_SEGMENT};objects.Append(cs20);
    OBJECT cs21={1,5,-1,-1,OBJECT::CUBE_SEGMENT};objects.Append(cs21);
    OBJECT cs22={2,6,-1,-1,OBJECT::CUBE_SEGMENT};objects.Append(cs22);
    OBJECT cs23={3,7,-1,-1,OBJECT::CUBE_SEGMENT};objects.Append(cs23);

    OBJECT cs0={12,13,-1,-1,OBJECT::CENTER_SEGMENT};objects.Append(cs0);
    OBJECT cs1={14,15,-1,-1,OBJECT::CENTER_SEGMENT};objects.Append(cs1);
    OBJECT cs2={16,17,-1,-1,OBJECT::CENTER_SEGMENT};objects.Append(cs2);

    objects.Sort();
    
    fout<<std::fixed;

    fout<<"\\documentclass{article}\n";
    fout<<"\\usepackage{pstricks}\n";
    fout<<"\n";
    fout<<"\\begin{document}\n";
    fout<<"\\pagestyle{empty}\n";
    fout<<"\n";
    fout<<"\\begin{pspicture}(6,6)\n";
    fout<<"\\psset{unit=4cm}\n";
    fout<<"\\newrgbcolor{color0}{0.95 0 0}";
    fout<<"\\newrgbcolor{color1}{0 0.75 0}";
    fout<<"\\newrgbcolor{color2}{0 0.25 1}";
    fout<<"\n";
    
    for(int i=0;i<objects.m;i++){
        if(i>0 && objects(i)==objects(i-1)) continue;
        objects(i).Draw(fout);}
    
    fout<<"\n";
    fout<<"\\end{pspicture}\n";
    fout<<"\\end{document}\n";

    fin.close();
    fout.close();

    return 0;
}
