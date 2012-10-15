//#####################################################################
// Copyright 2012, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Log/LOG.h>
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

T arrow_delay=0;

TV2 Project(const TV3& X)
{return TV2(X.y-0.4*X.x,X.z-0.3*X.x);}

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

    /* 04 */ TV3( 0.0, 0.525, 0.0 ),
    /* 05 */ TV3( 1.0, 0.475, 0.0 ),
    /* 06 */ TV3( 0.0, 0.475, 1.0 ),
    /* 07 */ TV3( 1.0, 0.525, 1.0 ),

    /* 08 */ TV3( 0.0, 0.0, 0.525 ),
    /* 09 */ TV3( 1.0, 0.0, 0.475 ),
    /* 10 */ TV3( 0.0, 1.0, 0.475 ),
    /* 11 */ TV3( 1.0, 1.0, 0.525 ),
    
    /* 12 */ TV3( 0.0, 0.5, 0.5 ),
    /* 13 */ TV3( 1.0, 0.5, 0.5 ),

    /* 14 */ TV3( 0.5, 0.0, 0.5 ),
    /* 15 */ TV3( 0.5, 1.0, 0.5 ),

    /* 16 */ TV3( 0.5, 0.5, 0.0 ),
    /* 17 */ TV3( 0.5, 0.5, 1.0 ),
    
    /* 18 */ TV3( 0.5, 0.5, 0.5 )};

ARRAY<TV2> V;

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

void Emit_Begin(std::ofstream& fout,T x0,T y0,T x1,T y1,T unit)
{
    fout<<std::fixed;
    fout<<"\\documentclass{article}\n";
    fout<<"\\usepackage{pstricks}\n";
    fout<<"\\usepackage[margin=0cm,papersize={"<<(x1-x0)*unit*1.001<<"cm,"<<(y1-y0)*unit*1.001<<"cm}]{geometry}\n";
    fout<<"\n";
    fout<<"\\begin{document}%\n";
    fout<<"\\psset{unit="<<unit<<"cm}%\n";
    fout<<"\\noindent%\n";
    fout<<"\\begin{pspicture}("<<x0<<","<<y0<<")("<<x1<<","<<y1<<")\n";
    fout<<"\\newrgbcolor{color0}{0.95 0 0}";
    fout<<"\\newrgbcolor{color1}{0 0.75 0}";
    fout<<"\\newrgbcolor{color2}{0 0.25 1}";
    fout<<"\\newrgbcolor{color3}{.7 .7 0}";
    fout<<"\n";
}    

void Emit_End(std::ofstream& fout)
{    
    fout<<"\n";
    fout<<"\\end{pspicture}\n";
    fout<<"\\end{document}\n";
}

void Emit_Arrow(const TV2& xp,T linewidth,std::ofstream& fout)
{    
    T length=0.15;

    fout<<"\\psline[linecolor=black,";
    fout<<"linewidth="<<linewidth;
    fout<<",arrowsize=0.06";
    fout<<",arrowlength=1";
    fout<<",arrowinset=.4";
    fout<<"]{c->}("<<xp.x-length<<","<<xp.y<<")("<<xp.x+length<<","<<xp.y<<")\n";
}

void Emit_Edge_Color(TV2 xp0,TV2 xp1,int color,bool visible,std::ofstream& fout)
{
//    TV2 c[2]={xp1,xp0};
    TV2 t=(xp1-xp0).Normalized();
    TV2 n=t.Rotate_Clockwise_90();
    T shift=0.0115;
    if(visible) n=-n;
    xp0+=shift*n;
    xp1+=shift*n;

    T radius=0.03;
    T linewidth=0.0175;
    std::string dash=",dash=0.015 0.01";

    T scale=1.4;
    T space=-0.015*scale;
    T length=0.04*scale;
    T width=0.025*scale;
    T back=0.005*scale;

    int spikes=1;

    for(int i=0;i<spikes;i++){
        TV2 a=xp1-(T)i*t*(length+space)-(radius+arrow_delay)*t;
        TV2 b=a-t*(length-back);
        TV2 c=a-t*length+n*width;
        fout<<"\\pspolygon[fillstyle=solid,linestyle=none,fillcolor=color"<<color<<",linecolor=color"<<color<<"]("<<a.x<<","<<a.y<<")("<<b.x<<","<<b.y<<")("<<c.x<<","<<c.y<<")\n";}

    fout<<"\\psline[linecolor=color"<<color<<",";
    if(visible) fout<<"linewidth="<<linewidth;
    else fout<<"linewidth="<<linewidth<<",linestyle=dashed"<<dash;
    fout<<"]{c-c}("<<xp0.x<<","<<xp0.y<<")("<<xp1.x<<","<<xp1.y<<")\n";
}

void Emit_Edge(TV2& xp0,TV2& xp1,int* color,bool visible,std::ofstream& fout)
{
    Emit_Edge_Color(xp1,xp0,color[0],visible,fout);
    Emit_Edge_Color(xp0,xp1,color[1],visible,fout);
}

void Emit_Point(TV2& xp,std::ofstream& fout,bool black=true)
{
    T r=0.03;

    if(black) fout<<"\\qdisk("<<xp.x<<","<<xp.y<<"){"<<r<<"}\n";
    else{
        fout<<"\\psset{linecolor=white}";
        fout<<"\\qdisk("<<xp.x<<","<<xp.y<<"){"<<r<<"}\n";
        fout<<"\\psset{linecolor=black}";
        fout<<"\\pscircle[linewidth=0.0075]("<<xp.x<<","<<xp.y<<"){"<<r<<"}\n";}
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
                Emit_Point(xp,fout);
                break;}
            case SEGMENT:{
                TV3 x0=X[v0+8];
                TV3 x1=X[v1+8];
                bool visible=Visible(x0,x1);
                TV2 xp0=Project(x0);
                TV2 xp1=Project(x1);

                int color[2];
                color[0]=c0;
                color[1]=c1;

                Emit_Edge(xp0,xp1,color,visible,fout);
                break;}
            case CUBE_POINT:{
                TV2 xp=Project(X[v0]);
                fout<<"\\psset{linecolor=color"<<c0<<"}";
                fout<<"\\qdisk("<<xp.x<<","<<xp.y<<"){0.035}\n";
                fout<<"\\psset{linecolor=black}";

                T alpha=M_PI/6;
                T start=0;//(alpha*c0)/2;
                T r=0.04625;
                int sectors=0;
                T offset=0;//-M_PI/2;

                for(int i=0;i<sectors;i++){
                    T phi=2*M_PI*i/sectors-start+offset;
                    for(int j=0;j<c0+1;j++){
                        fout<<"\\psline[linecolor=color"<<c0<<",";
                        fout<<"linewidth=0.0125";
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

    void Draw_Flat(std::ofstream& fout,bool black=true)
    {
        TV2 xp0=V(v0);
        TV2 xp1=V(v1);
        int color[2];
        color[0]=c0;
        color[1]=c1;
        Emit_Edge(xp0,xp1,color,true,fout);
        Emit_Point(xp0,fout,black);
        Emit_Point(xp1,fout,black);
    }

    bool Flip()
    {
        int tmp;
        tmp=v0;v0=v1;v1=tmp;
        tmp=c0;c0=c1;c1=tmp;
        return true;
    }
};

ARRAY<OBJECT> objects;

int main(int argc, char* argv[])
{
    PARSE_ARGS parse_args(argc,argv);
    std::string input_file_curves="in_curves.txt";
    std::string input_file_reduction="in_reduction.txt";
    std::string output_file="out";
    T scale=1,unit=1;
    bool EMIT_REDUCTION=false;
    bool EMIT_CURVES=false;
    bool EMIT_RULES=false;
    parse_args.Add("-ic",&input_file_curves,"file","input file with curves");
    parse_args.Add("-ir",&input_file_reduction,"file","input file with reduction");
    parse_args.Add("-o",&output_file,"file","output file");
    parse_args.Add("-scale",&scale,"scale","scale for graph reduction");
    parse_args.Add("-unit",&unit,"unit","choose units");
    parse_args.Add("-emit_curves",&EMIT_CURVES,"create output tex file with cube and curves");
    parse_args.Add("-emit_reduction",&EMIT_REDUCTION,"create output tex files with graph reduction process");
    parse_args.Add("-emit_rules",&EMIT_RULES,"create output tex files with reduction rules");
    parse_args.Parse();

    std::ifstream fin;
    std::ofstream fout;
    int v0,v1,v2,c0,c1;
    T x0,x1,y0,y1,z0,z1;
 
    // ##### CURVES ON CUBE #############################################################

    arrow_delay=.01*0;
    if(EMIT_CURVES||EMIT_REDUCTION)
    {
        fin.open(input_file_curves.c_str());
        fin>>x0>>y0>>x1>>y1;

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
        
        if(EMIT_CURVES)
        {
            LOG::cout<<"EMITTING CURVES"<<std::endl;
        
            fout.open((output_file+".tex").c_str());
            LOG::cout<<output_file+".tex"<<std::endl;
            Emit_Begin(fout,x0,y0,x1,y1,unit);
            for(int i=0;i<objects.m;i++){
                if(i>0 && objects(i)==objects(i-1)) continue;
                objects(i).Draw(fout);}
            Emit_End(fout);
            fout.close();
        }
        
        fin.close();
    }

    // ##### GRAPH REDUCTION #############################################################

    arrow_delay=.03*0;
    if(EMIT_REDUCTION)
    {
        LOG::cout<<"EMITTING GRAPH REDUCTION"<<std::endl;

        fin.open(input_file_reduction.c_str());
        fin>>x0>>y0>>x1>>y1;
        x0*=scale;
        y0*=scale;
        x1*=scale;
        y1*=scale;
    
        for(int i=0;i<18 && fin>>c0>>z0>>z1;i++){
            z0*=scale;
            z1*=scale;
            V.Append(TV2(z0,z1));}
        
        ARRAY<OBJECT> edges;
        for(int i=0;i<objects.m;i++)
            if(objects(i).type==OBJECT::SEGMENT)
                edges.Append(objects(i));
        
        fout.open((output_file+"_0.tex").c_str());
        LOG::cout<<output_file+"_0.tex"<<std::endl;
        Emit_Begin(fout,x0,y0,x1,y1,unit);
        for(int i=0;i<edges.m;i++) edges(i).Draw_Flat(fout);
        Emit_End(fout);
        fout.close();
        
        int count=0;
        while(fin>>v0>>v1>>v2>>c0>>c1){
            ARRAY<int> index;
            for(int i=0;i<edges.m;i++)
                if((edges(i).c0==c0 && edges(i).c1==c1) &&
                    (edges(i).v0==v0 || edges(i).v0==v1 || edges(i).v0==v2) &&
                    (edges(i).v1==v0 || edges(i).v1==v1 || edges(i).v1==v2)){
                    index.Append(i);
                    if(index.m==2) break;}
            assert(index.m==2);
            int a,b;
            if(edges(index(0)).v1==edges(index(1)).v0) a=index(0),b=index(1);
            else a=index(1),b=index(0);
            assert(edges(a).v1==edges(b).v0);
            // LOG::cout<<"remove "<<edges(a).v0<<" "<<edges(a).v1<<" "<<edges(a).c0<<" "<<edges(a).c1<<std::endl;
            // LOG::cout<<"remove "<<edges(b).v0<<" "<<edges(b).v1<<" "<<edges(b).c0<<" "<<edges(b).c1<<std::endl<<std::endl;
            edges(a).v1=edges(b).v1;
            edges.Remove_Index(b);
            if(a>b) a--;
            // LOG::cout<<"add "<<edges(a).v0<<" "<<edges(a).v1<<" "<<edges(a).c0<<" "<<edges(a).c1<<std::endl<<std::endl;
            for(b=0;b<edges.m;b++){
                if(a==b) continue;
                if((edges(a).v0==edges(b).v0 && edges(a).v1==edges(b).v1) ||
                    (edges(a).v0==edges(b).v1 && edges(a).v1==edges(b).v0 && edges(a).Flip())){
                    // LOG::cout<<"merge a "<<edges(a).v0<<" "<<edges(a).v1<<" "<<edges(a).c0<<" "<<edges(a).c1<<std::endl;
                    // LOG::cout<<"merge b "<<edges(b).v0<<" "<<edges(b).v1<<" "<<edges(b).c0<<" "<<edges(b).c1<<std::endl<<std::endl;
                    if(edges(a).c1==edges(b).c0) edges(a).c1=edges(b).c1;
                    else if(edges(a).c0==edges(b).c1) edges(a).c0=edges(b).c0;
                    edges.Remove_Index(b);
                    if(a>b) a--;
                    // LOG::cout<<"merge complete "<<edges(a).v0<<" "<<edges(a).v1<<" "<<edges(a).c0<<" "<<edges(a).c1<<std::endl<<std::endl;;                
                    if(edges(a).c0>edges(a).c1) edges(a).Flip();
                    if(edges(a).c0==edges(a).c1) edges.Remove_Index(a);
                    break;}}
            count++;
            if(edges.m){
                char buff[10];
                sprintf(buff,"%d",count);
                fout.open((output_file+"_"+buff+".tex").c_str());
                LOG::cout<<output_file+"_"+buff+".tex"<<std::endl;
                Emit_Begin(fout,x0,y0,x1,y1,unit);
                for(int i=0;i<edges.m;i++) edges(i).Draw_Flat(fout);
                Emit_End(fout);
                fout.close();}}
        fin.close();
    }

    // ##### REDUCTION RULES #############################################################

    if(EMIT_RULES)
    {
        LOG::cout<<"EMITTING REDUCTION RULES"<<std::endl;

        T height=0.2;
        T spacing_arrow=0.4;
        T spacing_arrow_tri=0.35;
        T spacing_lines=height/2;
        T ellipse_horizontal=0.12;
        T ellipse_vertical=0.07;
        T linewidth=0.01;
        
        {
            LOG::cout<<output_file+"_rule_A.tex"<<std::endl;
            fout.open((output_file+"_rule_A.tex").c_str());
            Emit_Begin(fout,x0,y0,x1,y1,unit);

            V.Remove_All();
            V.Append(TV2(0,-height));
            V.Append(TV2(0,height));
            OBJECT input={0,1,0,1,OBJECT::SEGMENT};
            input.Draw_Flat(fout);

            V.Remove_All();
            V.Append(TV2(2*spacing_arrow,-height));
            V.Append(TV2(2*spacing_arrow,height));
            OBJECT output={1,0,1,0,OBJECT::SEGMENT};
            output.Draw_Flat(fout);

            Emit_Arrow(TV2(spacing_arrow,0),linewidth,fout);

            Emit_End(fout);
            fout.close();
        }

        {
            LOG::cout<<output_file+"_rule_B.tex"<<std::endl;
            fout.open((output_file+"_rule_B.tex").c_str());
            Emit_Begin(fout,x0,y0,x1,y1,unit);

            V.Remove_All();
            V.Append(TV2(0,-height));
            V.Append(TV2(height*2/sqrt(3),height));
            V.Append(TV2(height*4/sqrt(3),-height));
            OBJECT input0={0,1,0,1,OBJECT::SEGMENT};
            OBJECT input1={1,2,0,1,OBJECT::SEGMENT};
            input0.Draw_Flat(fout);
            input1.Draw_Flat(fout);

            Emit_Arrow(TV2(height*4/sqrt(3)+spacing_arrow_tri,0),linewidth,fout);

            V.Remove_All();
            V.Append(TV2(height*4/sqrt(3)+2*spacing_arrow_tri,-height+spacing_lines));
            V.Append(TV2(height*4/sqrt(3)+2*spacing_arrow_tri+height*2/sqrt(3),height));
            V.Append(TV2(height*4/sqrt(3)+2*spacing_arrow_tri+height*4/sqrt(3),-height+spacing_lines));
            OBJECT output0={0,1,0,1,OBJECT::SEGMENT};
            OBJECT output1={1,2,0,1,OBJECT::SEGMENT};
            OBJECT output2={2,0,0,1,OBJECT::SEGMENT};
            output0.Draw_Flat(fout,false);
            output1.Draw_Flat(fout,false);
            output2.Draw_Flat(fout,false);

            V.Remove_All();
            V.Append(TV2(height*4/sqrt(3)+2*spacing_arrow_tri,-height));
            V.Append(TV2(height*4/sqrt(3)+2*spacing_arrow_tri+height*4/sqrt(3),-height));
            OBJECT output3={0,1,0,1,OBJECT::SEGMENT};
            output3.Draw_Flat(fout);

            Emit_End(fout);
            fout.close();
        }

        {
            LOG::cout<<output_file+"_rule_C.tex"<<std::endl;
            fout.open((output_file+"_rule_C.tex").c_str());
            Emit_Begin(fout,x0,y0,x1,y1,unit);

            V.Remove_All();
            V.Append(TV2(0,-height));
            V.Append(TV2(0,height));
            OBJECT input0={0,1,0,1,OBJECT::SEGMENT};
            input0.Draw_Flat(fout);

            V.Remove_All();
            V.Append(TV2(spacing_lines,-height));
            V.Append(TV2(spacing_lines,height));
            OBJECT input1={0,1,1,2,OBJECT::SEGMENT};
            input1.Draw_Flat(fout);

            fout<<"\\psellipse[bordercolor=white,border=0.0075,linewidth="<<linewidth<<"]("<<spacing_lines/2<<","<<height<<")("<<ellipse_horizontal<<","<<ellipse_vertical<<")\n";
            fout<<"\\psellipse[bordercolor=white,border=0.0075,linewidth="<<linewidth<<"]("<<spacing_lines/2<<","<<-height<<")("<<ellipse_horizontal<<","<<ellipse_vertical<<")\n";

            Emit_Arrow(TV2(spacing_lines+spacing_arrow,0),linewidth,fout);

            V.Remove_All();
            V.Append(TV2(2*spacing_arrow+spacing_lines,-height));
            V.Append(TV2(2*spacing_arrow+spacing_lines,height));
            OBJECT output={0,1,0,2,OBJECT::SEGMENT};
            output.Draw_Flat(fout);

            Emit_End(fout);
            fout.close();
        }

        {
            LOG::cout<<output_file+"_rule_D.tex"<<std::endl;
            fout.open((output_file+"_rule_D.tex").c_str());
            Emit_Begin(fout,x0,y0,x1,y1,unit);

            V.Remove_All();
            V.Append(TV2(0,-height));
            V.Append(TV2(0,height));
            OBJECT input={0,1,0,0,OBJECT::SEGMENT};
            input.Draw_Flat(fout);

            Emit_Arrow(TV2(spacing_arrow,0),linewidth,fout);

            T r=0.1;
            fout<<"\\pscircle[bordercolor=white,linewidth="<<linewidth<<"]("<<spacing_arrow+spacing_arrow_tri+r<<","<<0<<"){"<<r<<"}\n";
            fout<<"\\psline[linewidth="<<linewidth<<"]{c-c}("<<spacing_arrow+spacing_arrow_tri<<","<<-r<<")("<<spacing_arrow+spacing_arrow_tri+2*r<<","<<r<<")\n";

            Emit_End(fout);
            fout.close();
        }

        {
            LOG::cout<<output_file+"_rule_E.tex"<<std::endl;
            fout.open((output_file+"_rule_E.tex").c_str());
            Emit_Begin(fout,x0,y0,x1,y1,unit);

            V.Remove_All();
            V.Append(TV2(0,-height));
            V.Append(TV2(0,height));
            OBJECT input={0,1,0,1,OBJECT::SEGMENT};
            input.Draw_Flat(fout);

            V.Remove_All();
            V.Append(TV2(2*spacing_arrow,-height));
            V.Append(TV2(2*spacing_arrow,height));
            V.Append(TV2(2*spacing_arrow+height*sqrt(3),0));
            OBJECT output0={0,1,0,1,OBJECT::SEGMENT};
            OBJECT output1={1,2,0,1,OBJECT::SEGMENT};
            OBJECT output2={2,0,0,1,OBJECT::SEGMENT};
            output0.Draw_Flat(fout,false);
            output1.Draw_Flat(fout,false);
            output2.Draw_Flat(fout,false);
            fout<<"\\qdisk("<<V(2).x<<","<<V(2).y<<"){0.0125}\n";

            Emit_Arrow(TV2(spacing_arrow,0),linewidth,fout);

            Emit_End(fout);
            fout.close();
        }
    }
    
    return 0;
}
