#include "dk.h"


VECTOR_3D<double> DK::GetPosition(int p0, CYCLIC_ARRAY <DK_VERTEX <double> *> &maps, int target)
{
    if( target==p0)        return VECTOR_3D<double> (0, 0, 0);
    for(int i=1; i<=maps.m; i++)
    {
        if( maps(i)->index==target)    return maps(i)->temp_P;
    }
    printf("ERROR: cannot find\n");
    return VECTOR_3D<double> (0, 0, 0);
}
double DK::GetAngle(int center, int i, int j, CYCLIC_ARRAY<DK_VERTEX <double> *> &maps, int p0)
{
    VECTOR_3D<double> e1= maps(i)->temp_P-maps(center)->temp_P;
    VECTOR_3D<double> e2= maps(j)->temp_P-maps(center)->temp_P;
    VECTOR_3D<double> e3= maps(j)->temp_P-maps(i)->temp_P;
    VECTOR_3D<double> direction=VECTOR_3D<double>::Cross_Product(e1, e2);
    if( direction.z*100 > 0 ) return -1;

    double angle1=2-VECTOR_3D<double>::Dot_Product(e1, e2);
    double angle2=2-VECTOR_3D<double>::Dot_Product(-e1, e3);
    double angle3=2-VECTOR_3D<double>::Dot_Product(-e2, -e3);
    double maxangle=angle1;
    if(maxangle<angle2)    maxangle=angle2;
    if(maxangle<angle3) maxangle=angle3;
    return maxangle;
    //    return 2-VECTOR::Dot_Product(e1, e2);
}
bool DK::Is_Covered( CYCLIC_ARRAY <DK_VERTEX <double> *> &maps, int except)
{
    int test=0;
    except=maps(except)->index;
    for(int i=1; i<=maps.m; i++)
    {
        int index=maps(i)->index;
        if(index==except)    continue;
        for(int j=1; j<=(*currentmesh->neighbor_nodes)(index).m; j++)
        {
            int neighbor=(*currentmesh->neighbor_nodes)(index)(j);
            if( neighbor==except)    continue;
            for(int k=1; k<=maps.m; k++)
                if( maps(k)->index == neighbor)
                    test++;
        }
    }
    return test==2*(2*maps.m-5);
}

//there are at most 12 neighbors. indices[0] is just the center.
void DK::Build(int p0, DK_VERTEX_LIST<float, double> &vlist)
{
    ARRAY <int>  indices=(*currentmesh->neighbor_nodes)(p0);

    double angles[13];
    double sumofangles=0;
    for(int i=1; i<=indices.m; i++)
    {
        angles[i]=sumofangles;
        VECTOR_3D<double> e1, e2; 
        e1=emb_boundary_particles.X(indices(i))-emb_boundary_particles.X(p0);
        if(i!=indices.m)
            e2=emb_boundary_particles.X(indices(i+1))-emb_boundary_particles.X(p0);
        else
            e2=emb_boundary_particles.X(indices(1))-emb_boundary_particles.X(p0);
        sumofangles+=VECTOR_3D<double>::Angle_Between(e1, e2);
    }
    //set the position.
    CYCLIC_ARRAY <DK_VERTEX <double> *> maps;
    for( i=1; i<=indices.m; i++)
    {
        double r=(emb_boundary_particles.X(indices(i))-emb_boundary_particles.X(p0)).Magnitude();
        vlist(indices(i))->temp_P=VECTOR_3D<double>(r*cos(angles[i]*2*Pi/sumofangles), 
                            r*sin(angles[i]*2*Pi/sumofangles), 0);
        maps.Append(vlist(indices(i)));
    }
    //now, let's set the angles.
    for(i=1; i<=maps.m; i++)
        maps(i)->weight=GetAngle(i, i-1, i+1, maps, p0);
    
    ARRAY <DK_VERTEX <double> *> free_vlist;
    vlist(p0)->temp_P=VECTOR_3D<double>(0, 0, 0);
    vlist(p0)->tag=true;
    free_vlist.Append( vlist(p0) );
    for( int j=1; j<= (*currentmesh->incident_triangles)(p0).m; j++)
    {
        int t=(*currentmesh->incident_triangles)(p0)(j);
        int ti, tj, tk;
        ARRAY <int> *insides=tlist->remove(t, ti, tj, tk);
        VECTOR_3D<double> p1=GetPosition(p0, maps, ti);
        VECTOR_3D<double> p2=GetPosition(p0, maps, tj);
        VECTOR_3D<double> p3=GetPosition(p0, maps, tk);
        for( int k=1; k<=insides->m; k++)
        {
            double u=vlist((*insides)(k))->u;
            double v=vlist((*insides)(k))->v;
            vlist((*insides)(k))->temp_P=p1*u+p2*v+p3*(1-u-v);
            vlist((*insides)(k))->tag=true;
            free_vlist.Append( vlist( (*insides)(k) ) );
        }
    }

    while(1)
    {
        double minangle=(double)MAX_VALUE;
        int minp;
        //the ending conditions.
        if( maps.m<3)    
        {
            printf("ERROR: Irregular triangle when removing %d, %d indices\n", p0, maps.m);
            break;
        }
        else if(maps.m==3)    
        {
            if( tlist->ExistTriangle(1, 2, 3, maps)==false)
                tlist->Add(1, 2, 3,&free_vlist, maps, p0);
            else
                printf("3, Error: failed to add the triangle %d, %d, %d when removing \
                     %d.\n",maps(1)->index, maps(2)->index, maps(3)->index, p0);
            break;
        }
        else if( maps.m==4)
        {
            // the changed implications.
            if( maps(1)->weight+maps(3)->weight < maps(2)->weight+maps(4)->weight)
                minp=1;
            else minp=2;

            if( tlist->InvalidTriangle(minp-1, minp, minp+1, maps)==false &&
                tlist->InvalidTriangle(minp+1, minp+2, minp+3, maps)==false)
            {
                tlist->Add(minp-1, minp, minp+1, &free_vlist, maps, p0);
                tlist->Add(minp+1, minp+2, minp+3, &free_vlist, maps, p0);
            }
            else if ( tlist->InvalidTriangle(minp-2, minp-1, minp, maps)==false &&
                tlist->InvalidTriangle(minp, minp+1, minp+2, maps)==false)
            {
                tlist->Add(minp-2, minp-1, minp, &free_vlist, maps, p0);
                tlist->Add(minp, minp+1, minp+2, &free_vlist, maps, p0);
            }
            else 
                printf("Error: failed to add the triangle 4.\n");            
//            if( p0==14029)
//                printf("%d, %d, %d, %d index 4\n", (*maps)(1).index, 
//                (*maps)(2).index, (*maps)(3).index, (*maps)(4).index);
            break;
        }

        //find the smallest angle.
        for(int i=1; i<=maps.m; i++)
        {
            if(minangle > maps(i)->weight && maps(i)->weight>=0)
            {
                if( tlist->InvalidTriangle(i-1, i, i+1, maps)==false)
                {
                    if( (*currentmesh->neighbor_nodes)(maps(i)->index).m>3)
                    if( Is_Covered(maps, i)==false) 
                    {
                        minangle=maps(i)->weight;
                        minp=i;
                    }    //else    printf("Throw away because of covering.\n");
                }    //else    printf("Throw away because of exisiting\n");
            }
        }
        if( minangle > Pi)    printf("ERROR: cannot find the angle -> deadlock!\n");
        //remove the vertex minp.
        if( tlist->Add(minp-1, minp, minp+1, &free_vlist, maps, p0)==false)
            printf("ERROR: failed to add the triangle.\n");
        
        maps.Remove_Index(minp);
        maps(minp-1)->weight=GetAngle(minp-1, minp-2, minp, maps, p0);
        maps(minp)->weight=GetAngle(minp, minp-1, minp+1, maps, p0);
    }

    //clearance, make sure all free vertices are assigned.
    for( i=1; i<=free_vlist.m; i++)
        if( free_vlist(i)->tag==true)
        {
            printf("ERROR: unassigned %d when removing %d\n", free_vlist(i)->index, p0);
            exit(0);
        }
}

bool DK::Output(char *filename)
{
    FILE *f=fopen(filename, "w+");
    if(!f)    return false;


    for(int i=1; i<=vlist->vertex_sequence->m; i++)
    {
        if( (*vlist->vertex_sequence)(i)->real_u<0 ||
            (*vlist->vertex_sequence)(i)->real_v<0)
        {
            printf("ERROR: vertex %d not assigned.\n", i);
            return false;
        }
    }

    for(i=1; i<=vertices_info.m; i++)
    {
        int ti=currentmesh->triangles(1, i);
        int tj=currentmesh->triangles(2, i);
        int tk=currentmesh->triangles(3, i);
        double ti_realu, ti_realv, tj_realu, tj_realv, tk_realu, tk_realv;
        ti_realu=(*vlist)(ti)->real_u;
        ti_realv=(*vlist)(ti)->real_v;
        tj_realu=(*vlist)(tj)->real_u;
        tj_realv=(*vlist)(tj)->real_v;
        tk_realu=(*vlist)(tk)->real_u;
        tk_realv=(*vlist)(tk)->real_v;

        for( int j=1; j<=vertices_info(i).m; j++)
        {
            int index=vertices_info(i)(j);
            if( (*vlist)(index)->u < 0 || (*vlist)(index)->v < 0)
            {
                printf("ERROR: not defined\n");
                return false;
            }
            (*vlist)(index)->real_u= (*vlist)(index)->u * ti_realu +
                                    (*vlist)(index)->v * tj_realu +
                            (1-(*vlist)(index)->u-(*vlist)(index)->v) * tk_realu;
            (*vlist)(index)->real_v= (*vlist)(index)->u * ti_realv +
                                    (*vlist)(index)->v * tj_realv +
                            (1-(*vlist)(index)->u-(*vlist)(index)->v) * tk_realv;
        }
    }

    for( i=1; i<=vlist->nVertices; i++)
    {
        if( (*vlist)(i)->real_u < 0 || (*vlist)(i)->real_v < 0)
        {
            printf("Error: not assigned\n");
            return false;
        }
        fprintf(f, "%f %f\n", (*vlist)(i)->real_u, (*vlist)(i)->real_v);
    }
    
    fclose(f);
    return true;
}