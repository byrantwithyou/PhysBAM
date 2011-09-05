#include "dk.h"
#include "../../Public_Library/Matrices_And_Vectors/MATRIX_MXN.h"

/*
void LSCM ( DK &dk)
{
	LIST_ARRAY<int> &triangles=dk.currentmesh->triangles;
	LIST_ARRAY <DK_VERTEX <double> *> &vertex_sequence=*dk.vlist->vertex_sequence;

	int l=0;
	for(int i=1; i<=vertex_sequence.m; i++)
		if ( vertex_sequence(i)->tag==true)
			l++;

	VECTOR_ND L(l);
	VECTOR_ND U(2*l);
	int c=1;
	for( i=1; i<=vertex_sequence.m; i++)
		if ( vertex_sequence(i)->tag==true)
		{
			L(c)=i;
			U(c)=vertex_sequence(i)->real_u;
			U(c+l)=vertex_sequence(i)->real_v;
			//CString word;
			//word.Format("L %d: %d\n", c, i);
			//AfxMessageBox(word);
			c++;
		}


	int m=triangles.m;
	int n=vertex_sequence.m;

	MATRIX_MXN A(2*m, 2*n);

	for( i=1; i<=m; i++)
	{
		VECTOR_3D<double> p1=dk.emb_boundary_particles.X( triangles(1, i));
		VECTOR_3D<double> p2=dk.emb_boundary_particles.X( triangles(2, i));
		VECTOR_3D<double> p3=dk.emb_boundary_particles.X( triangles(3, i));
		
		VECTOR_3D<double> e1=p2-p1;
		VECTOR_3D<double> e2=p3-p1;
		
		double x1=0;
		double y1=0;
		double x2=e1.Magnitude();
		double y2=0;
		double x3=VECTOR_3D<double>::Dot_Product(e2, e1)/x2;
		double y3=sqrt( e2.Magnitude_Squared()-x3*x3);
	
		for(int j=1; j<=n; j++)
		{
			int index= vertex_sequence(j)->index;
//			printf("%d, %d, %d: %d\n", triangles(1, i), triangles(2, i), triangles(3,i),
//				index);
			if( triangles(1, i)==index)
			{
				A(i, j)=x3-x2;
				A(i, j+n)=-(y3-y2);
				A(i+m, j)=y3-y2;
				A(i+m, j+n)=x3-x2;
			}
			else if( triangles(2, i)==index)
			{
				A(i, j)=x1-x3;
				A(i, j+n)=-(y1-y3);
				A(i+m, j)=y1-y3;
				A(i+m, j+n)=x1-x3;
			}
			else if( triangles(3, i)==index)
			{
				A(i, j)=x2-x1;
				A(i, j+n)=-(y2-y1);
				A(i+m, j)=y2-y1;
				A(i+m, j+n)=x2-x1;
			}
		}
	}

	MATRIX_MXN Af(2*m, 2*n-2*l);
	MATRIX_MXN Ap(2*m, 2*l);


	for(i=1; i<=2*m; i++)
	{
		int ptr_Af=1;
		int ptr_Ap=1;
		for(int j=1; j<=n; j++)
		{
			if( ptr_Ap <=L.n && j==L(ptr_Ap))
			{
				Ap(i, ptr_Ap)=A(i, j);
				Ap(i, ptr_Ap+l)=A(i, j+n);
				ptr_Ap++;          
			}
			else
            {
				Af(i, ptr_Af)=A(i, j);
				Af(i, ptr_Af+n-l)=A(i, j+n);
				ptr_Af++;
			}
		}
	}


	MATRIX_MXN ata=Af.Transpose()*Af;
	MATRIX_NXN ATA(ata.n);
	for( i=1; i<=ata.n; i++)	for(int j=1; j<=ata.n; j++)
		ATA(i, j)=ata(i, j);
	ATA.LU_Inverse();
	
	for( i=1; i<=ata.n; i++)	for(int j=1; j<=ata.n; j++)
		ata(i, j)=(*ATA.inverse)(i, j);
	MATRIX_MXN B=Af.Transpose()*Ap;
	MATRIX_MXN X=ata*B;

/*	for(i=1; i<=X.m; i++)
	{
		float sum=0;
		for(int j=1; j<=X.n; j++)
			sum+=X(i, j);
	}
*/
/*	c=1;
	VECTOR_ND x=X*U;
//	CString word;
//	word.Format("%d, %d, %d\n", m, n, x.n);
//	AfxMessageBox(word);
	for( i=1; i<=vertex_sequence.m; i++)
		if ( vertex_sequence(i)->tag==false)
		{
			vertex_sequence(i)->real_u=-x(c);
			vertex_sequence(i)->real_v=-x(c+n-l);
	//		CString word;
	//		word.Format("L %f: %f\n", -x(c), -x(c+n-l));
	//		AfxMessageBox(word);
			c++;
		}
}
*/