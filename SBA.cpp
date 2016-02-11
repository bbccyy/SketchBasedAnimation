#include "stdafx.h"
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include<fstream>
#include<string>
#include<stdio.h>
#include <GL/glut.h>

using namespace std;


 struct mypoint{
	float x;
	float y;
	int x2;
	int y2;
	float xf;
	float yf;
	int num;
	int index;	
	struct mypoint* nextpoint;
	struct mypoint* frontpoint;
	struct mypoint* store;
    };


 struct trimesh{                      // store output tri_mesh which could be used to show pic in 3D.
	 float a1;
	 float a2;
	 float a3;
	 float b1;
	 float b2;
	 float b3;
	 float c1;
	 float c2;
	 float c3;
	 float d1;
	 float d2;
	 float d3;
	 int X;
	 struct trimesh* nextTri;
	 struct trimesh* frontTri;
	 struct point3D* link3D;
 };


 class scrPt
{
public:
	GLint x,y;
};

 struct point3D{
	 float a1;
	 float a2;
	 float a3;
	 int X;
	 int CheckifonlyOnestroke;
	 struct point3D* next3D;
	 struct point3D* front3D;
 };
 
static struct InRegion{   //struct that stores possible vertex for the given line segment.
	int x;
	int y;
	struct InRegion * nextIn;
	struct InRegion * frontIn;
    } *headerI ;


struct output {      // struct for triangle output.
	int v1[2];   // v1 store center p.
	float v2[2];
	int v3[2];     // v2 and v3 store line.
	float v4[2];    // float v2 and v4 in spineOut are store line segment.
	bool TwoInnerline;
	int X;
	struct output * nextOut;
	struct output * front;
	struct output * store;
	struct mypoint * storemy ;
    };



struct Center {          // store the center of circumcircle and its radius.
	double center[2];
	double distance;
};

GLfloat xRot = 0.0f;
GLfloat yRot = 0.0f;
GLfloat xTrans = 0;
GLfloat yTrans = 0;
GLfloat zTrans = -32.0;

int ox;
int oy;
int buttonState;
GLfloat xRotLength = 0.0f;
GLfloat yRotLength = 0.0f;

static scrPt endPt1, endPt2;
	   int index_all=1;
       int num_all=1;
	   GLsizei winWidth=600, winHeight =600; //Initial display window size.
       mypoint* header = NULL ;
	   mypoint* current= header;
	   static  mypoint * end_mypoint ;
	   bool seeround = true;
       mypoint * tempstore;
	   mypoint * baseringP = NULL;
static trimesh * triEx = NULL;
static output * spineOut = NULL;
static output * spineTri = NULL;  // used to store final triangle.
static mypoint * spine2p = NULL;  // used to store lines connect spine point and point on polygon.
static trimesh * tri3D = NULL;    // used to store tri mesh in 3D space and show up.
static int indexA = 1;   // for count cutting.
       int indexI = 2713;   // for painting the base ring.
       point3D * pointOut = NULL;     //  store the Painting strokes.
	   point3D * tempstrokes = NULL;
	   int strokecount = 1;   // count painting strokes.  dont delete.
	   int checknum=1;
	   float ur[7][2], dr[5][2];


	   GLfloat LightAmbient[]= { 0.5f, 0.5f, 0.5f, 1.0f }; 				// enviroment parameter.
       GLfloat LightDiffuse[]= { 1.0f, 1.0f, 1.0f, 1.0f };				 // diffuse light parameter.
       GLfloat LightPosition1[]= {0.0f, 100.0f, 950.0f, 1.0f };	
	   GLfloat LightPosition2[]= { 0.0f, -100.0f, -950.0f, 1.0f };
	   GLfloat LightPosition3[]= { 0.0f, 950.0f, 150.0f, 1.0f };	
	   GLfloat LightPosition4[]= { 0.0f, -950.0f, -150.0f, 1.0f };	
	   GLfloat LightPosition5[]= { 950.0f, 100.0f, 0.0f, 1.0f };	
	   GLfloat LightPosition6[]= { -950.0f, -100.0f, 0.0f, 1.0f };	

	   bool switchmesh = true;


struct mypoint *  readdata ()     //read points data in file.
{
 mypoint * P = new mypoint ,* Q = P ; 
 P->nextpoint =NULL;
 
 struct mypoint news1;
 fstream file("pointdata.txt",ios::in);
 while(1)
 {
  file.read((char *)& news1,sizeof(news1));  
  if(!file)break;
  else
  { 
   struct mypoint * V = new mypoint;

   V->num = news1.num;
   V->nextpoint=news1.nextpoint;
   V->x = news1.x;
   V->y = news1.y;   

   Q->nextpoint = V;   
   Q->nextpoint ->nextpoint = NULL;
   Q = Q->nextpoint;    
  }
 }
 file.close();
 return (P->nextpoint);   
}

void savedata(struct mypoint* head)    // save data.
{
	 mypoint head1;
	 fstream file("pointdata.txt",ios::out);

	 if(head != NULL)
	 {
		 while(head)
		 {
			 file.write((char*)head, sizeof(head1));

			 head = head->nextpoint;
		 }
	 }
	 file.close();
	 cout<< " SAVE COMPLETE!";
}

struct mypoint * changelength ( mypoint * header )     // divide chain header to make line segments.
{
	mypoint * current = header ;
	mypoint * header2 = new mypoint ;
	mypoint * currentNew = header2 ;
	mypoint * trans ;
	int numNew = 1;
	int innerNum =0 ;
	while ( current->nextpoint != NULL /*&& current->index == 1 */)  //in this ver. only the first stroke could be calculated.
	{
		if ( current->num%20==0 || current->num ==1 )
		{
			innerNum = 0;
			currentNew->x = current->x ;
			currentNew->y = current->y ;
			currentNew->num = numNew ;
			cout<<currentNew->num<<endl;
			numNew++;
			trans = currentNew ;
			currentNew->nextpoint = new mypoint ;
			currentNew = currentNew->nextpoint ;
		}
		current = current->nextpoint ;
		innerNum = innerNum + 1 ;
	}
	if (innerNum>=5)
	{
		currentNew->x = current->x ;
	    currentNew->y = current->y ;
	    currentNew->num = numNew ;
	    cout<<currentNew->num<<endl;
	    currentNew->nextpoint = NULL ;
	
	}
	else
	{
		trans->nextpoint = NULL ;
		delete currentNew ;
	}
	return header2;
}


//=========================S T A R T <C D T>=========================


bool Equal(float f1, float f2) {
    return (abs(f1 - f2) < 1e-4f);
}
//judge if 2 points are the same.
bool operator==(mypoint &p1, mypoint &p2) {
    return (Equal(p1.x, p2.x) && Equal(p1.y, p2.y));
}
// compare those 2 points, frist by x, if equal than compare y.
bool operator>(mypoint &p1, mypoint &p2) {
    return (p1.x > p2.x || (Equal(p1.x, p2.x) && p1.y > p2.y));
}
//out product of 2 vectors
float operator^(mypoint &p1, mypoint &p2) {
    return (p1.x * p2.y - p1.y * p2.x);
}
//determine the relationship between two line segments, than return the point of intersection:
//[do overlap] totally overlap[6]  both in one line and share one end point[5]  partially overlap[4]
//[no overlap] intersected at the end point[3]   intersect on line[2]  orthogonality[1]  no intersect[0]   wrong[-1]
int Intersection(mypoint p1, mypoint p2, mypoint p3, mypoint p4, mypoint &Int) {
    //make sure that p1!=p2，p3!=p4
    if (p1 == p2 || p3 == p4) {
        return -1; // -1 represents that at least one line segment overlapped from head to tail, and cannot construct a line segment.
    }
    //make sure start point alway at first.
    if (p1 > p2) {
        swap(p1, p2);
    }
    if (p3 > p4) {
        swap(p3, p4);
    }
    //judge if the 2 line segments are totally overlapped.
    if (p1 == p3 && p2 == p4) {
        return 6;
    }
    //calculate the line vectors.
    mypoint v1 = {p2.x - p1.x, p2.y - p1.y}, v2 = {p4.x - p3.x, p4.y - p3.y};
    //out product, if parallel the result is 0
    float Corss = v1 ^ v2;
    //if the start points are overlapped.
    if (p1 == p3) {
        Int = p1;
        //if the starting points are overlapped (as well parallel) return [5], if not parallel, but certainlly intersects at end points, return [3].
        return (Equal(Corss, 0) ? 5 : 3);
    }
    //if tail points are overlapped.
    if (p2 == p4) {
        Int = p2;
        //if the tail points are overlapped (as well parallel) return [5], if not parallel, but certainlly intersects at end points, return [3].
        return (Equal(Corss, 0) ? 5 : 3);
    }
    //if the 2 line segments connect with end points.
    if (p1 == p4) {
        Int = p1;
        return 3;
    }
    if (p2 == p3) {
        Int = p2;
        return 3;
    }
    //rearrage points from the starting point. if line 1 bigger, swap.
    if (p1 > p3) {
        swap(p1, p3);
        swap(p2, p4);
        //refresh out product.
        swap(v1, v2);
        Corss = v1 ^ v2;
    }
    //if parallel.
    if (Equal(Corss, 0)) {
        //out product of vector v1(p1,p2)and v2(p3,p4),judge if on the same line.
        mypoint vs = {p3.x - p1.x, p3.y - p1.y};
        //out product=0 than on same line,then determine if overlap.
        if (Equal(v1 ^ vs, 0)) {
            //the end point of form line>the start point of current line, then they overlap.
            if (p2 > p3) {
                Int = p3;
                return 4; 
            }
        }//if 3 points not in same line, the 2 parallel line will not on the same line.
        //no intersection.
        return 0;
    } //the rest are unparallel case.

    float ymax1 = p1.y, ymin1 = p2.y, ymax2 = p3.y, ymin2 = p4.y;
    if (ymax1 < ymin1) {
        swap(ymax1, ymin1);
    }
    if (ymax2 < ymin2) {
        swap(ymax2, ymin2);
    }

    if (p1.x > p4.x || p2.x < p3.x || ymax1 < ymin2 || ymin1 > ymax2) {
        return 0;
    }
    mypoint vs1 = {p1.x - p3.x, p1.y - p3.y}, vs2 = {p2.x - p3.x, p2.y - p3.y};
    mypoint vt1 = {p3.x - p1.x, p3.y - p1.y}, vt2 = {p4.x - p1.x, p4.y - p1.y};
    float s1v2, s2v2, t1v1, t2v1;

    if (Equal(s1v2 = vs1 ^ v2, 0) && p4 > p1 && p1 > p3) {
        Int = p1;
        return 2;
    }
    if (Equal(s2v2 = vs2 ^ v2, 0) && p4 > p2 && p2 > p3) {
        Int = p2;
        return 2;
    }
    if (Equal(t1v1 = vt1 ^ v1, 0) && p2 > p3 && p3 > p1) {
        Int = p3;
        return 2;
    }
    if (Equal(t2v1 = vt2 ^ v1, 0) && p2 > p4 && p4 > p1) {
        Int = p4;
        return 2;
    }
    if(s1v2 * s2v2 > 0 || t1v1 * t2v1 > 0) {
        return 0;
    } 

    float ConA = p1.x * v1.y - p1.y * v1.x;
    float ConB = p3.x * v2.y - p3.y * v2.x;

    Int.x = (ConB * v1.x - ConA * v2.x) / Corss;
    Int.y = (ConB * v1.y - ConA * v2.y) / Corss;

    return 1;
}

// find the end point of struct mypoint.
struct mypoint * Fin_end (mypoint *header, mypoint * end_mypoint)    // Find out the end of the chain.
{
	mypoint * current=header;
	while(current->nextpoint != NULL)
	{
	    current = current->nextpoint;
		cout<< current->num<<endl;
	}
	return end_mypoint=current;
}

// get the Angle from positive x-axis to a vector (x,y) with angle between 0 to 360°.
// used in function: Generate_Inregion for determing wether a point is inside the polygon.
double getRotateAngle(double x1, double y1, double x2, double y2)  
{ 
 const double epsilon = 1.0e-6;
 const double nyPI = acos(-1.0);
 double dist, dot, degree, angle;
 
 // normalize
 dist = sqrt( x1 * x1 + y1 * y1 );
 x1 /= dist;
 y1 /= dist;
 dist = sqrt( x2 * x2 + y2 * y2 );
 x2 /= dist;
 y2 /= dist;
 // dot product
 dot = x1 * x2 + y1 * y2;
 if ( fabs(dot-1.0) <= epsilon ) 
  angle = 0.0;
 else if ( fabs(dot+1.0) <= epsilon ) 
  angle = nyPI;
 else {
  double cross;
  
  angle = acos(dot);
  //cross product
  cross = x1 * y2 - x2 * y1;
  // vector p2 is clockwise from vector p1 
  // with respect to the origin (0.0)
  if (cross < 0 ) { 
   angle = 2 * nyPI - angle;
  }    
 }
 degree = angle *  180.0 / nyPI;
 return degree;
}

// delete struct InRegion as whole.
struct InRegion * deleteInre_whole ( InRegion * headerI )
{
	InRegion * currentIn = headerI ;
	while ( (currentIn=headerI)!= NULL )
	{
		headerI = currentIn->nextIn ;
		delete currentIn ;
	}
	headerI = new InRegion ;
	return headerI;
}

// only delete the given struct InRegion.
struct InRegion * Delete_Inregion ( InRegion * deleteable, InRegion * headerI )
{
	if ( deleteable==headerI )
	{
		headerI = headerI->nextIn;
		delete deleteable;
		headerI->frontIn=NULL;
		return headerI;
	}
	else if (deleteable->nextIn==NULL)
	{
		InRegion * currentIn = deleteable->frontIn;
		currentIn->nextIn = NULL;
		delete deleteable;
		return headerI;
	}
	else
	{
		InRegion * currentIn = deleteable->frontIn;
		currentIn->nextIn = deleteable->nextIn;
		currentIn = deleteable->nextIn;
		currentIn->frontIn = deleteable->frontIn;
		delete deleteable;
		return headerI;
	}
}  

// for every line segments in the polygon, find its corresponding line: Ax+By+C=0,
// then find the angle between [1,0] to that line(with anticlock direction).
// according to the angle and line function, determine points that locate on the right side of the line.
// for the given line segment, calculate the possible points and collect then in InRegion struct.
struct InRegion * Generate_Inregion ( mypoint * step1_current1, mypoint * step1_current2, mypoint * header, InRegion * headerI)
{
	//headerI = new InRegion;
	int x1 = step1_current1->x ;
	int y1 = step1_current1->y ;
	int x2 = step1_current2->x ;
	int y2 = step1_current2->y ;
	int X = x2-x1 ;
	int Y = y2-y1 ;

	double Angle = getRotateAngle(1,0,X,Y);

	//headerI->frontIn=NULL;
	mypoint * current = header;
	InRegion * currentIn = headerI ;
	InRegion * trans;

	if ( (Angle>0&&Angle<90) || (Angle>90&&Angle<180) )
	{
		while (current!=NULL)
		{
			trans = currentIn ;
			int yi = current->y ;
			int xi = current->x ;
			double xii = X*(yi-y1)/Y+x1 ;
			if ( xi < xii )
			{
				currentIn->x = xi;
				currentIn->y = yi;
				currentIn->nextIn = new InRegion ;
				currentIn = currentIn->nextIn; 
				currentIn->frontIn = trans;
			}
			current = current->nextpoint; 
		}
		trans=currentIn->frontIn;
		trans->nextIn = NULL;
		delete currentIn ;
	}

	else if ( Angle == 90 )
	{
		while (current!=NULL)
		{
			trans = currentIn ;
			int yi = current->y ;
			int xi = current->x ;  
			if ( xi < x1 )
			{
				currentIn->x = xi;
				currentIn->y = yi;
				currentIn->nextIn = new InRegion ;
				currentIn = currentIn->nextIn; 
				currentIn->frontIn = trans;
			}
			current = current->nextpoint; 
		}
		trans=currentIn->frontIn;
		trans->nextIn = NULL;
		delete currentIn ;
	}

	else if ( Angle == 180 )
	{
		while (current!=NULL)
		{
			trans = currentIn ;
			int yi = current->y ;
			int xi = current->x ;  
			if ( yi < y1 )
			{
				currentIn->x = xi;
				currentIn->y = yi;
				currentIn->nextIn = new InRegion ;
				currentIn = currentIn->nextIn; 
				currentIn->frontIn = trans;
			}
			current = current->nextpoint; 
		}
		trans=currentIn->frontIn;
		trans->nextIn = NULL;
		delete currentIn ;
	}

		else if ( Angle == 270 )
	{
		while (current!=NULL)
		{
			trans = currentIn ;
			int yi = current->y ;
			int xi = current->x ;  
			if ( xi > x1 )
			{
				currentIn->x = xi;
				currentIn->y = yi;
				currentIn->nextIn = new InRegion ;
				currentIn = currentIn->nextIn; 
				currentIn->frontIn = trans;
			}
			current = current->nextpoint; 
		}
		trans=currentIn->frontIn;
		trans->nextIn = NULL;
		delete currentIn ;
	}

		else if ( Angle == 0 || Angle == 360 )
	{
		while (current!=NULL)
		{
			trans = currentIn ;
			int yi = current->y ;
			int xi = current->x ;  
			if ( yi > y1 )
			{
				currentIn->x = xi;
				currentIn->y = yi;
				currentIn->nextIn = new InRegion ;
				currentIn = currentIn->nextIn; 
				currentIn->frontIn = trans;
			}
			current = current->nextpoint; 
		}
		trans=currentIn->frontIn;
		trans->nextIn = NULL;
		delete currentIn ;
	}

	else if ( (Angle>180&&Angle<270) || (Angle>270&&Angle<360) )
	{
		while (current!=NULL)
		{
			trans = currentIn ;
			int yi = current->y ;
			int xi = current->x ;
			double xii = X*(yi-y1)/Y+x1 ;
			if ( xi > xii )
			{
				currentIn->x = xi;
				currentIn->y = yi;
				currentIn->nextIn = new InRegion ;
				currentIn = currentIn->nextIn; 
				currentIn->frontIn = trans;
			}
			current = current->nextpoint; 
		}
		trans=currentIn->frontIn;
		trans->nextIn = NULL;
		delete currentIn ;
	}
	return headerI;

}

// use the output from functiong: Generate_Inregion, then refine the them in such way:
// first find a point in InRegion;
// secondly, according to this point, go through every edges on the polygon to see as if any intersection occurs.
// remove the point if it intersect with any polylines.
struct InRegion * Refine_Inregion (mypoint * Onepoint, InRegion * headerI, mypoint * header)  // Onepoint could be step1_current1 or step1_current2
{
	InRegion * currentIn = headerI;
	//int x1= Onepoint->x ;
	//int y1= Onepoint->y ;
	while ( currentIn != NULL )   // first find a point in InRegion
	{
		mypoint B;
		mypoint Inter;
		B.x = currentIn->x ;
		B.y = currentIn->y ;
		//int x2=currentIn->x ;
		//int y2=currentIn->y ;
		InRegion * deleteable = currentIn;
		currentIn=currentIn->nextIn ;
		mypoint * current = header;
		mypoint * trans = current;
		current = current->nextpoint ;
		do                        // secondly, according to this point, go through every edges on the polygon to see as if any intersection occurs.
		{
			if (current==NULL)
			{
				//int x4 = header->x ;
				//int y4 = header->y ;
				//int x3 = trans->x ;
				//int y3 = trans->y ;
				int Judge = Intersection (*Onepoint,B,*trans,*header,Inter);
				if (Judge>0&&Judge!=3&&Judge!=6)                              //Do intersect!(excluding totally overlap and connect with end both points)
				{
					headerI = Delete_Inregion(deleteable,headerI);
					break;
				}
				break;
			}
			else
			{
				//int x4 = current->x ;
				//int y4 = current->y ;
				//int x3 = trans->x ;
				//int y3 = trans->y ; 
				int Judge = Intersection (*Onepoint,B,*trans,*current,Inter);
				if (Judge>0&&Judge!=3&&Judge!=6)
				{
					headerI = Delete_Inregion(deleteable,headerI);
					break;
				}
				trans = current;
				current = current->nextpoint ;
			}
		}while(1);    //end of inner loop.
	}   // end of out loop.
	return headerI ;
}

// calculate the center and radius of a given triangle's circumcircle.
void calTriangleCircum(mypoint * pt1, mypoint * pt2, InRegion * pt3, Center * Center_R)
{
	double xa = pt1->x;
 	double ya = pt1->y;
 	double xb = pt2->x;
 	double yb = pt2->y;
 	double xc = pt3->x;
 	double yc = pt3->y;
 	//double x = ((y2-y1)*(y3*y3-y1*y1+x3*x3-x1*x1)-(y3-y1)*(y2*y2-y1*y1+x2*x2-x1*x1))/(2*(x3-x1)*(y2-y1)-2*((x2-x1)*(y3-y1))); 
	//double y = ((x2-x1)*(x3*x3-x1*x1+y3*y3-y1*y1)-(x3-x1)*(x2*x2-x1*x1+y2*y2-y1*y1))/(2*(y3-y1)*(x2-x1)-2*((y2-y1)*(x3-x1))); 
	double x = ((xb*xb+yb*yb-xa*xa-ya*ya)/2*(yc-ya)-(xc*xc+yc*yc-xa*xa-ya*ya)/2*(yb-ya))/((xb-xa)*(yc-ya)-(yb-ya)*(xc-xa));
	double y = ((xb*xb+yb*yb-xa*xa-ya*ya)/2*(xc-xa)-(xc*xc+yc*yc-xa*xa-ya*ya)/2*(xb-xa))/((yb-ya)*(xc-xa)-(xb-xa)*(yc-ya));

 	// center of circumcircle
 	Center_R->center[0]=x;
 	Center_R->center[1]=y;
 	// radius
 	double R = sqrt(pow(x-xa,2)+pow(y-ya,2));
 	Center_R->distance = R;
}



// further refine the data of struct Inregion.
// given a chain of Inregion, find a point, with two end points on the corresponding line segments,to form a triangle.
// computer the triangle's circumcircle, then see if any other points in the Inregion inside the region of circle.
// after several loop, the system will find a triangle with no other points inside its circumcircle.
// then save it in the struct output currentOut.
// additional information should add in currentOut as well, such as: a terminal triangle(with two edge lines), or a sleeve triangle(with only one edge line).
void Circle_Inregion ( mypoint * step1_current1, mypoint * step1_current2, InRegion * headerI, mypoint * header, output * currentOut,int* X, int* Y)
{
	InRegion * currentIn = headerI ;
	while ( currentIn != NULL )
	{
		InRegion * current_check = headerI ;
		//mypoint * current = header ;
		Center D;
		bool pointOK = true ;
		calTriangleCircum( step1_current1, step1_current2,currentIn, &D);

		while ( current_check != NULL )
		{
			int x1 = current_check->x ;
			int y1 = current_check->y ;
			double Radius = sqrt(pow(x1-D.center[0],2)+pow(y1-D.center[1],2));

			if ( (Radius+0.0001)<D.distance)  //compare the distances, and make sure that same point will not interrupt the process.
			{
				pointOK = false;
				break;
			}
			current_check = current_check->nextIn ;
		}  // end inner loop.

		mypoint * next2 = new mypoint;
		next2->x = -1;
		next2->y = -1;
		mypoint * trans_next2 = next2 ;

		if ( pointOK == true )
		{
			if (step1_current2->nextpoint!=NULL)
            next2 = step1_current2->nextpoint ;

			if ( step1_current1->num==1 && ( end_mypoint->x == currentIn->x && end_mypoint->y == currentIn->y))
			{
				currentOut->v1[0] = step1_current1->x ;
				currentOut->v1[1] = step1_current1->y ;
				currentOut->v2[0] = step1_current2->x ;
				currentOut->v2[1] = step1_current2->y ;
				currentOut->v3[0] = currentIn->x ;
				currentOut->v3[1] = currentIn->y ;
				*X=currentIn->x ;
				*Y=currentIn->y ;
				currentOut->TwoInnerline = false ;
				currentOut->nextOut = new output ;
			//	currentOut = currentOut->nextOut ;
			}
			else if ((next2->x == currentIn->x && next2->y == currentIn->y) ||
				( step1_current2->num == end_mypoint->num && (currentIn->x == header->x && currentIn->y == header->y)))
			{
				currentOut->v1[0] = step1_current2->x ;
				currentOut->v1[1] = step1_current2->y ;
				currentOut->v2[0] = currentIn->x ;
				currentOut->v2[1] = currentIn->y ;
				currentOut->v3[0] = step1_current1->x ;
				currentOut->v3[1] = step1_current1->y ;
				*X=currentIn->x;
				*Y=currentIn->y;
				currentOut->TwoInnerline = false ;
				currentOut->nextOut = new output ;
			//	currentOut = currentOut->nextOut ;
			}
			else
			{
				currentOut->v1[0] = step1_current1->x ;
				currentOut->v1[1] = step1_current1->y ;
				currentOut->v2[0] = currentIn->x ;
				currentOut->v2[1] = currentIn->y ;
				currentOut->v3[0] = step1_current2->x ;
				currentOut->v3[1] = step1_current2->y ;
				*X=currentIn->x;
				*Y=currentIn->y;
				currentOut->TwoInnerline = true ;
				currentOut->nextOut = new output ;
				//currentOut = currentOut->nextOut ;
			}  //end inner if.
			delete trans_next2 ;
			return ;
		}  //end out if.
		delete trans_next2 ;
		currentIn = currentIn->nextIn ;
	} // end out loop.
	cout<< "error: NO point available to form CDT."<<endl ;
}

// using all the above function to complete the CDT process.
// make sure also link the start point and the end point.
struct output * CDT( mypoint * header )
{
	//Fin_end( header, end_mypoint );  headerOut
	output * headerOut=new output, *currentOut;

	currentOut = headerOut ;
	mypoint * step1_current1;
	mypoint * step1_current2;
	step1_current1 = header ;
	step1_current2 = step1_current1->nextpoint ;
	bool type = false;
	headerI = new InRegion ;

	do
	{
		mypoint * trans;
		int X , Y ;

		if ( step1_current1->num == 1 )
		{
			headerI = Generate_Inregion (step1_current1,step1_current2,header,headerI);
			headerI = Refine_Inregion (step1_current1,headerI,header);
			headerI = Refine_Inregion (step1_current2,headerI,header);
			Circle_Inregion (step1_current1,step1_current2,headerI,header,currentOut,&X,&Y);
			currentOut=currentOut->nextOut;
			headerI = deleteInre_whole(headerI);
			trans = step1_current2->nextpoint ;
			if ( trans->x==X && trans->y==Y )     // if the finded point = the point next to step1_current2.
			{
				step1_current1=trans ;
				step1_current2=trans->nextpoint ;    // move one more step.
				continue ;
			}
			else if ( end_mypoint->x==X && end_mypoint->y==Y )   // if for the first 2 step(#1 & #2)  the finded point is the end point.
			{
				type = true;
			}
			step1_current1 = step1_current2 ;
			step1_current2 = step1_current2->nextpoint ;
			continue ;
		}

		if ( step1_current1==end_mypoint && type == true )
		{
			currentOut->nextOut  = NULL ;
			break;
		}

		headerI = Generate_Inregion (step1_current1,step1_current2,header,headerI);
		headerI = Refine_Inregion (step1_current1,headerI,header);
		headerI = Refine_Inregion (step1_current2,headerI,header);
		Circle_Inregion (step1_current1,step1_current2,headerI,header,currentOut,&X,&Y);
		currentOut=currentOut->nextOut;
		headerI = deleteInre_whole(headerI);
		trans = step1_current2->nextpoint ;
		if ( trans == NULL )                     // if step1_current2 travels to the end.
		{
			if ( header->x==X && header->y==Y )    // check if the finded point is the frist point.
			{
				currentOut->nextOut  = NULL ;
				break;
			}
            step1_current1 = step1_current2 ;      //step1_current2 will move to the head, step1_current1 in this case is the end.
			step1_current2 = header ;
			continue ;
		}
		if ( step1_current2==header)              // if the step1_c2 already at header, close the chain of currentOut, break!
		{
			currentOut->nextOut  = NULL ;
			break ;
		}
		if ( trans->x==X && trans->y==Y )         // if the finded point is the next point of step1_c2, move them one more step.
		{
			if (trans->nextpoint == NULL){        // this is a very rare problem, must fix it to avoid bug. (how about step1_c1 and s1_c2 are #end-2 and #end-1??, 
				step1_current1=trans ;            // if both move one more step at this situation, the s1_c2 will meets nothing!)
				step1_current2=header ;}
			else
			{
				step1_current1=trans ;
				step1_current2=trans->nextpoint ;
			}
		}
		else                                     // normal case, move one step ahead.
		{
			step1_current1 = step1_current2 ;
			step1_current2 = step1_current2->nextpoint ;
		}
	} while ( step1_current2 != header->nextpoint );

    currentOut = headerOut ;
	output * transO ;
	while (1)
	{
		transO = currentOut->nextOut  ;
		if ( transO->nextOut == NULL )
		{
			currentOut->nextOut  = NULL ;
			delete transO ;
			break ;
		}
		currentOut = currentOut->nextOut ;
	}
	return headerOut;
} 


//============ Find out Junctriangle=============



struct link{
	int va[2];
	int vb[2];
	struct link* nextlink;
	struct link* frontlink;
};


struct link * deleteLink_whole ( link * headerLink )
{
	if (headerLink->nextlink!=NULL){
	link * currentLink = headerLink ;
	while ( (currentLink=headerLink)!= NULL )
	{
		headerLink = currentLink->nextlink ;
		delete currentLink ;
	}
	}
	return NULL;
}


/* if possible, this function will return a pointer to an output struct which contains the three vetexs of junction triangle
   but, if there are no junction triangle to the given line segment, function will return NULL.*/
struct output * Junctri ( link * headerlink, link * currentlink ,output * junctriangle)
{
	link * va_sub = new link ,* currentva=va_sub;
	link * vb_sub = new link ,* currentvb=vb_sub;
	va_sub->va[1]=NULL;        //initialize sub chain elements with NULL.
	vb_sub->va[1]=NULL;
	link * current = headerlink ;
	int va1 = currentlink->va[0] ;
	int va2 = currentlink->va[1] ;  // line segment va---vb   point va
	int vb1 = currentlink->vb[0] ; 
	int vb2 = currentlink->vb[1] ;  // point vb
	link * trans;

	while ( current!=NULL )
	{
		if ((current->va[0]==va1&&current->va[1]==va2) && (current->vb[0]!=vb1 || current->vb[1]!=vb2))    // for given endpoint (va1 va2).
		{
			currentva->va[0]=current->vb[0];
			currentva->va[1]=current->vb[1];
			trans = currentva;
			currentva->nextlink = new link ;
			currentva=currentva->nextlink ;
			currentva->frontlink = trans ;
		}
		else if ((current->vb[0]==va1&&current->vb[1]==va2) && (current->va[0]!=vb1 || current->va[1]!=vb2) ) 
		{
			currentva->va[0]=current->va[0];
			currentva->va[1]=current->va[1];
			trans = currentva;
			currentva->nextlink = new link ;
			currentva=currentva->nextlink ;
			currentva->frontlink = trans ;
		}
		current = current->nextlink ;
	}
	if ( va_sub->va[1]!= NULL )
	{
		trans = currentva->frontlink ;
		trans->nextlink = NULL;
		delete currentva ;
	}

	current = headerlink ;

	while ( current!=NULL )
	{
		if ((current->va[0]==vb1&&current->va[1]==vb2) && (current->vb[0]!=va1 || current->vb[1]!=va2))    // for given endpoint (vb1 vb2).
		{
			currentvb->va[0]=current->vb[0];
			currentvb->va[1]=current->vb[1];
			trans = currentvb;
			currentvb->nextlink = new link ;
			currentvb=currentvb->nextlink ;
			currentvb->frontlink = trans ;
		}
		else if ((current->vb[0]==vb1&&current->vb[1]==vb2) && (current->va[0]!=va1 || current->va[1]!=va2) ) 
		{
			currentvb->va[0]=current->va[0];
			currentvb->va[1]=current->va[1];
			trans = currentvb;
			currentvb->nextlink = new link ;
			currentvb=currentvb->nextlink ;
			currentvb->frontlink = trans ;
		}
		current = current->nextlink ;
	}
	if ( vb_sub->va[1]!= NULL )
	{
		trans = currentvb->frontlink ;
		trans->nextlink = NULL;
		delete currentvb ;
	}

    // now we have va_sub and vb_sub, time to find there junction!

	if ( va_sub->va[1]==NULL || vb_sub->va[1]==NULL ) return NULL ;

	currentva = va_sub;
	while ( currentva!=NULL )    // for every element in va_sub.
	{
		currentvb = vb_sub ;
		while ( currentvb != NULL )   // search every points in vb_sub.
		{
			if ( currentva->va[0]==currentvb->va[0] && currentva->va[1]==currentvb->va[1] )  //if find junction.
			{
				output * currentjunc = junctriangle ;
				while ( currentjunc != NULL )           // search to make sure that no redundent values in chain junctriangle (already have the same values).
				{
					if ( (currentjunc->v1[0]==va1 && currentjunc->v1[1]==va2 && currentjunc->v2[0]==vb1 && currentjunc->v2[1]==vb2
						&& currentjunc->v3[0]==currentva->va[0] && currentjunc->v3[1]==currentva->va[1] )     
						||  (currentjunc->v1[0]==va1 && currentjunc->v1[1]==va2 && currentjunc->v2[0]==currentva->va[0] && currentjunc->v2[1]==currentva->va[1]
						&& currentjunc->v3[0]==vb1 && currentjunc->v3[1]==vb2 )   
					    || (currentjunc->v1[0]==vb1 && currentjunc->v1[1]==vb2 && currentjunc->v2[0]==va1 && currentjunc->v2[1]==va2
						&& currentjunc->v3[0]==currentva->va[0] && currentjunc->v3[1]==currentva->va[1] ) 
						|| (currentjunc->v1[0]==vb1 && currentjunc->v1[1]==vb2 && currentjunc->v2[0]==currentva->va[0]&& currentjunc->v2[1]==currentva->va[1]
						&& currentjunc->v3[0]==va1 && currentjunc->v3[1]==va2 ) 
						|| (currentjunc->v1[0]==currentva->va[0] && currentjunc->v1[1]==currentva->va[1] && currentjunc->v2[0]==va1 && currentjunc->v2[1]==va2
						&& currentjunc->v3[0]==vb1 && currentjunc->v3[1]==vb2 ) 
						|| (currentjunc->v1[0]==currentva->va[0] && currentjunc->v1[1]==currentva->va[1] && currentjunc->v2[0]==vb1 && currentjunc->v2[1]==vb2
						&& currentjunc->v3[0]==va1 && currentjunc->v3[1]==va2 ) )
						{
							va_sub = deleteLink_whole(va_sub);
							vb_sub = deleteLink_whole(vb_sub);
							return NULL ;
						}

					currentjunc = currentjunc->nextOut ;
				}
				output * junctemp = new output ;
                junctemp->v1[0]=va1;
			    junctemp->v1[1]=va2;
				junctemp->v2[0]=vb1;
				junctemp->v2[1]=vb2;
				junctemp->v3[0]=currentva->va[0];
				junctemp->v3[1]=currentva->va[1];
				va_sub = deleteLink_whole(va_sub);
				vb_sub = deleteLink_whole(vb_sub);
				return junctemp ;
			} // end if find junction

			// if no junctriangle here, try to find junc quad.
			int qa0 = currentva->va[0];
			int qa1 = currentva->va[1];
			int qb0 = currentvb->va[0];
			int qb1 = currentvb->va[1];
			current = headerlink ;

			while (current!=NULL)
			{
				if ( (current->va[0]==qa0 && current->va[1]==qa1 && current->vb[0]==qb0 && current->vb[1]==qb1)
					||(current->va[0]==qb0 && current->va[1]==qb1 && current->vb[0]==qa0 && current->vb[1]==qa1))  // if find inner quad.
				{
					trans = headerlink->nextlink ;
					link * insert = new link;
					insert->va[0] = qa0 ;
					insert->va[1] = qa1 ;
					insert->vb[0] = vb1 ;
					insert->vb[1] = vb2 ;
					insert->frontlink = headerlink ;
					insert->nextlink = trans ;
					trans->frontlink = insert ;
					break ;
				}
				current = current->nextlink ;
			}                                       // just find out the quad and divide it into 2 pieces, store and refresh the new inner edge into 'headerlink'

			currentvb = currentvb->nextlink ;
		} // end inner while about vb_sub.
		
		currentva = currentva->nextlink ;
	} // end out while for va_sub.

	va_sub = deleteLink_whole(va_sub);
	vb_sub = deleteLink_whole(vb_sub);
	return NULL;
}

// return junctriangle in struct output.
// function will find out all triangle with three inner lines(junctriangle).
struct output * adjust_chain( output * headerOut )
{
	output * junctriangle = NULL ;  //used to store inner triangle without edges on the polygon.
	link * headerlink = new link ;
    link * currentLink = headerlink ;
	output * currentOut = headerOut ;
	output * Junc, * currentJunc = NULL ;
	link * trans ;

	do
	{
		if (currentOut->TwoInnerline ==true )
		{
			currentLink->va[0]=currentOut->v1[0];
			currentLink->va[1]=currentOut->v1[1];
			currentLink->vb[0]=currentOut->v2[0];
			currentLink->vb[1]=currentOut->v2[1];
			trans = currentLink ;
			currentLink->nextlink = new link ;
			currentLink = currentLink->nextlink ;
			currentLink->frontlink = trans ;

			currentLink->va[0]=currentOut->v2[0];
			currentLink->va[1]=currentOut->v2[1];
			currentLink->vb[0]=currentOut->v3[0];
			currentLink->vb[1]=currentOut->v3[1];
		}
		else if ( currentOut->TwoInnerline == false )
		{
			currentLink->va[0]=currentOut->v2[0];
			currentLink->va[1]=currentOut->v2[1];
			currentLink->vb[0]=currentOut->v3[0];
			currentLink->vb[1]=currentOut->v3[1];
		}
		trans = currentLink ;
		currentLink->nextlink = new link ;
		currentLink = currentLink->nextlink ;
		currentLink->frontlink = trans ;

		currentOut = currentOut->nextOut ;
	}while (currentOut != NULL );

	trans = currentLink ;
	currentLink = currentLink->frontlink ;	
	currentLink->nextlink  = NULL ;
	delete trans ;
	currentLink = headerlink ;       // complete link: headerlink (store all inner edges in link)

	while ( currentLink != NULL )
	{
		Junc = Junctri( headerlink, currentLink, junctriangle );    // find the junctriangle.
		if (Junc != NULL )
		{
			if ( junctriangle == NULL )
			{
				junctriangle = Junc ;
				junctriangle->nextOut = NULL ;
				currentLink=currentLink->nextlink ;
				currentJunc = junctriangle ; 
				continue ;
			}
			currentJunc->nextOut = Junc ;
			currentJunc = currentJunc->nextOut ;
			currentJunc->nextOut = NULL ;
		}
		currentLink=currentLink->nextlink ;
	}
	if (currentJunc)
	currentJunc->nextOut = NULL ;

	headerlink = deleteLink_whole(headerlink);  //delete headerlink.
	return junctriangle ;

}

//============ Finish Junctriangle===============

//============ Spine ==============

 //reset all output_X to 0
 //need to do this against headerOut and junctriangle.
struct output * reset_X(output* headerOut)  
{    
	output *current = headerOut;
	while(current!=NULL)
	{
		current->X = 0;
		current = current->nextOut ;
	}
	return headerOut;
}


//[0]no found  [1]sleeve tri found  [2]new junctri found  [3]old junctri. 
//[va]and [vb] stores the inner edge needed to determine. 
int determine_if_edge_useable(int va1,int va2,int vb1, int vb2,output *headerOut,output *junctriangle)
{

	output * currentIn = headerOut;     // in headerOut.
	while (currentIn!=NULL)
	{
		if (currentIn->TwoInnerline==true && currentIn->X==0 )  //the sleeve tri should not be used and deleted.
		{
			int la1=currentIn->v1[0];
			int la2=currentIn->v1[1];
			int lb1=currentIn->v2[0];
			int lb2=currentIn->v2[1];
			int lc1=currentIn->v3[0];
			int lc2=currentIn->v3[1];

			if ( (va1==la1&&va2==la2&&vb1==lb1&&vb2==lb2)||(va1==lb1&&va2==lb2&&vb1==la1&&vb2==la2)  //determine wether current sleeve tri share the same edge[va vb].
				||(va1==lc1&&va2==lc2&&vb1==lb1&&vb2==lb2)||(va1==lb1&&va2==lb2&&vb1==lc1&&vb2==lc2))
				return 1;
		}
		currentIn=currentIn->nextOut ;
	}
	currentIn = junctriangle ;       // in junctriangle.
	while (currentIn!=NULL)
	{
		int la1=currentIn->v1[0];
		int la2=currentIn->v1[1];
		int lb1=currentIn->v2[0];
		int lb2=currentIn->v2[1];
		int lc1=currentIn->v3[0];
		int lc2=currentIn->v3[1];

		if ((va1==la1&&va2==la2&&vb1==lb1&&vb2==lb2)||(va1==lb1&&va2==lb2&&vb1==la1&&vb2==la2)  //determine wehter current sleeve tri share the same edge[va vb].
			||(va1==lc1&&va2==lc2&&vb1==lb1&&vb2==lb2)||(va1==lb1&&va2==lb2&&vb1==lc1&&vb2==lc2)
			||(va1==la1&&va2==la2&&vb1==lc1&&vb2==lc2)||(va1==lc1&&va2==lc2&&vb1==la1&&vb2==la2))
		{
			if (currentIn->X!=0)
				return 3;           // the finded junctriangle already be used.
			else
			{
			//	currentIn->X=1;     // label that this junctriangle is used.
				return 2;           // return 2 to demonstriate that 
			}
		}
		currentIn=currentIn->nextOut ;
	}
	return 0;  //there no sutable triangle here, this edge must be an edge of inner polygon.
}

// Input a line with 2 end points, Output a sturct output, with disposalable point at frist, new edge at second and third.
//[0] No sleeve found,  [1] sleeve found,  
struct output * find_next(mypoint* header,output* headerOut, int a1,int a2,int b1,int b2)
{
	 mypoint *current = header;
	 int num_start;
	 int num_end;
	 int num_total;
	 while(current!=NULL)
	 {
		 if (current->x==a1&&current->y==a2)
			 num_start=current->num;
		 if (current->x==b1&&current->y==b2)
			 num_end=current->num;
		 num_total=current->num ;
		 current=current->nextpoint;
	 }
	 if(num_start>num_end)
	 {int trans = num_start;
	 num_start=num_end;
	 num_end=trans;
	 }

	 if (abs(num_start-num_end)>(num_total-num_end+num_start))   //one point at the end  and the other at the beginning.
		// swap(num_start,num_end);
		{int trans = num_start;
	     num_start=num_end;
	     num_end=trans; }

	 current = header;
	 while(current->num!=num_start){     //redefine input start. 
		 current=current->nextpoint;}
	 a1 = current->x ;
	 a2 = current->y ;

	 current = header;
	 while(current->num!=num_end){     //redefine input end. 
		 current=current->nextpoint;}
	 b1 = current->x ;
	 b2 = current->y ;

	 //case 1,start-1, end+0:
	 int start = num_start;
	 int end= num_end;
	 start-=1;

	 if (start==0)
		 start=num_total;

	 current = header;
	 while(current->num!=start){     //find out start value in header.
		 current=current->nextpoint;
	 }
	 int s1=current->x ;
	 int s2=current->y ;

	 output *currentOut = headerOut ;
	 while (currentOut!=NULL)
	 {   // must be a sleeve tri as well as not used.
		 if (currentOut->TwoInnerline==true && currentOut->X==0)  
		 {
			int la1=currentOut->v1[0];
			int la2=currentOut->v1[1];
			int lb1=currentOut->v2[0];
			int lb2=currentOut->v2[1];
			int lc1=currentOut->v3[0];
			int lc2=currentOut->v3[1];

			if ( (s1==la1&&s2==la2&&b1==lb1&&b2==lb2)||(s1==lb1&&s2==lb2&&b1==la1&&b2==la2)  //determine wether current sleeve tri share the same edge[va vb].
				||(s1==lc1&&s2==lc2&&b1==lb1&&b2==lb2)||(s1==lb1&&s2==lb2&&b1==lc1&&b2==lc2))
			{
				output *store = new output;
				store->v1[0]=a1;
				store->v1[1]=a2;     //point need to add to disposal cluster.
				store->v2[0]=s1;
				store->v2[1]=s2;     //new line from num_start-1.
				store->v3[0]=b1;
				store->v3[1]=b2;     // to num_end.
				store->X=1;          // sleeve triangle found.
				return store;
			}
		 }
		 currentOut=currentOut->nextOut ;
	 }// finish start point.


	 //case 2,start-0, end+1:
	 start = num_start;
	 end= num_end;
	 end+=1;

	 if (end>num_total)
		 end=1;

	 current = header;
	 while(current->num!=end){     //find out start value in header.
		 current=current->nextpoint;
	 }
	 int e1=current->x ;
	 int e2=current->y ;

	 currentOut = headerOut ;
	 while (currentOut!=NULL)
	 {
		 if (currentOut->TwoInnerline==true && currentOut->X==0)
		 {
			int la1=currentOut->v1[0];
			int la2=currentOut->v1[1];
			int lb1=currentOut->v2[0];
			int lb2=currentOut->v2[1];
			int lc1=currentOut->v3[0];
			int lc2=currentOut->v3[1];

			if ( (e1==la1&&e2==la2&&a1==lb1&&a2==lb2)||(e1==lb1&&e2==lb2&&a1==la1&&a2==la2)  //determine wether current sleeve tri share the same edge[va vb].
				||(e1==lc1&&e2==lc2&&a1==lb1&&a2==lb2)||(e1==lb1&&e2==lb2&&a1==lc1&&a2==lc2))
			{
				output *store = new output;
				store->v1[0]=b1;
				store->v1[1]=b2;     //point needed to add to disposal cluster.
				store->v2[0]=e1;
				store->v2[1]=e2;     //new line from num_start-1.
				store->v3[0]=a1;
				store->v3[1]=a2;     // to num_start.
				store->X=1;          // sleeve triangle found.
				return store;
			}
		 }
		 currentOut=currentOut->nextOut ;
	 }// finish end point.

	 output *store = new output;
	 store->v1[0]=NULL;
	 store->v1[1]=NULL;    
	 store->v2[0]=NULL;
	 store->v2[1]=NULL;    
	 store->v3[0]=NULL;
	 store->v3[1]=NULL;    
	 store->X=0;          // No sleeve triangle found.
	 return store;
}

// Delete chain of mypoint as whole.
void Delete_mypoint (mypoint * headermy)
{
	mypoint * current = headermy ;
	while ((current = headermy)!=NULL)
	{
		headermy = current->nextpoint ;
		delete current ;
	}
}

//used to add struct into ouput headerOut.
//[8] denotes only line segment stores here.
struct output * Add_inside(output *headerOut,int a1, int a2, float f1, float f2)
{
	output * trans = new output ;
	output * temp ;
	trans->v2[0] = f1;
	trans->v2[1] = f2;
	trans->v1[0] = a1;
	trans->v1[1] = a2;
	trans->X = 8 ;      //[8] denotes only line segment stores here.
	temp = headerOut->nextOut ;
	temp->front = trans ;
	trans->nextOut = temp ;
	trans->front = headerOut ;
	headerOut->nextOut = trans ;
	return headerOut ;
}


//Add [a] and [f] into mypoint * spine2p. 
struct mypoint * Add_mypoint(mypoint *spine2p ,int a1, int a2, float f1, float f2)
{
	if (spine2p!=NULL)
	{
	mypoint * trans = new mypoint ;
	mypoint * temp ;
	trans->xf = f1;
	trans->yf = f2;
	trans->x = a1;
	trans->y = a2;
	temp = spine2p->nextpoint ;
	trans->nextpoint = temp ;
	trans->frontpoint = spine2p ;
	spine2p->nextpoint = trans ;
	if (temp!=NULL)
		temp->frontpoint = trans ;
	return spine2p ;
	}
	if (spine2p==NULL)
	{
		spine2p = new mypoint ;
		spine2p->xf = f1;
		spine2p->yf = f2;
		spine2p->x = a1;
		spine2p->y = a2;
		spine2p->nextpoint = NULL ;
		return spine2p ;
	}
}



// result = the finial no merged edge, with the center points.
// Add_tri just add line segments and disable the used sleeve triangles in headerOut.
// [5] reperesnts disable the very triangle in headerOut.
struct output * Add_tri(mypoint *tempstore,output *headerOut,output * result)
{

	mypoint * tempmy = tempstore ;
	while (tempmy!=NULL)
	{
		output * currentO = headerOut ;
		while (currentO!=NULL)
		{
			if(currentO->TwoInnerline==true && ( (tempmy->x==currentO->v1[0]&&tempmy->y==currentO->v1[1])|| 
				(tempmy->x==currentO->v2[0]&&tempmy->y==currentO->v2[1])||(tempmy->x==currentO->v3[0]&&tempmy->y==currentO->v3[1])))
				
				currentO->X = 5 ;            // [5]    disable this triangle.

			currentO=currentO->nextOut ;
		}
		tempmy = tempmy->nextpoint ;
	}

	//headerOut = Add_inside(headerOut ,result->v1[0],result->v1[1],result->v2[0],result->v2[1]) ;
	//headerOut = Add_inside(headerOut ,result->v3[0],result->v3[1],result->v2[0],result->v2[1]) ;

	spine2p = Add_mypoint(spine2p ,result->v1[0],result->v1[1],result->v2[0],result->v2[1]) ;
	spine2p = Add_mypoint(spine2p ,result->v3[0],result->v3[1],result->v2[0],result->v2[1]) ;


	tempmy = tempstore ;
	while (tempmy!=NULL)
	{
		//headerOut = Add_inside(headerOut ,tempmy->x,tempmy->y ,result->v2[0],result->v2[1]) ;
		spine2p = Add_mypoint(spine2p ,tempmy->x,tempmy->y ,result->v2[0],result->v2[1]) ;
		tempmy = tempmy->nextpoint ;
	}
	return headerOut;
}

// Only disable those merged sleeve triangles. if delete, set [5] to sleeve in headerOut.
struct output * Disable_tri(mypoint *tempstore,output *headerOut)
{
	mypoint * tempmy = tempstore ;
	while (tempmy!=NULL)
	{
		output * currentO = headerOut ;
		while (currentO!=NULL)
		{   // if the triangle is sleeve, mean while the triangle contains points in tempstore. Thus, we know it should be disable.
			if(currentO->TwoInnerline==true && ( (tempmy->x==currentO->v1[0]&&tempmy->y==currentO->v1[1])|| 
				(tempmy->x==currentO->v2[0]&&tempmy->y==currentO->v2[1])||(tempmy->x==currentO->v3[0]&&tempmy->y==currentO->v3[1])))
				
				currentO->X = 5 ;            // [5]    disable this triangle.

			currentO=currentO->nextOut ;
		}
		tempmy = tempmy->nextpoint ;
	}
	return headerOut ;
}


//Used to store those FAN triangles on the terminal region.
//X = [1] denote Only ONE spine point here.
//spineTri-> [v1]  [v3] is edge points  [v4] is the spine point or so called center point.
void Save2spinetri(output * result , mypoint * header , float c1, float c2)
{
	int a1 = result->v1[0];
	int a2 = result->v1[1];
	int b1 = result->v3[0];
	int b2 = result->v3[1];
	int num_v1;
	int num_v3;
	int NUM_end = end_mypoint->num ;
	int temp;
	mypoint * current = header ;

	while ( current!=NULL )   // find v1 and v3's corresponding NUM.
	{
		if (a1==current->x && a2==current->y)
			num_v1 = current->num;
		if (b1==current->x && b2==current->y)
			num_v3 = current->num;
		current = current->nextpoint ;
	}
	if (num_v1>num_v3)
		swap(num_v1,num_v3);
	if ( abs(num_v3-num_v1)>(NUM_end-num_v3+num_v1))   // so that make sure num_v1 always be the START point and num_v3 be the END.
		swap(num_v1,num_v3);
	while(num_v1!=num_v3)
	{
		current = header ;
		temp = num_v1;
		num_v1++;
		if (num_v1>NUM_end)
			num_v1 = 1;
		while (current!=NULL)
		{
			if (current->num == temp )
			{
				a1 = current->x ;
				a2 = current->y ;
			}
			if (current->num == num_v1)
			{
				b1 = current->x ;
				b2 = current->y ;
			}
			current = current->nextpoint ;
		}
		if (spineTri!=NULL)
		{
			output * trans = new output ;
	        output * tempOut ;
	        trans->v1[0] = a1;
			trans->v1[1] = a2;
			trans->v3[0] = b1;
			trans->v3[1] = b2;
			trans->v4[0] = c1;
			trans->v4[1] = c2;
			trans->X = 1 ;
	        tempOut = spineTri->nextOut ;
	        trans->nextOut = tempOut ;
	        trans->front = spineTri ;
	        spineTri->nextOut = trans ;
	        if (tempOut!=NULL)
				tempOut->front = trans ;
		}

		if (spineTri==NULL)
		{
			spineTri = new output ;
			spineTri->v1[0] = a1;
			spineTri->v1[1] = a2;
			spineTri->v3[0] = b1;
			spineTri->v3[1] = b2;
			spineTri->v4[0] = c1;
			spineTri->v4[1] = c2;
			spineTri->X = 1 ; // [1]denotes only ONE spine point in this tri.
			spineTri->nextOut = NULL ;
		}
	
	}
}



// Store 3 points in spineTri, 2 must be float points(points on spine [spineOut])
// X = [2] denotes that this tri contains 2 spine points. 
// TWO spine Point in [v2] and [v4] ,   ONE edge point stores in [v1] .
void Save2spinetri_2p(int c1 , int c2, float b1 , float b2 , float a1, float a2)
{
	if (spineTri!=NULL)
		{
			output * trans = new output ;
	        output * tempOut ;
	        trans->v2[0] = a1;
			trans->v2[1] = a2;
			trans->v4[0] = b1;
			trans->v4[1] = b2;
			trans->v1[0] = c1;
			trans->v1[1] = c2;
			trans->X = 2 ;                       // X = [2] denotes that this tri contains 2 spine points
	        tempOut = spineTri->nextOut ;
	        trans->nextOut = tempOut ;
	        trans->front = spineTri ;
	        spineTri->nextOut = trans ;
	        if (tempOut!=NULL)
				tempOut->front = trans ;
		}

		if (spineTri==NULL)
		{
			spineTri = new output ;
			spineTri->v2[0] = a1;
			spineTri->v2[1] = a2;
			spineTri->v4[0] = b1;
			spineTri->v4[1] = b2;
			spineTri->v1[0] = c1;
			spineTri->v1[1] = c2;
			spineTri->X = 2 ; 
			spineTri->nextOut = NULL ;
		}
}



// The same as Save2spinetri_2p exclude only contains 1 spine point.
// X = [1] denotes ONE spine point.
// [v1] and [v3] are edge points, [v4] is float spine point. 
void Save2spinetri_1p(int a1, int a2, int b1 , int b2 , float c1 , float c2)
{
	if (spineTri!=NULL)
		{
			output * trans = new output ;
	        output * tempOut ;
	        trans->v1[0] = a1;
			trans->v1[1] = a2;
			trans->v3[0] = b1;
			trans->v3[1] = b2;
			trans->v4[0] = c1;
			trans->v4[1] = c2;
			trans->X = 1 ;                       // X = [1] denotes that this tri contains 1 spine points
	        tempOut = spineTri->nextOut ;
	        trans->nextOut = tempOut ;
	        trans->front = spineTri ;
	        spineTri->nextOut = trans ;
	        if (tempOut!=NULL)
				tempOut->front = trans ;
		}

		if (spineTri==NULL)
		{
			spineTri = new output ;
			spineTri->v1[0] = a1;
			spineTri->v1[1] = a2;
			spineTri->v3[0] = b1;
			spineTri->v3[1] = b2;
			spineTri->v4[0] = c1;
			spineTri->v4[1] = c2;
			spineTri->X = 1 ; 
			spineTri->nextOut = NULL ;
		}
}



//for a given terminal triangle, do it!
//return an output struct contains information necessary. v1 and v3 are edge, v2 are centerP, if[9],then v2=NULL.
//[9]  denotes inner Ploygon.{with tempstore combined}  [0] denots OK   [2] is about merge with junctri
//[3] merge with old junctri~ 
struct output * refine_terminal(output* headerOut,output* junctriangle,mypoint* header,output *terminal)
{
	//Preprocessing.
	output * current = headerOut ;
	mypoint * tempstore = new mypoint;  //store points between va and vb.
	output *result = new output ;    // store a center point and its edge as a result.

	int d1 = terminal->v1[0];
	int d2 = terminal->v1[1];     //dispossal point, should store in tempstore.
	int a1 = terminal->v2[0];
	int a2 = terminal->v2[1]; 
	int b1 = terminal->v3[0];
	int b2 = terminal->v3[1];     // end points of the inner edge.  [a] and [b]
	tempstore->x = d1;
	tempstore->y = d2;
	tempstore->nextpoint = NULL ;
	mypoint *currentT = tempstore ;
	int count = 0;

	//Main loop.
	while (1)
	{
		float centerP1 = ((float)a1+(float)b1)/2;    //have to reflash a and b at the end.
		float centerP2 = ((float)a2+(float)b2)/2;

		currentT = tempstore ;
		bool need_to_merge=true ;
		while (currentT!=NULL)      //search points in tempstore to see if need to merge this edge.
		{
			if (
			sqrt(pow((currentT->x - centerP1),2)+ pow( (currentT->y - centerP2),2))>
			(sqrt(((float)a1-(float)b1)*((float)a1-(float)b1)+((float)a2-(float)b2)*((float)a2-(float)b2))/2))
			{
			need_to_merge = false ;
			break;}
			currentT = currentT->nextpoint ;
		}

		if (need_to_merge == false)
		{
			result->v1[0]= a1;
			result->v1[1]= a2;
			result->v2[0]= centerP1;
			result->v2[1]= centerP2;
			result->v3[0]= b1;
			result->v3[1]= b2 ;
			result->storemy = tempstore ;
			result->nextOut = NULL ;
			result->X = 0 ;            // every think OK.
			headerOut = Add_tri(tempstore,headerOut,result) ;
			Save2spinetri(result, header, centerP1 , centerP2);   // save to spineTri as final tri mesh in 2D.

			if (spineTri==NULL)
				spineTri=NULL;


			//Delete_mypoint(tempstore);    
			return result ;
		}

		headerOut = Disable_tri(tempstore, headerOut);  //reflesh headerOut.

		//[0]no found  [1]sleeve tri found  [2]new junctri found  [3]old junctri. 
		int Type = determine_if_edge_useable(a1,a2,b1,b2,headerOut,junctriangle);

		if (Type==0)           //[0]no found
		{
			result->v1[0]=a1 ;
			result->v1[1]= a2 ;
			result->v2[0]= NULL;
			result->v2[1]= NULL;
			result->v3[0]= b1 ;
			result->v3[1]= b2 ;
			result->X = 9 ;    //  [9]  denotes that this result contains an edge of inner Ploygon.
			result->storemy = tempstore ;     // tempstore will be used in future.
			result->nextOut = NULL ;
			//Delete_mypoint(tempstore);    // No need to delete tempstore here,cause will use in future.
			return result ;
		}

		if (Type==2)
		{
			output * currentJunc = junctriangle;
			while (currentJunc!=NULL)
			{
				int la1=currentJunc->v1[0];
		        int la2=currentJunc->v1[1];
		        int lb1=currentJunc->v2[0];
		        int lb2=currentJunc->v2[1];
		        int lc1=currentJunc->v3[0];
		        int lc2=currentJunc->v3[1];
				//determine wether current edge [a b] share the same edge in Junctri.
		        if ((a1==la1&&a2==la2&&b1==lb1&&b2==lb2)||(a1==lb1&&a2==lb2&&b1==la1&&b2==la2)  
			        ||(a1==lc1&&a2==lc2&&b1==lb1&&b2==lb2)||(a1==lb1&&a2==lb2&&b1==lc1&&b2==lc2)
			        ||(a1==la1&&a2==la2&&b1==lc1&&b2==lc2)||(a1==lc1&&a2==lc2&&b1==la1&&b2==la2))
					break;
				currentJunc=currentJunc->nextOut ;
			}
			currentJunc->X = 1; //make sure other than 0 so that we know it's been used.
			int la1=currentJunc->v1[0];
		    int la2=currentJunc->v1[1];
		    int lb1=currentJunc->v2[0];
		    int lb2=currentJunc->v2[1];
		    int lc1=currentJunc->v3[0];
		    int lc2=currentJunc->v3[1];

			output * JuncEdges = new output;
			if ((a1==la1&&a2==la2&&b1==lb1&&b2==lb2)||(a1==lb1&&a2==lb2&&b1==la1&&b2==la2))
			{
				JuncEdges->v1[0]=la1;
				JuncEdges->v1[1]=la2;
				JuncEdges->v2[0]=lc1;
				JuncEdges->v2[1]=lc2;
				JuncEdges->v3[0]=lb1;
				JuncEdges->v3[1]=lb2;
			}
			if ((a1==lc1&&a2==lc2&&b1==lb1&&b2==lb2)||(a1==lb1&&a2==lb2&&b1==lc1&&b2==lc2))
			{
				JuncEdges->v1[0]=lb1;
				JuncEdges->v1[1]=lb2;
				JuncEdges->v2[0]=la1;
				JuncEdges->v2[1]=la2;
				JuncEdges->v3[0]=lc1;
				JuncEdges->v3[1]=lc2;
			}
			if ((a1==la1&&a2==la2&&b1==lc1&&b2==lc2)||(a1==lc1&&a2==lc2&&b1==la1&&b2==la2))
			{
				JuncEdges->v1[0]=la1;
				JuncEdges->v1[1]=la2;
				JuncEdges->v2[0]=lb1;
				JuncEdges->v2[1]=lb2;
				JuncEdges->v3[0]=lc1;
				JuncEdges->v3[1]=lc2;
			}
			float jc1 = (la1+lb1)/2;
			float jc2 = (la2+lb2)/2;

			float Jcenter1 = (2*jc1+lc1)/3;
			float Jcenter2 = (2*jc2+lc2)/3;    // calculate the center of Junctriangle.

			//mypoint * JuncCenter = new mypoint;
			//JuncCenter->xf = Jcenter1;
			//JuncCenter->yf = Jcenter2;

			result->v1[0]= a1;
			result->v1[1]= a2;
			result->v2[0]= Jcenter1;
			result->v2[1]= Jcenter2;
			result->v3[0]= b1;
			result->v3[1]= b2 ;
			result->storemy = tempstore ;
			result->X = 2 ;            // this is about merge with junctriangle.
			result->nextOut = NULL ;


			currentJunc->store = JuncEdges;  // the free edges are stored in junctriangle.
			currentJunc->v4[0]=Jcenter1;
			currentJunc->v4[1]=Jcenter2;       // the junc center is store here.

			Save2spinetri(result, header, Jcenter1 , Jcenter2); 
			headerOut = Add_tri(tempstore,headerOut,result) ;
			//Delete_mypoint(tempstore);
			return result ;
		}

		if (Type==3)
		{
			output * currentJunc = junctriangle;
			while (currentJunc!=NULL)
			{
				int la1=currentJunc->v1[0];
		        int la2=currentJunc->v1[1];
		        int lb1=currentJunc->v2[0];
		        int lb2=currentJunc->v2[1];
		        int lc1=currentJunc->v3[0];
		        int lc2=currentJunc->v3[1];
				//determine wether current edge [a b] share the same edge in Junctri.
		        if ((a1==la1&&a2==la2&&b1==lb1&&b2==lb2)||(a1==lb1&&a2==lb2&&b1==la1&&b2==la2)  
			        ||(a1==lc1&&a2==lc2&&b1==lb1&&b2==lb2)||(a1==lb1&&a2==lb2&&b1==lc1&&b2==lc2)
			        ||(a1==la1&&a2==la2&&b1==lc1&&b2==lc2)||(a1==lc1&&a2==lc2&&b1==la1&&b2==la2))
					break;
				currentJunc=currentJunc->nextOut ;
			}
			//mypoint * JuncCenter = currentJunc->store ; 
			currentJunc->X = 2;  // denotes that this junctri is used twice.
			output * juncStore = currentJunc->store ;

			int la1=juncStore->v1[0];
		    int la2=juncStore->v1[1];
		    int lb1=juncStore->v2[0];
		    int lb2=juncStore->v2[1];
		    int lc1=juncStore->v3[0];
		    int lc2=juncStore->v3[1];

			if( (a1==la1&&a2==la2&&b1==lb1&&b2==lb2)||(a1==lb1&&a2==lb2&&b1==la1&&b2==la2) )
			{
				juncStore->v1[0]= lb1;
				juncStore->v1[1]= lb2;
				juncStore->v2[0]= lc1;
				juncStore->v2[1]= lc2;
				juncStore->v3[0]= NULL;
				juncStore->v3[1]= NULL;
			}
			if ( (a1==lc1&&a2==lc2&&b1==lb1&&b2==lb2)||(a1==lb1&&a2==lb2&&b1==lc1&&b2==lc2) )
			{
				juncStore->v1[0]= la1;
				juncStore->v1[1]= la2;
				juncStore->v2[0]= lb1;
				juncStore->v2[1]= lb2;
				juncStore->v3[0]= NULL;
				juncStore->v3[1]= NULL;
			}

			result->v1[0]= a1;
			result->v1[1]= a2;
			result->v2[0]= currentJunc->v4[0] ;
			result->v2[1]= currentJunc->v4[1] ;
			result->v3[0]= b1;
			result->v3[1]= b2 ;
			result->storemy = tempstore ;
			result->X = 3 ;            // this is about merge with old junctriangle.
			result->nextOut = NULL ;

			Save2spinetri(result, header, currentJunc->v4[0] , currentJunc->v4[1]); 
			headerOut = Add_tri(tempstore,headerOut,result) ;
			//Delete_mypoint(tempstore);
			return result ;
		}

		if ( Type==1 )
		{
			output * Newfind = find_next(header, headerOut, a1, a2, b1, b2);

			mypoint * tempstoreNew = new mypoint ;
			tempstoreNew->x = Newfind->v1[0];
			tempstoreNew->y = Newfind->v1[1];

			tempstoreNew->nextpoint = tempstore ;
			tempstore = tempstoreNew ;             // increase tempstore(points between [a] and [b] +1)

			a1 = Newfind->v2[0];
			a2 = Newfind->v2[1];
			b1 = Newfind->v3[0];
			b2 = Newfind->v3[1];
		}

		count +=1;             // see if infinit loop.
		cout << count <<endl;
  	} // end while(1).

}

// go through headerOut to search for terminal tri, calculate its corresponding result and link them as a chain.
struct output * termianl_chain(output* headerOut,output* junctriangle,mypoint* header)
{
	output * terminal = headerOut ;
	output *result = NULL;
	output *headerresult = NULL;
	while (terminal!=NULL)
	{
		if (terminal->TwoInnerline==false)
		{
			result = refine_terminal(headerOut,junctriangle,header,terminal);
		}
		if (headerresult==NULL && result!=NULL)
		{
			headerresult = result;
			result = NULL ;
			terminal=terminal->nextOut ;
			continue ;
		}
		if (headerresult!=NULL && result!=NULL)
		{
			result->nextOut = headerresult ;
			headerresult->front = result ;
			headerresult = result ;
			result = NULL ;
		}
		terminal = terminal->nextOut ;
	}
	if (spine2p==NULL)
		spine2p=NULL ;
	return headerresult;
}


// get spine from junctriangle only.  spineOut{}
// Do this only if their values in junctriangle.
// Also save   Spine2p{}   and   SpineTri{}   about Junctri. 
struct output * refine_junctriangle(output *junctriangle, output * headerOut)
{
	if (junctriangle==NULL)
		return NULL;
	output * current = junctriangle ;
	output * spineOut = new output ;
	output * currentS = spineOut ;
	output * trans ;

	while (current!=NULL)
	{
		if (current->X ==0 )   // no use junctri.
		{
			int la1=current->v1[0];
		    int la2=current->v1[1];
		    int lb1=current->v2[0];
		    int lb2=current->v2[1];
		    int lc1=current->v3[0];
		    int lc2=current->v3[1];

			float jc1 = ((float)la1+(float)lb1)/2;
			float jc2 = ((float)la2+(float)lb2)/2;

			float ja1 = ((float)lb1+(float)lc1)/2;
			float ja2 = ((float)lb2+(float)lc2)/2;

			float jb1 = ((float)la1+(float)lc1)/2;
			float jb2 = ((float)la2+(float)lc2)/2;

			//float Jcenter1 = current->v4[0] ;
			//float Jcenter2 = current->v4[1] ;

			float Jcenter1 = (2*jc1+lc1)/3;
			float Jcenter2 = (2*jc2+lc2)/3;    // calculate the center of Junctriangle.

			//headerOut = Add_inside(headerOut ,la1,la2,Jcenter1,Jcenter2) ;
			//headerOut = Add_inside(headerOut ,lb1,lb2,Jcenter1,Jcenter2) ;
			//headerOut = Add_inside(headerOut ,lc1,lc2,Jcenter1,Jcenter2) ;
			spine2p = Add_mypoint(spine2p ,la1,la2,Jcenter1,Jcenter2) ;
			spine2p = Add_mypoint(spine2p ,lb1,lb2,Jcenter1,Jcenter2) ;
			spine2p = Add_mypoint(spine2p ,lc1,lc2,Jcenter1,Jcenter2) ;     // store 2 points spine2edge line.

			Save2spinetri_2p(la1, la2, jc1 , jc2 , Jcenter1 , Jcenter2 );
			Save2spinetri_2p(la1, la2, jb1 , jb2 , Jcenter1 , Jcenter2 );
			Save2spinetri_2p(lb1, lb2, jc1 , jc2 , Jcenter1 , Jcenter2 );
			Save2spinetri_2p(lb1, lb2, ja1 , ja2 , Jcenter1 , Jcenter2 );
			Save2spinetri_2p(lc1, lc2, ja1 , ja2 , Jcenter1 , Jcenter2 );
			Save2spinetri_2p(lc1, lc2, jb1 , jb2 , Jcenter1 , Jcenter2 );   // store final 2D spineTri.


			currentS->v2[0]=jc1; 
			currentS->v2[1]=jc2; 
			currentS->v4[0]=Jcenter1; 
			currentS->v4[1]=Jcenter2;
			trans = currentS;
			currentS->nextOut = new output ;
			currentS=currentS->nextOut ;
			currentS->front=trans;

			currentS->v2[0]=ja1; 
			currentS->v2[1]=ja2; 
			currentS->v4[0]=Jcenter1; 
			currentS->v4[1]=Jcenter2;
			trans = currentS;
			currentS->nextOut = new output ;
			currentS=currentS->nextOut ;
			currentS->front=trans;

			currentS->v2[0]=jb1; 
			currentS->v2[1]=jb2; 
			currentS->v4[0]=Jcenter1; 
			currentS->v4[1]=Jcenter2;
			trans = currentS;
			currentS->nextOut = new output ;
			currentS=currentS->nextOut ;
			currentS->front=trans;
		}
		if ( current->X ==1 )   //use once junctri.
		{
			output * JuncEdges = current->store ;

			int la1=JuncEdges->v1[0];
		    int la2=JuncEdges->v1[1];
		    int lb1=JuncEdges->v2[0];
		    int lb2=JuncEdges->v2[1];
		    int lc1=JuncEdges->v3[0];
		    int lc2=JuncEdges->v3[1];

			float jc1 = ((float)la1+(float)lb1)/2;
			float jc2 = ((float)la2+(float)lb2)/2;

			float ja1 = ((float)lb1+(float)lc1)/2;
			float ja2 = ((float)lb2+(float)lc2)/2;

			float Jcenter1 = current->v4[0] ;
			float Jcenter2 = current->v4[1] ;

			//headerOut = Add_inside(headerOut ,lb1,lb2,Jcenter1,Jcenter2) ;
			spine2p = Add_mypoint(spine2p ,lb1,lb2,Jcenter1,Jcenter2) ;

			Save2spinetri_2p(la1, la2, jc1 , jc2 , Jcenter1 , Jcenter2 );
			Save2spinetri_2p(lb1, lb2, jc1 , jc2 , Jcenter1 , Jcenter2 );
			Save2spinetri_2p(lb1, lb2, ja1 , ja2 , Jcenter1 , Jcenter2 );
			Save2spinetri_2p(lc1, lc2, ja1 , ja2 , Jcenter1 , Jcenter2 );


			currentS->v2[0]=jc1; 
			currentS->v2[1]=jc2; 
			currentS->v4[0]=Jcenter1; 
			currentS->v4[1]=Jcenter2;
			trans = currentS;
			currentS->nextOut = new output ;
			currentS=currentS->nextOut ;
			currentS->front=trans;

			currentS->v2[0]=ja1; 
			currentS->v2[1]=ja2; 
			currentS->v4[0]=Jcenter1; 
			currentS->v4[1]=Jcenter2;
			trans = currentS;
			currentS->nextOut = new output ;
			currentS=currentS->nextOut ;
			currentS->front=trans;
		}
		if ( current->X ==2 )   // use twice junctri.
		{
			output * JuncEdges = current->store ;

			int la1=JuncEdges->v1[0];
		    int la2=JuncEdges->v1[1];
		    int lb1=JuncEdges->v2[0];
		    int lb2=JuncEdges->v2[1];

			float jc1 = ((float)la1+(float)lb1)/2;
			float jc2 = ((float)la2+(float)lb2)/2;

			float Jcenter1 = current->v4[0] ;
			float Jcenter2 = current->v4[1] ;

			Save2spinetri_2p(la1, la2, jc1 , jc2 , Jcenter1 , Jcenter2 );
			Save2spinetri_2p(lb1, lb2, jc1 , jc2 , Jcenter1 , Jcenter2 );

			currentS->v2[0]=jc1; 
			currentS->v2[1]=jc2; 
			currentS->v4[0]=Jcenter1; 
			currentS->v4[1]=Jcenter2;
			trans = currentS;
			currentS->nextOut = new output ;
			currentS=currentS->nextOut ;
			currentS->front=trans;
		}
		
		current = current->nextOut ;
		}
		trans->nextOut = NULL;
		delete currentS ;

		return spineOut ;
}



// Input 2 value (center of 2 lines), this func will store them in spineOut.
// If first use, make sure spineOut = NULL.
struct output * Add2spine(float x1, float x2,float y1, float y2, output * spineOut)
{
	if (spineOut!=NULL)
	{
		output * currentS = new output ;
	    output * trans = spineOut->nextOut;

	        currentS->v2[0]=x1; 
			currentS->v2[1]=x2; 
			currentS->v4[0]=y1; 
			currentS->v4[1]=y2;
			currentS->nextOut = trans;
			currentS->front = spineOut;
		    spineOut->nextOut = currentS ;
			if (trans!=NULL)
				trans->front =currentS;
			return spineOut ;
	}
	if (spineOut==NULL)
	{
		output * spine =new output;
		spine->v2[0] = x1;
		spine->v2[1] = x2;
		spine->v4[0] = y1;
		spine->v4[1] = y2;
		spine->nextOut = NULL;
		return spine ;
	}
}

// store inner polygon line to mypoint.
struct mypoint * Add2inploy(int x1, int x2,int y1, int y2, mypoint * spineOut)
{
	if (spineOut!=NULL)
	{
		mypoint * currentS = new mypoint ;
	    mypoint * trans = spineOut->nextpoint ;

	        currentS->x =x1; 
			currentS->y =x2; 
			currentS->x2=y1; 
			currentS->y2=y2;
			currentS->index=0;
			currentS->nextpoint = trans;
			currentS->frontpoint = spineOut;
		    spineOut->nextpoint = currentS ;
			if (trans!=NULL)
				trans->frontpoint =currentS;
			return spineOut ;
	}
	if (spineOut==NULL)
	{
		mypoint * spine = new mypoint;
		spine->x = x1;
		spine->y = x2;
		spine->x2 = y1;
		spine->y2 = y2;
		spine->index=0;
		spine->nextpoint = NULL;
		return spine ;
	}
}


// just see if the edge[a b] can find its companion.
bool iflineok(int a1, int a2, int b1, int b2, output * headerOut, output *junctriangle)
{
	output * current =headerOut ;
	while (current!=NULL)
	{
		int la1=current->v1[0];
		int la2=current->v1[1];
		int lb1=current->v2[0];
		int lb2=current->v2[1];
		int lc1=current->v3[0];
		int lc2=current->v3[1];
		                                       //I make the be checked tri'X = [7] so that no need to check itself.
		if (current->TwoInnerline==true && current->X!=7)  //see if can find sleeve.  {also may junc or terminal}
		{
			//determine wether current edge [a b] share the same edge in Junctri.{for sleeve, 2 edge is enough~}
		    if ((a1==la1&&a2==la2&&b1==lb1&&b2==lb2)||(a1==lb1&&a2==lb2&&b1==la1&&b2==la2)  
			    ||(a1==lc1&&a2==lc2&&b1==lb1&&b2==lb2)||(a1==lb1&&a2==lb2&&b1==lc1&&b2==lc2)
			    ||(a1==la1&&a2==la2&&b1==lc1&&b2==lc2)||(a1==lc1&&a2==lc2&&b1==la1&&b2==la2))
				return true;
		}
		if (current->TwoInnerline==false)  //see if its a terminal edge.
		{
			if ((a1==lc1&&a2==lc2&&b1==lb1&&b2==lb2)||(a1==lb1&&a2==lb2&&b1==lc1&&b2==lc2))
				return true;
		}
		current = current->nextOut ;
	}

	current = junctriangle;
	while( current!=NULL )    //see if in the junctirangle.
	{
		int la1=current->v1[0];
		int la2=current->v1[1];
		int lb1=current->v2[0];
		int lb2=current->v2[1];
		int lc1=current->v3[0];
		int lc2=current->v3[1];

		if ((a1==la1&&a2==la2&&b1==lb1&&b2==lb2)||(a1==lb1&&a2==lb2&&b1==la1&&b2==la2)  
			||(a1==lc1&&a2==lc2&&b1==lb1&&b2==lb2)||(a1==lb1&&a2==lb2&&b1==lc1&&b2==lc2)
			||(a1==la1&&a2==la2&&b1==lc1&&b2==lc2)||(a1==lc1&&a2==lc2&&b1==la1&&b2==la2))
			return true;
		current = current->nextOut ;
	}
	return false ;
}



// just see if the edge[a b] on terminal tri can find its companion.
bool ifterminalok(int a1, int a2, int b1, int b2, output * headerOut, output *junctriangle)
{
	output * current =headerOut ;
	while (current!=NULL)
	{
		int la1=current->v1[0];
		int la2=current->v1[1];
		int lb1=current->v2[0];
		int lb2=current->v2[1];
		int lc1=current->v3[0];
		int lc2=current->v3[1];
		
		if (current->TwoInnerline==true && current->X!=7)  //see if can find sleeve.  {also may junc or terminal}
		{
			//determine wether current edge [a b] share the same edge in Junctri.{for sleeve, 2 edge is enough~}
		    if ((a1==la1&&a2==la2&&b1==lb1&&b2==lb2)||(a1==lb1&&a2==lb2&&b1==la1&&b2==la2)  
			    ||(a1==lc1&&a2==lc2&&b1==lb1&&b2==lb2)||(a1==lb1&&a2==lb2&&b1==lc1&&b2==lc2)
			    ||(a1==la1&&a2==la2&&b1==lc1&&b2==lc2)||(a1==lc1&&a2==lc2&&b1==la1&&b2==la2))
				return true;
		}
		current = current->nextOut ;
	}

	current = junctriangle;
	while( current!=NULL )    //see if in the junctirangle.
	{
		int la1=current->v1[0];
		int la2=current->v1[1];
		int lb1=current->v2[0];
		int lb2=current->v2[1];
		int lc1=current->v3[0];
		int lc2=current->v3[1];

		if ((a1==la1&&a2==la2&&b1==lb1&&b2==lb2)||(a1==lb1&&a2==lb2&&b1==la1&&b2==la2)  
			||(a1==lc1&&a2==lc2&&b1==lb1&&b2==lb2)||(a1==lb1&&a2==lb2&&b1==lc1&&b2==lc2)
			||(a1==la1&&a2==la2&&b1==lc1&&b2==lc2)||(a1==lc1&&a2==lc2&&b1==la1&&b2==la2))
			return true;
		current = current->nextOut ;
	}
	return false ;
}


//search whole headerOut to find out inside polygon edges.
//Polygon Edges store in header  InPoly.
//If no founded at all, will return directly InPoly = NULL.
struct mypoint * findinpoly (output * headerOut,output * junctriangle)
{
	output * currentO = headerOut ;
	mypoint * InPoly = NULL;
	bool type =true;
	while (currentO!=NULL)
	{
		type = true ;                         // X!=8 because [8] is not tri but line seg!
		if ( currentO->TwoInnerline==true && currentO->X!=8 ) // sleeve that haven't been delete.
		{
			if (currentO->X ==5)
				type = false;
			currentO->X = 7 ;  // set to [7], but need to retrun to 0 at the end.
			int a1 = currentO->v1[0];
			int a2 = currentO->v1[1];
			int b1 = currentO->v2[0];
			int b2 = currentO->v2[1];
			if (!iflineok(a1,a2,b1,b2,headerOut,junctriangle))   // this is the inner Polygon line!
				InPoly = Add2inploy(a1,a2,b1,b2,InPoly);  // store inner polygon line to mypoint, return the header.
			 a1 = currentO->v2[0];
			 a2 = currentO->v2[1];
			 b1 = currentO->v3[0];
			 b2 = currentO->v3[1];
			if (!iflineok(a1,a2,b1,b2,headerOut,junctriangle))   // check another line of this sleeve tri.
				InPoly = Add2inploy(a1,a2,b1,b2,InPoly);
			if (type==false){
				currentO->X = 5;}
			else{
				currentO->X = 0 ;}
		}                                     // X!=8 because [8] is not tri but line seg!
		if ( currentO->TwoInnerline==false && currentO->X!=8)
		{
			int a1 = currentO->v2[0];
			int a2 = currentO->v2[1];
			int b1 = currentO->v3[0];
			int b2 = currentO->v3[1];
			if (!ifterminalok(a1,a2,b1,b2,headerOut,junctriangle))   // check another line of this terminal tri.
				InPoly = Add2inploy(a1,a2,b1,b2,InPoly);
		}
		currentO = currentO->nextOut ;
	}
	return InPoly ;
}


// mypoint * InPoly = findinpoly(headerOut, junctriangle);
// assuming InPoly is generated.
// if found, return a chain.
// if not, return NULL.
struct mypoint * OneChainPoly(output * headerOut, output * junctriangle, mypoint * InPoly)
{
	mypoint * current = InPoly ;
	mypoint * PolyChain = NULL ;
	while (current!=NULL)
	{
		if (current->index==0)    // if available , index will remain 0.
			break;
		current=current->nextpoint;
	}
	if (current==NULL)    //if no value could use, return NULL. 
		return NULL ;

	int A1=current->x;
	int A2=current->y;
	int B1=current->x2;
	int B2=current->y2;
	current->index=1;  //remove from searching base itself.

	PolyChain = Add2inploy(B1,B2,A1,A2,PolyChain);  // first element of polygon store here.
	while (1)
	{
		mypoint * currentIn = InPoly ;
		while (currentIn!=NULL)
		{
			if ((currentIn->x==A1&&currentIn->y==A2)&&currentIn->index==0)
			{
				PolyChain = Add2inploy(currentIn->x,currentIn->y,currentIn->x2,currentIn->y2,PolyChain);
				A1 = currentIn->x2 ;
				A2 = currentIn->y2 ;
				currentIn->index = 1; // if fit, then remove this point.
				break ;
			}
			if ((currentIn->x2==A1&&currentIn->y2==A2)&&currentIn->index==0)
			{
				PolyChain = Add2inploy(currentIn->x2,currentIn->y2,currentIn->x,currentIn->y,PolyChain);
				A1 = currentIn->x ;
				A2 = currentIn->y ;
				currentIn->index = 1; // if fit, then remove this point.
				break ;
			}
			currentIn = currentIn->nextpoint ;
		}
		if(A1==B1 && A2==B2)
			break;
	}
	return PolyChain ;
}

// return a chain of polygon, if none , return NULL .
struct mypoint * ChainofPolyChain(output * headerOut, output * junctriangle )
{
	mypoint * InPoly = findinpoly(headerOut, junctriangle );
	if ( InPoly==NULL)
		return NULL ;
	mypoint * Chain_Polychain = new mypoint ;
	mypoint * current = Chain_Polychain;
	mypoint * trans ;
	int num=0;
	while (1)
	{
		mypoint * OneChain =  OneChainPoly(headerOut, junctriangle, InPoly );
		if (OneChain==NULL)
			break;
		num++;
		current->store = OneChain;
		current->num =num;
		current->nextpoint= new mypoint ;
		trans = current ;
		current=current->nextpoint ;
		current->frontpoint = trans ;
		current->nextpoint = NULL ;
	}
	trans = current->frontpoint ;
	trans->nextpoint = NULL ;
	delete current ;
	return Chain_Polychain ; 
}


//For a given PolyChain, find its center, search through result{} to disable merged line[9]
//store outside point to center in headerOut [8] -> line segment ;
//then connect center to each other lines' center in PolyChian, store them to spineOut.
//return spineOut. but must have value first.
struct output * calculate_Inpoly(mypoint * PolyChain , output * headerresult, output * headerOut,output * spineOut)
{
	mypoint * current = PolyChain ;  //PolyChain: [x y] [x2 y2]  index=0 ;
	float a1 = current->x ;
	float a2 = current->y ;
	float b1 = current->x2 ;
	float b2 = current->y2 ;

	current=current->nextpoint->nextpoint;

	float c1 = current->x ;
	float c2 = current->y ;
	float d1 = current->x2 ;
	float d2 = current->y2 ;

	float Cab1 = (a1+b1)/2;
	float Cab2 = (a2+b2)/2;

	float Ccd1 = (c1+d1)/2;
	float Ccd2 = (c2+d2)/2;

	float Cbc1 = (b1+c1)/2;
	float Cbc2 = (b2+c2)/2;

	float Cad1 = (a1+d1)/2;
	float Cad2 = (a2+d2)/2;

	float M = (Cab2-Ccd2)/(Cab1-Ccd1);
	float N = (Cad2-Cbc2)/(Cad1-Cbc1);

	float CenterX = (M*Ccd1-N*Cbc1+Cbc2-Ccd2)/(M-N);
	float CenterY = ((1/M)*Ccd2-(1/N)*Cbc2+Cbc1-Ccd1)/(1/M-1/N);  // find the center of Polygon.


	output * currentR = headerresult ;
	mypoint * tempmy ;
	while (currentR!=NULL)
	{
		if (currentR->X ==9)
		{
			tempmy = currentR->storemy ;    // Add line seg to headerOut.
          //  headerOut = Add_inside(headerOut ,currentR->v1[0],currentR->v1[1],CenterX,CenterY) ;
	      //  headerOut = Add_inside(headerOut ,currentR->v3[0],currentR->v3[1],CenterX,CenterY) ;
			spine2p = Add_mypoint (spine2p ,currentR->v1[0],currentR->v1[1],CenterX,CenterY) ;
			spine2p = Add_mypoint (spine2p ,currentR->v3[0],currentR->v3[1],CenterX,CenterY) ;
			Save2spinetri(currentR, header, CenterX,CenterY);

            while (tempmy!=NULL)
			{
                 //headerOut = Add_inside(headerOut ,tempmy->x,tempmy->y ,CenterX,CenterY) ;
				 spine2p = Add_mypoint (spine2p ,tempmy->x,tempmy->y ,CenterX,CenterY) ;
		         tempmy = tempmy->nextpoint ;
	        }
			current = PolyChain ;
			while (current!=NULL)
			{
				int x1 = current->x ;
			    int x2 = current->y ;
			    int y1 = current->x2 ;
			    int y2 = current->y2 ;

				if ( (currentR->v1[0]==x1&&currentR->v1[1]==x2&&currentR->v3[0]==y1&&currentR->v3[1]==y2 )     
					|| (currentR->v1[0]==y1&&currentR->v1[1]==y2&&currentR->v3[0]==x1&&currentR->v3[1]==x2 ))
				{
					current->index = 1 ;   // disable this merged edge in inner polygon.
					break;
				}
				current = current->nextpoint ;
			}
		}
		currentR = currentR->nextOut ;
	  }

	  current = PolyChain ;

	  output * currentS = spineOut;     // the spineOut must have value first! so recommand to do this aftrer connect sleeve tri.
	  while (currentS->nextOut!=NULL)
	  {
		  currentS = currentS->nextOut;
	  }

	  output * trans ;
	  while (current!=NULL)
	  {
		  if (current->index==0)
		  {
			  int x1 = current->x ;
			  int x2 = current->y ;
			  int y1 = current->x2 ;
			  int y2 = current->y2 ;
			  float C1 = ((float)x1+(float)y1)/2;
			  float C2 = ((float)x2+(float)y2)/2;

			  //headerOut = Add_inside(headerOut ,x1,x2,CenterX,CenterY) ;
			  //headerOut = Add_inside(headerOut ,y1,y2,CenterX,CenterY) ;
			  spine2p = Add_mypoint (spine2p ,x1,x2,CenterX,CenterY) ;
			  spine2p = Add_mypoint (spine2p ,y1,y2,CenterX,CenterY) ;

			  Save2spinetri_2p(x1, x2, C1 , C2 , CenterX,CenterY);
			  Save2spinetri_2p(y1, y2, C1 , C2 , CenterX,CenterY);      // Save the final spineTri.

			  trans = currentS;
			  currentS->nextOut = new output ;
			  currentS=currentS->nextOut ;
			  currentS->v2[0]=C1; 
			  currentS->v2[1]=C2; 
			  currentS->v4[0]=CenterX; 
			  currentS->v4[1]=CenterY;
			  currentS->nextOut = NULL ;
			  currentS->front=trans;
		  }
		  current = current->nextpoint ;
	  }

	  return spineOut;

 }


 // Connect all available sleeve tri. Output in spineOut.
 output * connectsleeve (output * headerOut, output * spineOut)
 {
	 output * currentO = headerOut ;
	 output * currentS;
	 output * trans ;
	 if (spineOut==NULL)
	 {
		 spineOut = new output;
		 currentS = spineOut;
	 }
	 else
	 {
		 currentS = spineOut;
	      while (currentS->nextOut!=NULL)
	     {
		  currentS = currentS->nextOut;
	     }
		  trans = currentS ;
		  currentS->nextOut = new output ;
		  currentS = currentS->nextOut ;
		  currentS->front = trans ;
	 }

	 while (currentO!=NULL)
	 {
		 if (currentO->TwoInnerline==true && currentO->X!=5&&currentO->X!=8)
		 {
			 int x1 = currentO->v1[0] ;
			 int x2 = currentO->v1[1] ;
			 int y1 = currentO->v2[0] ;
			 int y2 = currentO->v2[1] ;
			 float C1 = ((float)x1+(float)y1)/2;
			 float C2 = ((float)x2+(float)y2)/2;

			 int z1 = currentO->v3[0] ;
			 int z2 = currentO->v3[1] ;

			 float D1 = ((float)y1+(float)z1)/2;
			 float D2 = ((float)y2+(float)z2)/2;

			
			 currentS->v2[0]=C1; 
			 currentS->v2[1]=C2; 
			 currentS->v4[0]=D1; 
			 currentS->v4[1]=D2;
			 currentS->nextOut =  new output ;
			 trans = currentS;
			 currentS=currentS->nextOut ;
             currentS->front=trans;
			 currentS->nextOut = NULL ;
		 }
		 currentO = currentO->nextOut ;
	 }
	 trans = currentS->front ;
	 trans->nextOut = NULL ;
	 delete currentS ;
	 return spineOut ;
 }


 // Go through every sub_Polychain and calculate the spineOut, return it.
 struct output * connectPolychain(mypoint* Chain_Polychain, output * headerresult, output * headerOut, output * spineOut)
 {
	 if (Chain_Polychain==NULL)
		 return spineOut;
	 
	 mypoint * current = Chain_Polychain ;
	 mypoint * store ;
	 while (current!=NULL)
	 {
		 store = current->store ;
		 spineOut = calculate_Inpoly(store , headerresult , headerOut , spineOut);  //reflesh spineOut by checking the given PolyChain: store.
		 current = current->nextpoint ;
	 }
	 return spineOut ;
 }


 // Must put after <calculate_Inpoly> 
void dividesleeve( output * headerOut )
 {
	 output * current = headerOut ;
	 while ( current!= NULL )
	 {
		 if (current->X==0 && current->TwoInnerline==true )
		 {
			int a1 = current->v1[0];
			int a2 = current->v1[1];
			int b1 = current->v2[0];
			int b2 = current->v2[1];
			int c1 = current->v3[0];
			int c2 = current->v3[1];

			 float Cab1 = ((float)a1+(float)b1)/2;
	         float Cab2 = ((float)a2+(float)b2)/2;

	         float Cbc1 = ((float)c1+(float)b1)/2;
	         float Cbc2 = ((float)c2+(float)b2)/2;

			 spine2p = Add_mypoint (spine2p ,c1,c2,Cab1,Cab2) ;
			 spine2p = Add_mypoint (spine2p ,a1,a2,Cab1,Cab2) ;
			 spine2p = Add_mypoint (spine2p ,b1,b2,Cab1,Cab2) ;
			 spine2p = Add_mypoint (spine2p ,c1,c2,Cbc1,Cbc2) ;
			 spine2p = Add_mypoint (spine2p ,b1,b2,Cbc1,Cbc2) ;

			 Save2spinetri_1p(a1, a2, c1, c2, Cab1, Cab2);

			 Save2spinetri_2p(c1, c2, Cab1, Cab2, Cbc1, Cbc2);
			 Save2spinetri_2p(b1, b2, Cab1, Cab2, Cbc1, Cbc2);

		 }
		 current = current->nextOut ; 
	 }
 }


 // used to remove same elements in spine2p;
 void adjust_spine2p(void)
 {
	 mypoint * current = spine2p;
	 mypoint * currentIn;
	 mypoint * trans1;
	 mypoint * trans2;
	 int x; 
	 int y;
	 float xf;
	 float yf;

	 while ( current!=NULL)
	 {
		 x=current->x;
		 y=current->y;
		 xf=current->xf;
		 yf=current->yf;

		 currentIn = current->nextpoint ;
		 while (currentIn!=NULL)
		 {    // if two elements in spine2p are the same. delete the latter one.
			 if (currentIn->x==x&&currentIn->y==y&&currentIn->xf==xf&&currentIn->yf==yf)
			 {
				 trans1 = currentIn->nextpoint ;
				 trans2 = currentIn;
				 currentIn = currentIn->frontpoint ;
				 currentIn->nextpoint = trans1 ;
				 delete trans2;
				 if (trans1!=NULL)
					 trans1->frontpoint = currentIn ;
			 }
			 currentIn = currentIn->nextpoint ;
		 }
		 current = current->nextpoint ;
	 }
 }




 //============ End Spine ==============

 //============ Build up in 3D START ==============



 // Input spine Center point [a] Output its hight in Z-axis.
float Hight( float a1, float a2 )
 {
	 mypoint * current = spine2p;
	 float b1;
	 float b2;
	 float m1, m2;
	 float Length = 0, Length1=0;
	 int num = 0, num1 = 0;
	 while (current!=NULL)
	 {
		 if (a1==current->xf && a2==current->yf)
		 {
			 b1 = current->x ;
			 b2 = current->y ;
			 num++;
			 Length = Length + sqrt(pow((a1-b1),2)+pow((a2-b2),2)) ; 
		 }
		 current=current->nextpoint ;
	 }

	 output * currentO = spineOut;
	 while (currentO!=NULL)
	 {
		 if (abs(currentO->v2[0]-a1)<1.1 && abs(currentO->v2[1]-a2)<1.1)
		 {
			 m1 = currentO->v4[0];
			 m2 = currentO->v4[1];
			 break;
		 }
		 if (abs(currentO->v4[0]-a1)<1.1 && abs(currentO->v4[1]-a2)<1.1)
		 {
			 m1 = currentO->v2[0];
			 m2 = currentO->v2[1];
			 break;
		 }
		 currentO = currentO->nextOut;
	 }
	 current = spine2p;
	 while (current!=NULL)
	 {
		 if (m1==current->xf && m2==current->yf)
		 {
			 b1 = current->x ;
			 b2 = current->y ;
			 num1++;
			 Length1 = Length1 + sqrt(pow((m1-b1),2)+pow((m2-b2),2)) ; 
		 }
		 current=current->nextpoint ;
	 }
	 float count = (float)num / (float)num1;
	 if ( count>1.9 )
	 {
		 Length = Length1/num1;
		 return Length*0.7*0.97;
	 }
	 Length = Length/num ;
	 return Length*0.7 ;
 }


 // Input an array[5][3]  2D line and hight, output the 3D coordinate in this array[][].
void Hight_sub(float Hight , float p[5][3] , float c1 , float c2 , float e1 , float e2)
{

	p[0][0] = c1;
	p[0][1] = c2;
	p[0][2] = Hight;
	p[4][0] = e1;
	p[4][1] = e2;
	p[4][2] = 0;

	p[2][0] = (2*e1+c1)/3 ;
	p[2][1] = (2*e2+c2)/3 ;
	p[2][2] = sqrt(5.0)/3*Hight ;

	p[1][0] = (e1+2*c1)/3 ;
	p[1][1] = (e2+2*c2)/3 ;
	p[1][2] = sqrt(8.0)/3*Hight ;

	p[3][0] = (5*e1+c1)/6 ;
	p[3][1] = (5*e2+c2)/6 ;
	p[3][2] = sqrt(11.0)/6*Hight ;

}


// Input 3 space points, store in static tri3D.
void Save2trimesh(float a1, float a2, float a3, float b1 , float b2 , float b3, float c1 , float c2, float c3)
{
	if (tri3D!=NULL)
		{
			trimesh * trans = new trimesh ;
	        trimesh * tempOut ;
	        trans->a1 = a1;
			trans->a2 = a2;
			trans->a3 = a3;
			trans->b1 = b1;
			trans->b2 = b2;
			trans->b3 = b3;
			trans->c1 = c1;
			trans->c2 = c2;
			trans->c3 = c3;
			trans->X = 0 ;   

			trimesh * trans2 = new trimesh ;
			trans2->a1 = a1;
			trans2->a2 = a2;
			trans2->a3 = -a3;
			trans2->b1 = b1;
			trans2->b2 = b2;
			trans2->b3 = -b3;
			trans2->c1 = c1;
			trans2->c2 = c2;
			trans2->c3 = -c3;
			trans2->X = 0 ; 

			tempOut = tri3D->nextTri ;
			trans->nextTri = trans2 ;
			trans2->frontTri = trans ;
			trans2->nextTri = tempOut ;
			trans->frontTri = tri3D ;
	        tri3D->nextTri = trans ;
	        if (tempOut!=NULL)
				tempOut->frontTri = trans2 ;
		}

		if (tri3D==NULL)
		{
			trimesh * tempOut ;
			tri3D = new trimesh ;
			tri3D->a1 = a1;
			tri3D->a2 = a2;
			tri3D->a3 = a3;
			tri3D->b1 = b1;
			tri3D->b2 = b2;
			tri3D->b3 = b3;
			tri3D->c1 = c1;
			tri3D->c2 = c2;
			tri3D->c3 = c3;
			tri3D->X = 0 ; 
			tri3D->nextTri = new trimesh ;

			tempOut = tri3D->nextTri ;

			tempOut->a1 = a1;
			tempOut->a2 = a2;
			tempOut->a3 = a3;
			tempOut->b1 = b1;
			tempOut->b2 = b2;
			tempOut->b3 = b3;
			tempOut->c1 = c1;
			tempOut->c2 = c2;
			tempOut->c3 = c3;
			tempOut->X = 0 ; 

			tempOut->nextTri = NULL ;
		}
}


void makefull3Dtri(void)
{
	output * current = spineTri ;

	float a[5][3];
	float b[5][3];
	float e1;
	float e2;
	float c1;
	float c2;
	float m1;
	float m2;
	float H;
	float H2;

	while (current!=NULL)
	{
		if (current->X==1)
		{
			e1 = current->v1[0];
			e2 = current->v1[1];
			m1 = current->v3[0];
			m2 = current->v3[1];    // [m] and [e] are all edge points.
			c1 = current->v4[0];
			c2 = current->v4[1];

			H = Hight(c1, c2);
			Hight_sub(H, a, c1, c2, e1, e2);
			Hight_sub(H, b, c1, c2, m1, m2);   // now there are 2 arrays storing spacial points in order.   a[][] and b[][].

			Save2trimesh(a[4][0],a[4][1],a[4][2],b[4][0],b[4][1],b[4][2],b[3][0],b[3][1],b[3][2]);

			Save2trimesh(a[4][0],a[4][1],a[4][2],a[3][0],a[3][1],a[3][2],b[3][0],b[3][1],b[3][2]);

			Save2trimesh(a[3][0],a[3][1],a[3][2],b[3][0],b[3][1],b[3][2],b[2][0],b[2][1],b[2][2]);

			Save2trimesh(a[3][0],a[3][1],a[3][2],a[2][0],a[2][1],a[2][2],b[2][0],b[2][1],b[2][2]);

			Save2trimesh(a[2][0],a[2][1],a[2][2],b[2][0],b[2][1],b[2][2],b[1][0],b[1][1],b[1][2]);

			Save2trimesh(a[2][0],a[2][1],a[2][2],a[1][0],a[1][1],a[1][2],b[1][0],b[1][1],b[1][2]);

			Save2trimesh(a[1][0],a[1][1],a[1][2],b[1][0],b[1][1],b[1][2],a[0][0],a[0][1],a[0][2]);
		}

		if (current->X==2)
		{
			e1 = current->v1[0];
			e2 = current->v1[1];
			m1 = current->v2[0];
			m2 = current->v2[1];    // [m] and [c] are all Center points.
			c1 = current->v4[0];
			c2 = current->v4[1];

			H = Hight(c1, c2);
			H2 = Hight(m1, m2);

			if (H==0)
				H = H2;

			Hight_sub(H, a, c1, c2, e1, e2);
			Hight_sub(H2, b, m1, m2, e1, e2); 

			Save2trimesh(a[0][0],a[0][1],a[0][2],b[0][0],b[0][1],b[0][2],b[1][0],b[1][1],b[1][2]);

			Save2trimesh(a[0][0],a[0][1],a[0][2],a[1][0],a[1][1],a[1][2],b[1][0],b[1][1],b[1][2]);

			Save2trimesh(a[1][0],a[1][1],a[1][2],b[1][0],b[1][1],b[1][2],b[2][0],b[2][1],b[2][2]);

			Save2trimesh(a[1][0],a[1][1],a[1][2],a[2][0],a[2][1],a[2][2],b[2][0],b[2][1],b[2][2]);

			Save2trimesh(a[2][0],a[2][1],a[2][2],b[2][0],b[2][1],b[2][2],b[3][0],b[3][1],b[3][2]);

			Save2trimesh(a[2][0],a[2][1],a[2][2],a[3][0],a[3][1],a[3][2],b[3][0],b[3][1],b[3][2]);

			Save2trimesh(a[3][0],a[3][1],a[3][2],b[3][0],b[3][1],b[3][2],a[4][0],a[4][1],a[4][2]);
		}

		current=current->nextOut ; 
	}
}


void drawtri3D( void )
{
	glLineWidth(1.0f);
	trimesh * current = tri3D ;
	//int a = winWidth/2;
	//int b = winHeight/2;
	while (current!=NULL)
	{
		glBegin(GL_LINE_LOOP);
		glColor3f(1.0f, 0.0f, 0.0f);
		glVertex3f(current->a1,current->a2,current->a3);
		glVertex3f(current->b1,current->b2,current->b3);
		glVertex3f(current->c1,current->c2,current->c3);
		glEnd();
		current = current->nextTri ;
	}
}


void drawcuttingedge( void )
{
	glLineWidth(2.0f);
	trimesh * current = tri3D ;
	//int a = winWidth/2;
	//int b = winHeight/2;
	while (current!=NULL)
	{
		if (current->X==54)
		{
		    glBegin(GL_LINES);
		    glColor3f(1.0f, 1.0f, 1.0f);
		    glVertex3f(current->b1,current->b2,current->b3);
			glVertex3f(current->c1,current->c2,current->c3);
		    glEnd();
		}
		current = current->nextTri ;
	}
}

void drawtri3Dreal( void )
{
	//glLineWidth(1.0f);
	trimesh * current = tri3D ;
	while (current!=NULL)
	{
		glBegin(GL_TRIANGLES);
		glColor3f(0.7f, 0.7f, 0.7f);
		glVertex3f(current->a1,current->a2,current->a3);
		glVertex3f(current->b1,current->b2,current->b3);
		glVertex3f(current->c1,current->c2,current->c3);
		glEnd();
		current = current->nextTri ;
	}
}



// CUBE for TEST
void draw_box(float ox, float oy, float oz, float width, float height, float length)
{
    glLineWidth(1.0f);

    glBegin(GL_LINES);   
        glColor3f(1.0f, 0.0f, 0.0f);
		
		//1
        glVertex3f(ox, oy, oz);
        glVertex3f(ox+width, oy, oz);

		//2
        glVertex3f(ox, oy, oz);
        glVertex3f(ox, oy+height, oz);

		//3
        glVertex3f(ox, oy, oz);
        glVertex3f(ox, oy, oz+length);

		//4
        glVertex3f(ox+width, oy, oz);
        glVertex3f(ox+width, oy+height, oz);

		//5
        glVertex3f(ox+width, oy+height, oz);
        glVertex3f(ox, oy+height, oz);

		//6
        glVertex3f(ox, oy+height, oz+length);
        glVertex3f(ox, oy, oz+length);

		//7
        glVertex3f(ox, oy+height, oz+length);
        glVertex3f(ox, oy+height, oz);

		//8
        glVertex3f(ox+width, oy, oz);
        glVertex3f(ox+width, oy, oz+length);

		//9
        glVertex3f(ox, oy, oz+length);
        glVertex3f(ox+width, oy, oz+length);

		//10
        glVertex3f(ox+width, oy+height, oz);
        glVertex3f(ox+width, oy+height, oz+length);

		//11
        glVertex3f(ox+width, oy+height, oz+length);
        glVertex3f(ox+width, oy, oz+length);

		//12
        glVertex3f(ox, oy+height, oz+length);
        glVertex3f(ox+width, oy+height, oz+length);
    glEnd();
}



 //=============== Build up in 3D END ================

 //======== refine tri3D ========

// Input 3 space points ans store in tri3D without [-Z] value( unlike Save2trimesh)
void Save2trimesh3D(float a1, float a2, float a3, float b1 , float b2 , float b3, float c1 , float c2, float c3, int X)
{
	if (tri3D!=NULL)
		{
			trimesh * trans = new trimesh ;
	        trimesh * tempOut ;
	        trans->a1 = a1;
			trans->a2 = a2;
			trans->a3 = a3;
			trans->b1 = b1;
			trans->b2 = b2;
			trans->b3 = b3;
			trans->c1 = c1;
			trans->c2 = c2;
			trans->c3 = c3;
			trans->X = X ;   

			tempOut = tri3D->nextTri ;
			trans->frontTri = tri3D ;
	        tri3D->nextTri = trans ;
			trans->nextTri = tempOut ;
	        if (tempOut!=NULL)
				tempOut->frontTri = trans ;
		}

		if (tri3D==NULL)
		{
			trimesh * tempOut ;
			tri3D = new trimesh ;
			tri3D->a1 = a1;
			tri3D->a2 = a2;
			tri3D->a3 = a3;
			tri3D->b1 = b1;
			tri3D->b2 = b2;
			tri3D->b3 = b3;
			tri3D->c1 = c1;
			tri3D->c2 = c2;
			tri3D->c3 = c3;
			tri3D->X = X ; 

			tri3D->nextTri = NULL ;
		}
}


// refine and subdivide large triangles in tri3D.
// need the help of func: Save2trimesh3D.
void refinewholetri3D ( void )
{
	trimesh * current = tri3D;
	trimesh * trans1, * trans2 ;
	float a1,a2,a3,b1,b2,b3,c1,c2,c3,Dab,Dbc,Dca,D1,D2,D3,v1,v2,v3,
		 m1,m2,m3, p1,p2,p3,q1,q2,q3;

	current = tri3D;

	while (current!=NULL)
	{
		a1 = current->a1;
		a2 = current->a2;
		a3 = current->a3;
		b1 = current->b1;
		b2 = current->b2;
		b3 = current->b3;
		c1 = current->c1;
		c2 = current->c2;
		c3 = current->c3;

		Dab = sqrt((a1-b1)*(a1-b1)+(a2-b2)*(a2-b2)+(a3-b3)*(a3-b3));
		Dbc = sqrt((b1-c1)*(b1-c1)+(b2-c2)*(b2-c2)+(b3-c3)*(b3-c3));
		Dca = sqrt((c1-a1)*(c1-a1)+(c2-a2)*(c2-a2)+(c3-a3)*(c3-a3));

		D1 = Dab;
		D2 = Dbc;
		D3 = Dca;

		if (D1>D2)
			swap(D1,D2);
		if (D2>D3)
			swap(D2,D3);   // D3 max.
		if (D1*D2/2 > 127.0)
		{
			if (current == tri3D )    // delete current first.
			{
				tri3D = tri3D->nextTri;
				delete current;
				current = tri3D;
			}
			else
			{
				trans1 = current->frontTri;
				trans2 = current->nextTri;
				delete current;
				current = trans2;
				trans1->nextTri = current;
				if (current!=NULL)
					current->frontTri = trans1;
			}

			p1 = (a1+b1)/2;
			p2 = (a2+b2)/2;
			p3 = (a3+b3)/2;

			q1 = (b1+c1)/2;
			q2 = (b2+c2)/2;
			q3 = (b3+c3)/2;

			m1 = (c1+a1)/2;
			m2 = (c2+a2)/2;
			m3 = (c3+a3)/2;

			Save2trimesh3D(p1,p2,p3,q1,q2,q3,b1,b2,b3,0);
			Save2trimesh3D(p1,p2,p3,m1,m2,m3,q1,q2,q3,0);
			Save2trimesh3D(a1,a2,a3,m1,m2,m3,p1,p2,p3,0);
			Save2trimesh3D(c1,c2,c3,m1,m2,m3,q1,q2,q3,0);

			continue;
		}
		current = current->nextTri ;
	}
}


// refine and subdivide large triangles in tri3D.
// need the help of func: Save2trimesh3D.

//=========== End refine tri3D ============


 //==================== MakeSpine is to creat oubject in tri3D=================

 void drawspine(output * spineOut)
 {
	 output * currentS = spineOut ;
		while (currentS!=NULL)
		{
			glBegin(GL_LINES);
		    glColor3f(1.0,1.0,0.0);
			glVertex2f((GLfloat)currentS->v2[0],(GLfloat)currentS->v2[1]);
			glVertex2f((GLfloat)currentS->v4[0],(GLfloat)currentS->v4[1]);
			glEnd () ;
			glFlush () ;
			currentS = currentS->nextOut ;
		}
 }


  void drawspine2P(void)
 {
	 mypoint * currentS = spine2p ;
		while (currentS!=NULL)
		{
			glBegin(GL_LINES);
		    glColor3f(0.0,1.0,1.0);
			glVertex2f((GLfloat)currentS->x,(GLfloat)currentS->y);
			glVertex2f((GLfloat)currentS->xf,(GLfloat)currentS->yf);
			glEnd () ;
			glFlush () ;
			currentS = currentS->nextpoint ;
		}
 }


   void drawspineTri(void)
 {
	 output * currentT = spineTri ;
		while (currentT!=NULL)
		{
			glBegin(GL_LINE_LOOP);
		    glColor3f(0.0,1.0,1.0);
			if (currentT->X==1)
			{
			glVertex2f((GLfloat)currentT->v1[0],(GLfloat)currentT->v1[1]);
			glVertex2f((GLfloat)currentT->v3[0],(GLfloat)currentT->v3[1]);
			glVertex2f((GLfloat)currentT->v4[0],(GLfloat)currentT->v4[1]);
			}
			if (currentT->X==2)
			{
				glVertex2f((GLfloat)currentT->v1[0],(GLfloat)currentT->v1[1]);
			    glVertex2f((GLfloat)currentT->v2[0],(GLfloat)currentT->v2[1]);
			    glVertex2f((GLfloat)currentT->v4[0],(GLfloat)currentT->v4[1]);
			}
			glEnd () ;
			glFlush () ;
			currentT = currentT->nextOut ;
		}
 }


 struct output * MakeSpine( output * headerOut , output * junctriangle , output * headerresult )
 {
	 spineOut = refine_junctriangle(junctriangle, headerOut);

	 mypoint * Chain_Polychain = ChainofPolyChain(headerOut, junctriangle);

	 spineOut = connectsleeve(headerOut , spineOut);

	 spineOut = connectPolychain(Chain_Polychain , headerresult , headerOut , spineOut);

	 //drawspine(spineOut);

	 dividesleeve( headerOut );

	 adjust_spine2p();

	 glClear (GL_COLOR_BUFFER_BIT);

     //drawspine2P();

	 if (spineTri==NULL)
		 spineTri=NULL;
	
	 //glClear (GL_COLOR_BUFFER_BIT);

	 //drawspineTri();

	 makefull3Dtri();


	 if (tri3D!=NULL)
		 seeround=true ;

	 refinewholetri3D();

	 return spineOut ;
 }


//==================== NO USE, SIMPLY FOR TEST, END =================


 //============== from painting etc===============

 //============3D Painting============

 // Input 3 value (the position of allocated stroke point), this func will store them in pointOut IN ORDER.
// If first use, make sure pointOut = NULL.
struct point3D * Add2point3D(float x1, float x2,float x3, point3D * pointOut, int index, int check)
{
	if (pointOut!=NULL)
	{
		point3D * currentP = new point3D ;
		point3D * current = pointOut;
		currentP->a1=x1; 
		currentP->a2=x2; 
		currentP->a3=x3; 
		currentP->X=index;
		currentP->CheckifonlyOnestroke=check;
		while(current->next3D!=NULL)
		{
			current = current->next3D;
		}
		current->next3D = currentP;
		currentP->front3D = current;
		currentP->next3D = NULL;
        return pointOut ;
	}
	if (pointOut==NULL)
	{
		point3D * pointOut =new point3D;
		pointOut->a1 = x1;
		pointOut->a2 = x2;
		pointOut->a3 = x3;
		pointOut->X = index;
		pointOut->CheckifonlyOnestroke=check;
		pointOut->next3D = NULL;
		return pointOut ;
	}
}

 // Delete chain of point3D.
void Delete_point3D (point3D * headermy)
{
	point3D * current = headermy ;
	while ((current = headermy)!=NULL)
	{
		headermy = current->next3D ;
		delete current ;
	}
}


// Delete trimesh as whole.
void Delete_trimesh (trimesh * headermy)
{
	trimesh * current = headermy ;
	while ((current = headermy)!=NULL)
	{
		headermy = current->nextTri ;
		delete current ;
	}
}


 // Input the current spatial triangle, an orig point from drawing, a dir vector from camera to 0, and three parapmeters.
 // judge if the line intersects with given triangle.
 // The intersection is at [orig + dir * t].
 bool TriangleIntersect(trimesh * current , point3D orig, point3D dir, float *t, float *u, float *v)
 {
	 //E1 and E2.  E1=v1-v0;  E2=v2-v0.
	 point3D E1,E2,P ;
	 E1.a1 = current->b1 - current->a1 ;
	 E1.a2 = current->b2 - current->a2 ;
	 E1.a3 = current->b3 - current->a3 ;

	 E2.a1 = current->c1 - current->a1 ;
	 E2.a2 = current->c2 - current->a2 ;
	 E2.a3 = current->c3 - current->a3 ;

	 // P.  dir Cross E2.
	 P.a1 = dir.a2*E2.a3 - E2.a2*dir.a3;
	 P.a2 = dir.a3*E2.a1 - E2.a3*dir.a1;
	 P.a3 = dir.a1*E2.a2 - E2.a1*dir.a2;

	 // determinant.  E1 Dot P.
	 float det = E1.a1*P.a1 + E1.a2*P.a2 + E1.a3*P.a3 ;

	 // keep det > 0, modify T accordingly
	 point3D T;
	 if( det >0 )
    {
		T.a1 = orig.a1 - current->a1;
		T.a2 = orig.a2 - current->a2;
		T.a3 = orig.a3 - current->a3;
	 }
	 else
	 {
		T.a1 = current->a1 - orig.a1;
		T.a2 = current->a2 - orig.a2;
		T.a3 = current->a3 - orig.a3;
		det = -det ;
	 }

	  // If determinant is near zero, ray lies in plane of triangle.
	 if( det <0.0001f )
        return false;

	 // Calculate u and make sure u <= 1.  u = T Dot P.
	 *u = T.a1*P.a1+T.a2*P.a2+T.a3*P.a3 ;
	 if( *u <0.0f||*u > det )
        return false;

	 // Q. Q = T Cross E1.
	 point3D Q;
	 Q.a1 = T.a2*E1.a3 - E1.a2*T.a3;
	 Q.a2 = T.a3*E1.a1 - E1.a3*T.a1;
	 Q.a3 = T.a1*E1.a2 - E1.a1*T.a2;

	 // Calculate v and make sure u + v <= 1.  v = dir Dot Q.
	 *v = dir.a1*Q.a1+dir.a2*Q.a2+dir.a3*Q.a3 ;
	 if( *v <0.0f||*u +*v > det )
        return false;

	 // Calculate t, scale parameters, ray intersects triangle. t = E2 Dot Q.
	 *t = E2.a1*Q.a1+E2.a2*Q.a2+E2.a3*Q.a3 ;

	 float fInvDet =1.0f/ det;
     *t *= fInvDet;
     *u *= fInvDet;
     *v *= fInvDet;

	 return true;
 }


 // Put data in mypoint to point3D by adding Z-axis.
 struct point3D* tranStrokes( mypoint * header)
 {
	if (header==NULL)
		 return NULL;
	 mypoint * current = header ;
	 point3D * strokes3D = new point3D ;
	 point3D * trans,* current3D=strokes3D;
	 int IndexS = current->index;
	 while (current!=NULL)  
	 {
		 if (IndexS!=current->index)
		 {
			 strokecount++;
			 IndexS = current->index;
		 }

		 current3D->a1 = current->x ;
		 current3D->a2 = current->y ;
		 current3D->a3 = 1000.0;
		 current3D->X = strokecount ;

		 current = current->nextpoint ;
		 trans = current3D;
		 current3D->next3D = new point3D ;
		 current3D = current3D->next3D ;
		 current3D->front3D = trans ;
		 trans->next3D = current3D ;
	 }

	 trans = current3D->front3D ;
	 delete current3D;
	 trans->next3D = NULL ;
	 strokecount++;
	 return strokes3D;

 }


 // Return the Camera position after rotation.
struct point3D* rotateCamera (float Xroat, float Yroat)
{
	Xroat = Xroat*3.1416/180;
	Yroat = Yroat*3.1416/180;
	point3D * Camroat = new point3D ;
	float a1 = 0.0;
	float a2 = -1000*sin(-Xroat);
	float a3 = 1000*cos(-Xroat);

	Camroat->a1 = a3*sin(-Yroat)+a1*cos(-Yroat);
    Camroat->a2 = a2;
	Camroat->a3 = a3*cos(-Yroat)-a1*sin(-Yroat);

	return Camroat;
}


// Return all Strokes' position after rotation.
struct point3D* rotateStrokes (float Xroat, float Yroat, point3D * strokes3D)
{
	Xroat = Xroat*3.1416/180;
	Yroat = Yroat*3.1416/180;
	point3D * current = strokes3D ;
	float a1,a2,a3;
	while (current!=NULL)
	{
		a1 = current->a1;
		a2 = current->a2 * cos(-Xroat) - current->a3 * sin(-Xroat);
		a3 = current->a2 * sin(-Xroat) + current->a3 * cos(-Xroat);

		current->a1 = a3*sin(-Yroat)+a1*cos(-Yroat);
		current->a2 = a2;
		current->a3 = a3*cos(-Yroat)-a1*sin(-Yroat);

		current = current->next3D ;
	}
	return strokes3D;
}


//This function is preparing for show draw lines on screen.
//tempstrokes
//use a point of header, output its rotated position and store it in tempstrokes.
struct point3D* rotateOnepoint (float Xroat, float Yroat, mypoint * currentP, float sign)
{
	if (sign>0.0)
	{
		Xroat = -Xroat;
		Yroat = -Yroat;
	}
	Xroat = Xroat*3.1416/180;
	Yroat = Yroat*3.1416/180;
	float a1,a2,a3,x1,x2,x3;
	
	a1 = (float)currentP->x;
	a2 = (float)currentP->y * cos(-Xroat) - 1000.0 * sin(-Xroat);
	a3 = (float)currentP->y * sin(-Xroat) + 1000.0 * cos(-Xroat);

	x1 = a3*sin(-Yroat)+a1*cos(-Yroat);
	x2 = a2;
	x3 = a3*cos(-Yroat)-a1*sin(-Yroat);

	tempstrokes = Add2point3D(x1, x2, x3, tempstrokes, 0, 0);

	return tempstrokes;
}

// Get the spatial stroke points by two rotation angle, spatial triangles and some raw strokes data.
// pointOut stores the set of spatial vertexes of the given strokes.
struct point3D* locate3Dstrokes(trimesh * tri3D, float Xroat, float Yroat, mypoint * header, point3D * pointOut)
{
	//point3D * pointOut = NULL ;
	point3D * Camroat = rotateCamera (Xroat, Yroat);
	point3D * strokes3D = tranStrokes (header);
	strokes3D = rotateStrokes ( Xroat, Yroat, strokes3D);
	float t, u, v, a, b, c;
	point3D dir;
	dir.a1 = 0-Camroat->a1 ;
	dir.a2 = 0-Camroat->a2 ;
	dir.a3 = 0-Camroat->a3 ;
	point3D * current=strokes3D;
	point3D orig;
	trimesh * currentri=tri3D;
	bool type = false;
	float d = 99999.9;

	while (current!=NULL)
	{
		orig.a1=current->a1;
		orig.a2=current->a2;
		orig.a3=current->a3;
		currentri = tri3D;
		while (currentri!=NULL)
		{
			type = TriangleIntersect(currentri, orig, dir, &t, &u, &v);
			if (type==true)
			{
				if (d>t)
					d=t;
			}
			currentri = currentri->nextTri ;
			type = false ;
		}
		if ( 99999.9-d<1)
		{
			current = current->next3D ;
			continue;
		}
		d = d-0.0001;
		a = orig.a1+dir.a1*d;    // use orig + dir * d to calculate stroke point's location.
		b = orig.a2+dir.a2*d;
		c = orig.a3+dir.a3*d;
		pointOut = Add2point3D( a, b, c, pointOut, current->X, current->CheckifonlyOnestroke );
		d=99999.9;
		current = current->next3D ;
	}

	return pointOut ;
}

//============End 3D Painting=============



//============ For test Painting =============
// draw strokes(using pointOut).
void drawstrokes(void)
{
	glLineWidth(2.0f);
	point3D * current = pointOut ;
	int a;
	if (pointOut==NULL)
		return ;
	
	int IndexS = current->X ;

	while(1)
	{
		glBegin(GL_LINE_STRIP);
		glColor3f (0.0, 0.0, 1.0);
		while ( current!=NULL && IndexS == current->X )
		{
			if (current->CheckifonlyOnestroke!=97)
			{
			    glVertex3f(current->a1,current->a2,current->a3);
			}
			current = current->next3D ;
		}
		glEnd();
	    glFlush();

		if (current==NULL)
			break;
		IndexS = current->X ;
	}
	
}

void drawstrokes_basering(void)
{
	glLineWidth(3.0f);
	glColor3f(1.0f, 1.0f, 1.0f);
	point3D * current = pointOut ;
	if (pointOut==NULL)
		return ;
	
	int IndexS ;

	while(1)
	{
		if (current->CheckifonlyOnestroke==97)
		{
			IndexS = current->X;
			while ( current->next3D!=NULL && current->X == IndexS )
		    {
			    glBegin(GL_LINES);
			    glVertex3f(current->a1,current->a2,current->a3);
		        current = current->next3D ;
			    glVertex3f(current->a1,current->a2,current->a3);
			    glEnd();
	            glFlush();
		    }
		}

		if (current->next3D==NULL)
			break;
		current = current->next3D ;
	}
	
}


// draw hand-drawing strokes.
void drawLineSegment (void)
{
	point3D * currentP = tempstrokes ;
	glColor3f (0.0, 1.0, 0.0);
	while (currentP!=NULL)
	{
	   glBegin (GL_POINTS);
	   glVertex3f (currentP->a1,currentP->a2,currentP->a3);  // previous point in movement
       glEnd();
	   currentP = currentP->next3D;
	}
}


//==============   Extrusion    ============



// Get the gravity center of any polygon which stores in the header and index==1.
// output->a1 = CenterX ;
// output->a2 = CenterY .
struct point3D* getGravitycenter( mypoint * header)
{
	mypoint * current = header;
	float Distance=0 ;
	float Discheck=0 ;
	int num=0 ;
	float a1,a2,b1,b2,c1,c2;
	point3D *P=NULL, *Q=NULL;
	while(current!=NULL && current->index==1)
	{
		if (current->nextpoint==NULL || current->nextpoint->index!=1)
		{
			Distance += sqrt(pow((float)current->x-(float)header->x,2)+pow((float)current->y-(float)header->y,2)) ;
			break;
		}
		a1 = current->x;
		a2 = current->y;
		current = current->nextpoint;
		b1 = current->x;
		b2 = current->y;

		Distance += sqrt(pow(a1-b1,2)+pow(a2-b2,2));

	}
	Distance = Distance/6.0;

	current = header;
	while (num!=5)
	{
		a1 = current->x;
		a2 = current->y;
		current = current->nextpoint;
		b1 = current->x;
		b2 = current->y;

		Discheck += sqrt(pow(a1-b1,2)+pow(a2-b2,2));
		if ( Discheck > Distance )
		{
			Discheck = 0.0;
			P = Add2point3D( b1, b2, 0.0, P, 0, 0);
			num++;
		}
	}

	point3D *currentP = P;
	c1 = header->x ;
	c2 = header->y ;
	float lc1, lc2, Center1, Center2, Ca1, Ca2, Cb1, Cb2, Cc1, Cc2, Cd1, Cd2;
	while ( currentP->next3D!=NULL )
	{
		a1 = currentP->a1;
		a2 = currentP->a2;
		currentP = currentP->next3D ;
		b1 = currentP->a1;
		b2 = currentP->a2;

		lc1 = (a1+b1)/2;
		lc2 = (a2+b2)/2;

		Center1 = (2*lc1+c1)/3;
		Center2 = (2*lc2+c2)/3;

		Q = Add2point3D( Center1, Center2, 0.0, Q, 0, 0);
	}
	a1 = Q->a1;
	a2 = Q->a2;
	b1 = Q->next3D->a1;
	b2 = Q->next3D->a2;
	c1 = Q->next3D->next3D->a1;
	c2 = Q->next3D->next3D->a2;
	lc1 = Q->next3D->next3D->next3D->a1;
	lc2 = Q->next3D->next3D->next3D->a2;

	Ca1 = (a1+b1)/2;
	Ca2 = (a2+b2)/2;
	Cb1 = (c1+b1)/2;
	Cb2 = (c2+b2)/2;
	Cc1 = (c1+lc1)/2;
	Cc2 = (c2+lc2)/2;
	Cd1 = (lc1+a1)/2;
	Cd2 = (lc2+a2)/2;

	float M = (Ca2-Cc2)/(Ca1-Cc1);
	float N = (Cd2-Cb2)/(Cd1-Cb1);

	float CenterX = (M*Cc1-N*Cb1+Cb2-Cc2)/(M-N);
	float CenterY = ((1/M)*Cc2-(1/N)*Cb2+Cb1-Cc1)/(1/M-1/N);  // find the center of Polygon.

	Delete_point3D(P);
	Delete_point3D(Q);
	point3D * output = new point3D;
	output->a1 = CenterX ;
	output->a2 = CenterY ;

	return output;
}

// Add frontpoint of header.
void addfrontpoint(void)
{
	if (header==NULL || header->nextpoint==NULL)
		return;
	mypoint * current = header;
	mypoint * trans;
	while ( current->nextpoint!=NULL)
	{
		trans = current;
		current = current->nextpoint;
		current->frontpoint = trans;
	}
}

// Output the gravity center, and two intersected points on base ring.
// a & b -> intersections   c-> G center.
// a-> up intersection.
// Vector_Center->a == second stroke's center point of start and end.  
// Vector_Center->b == second stroke's feature vector.  (b1-0, b2-0)
struct trimesh * baseringLine( mypoint * header, trimesh * Vector_Center )
{
	mypoint * current = header ;

	while (current!=NULL)
	{
		if (current->index==2)
			break;
		current = current->nextpoint ;
	}
	float x1 = current->x;    // first point of stroke2.
	float x2 = current->y;

	int count = 0;
	float temp1;
	float temp2;

	while (current!=NULL)
	{
		if ( count == 10 )
		{
			temp1 = current->x;     // a random point in stroke2.
			temp2 = current->y;
		}
		if (current->nextpoint==NULL || current->nextpoint->index==3)
			break;
		current = current->nextpoint ;
		count++;
	}

	float y1 = current->x;     // end point of stroke2.
	float y2 = current->y;

	point3D *target = getGravitycenter(header);
	float t1 = target->a1;    // the gravity center of stroke1.
	float t2 = target->a2;

	float lc1 = (x1+y1)/2;    // mid-point of first and end points of stroke2.
	float lc2 = (x2+y2)/2;

	float v1, v2, v3, v4;

	v3 = temp1 - lc1;
	v4 = temp2 - lc2;

	v1 = y2 - x2;
	v2 = x1 - y1;

	if ( acos((v3*v1+v4*v2)/sqrt(v3*v3+v4*v4)/sqrt(v1*v1+v2*v2))>3.1416/2)
	{
		v1 = x2 - y2;
		v2 = y1 - x1;
	}

	Vector_Center->a1 = lc1;
	Vector_Center->a2 = lc2;
	Vector_Center->b1 = v1;
	Vector_Center->b2 = v2;
	Vector_Center->c1 = x1;
	Vector_Center->c2 = x2;
	Vector_Center->d1 = y1;
	Vector_Center->d2 = y2;

	t1 = t1-lc1;   // now t is a vector from center of line xy to target.
	t2 = t2-lc2;

	x1 = x1+t1;    // use the vector t to allocate new location of xy.
	x2 = x2+t2;
	y1 = y1+t1;
	y2 = y2+t2;

	float xy1 = x1-target->a1;     // a vector form G center to x(first point).
	float xy2 = x2-target->a2; 
	float yg1 = y1-target->a1;     // a vector form G to y(end point).
	float yg2 = y2-target->a2;

	x1 = x1+xy1;
	x2 = x2+xy2;
	y1 = y1+yg1;
	y2 = y2+yg2;        // enlarge xy one times.

	mypoint p1, p2, p3, p4, pg;
	p1.x = x1;
	p1.y = x2;
	p2.x = y1;
	p2.y = y2;
	pg.x = target->a1;
	pg.y = target->a2;

	current = header ;
	mypoint Int;
	int type;
	trimesh * A = new trimesh;

	while (current!=NULL && current->index==1)
	{
		p3.x = current->x;
		p3.y = current->y;
		current = current->nextpoint;
		if (current!=NULL && current->index==1)
		{
			p4.x = current->x;
		    p4.y = current->y;
		}
		else 
		{
			p4.x = header->x;
			p4.y = header->y;
		}

		type = Intersection( p1, pg, p3, p4, Int);

		if (type == 2 || type == 3 || type == 1)                  // If intersected between first point and G center.
		{
			A->a1 = p4.x;
			A->a2 = p4.y;
		}

		type = Intersection( p2, pg, p3, p4, Int);

		if (type == 2 || type == 3 || type == 1)                  // If intersected between end point and G center.
		{
			A->b1 = p4.x;
			A->b2 = p4.y;
		}
	}

	A->c1 = target->a1;
	A->c2 = target->a2;
	return A;
}

// for a given line segment, find the possible intersection with a stroke.
// the intersection point is store in the struct A.
struct trimesh * searchintersection (trimesh * A, mypoint * header, float x1, float x2, float y1, float y2)
{
	mypoint p1, p2, p3, p4;
	p1.x = x1;
	p1.y = x2;
	p2.x = y1;
	p2.y = y2;

	mypoint * current = header ;
	mypoint Int;

	int type;
	while (current!=NULL && current->index==1)
	{
		p3.x = current->x;
		p3.y = current->y;
		current = current->nextpoint;
		if (current!=NULL && current->index==1)
		{
			p4.x = current->x;
		    p4.y = current->y;
		}
		else 
		{
			p4.x = header->x;
			p4.y = header->y;
		}

		type = Intersection( p1, p2, p3, p4, Int);

		if (type == 3 || type == 2 || type == 1)                  // If intersected.
		{		
			A->a1 = p4.x;
			A->a2 = p4.y;
			break;
		}
	}
	return A;
}



// given center line of a ring stroke, divide the line by 6, find the corresponding Hight of each segment.
// output is for uper points-> u[5] and for done points-> d[5] in order.
// float x is first intersection,  y is the second.
void findhight (float x1, float x2, float y1, float y2, mypoint * header,float u[5],float d[5],float uu[5][2], float dd[5][2])
{
	float v1 = x1-y1;
	float v2 = x2-y2;

	float orth1 = -10*v2;
	float orth2 = 10*v1;

	float inter11 = y1+5*v1/6;
	float inter12 = y2+5*v2/6;
	float inter21 = y1+4*v1/6;
	float inter22 = y2+4*v2/6;
	float inter31 = y1+3*v1/6;
	float inter32 = y2+3*v2/6;
	float inter41 = y1+2*v1/6;
	float inter42 = y2+2*v2/6;
	float inter51 = y1+v1/6;
	float inter52 = y2+v2/6;

	float outer11 = inter11 + orth1;
	float outer12 = inter12 + orth2;
	float outer21 = inter21 + orth1;
	float outer22 = inter22 + orth2;
	float outer31 = inter31 + orth1;
	float outer32 = inter32 + orth2;
	float outer41 = inter41 + orth1;
	float outer42 = inter42 + orth2;
	float outer51 = inter51 + orth1;
	float outer52 = inter52 + orth2;

	trimesh * A = new trimesh;
	float D=0;
	A = searchintersection(A, header, inter11, inter12, outer11, outer12);
	D = sqrt(pow(inter11-A->a1,2)+pow(inter12-A->a2,2));
	u[0] = D;
	uu[0][0] = A->a1;
	uu[0][1] = A->a2;

	A = searchintersection(A, header, inter21, inter22, outer21, outer22);
	D = sqrt(pow(inter21-A->a1,2)+pow(inter22-A->a2,2));
	u[1] = D;
	uu[1][0] = A->a1;
	uu[1][1] = A->a2;

	A = searchintersection(A, header, inter31, inter32, outer31, outer32);
	D = sqrt(pow(inter31-A->a1,2)+pow(inter32-A->a2,2));
	u[2] = D;
	uu[2][0] = A->a1;
	uu[2][1] = A->a2;

	A = searchintersection(A, header, inter41, inter42, outer41, outer42);
	D = sqrt(pow(inter41-A->a1,2)+pow(inter42-A->a2,2));
	u[3] = D;
	uu[3][0] = A->a1;
	uu[3][1] = A->a2;

	A = searchintersection(A, header, inter51, inter52, outer51, outer52);
	D = sqrt(pow(inter51-A->a1,2)+pow(inter52-A->a2,2));
	u[4] = D;
	uu[4][0] = A->a1;
	uu[4][1] = A->a2;

	orth1 = 10*v2;
	orth2 = -10*v1;

	outer11 = inter11 + orth1;
	outer12 = inter12 + orth2;
	outer21 = inter21 + orth1;
	outer22 = inter22 + orth2;
	outer31 = inter31 + orth1;
	outer32 = inter32 + orth2;
	outer41 = inter41 + orth1;
	outer42 = inter42 + orth2;
	outer51 = inter51 + orth1;
	outer52 = inter52 + orth2;

	A = searchintersection(A, header, inter11, inter12, outer11, outer12);
	D = sqrt(pow(inter11-A->a1,2)+pow(inter12-A->a2,2));
	d[0] = D;
	dd[0][0] = A->a1;
	dd[0][1] = A->a2;

	A = searchintersection(A, header, inter21, inter22, outer21, outer22);
	D = sqrt(pow(inter21-A->a1,2)+pow(inter22-A->a2,2));
	d[1] = D;
	dd[1][0] = A->a1;
	dd[1][1] = A->a2;

	A = searchintersection(A, header, inter31, inter32, outer31, outer32);
	D = sqrt(pow(inter31-A->a1,2)+pow(inter32-A->a2,2));
	d[2] = D;
	dd[2][0] = A->a1;
	dd[2][1] = A->a2;

	A = searchintersection(A, header, inter41, inter42, outer41, outer42);
	D = sqrt(pow(inter41-A->a1,2)+pow(inter42-A->a2,2));
	d[3] = D;
	dd[3][0] = A->a1;
	dd[3][1] = A->a2;

	A = searchintersection(A, header, inter51, inter52, outer51, outer52);
	D = sqrt(pow(inter51-A->a1,2)+pow(inter52-A->a2,2));
	d[4] = D;
	dd[4][0] = A->a1;
	dd[4][1] = A->a2;

}




// store mypoint in order.
// use to store base ring.
struct mypoint * Add2mypoint(float x1, float x2, mypoint * pointOut)
{
	if (pointOut!=NULL)
	{
		mypoint * currentP = new mypoint ;
		mypoint * current = pointOut;
		currentP->x=x1; 
		currentP->y=x2; 

		while(current->nextpoint!=NULL)
		{
			current = current->nextpoint;
		}
		current->nextpoint = currentP;
		currentP->frontpoint = current;
		currentP->nextpoint = NULL;
        return pointOut ;
	}
	if (pointOut==NULL)
	{
		mypoint * pointOut =new mypoint;
		pointOut->x = x1;
		pointOut->y = x2;
		pointOut->nextpoint = NULL;
		return pointOut ;
	}
}


// store 2D refined base ring into mypoint * baseringP.
void findbaseringP (float x1, float x2, float y1, float y2, float u[5],float d[5])
{
	float v1 = x1-y1;
	float v2 = x2-y2;

	float orth1 = -10*v2;
	float orth2 = 10*v1;
	float sum = sqrt(orth1*orth1+orth2*orth2);
	orth1 = orth1/sum;
	orth2 = orth2/sum;   // normal vector.
	

	float inter11 = y1+5*v1/6;
	float inter12 = y2+5*v2/6;
	float inter21 = y1+4*v1/6;
	float inter22 = y2+4*v2/6;
	float inter31 = y1+3*v1/6;
	float inter32 = y2+3*v2/6;
	float inter41 = y1+2*v1/6;
	float inter42 = y2+2*v2/6;
	float inter51 = y1+v1/6;
	float inter52 = y2+v2/6;

	float outer11 = inter11 + orth1*u[0];
	float outer12 = inter12 + orth2*u[0];
	float outer21 = inter21 + orth1*u[1];
	float outer22 = inter22 + orth2*u[1];
	float outer31 = inter31 + orth1*u[2];
	float outer32 = inter32 + orth2*u[2];
	float outer41 = inter41 + orth1*u[3];
	float outer42 = inter42 + orth2*u[3];
	float outer51 = inter51 + orth1*u[4];
	float outer52 = inter52 + orth2*u[4];

	baseringP = Add2mypoint(x1,x2,baseringP);
	baseringP = Add2mypoint(outer11,outer12,baseringP);
	baseringP = Add2mypoint(outer21,outer22,baseringP);
	baseringP = Add2mypoint(outer31,outer32,baseringP);
	baseringP = Add2mypoint(outer41,outer42,baseringP);
	baseringP = Add2mypoint(outer51,outer52,baseringP);



	orth1 = -orth1;
	orth2 = -orth2;

	outer11 = inter11 + orth1*d[0];
	outer12 = inter12 + orth2*d[0];
	outer21 = inter21 + orth1*d[1];
	outer22 = inter22 + orth2*d[1];
	outer31 = inter31 + orth1*d[2];
	outer32 = inter32 + orth2*d[2];
	outer41 = inter41 + orth1*d[3];
	outer42 = inter42 + orth2*d[3];
	outer51 = inter51 + orth1*d[4];
	outer52 = inter52 + orth2*d[4];

	baseringP = Add2mypoint(y1,y2,baseringP);
	baseringP = Add2mypoint(outer51,outer52,baseringP);
	baseringP = Add2mypoint(outer41,outer42,baseringP);
	baseringP = Add2mypoint(outer31,outer32,baseringP);
	baseringP = Add2mypoint(outer21,outer22,baseringP);
	baseringP = Add2mypoint(outer11,outer12,baseringP);

}




//Add [a] and [f] into mypoint * spine2p IN ORDER.
struct mypoint * Add_mypoint_inorder(mypoint *spine2p ,float a1, float a2, float f1, float f2, int num)
{
	if (spine2p!=NULL)
	{
	mypoint * trans = new mypoint ;
	mypoint * current = spine2p;
	trans->xf = f1;
	trans->yf = f2;
	trans->x = a1;
	trans->y = a2;
	trans->num = num;

	while(current->nextpoint!=NULL)
		{
			current = current->nextpoint;
		}
	current->nextpoint = trans;
	trans->frontpoint = current;
	trans->nextpoint = NULL;
	return spine2p ;
	}
	if (spine2p==NULL)
	{
		spine2p = new mypoint ;
		spine2p->xf = f1;
		spine2p->yf = f2;
		spine2p->x = a1;
		spine2p->y = a2;
		spine2p->num = num;
		spine2p->nextpoint = NULL ;
		return spine2p ;
	}
}


// refine the current to end of index = indexN stroke in header so that the length of each line segment is the same.
// may be looks nice.
// DATA structure: stroke2->a1,a2 & num.
struct mypoint * refine2stroke( mypoint * current, int indexN, float length )
{
	mypoint * stroke2 = NULL;

	mypoint * header2 = current ;

	float Distance=0 ;
	float a1,a2,b1,b2;

	while(current!=NULL && current->index==indexN)
	{
		a1 = current->x;
		a2 = current->y;
		if (current->nextpoint==NULL || current->nextpoint->index!=indexN)
			break;
		current = current->nextpoint;
		b1 = current->x;
		b2 = current->y;

		Distance += sqrt(pow(a1-b1,2)+pow(a2-b2,2));

	}
	float tempb1 = b1;
	float tempb2 = b2;   // this is the end point.

	int segNum = floor(Distance/length);
	Distance = Distance / segNum;
	float distancecheck=0;
	int Num = 0;
	int count = 1;

	current = header2;
	stroke2 = Add_mypoint_inorder(stroke2, header2->x, header2->y, 0.0, 0.0, count);  // store the first point.


	while(current!=NULL && current->index==indexN)
	{
		a1 = current->x;
		a2 = current->y;
		if (current->nextpoint==NULL || current->nextpoint->index!=indexN)
			break;
		current = current->nextpoint;
		b1 = current->x;
		b2 = current->y;

		distancecheck += sqrt(pow(a1-b1,2)+pow(a2-b2,2));
		if (distancecheck>Distance)
		{
			distancecheck=0;
			count++;
			Num++;
			stroke2 = Add_mypoint_inorder(stroke2, b1, b2, 0.0, 0.0, count);
			
			if (Num==segNum-1)
				break;
		}
	}
	count++;
	stroke2 = Add_mypoint_inorder(stroke2, tempb1, tempb2, 0.0, 0.0, count);

	return stroke2;
}


// vector from  a to b   and  a  to  c.
// return abs(Angle-pi/2)
float seeangle( float a1, float a2, float b1, float b2, float c1, float c2)
{
	float v1 = b1 - a1;
	float v2 = b2 - a2;
	float vx1 = c1 - a1;
	float vx2 = c2 - a2;

	float Angle = acos((vx1*v1 + vx2*v2)/sqrt(vx1*vx1+vx2*vx2)/sqrt(v1*v1+v2*v2));
	return abs(Angle-3.1416/2);
}


// find the inline of stroke2 than store them in Inline IN ORDER.
// Inline->x,y always store the Start points, and Inline->xf,yf always stroe the End points.
// if num==0, this mypoint stores the converge point.
struct mypoint * findstr2inline(mypoint * stroke2)
{
	mypoint * Inline = NULL;

	mypoint * currentS = stroke2;
	while ( currentS->nextpoint!=NULL)
	{
		currentS = currentS->nextpoint;
	}
	mypoint * currentE = currentS;
	currentS = stroke2;

	mypoint * transS1, * transS2, * transE1, * transE2;
	float sum1, sum2, sum3, type1, type2, type3;
	int num = 1;
	Inline = Add_mypoint_inorder(Inline, currentS->x, currentS->y, currentE->x, currentE->y, num);   // the first Inline.

	while ( currentS->nextpoint != currentE->frontpoint )
	{

	// TYPE 1.
	transE1 = currentE->frontpoint;
	transE2 = transE1->frontpoint;
	transS1 = currentS->nextpoint;
	sum1 = seeangle(currentS->x, currentS->y, transS1->x, transS1->y, transE1->x, transE1->y);
	sum1 += seeangle(transE1->x, transE1->y, transE2->x, transE2->y, currentS->x, currentS->y);

	//TYPE 2.
	transS1 = currentS->nextpoint;
	transS2 = transS1->nextpoint;
	transE1 = currentE->frontpoint;
	sum2 = seeangle(transS1->x, transS1->y, transS2->x, transS2->y, currentE->x, currentE->y);
	sum2 += seeangle(currentE->x, currentE->y, transE1->x, transE1->y, transS1->x, transS1->y);

	//TYPE 3.
	transE1 = currentE->frontpoint;
	transE2 = transE1->frontpoint;
	transS1 = currentS->nextpoint;
	transS2 = transS1->nextpoint;
	sum3 = seeangle(transS1->x, transS1->y, transS2->x, transS2->y, transE1->x, transE1->y);
	sum3 += seeangle(transE1->x, transE1->y, transE2->x, transE2->y, transS1->x, transS1->y);

	type1 = sum1;
	type2 = sum2;
	type3 = sum3;

	if ( sum1< sum2 )
		swap(sum1, sum2);
	if ( sum2< sum3 )
		swap(sum2, sum3);
	
	if ( sum3==type1)
	{
		num++;
		Inline = Add_mypoint_inorder( Inline, currentS->x, currentS->y, currentE->frontpoint->x, currentE->frontpoint->y, num);
		currentE = currentE->frontpoint ;
	}

	if ( sum3==type2)
	{
		num++;
		Inline = Add_mypoint_inorder( Inline, currentS->nextpoint->x, currentS->nextpoint->y, currentE->x, currentE->y, num);
		currentS = currentS->nextpoint ;
	}

	if ( sum3==type3)
	{
		num++;
		Inline = Add_mypoint_inorder( Inline,currentS->nextpoint->x, currentS->nextpoint->y,currentE->frontpoint->x, currentE->frontpoint->y, num);
		currentE = currentE->frontpoint ;
		currentS = currentS->nextpoint ;
	}

	}

	// also store the converge point with num = 0 
	Inline = Add_mypoint_inorder( Inline,currentS->nextpoint->x, currentS->nextpoint->y,currentE->frontpoint->x, currentE->frontpoint->y, 0);

	return Inline;
}


// Input 3 space points, store in static triEx.
void Save2trimeshEX(float a1, float a2, float a3, float b1 , float b2 , float b3, float c1 , float c2, float c3, int X)
{
	if (triEx!=NULL)
		{
			trimesh * trans = new trimesh ;
	        trimesh * tempOut ;
	        trans->a1 = a1;
			trans->a2 = a2;
			trans->a3 = a3;
			trans->b1 = b1;
			trans->b2 = b2;
			trans->b3 = b3;
			trans->c1 = c1;
			trans->c2 = c2;
			trans->c3 = c3;
			trans->X = X ;   

			tempOut = triEx->nextTri ;
			trans->frontTri = triEx ;
	        triEx->nextTri = trans ;
			trans->nextTri = tempOut ;
	        if (tempOut!=NULL)
				tempOut->frontTri = trans ;
		}

		if (triEx==NULL)
		{
			trimesh * tempOut ;
			triEx = new trimesh ;
			triEx->a1 = a1;
			triEx->a2 = a2;
			triEx->a3 = a3;
			triEx->b1 = b1;
			triEx->b2 = b2;
			triEx->b3 = b3;
			triEx->c1 = c1;
			triEx->c2 = c2;
			triEx->c3 = c3;
			triEx->X = X ; 

			triEx->nextTri = NULL ;
		}
}




// given 2 adjacent inner lines of stroke2 with their end points(la-lb and lc-ld).
// given 'standard' length of base ring,
// given which side should be consider, istrightside==1, than input u[5], elseif==2 input d[5] elseif == 0 need to input both u and d since is the last point.
// output triangles stored in the triEX.
void Inline2triEX(float la1, float la2, float lb1, float lb2, float lc1, float lc2, float ld1, float ld2, float x[5], float standard, int isrightside)
{
	float u[5];
	u[0]=x[0];
	u[1]=x[1];
	u[2]=x[2];
	u[3]=x[3];
	u[4]=x[4];

	if (isrightside==2)  //convert the value to the other side of Z-axis.   2->down[]  1->up[]   0->the last point.
	{
		u[0]=-x[0];
		u[1]=-x[1];
		u[2]=-x[2];
		u[3]=-x[3];
		u[4]=-x[4];
	}
	int TYPE=0;
	if ( la1==lc1 && la2==lc2 )
		TYPE=1;
	else if (lb1==ld1 && lb2==ld2)
		TYPE=2;
	else
		TYPE=0;

	if (isrightside==3)   // if is the last point.  TYPE=3.
		TYPE=3;

	float D1 = sqrt(pow(la1-lb1,2)+pow(la2-lb2,2));
	float D2 = sqrt(pow(lc1-ld1,2)+pow(lc2-ld2,2));
	D1 = D1/standard;
	D2 = D2/standard;


	float v1 = la1-lb1;
	float v2 = la2-lb2;

	float inter11 = lb1+5*v1/6;
	float inter12 = lb2+5*v2/6;
	float inter21 = lb1+4*v1/6;
	float inter22 = lb2+4*v2/6;
	float inter31 = lb1+3*v1/6;
	float inter32 = lb2+3*v2/6;
	float inter41 = lb1+2*v1/6;
	float inter42 = lb2+2*v2/6;
	float inter51 = lb1+v1/6;
	float inter52 = lb2+v2/6;

	v1 = lc1-ld1;
	v2 = lc2-ld2;

	float outer11 = ld1+5*v1/6;
	float outer12 = ld2+5*v2/6;
	float outer21 = ld1+4*v1/6;
	float outer22 = ld2+4*v2/6;
	float outer31 = ld1+3*v1/6;
	float outer32 = ld2+3*v2/6;
	float outer41 = ld1+2*v1/6;
	float outer42 = ld2+2*v2/6;
	float outer51 = ld1+v1/6;
	float outer52 = ld2+v2/6;
	//float outer52 = ld2+1/6*v2;

	if (TYPE==0)
	{
		Save2trimeshEX(la1,la2,0.0,lc1,lc2,0.0,outer11,outer12,u[0]*D2,0);
		Save2trimeshEX(la1,la2,0.0,inter11,inter12,u[0]*D1,outer11,outer12,u[0]*D2,0);
		Save2trimeshEX(inter11,inter12,u[0]*D1,outer11,outer12,u[0]*D2,outer21,outer22,u[1]*D2,0);
		Save2trimeshEX(inter11,inter12,u[0]*D1,inter21,inter22,u[1]*D1,outer21,outer22,u[1]*D2,0);
		Save2trimeshEX(inter21,inter22,u[1]*D1,outer21,outer22,u[1]*D2,outer31,outer32,u[2]*D2,0);
		Save2trimeshEX(inter21,inter22,u[1]*D1,inter31,inter32,u[2]*D1,outer31,outer32,u[2]*D2,0);
		Save2trimeshEX(inter31,inter32,u[2]*D1,outer31,outer32,u[2]*D2,outer41,outer42,u[3]*D2,0);
		Save2trimeshEX(inter31,inter32,u[2]*D1,inter41,inter42,u[3]*D1,outer41,outer42,u[3]*D2,0);
		Save2trimeshEX(inter41,inter42,u[3]*D1,outer41,outer42,u[3]*D2,outer51,outer52,u[4]*D2,0);
		Save2trimeshEX(inter41,inter42,u[3]*D1,inter51,inter52,u[4]*D1,outer51,outer52,u[4]*D2,0);
		Save2trimeshEX(inter51,inter52,u[4]*D1,outer51,outer52,u[4]*D2,ld1,ld2,0.0,0);
		Save2trimeshEX(inter51,inter52,u[4]*D1,lb1,lb2,0.0,ld1,ld2,0.0,0);
	}
	else if (TYPE==1)
	{
		Save2trimeshEX(la1,la2,0.0,inter11,inter12,u[0]*D1,outer11,outer12,u[0]*D2,0);
		Save2trimeshEX(inter11,inter12,u[0]*D1,outer11,outer12,u[0]*D2,outer21,outer22,u[1]*D2,0);
		Save2trimeshEX(inter11,inter12,u[0]*D1,inter21,inter22,u[1]*D1,outer21,outer22,u[1]*D2,0);
		Save2trimeshEX(inter21,inter22,u[1]*D1,outer21,outer22,u[1]*D2,outer31,outer32,u[2]*D2,0);
		Save2trimeshEX(inter21,inter22,u[1]*D1,inter31,inter32,u[2]*D1,outer31,outer32,u[2]*D2,0);
		Save2trimeshEX(inter31,inter32,u[2]*D1,outer31,outer32,u[2]*D2,outer41,outer42,u[3]*D2,0);
		Save2trimeshEX(inter31,inter32,u[2]*D1,inter41,inter42,u[3]*D1,outer41,outer42,u[3]*D2,0);
		Save2trimeshEX(inter41,inter42,u[3]*D1,outer41,outer42,u[3]*D2,outer51,outer52,u[4]*D2,0);
		Save2trimeshEX(inter41,inter42,u[3]*D1,inter51,inter52,u[4]*D1,outer51,outer52,u[4]*D2,0);
		Save2trimeshEX(inter51,inter52,u[4]*D1,outer51,outer52,u[4]*D2,ld1,ld2,0.0,0);
		Save2trimeshEX(inter51,inter52,u[4]*D1,lb1,lb2,0.0,ld1,ld2,0.0,0);
	}
	else if (TYPE==2)
	{
		Save2trimeshEX(la1,la2,0.0,lc1,lc2,0.0,outer11,outer12,u[0]*D2,0);
		Save2trimeshEX(la1,la2,0.0,inter11,inter12,u[0]*D1,outer11,outer12,u[0]*D2,0);
		Save2trimeshEX(inter11,inter12,u[0]*D1,outer11,outer12,u[0]*D2,outer21,outer22,u[1]*D2,0);
		Save2trimeshEX(inter11,inter12,u[0]*D1,inter21,inter22,u[1]*D1,outer21,outer22,u[1]*D2,0);
		Save2trimeshEX(inter21,inter22,u[1]*D1,outer21,outer22,u[1]*D2,outer31,outer32,u[2]*D2,0);
		Save2trimeshEX(inter21,inter22,u[1]*D1,inter31,inter32,u[2]*D1,outer31,outer32,u[2]*D2,0);
		Save2trimeshEX(inter31,inter32,u[2]*D1,outer31,outer32,u[2]*D2,outer41,outer42,u[3]*D2,0);
		Save2trimeshEX(inter31,inter32,u[2]*D1,inter41,inter42,u[3]*D1,outer41,outer42,u[3]*D2,0);
		Save2trimeshEX(inter41,inter42,u[3]*D1,outer41,outer42,u[3]*D2,outer51,outer52,u[4]*D2,0);
		Save2trimeshEX(inter41,inter42,u[3]*D1,inter51,inter52,u[4]*D1,outer51,outer52,u[4]*D2,0);
		Save2trimeshEX(inter51,inter52,u[4]*D1,outer51,outer52,u[4]*D2,lb1,lb2,0.0,0);
	}
	else if (TYPE==3)
	{
		Save2trimeshEX(la1,la2,0.0,lc1,lc2,0.0,inter11,inter12,u[0]*D1,0);
		Save2trimeshEX(inter11,inter12,u[0]*D1,inter21,inter22,u[1]*D1,lc1,lc2,0.0,0);
		Save2trimeshEX(inter21,inter22,u[1]*D1,inter31,inter32,u[2]*D1,lc1,lc2,0.0,0);
		Save2trimeshEX(inter31,inter32,u[2]*D1,inter41,inter42,u[3]*D1,lc1,lc2,0.0,0);
		Save2trimeshEX(inter41,inter42,u[3]*D1,inter51,inter52,u[4]*D1,lc1,lc2,0.0,0);
		Save2trimeshEX(inter51,inter52,u[4]*D1,lb1,lb2,0.01,lc1,lc2,0.0,0);
	}
}


// output is the direction vector of stoke2
// The 3D stoke triangles will be stored in triEX temperaily.
// After rotation and translation, triEX will restore in tri3D and triEX will be removed.
struct trimesh * makefull3D(mypoint * header)
{
	trimesh * Vector_Center = new trimesh;               // Vector_Center,->a the center of Start & End point of stroke2,  and ->b the direction vector.
	                                                     // Vector_Center->c  and ->d  the start and end point of stroke2.
	trimesh * A = baseringLine(header, Vector_Center);   // A->a,b intersection,  ->c  gravity center.
	float x1,x2,y1,y2,v1,v2;
	float u[5], d[5] ,uu[7][2], dd[5][2];  // uu and dd store the intersection for base ring.
	v1 = A->c1 - Vector_Center->a1;
	v2 = A->c2 - Vector_Center->a2;

	x1 = v1 + Vector_Center->c1;
	x2 = v2 + Vector_Center->c2;
	y1 = v1 + Vector_Center->d1;
	y2 = v2 + Vector_Center->d2;

	findhight(A->a1, A->a2, A->b1, A->b2, header, u, d, uu, dd);

	findbaseringP(x1,x2,y1,y2,u,d);    // crate a chain of mypoint with refined base ring.

	float standard = sqrt(pow(A->a1-A->b1,2)+pow(A->a2-A->b2,2));   // standard distance corresponding to u[5] and d[5].

	mypoint * currentH = header;

	while (currentH!=NULL)
	{
		if (currentH->index==2)
			break;
		currentH = currentH->nextpoint ;
	}

	mypoint * stroke2 = refine2stroke( currentH , 2, 15.0 );

	mypoint * Inline = findstr2inline( stroke2 );

	mypoint * current = Inline;
	mypoint * temp;
	while (current->num!=0)
	{
		temp = current;
		current = current->nextpoint;
		if ( current->num==0)
			break;
		Inline2triEX(temp->x,temp->y,temp->xf,temp->yf,current->x,current->y,current->xf,current->yf,u,standard,1);
		Inline2triEX(temp->x,temp->y,temp->xf,temp->yf,current->x,current->y,current->xf,current->yf,d,standard,2);
	}
	temp = current->frontpoint;
	d[0]=-d[0];
	d[1]=-d[1];
	d[2]=-d[2];
	d[3]=-d[3];
	d[4]=-d[4];
	Inline2triEX(temp->x,temp->y,temp->xf,temp->yf,current->x,current->y,1+current->x,1+current->y,u,standard,3);
	Inline2triEX(temp->x,temp->y,temp->xf,temp->yf,current->x,current->y,1+current->x,1+current->y,d,standard,3);

	Vector_Center->nextTri = A;
	return Vector_Center;
}


//=====extrusion for object creating=====


//=====  for rotation and translating====

//only pick up the first stroke in header and store in header1.
struct mypoint *pick1stroke(mypoint* header)
{
	mypoint * current = header;
	mypoint * header1 = NULL;
	while (current!=NULL && current->index==1)
	{
		header1 = Add_mypoint_inorder( header1, current->x, current->y, 0.0, 0.0, current->num);
		current = current->nextpoint;
	}
	return header1;
}


// given angle and corss product vector, compute the 3by3 rotation matrix.
void rotatematrix(float Angle, float vx, float vy, float vz,float m[3][3])
{
	m[0][0] = vx*vx*(1-cos(Angle))+cos(Angle);
	m[0][1] = vx*vy*(1-cos(Angle))-vz*sin(Angle);
	m[0][2] = vx*vz*(1-cos(Angle))+vy*sin(Angle);
	m[1][0] = vx*vy*(1-cos(Angle))+vz*sin(Angle);
	m[1][1] = vy*vy*(1-cos(Angle))+cos(Angle);
	m[1][2] = vz*vy*(1-cos(Angle))-vx*sin(Angle);
	m[2][0] = vx*vz*(1-cos(Angle))-vy*sin(Angle);
	m[2][1] = vy*vz*(1-cos(Angle))+vx*sin(Angle);
	m[2][2] = vz*vz*(1-cos(Angle))+cos(Angle);
}


//  trimesh * current = one of triEx.
//  output tri3D with rotated and translated values of corresponding triEx.
void transtriEX(float t1, float t2, float t3, float m[3][3], float n[3][3], trimesh * current)
{
	float a1 = current->a1;
	float a2 = current->a2;
	float a3 = current->a3;
	float b1 = current->b1;
	float b2 = current->b2;
	float b3 = current->b3;
	float c1 = current->c1;
	float c2 = current->c2;
	float c3 = current->c3;

	current->a1 = m[0][0]*a1+m[0][1]*a2+m[0][2]*a3;
	current->a2 = m[1][0]*a1+m[1][1]*a2+m[1][2]*a3;
	current->a3 = m[2][0]*a1+m[2][1]*a2+m[2][2]*a3;

	current->b1 = m[0][0]*b1+m[0][1]*b2+m[0][2]*b3;
	current->b2 = m[1][0]*b1+m[1][1]*b2+m[1][2]*b3;
	current->b3 = m[2][0]*b1+m[2][1]*b2+m[2][2]*b3;

	current->c1 = m[0][0]*c1+m[0][1]*c2+m[0][2]*c3;
	current->c2 = m[1][0]*c1+m[1][1]*c2+m[1][2]*c3;
	current->c3 = m[2][0]*c1+m[2][1]*c2+m[2][2]*c3;
	
	a1 = current->a1;
	a2 = current->a2;
	a3 = current->a3;
	b1 = current->b1;
	b2 = current->b2;
	b3 = current->b3;
	c1 = current->c1;
	c2 = current->c2;
	c3 = current->c3;

	current->a1 = n[0][0]*a1+n[0][1]*a2+n[0][2]*a3;
	current->a2 = n[1][0]*a1+n[1][1]*a2+n[1][2]*a3;
	current->a3 = n[2][0]*a1+n[2][1]*a2+n[2][2]*a3;

	current->b1 = n[0][0]*b1+n[0][1]*b2+n[0][2]*b3;
	current->b2 = n[1][0]*b1+n[1][1]*b2+n[1][2]*b3;
	current->b3 = n[2][0]*b1+n[2][1]*b2+n[2][2]*b3;

	current->c1 = n[0][0]*c1+n[0][1]*c2+n[0][2]*c3;
	current->c2 = n[1][0]*c1+n[1][1]*c2+n[1][2]*c3;
	current->c3 = n[2][0]*c1+n[2][1]*c2+n[2][2]*c3;
	
	a1 = current->a1;
	a2 = current->a2;
	a3 = current->a3;
	b1 = current->b1;
	b2 = current->b2;
	b3 = current->b3;
	c1 = current->c1;
	c2 = current->c2;
	c3 = current->c3;

	
	current->a1 = a1 + t1;
	current->a2 = a2 + t2;
	current->a3 = a3 + t3;

	current->b1 = b1 + t1;
	current->b2 = b2 + t2;
	current->b3 = b3 + t3;

	current->c1 = c1 + t1;
	current->c2 = c2 + t2;
	current->c3 = c3 + t3;
	
	Save2trimesh3D(current->a1,current->a2,current->a3,current->b1,current->b2,current->b3,current->c1,current->c2,current->c3,0);
	
}


// rotate and tanslate triEx then store in tri3D.
void triEX2tri3D (mypoint * header,float Xroat, float Yroat)
{
	// first make the full 3D of stroke2 on X-Y plane.
	trimesh * localinfo = makefull3D(header);    // localinfo = Vector_Center + A.

	mypoint * header1 = pick1stroke(header);     // get the first stroke in header1.

	// need to combine first intersection and the G center with header1.
	mypoint * temp1 = new mypoint, * temp2 = new mypoint;     // add localinformation into header1(first two of header1).
	temp1->x = localinfo->nextTri->c1;
	temp1->y = localinfo->nextTri->c2;     // gravity center.
	temp1->index = 1;
	temp1->num = 1;
	temp1->nextpoint = temp2;
	temp2->x = localinfo->nextTri->a1;
	temp2->y = localinfo->nextTri->a2;     // up intersection.
	temp2->index = 1;
	temp2->num = 1;
	temp2->nextpoint = header1;
	header1->frontpoint = temp2;
	header1 = temp1;

	point3D * baseringOut = NULL;
	point3D * refined_basering = NULL;
	point3D * current, * temper1, * temper2;
	baseringOut = locate3Dstrokes(tri3D, Xroat, Yroat, header1, baseringOut);   // find the spatial location of basering.
	refined_basering = locate3Dstrokes(tri3D, Xroat, Yroat, baseringP, refined_basering);

	current = baseringOut->next3D->next3D;    // the actually first point of basering.


	point3D * currentA = refined_basering;
	while (currentA!=NULL)    // store in pointOut to show it up.
	{
		pointOut = Add2point3D(currentA->a1,currentA->a2,currentA->a3,pointOut,indexI,97);
		currentA = currentA->next3D;
	}
	indexI++;

	

	float sumXY=0.0, sumYZ=0.0, sumZX=0.0;

	while ( current!=NULL )
	{
		temper1 = current;
		current = current->next3D ;
		temper2 = current;
		if ( current==NULL)
			temper2 = baseringOut->next3D->next3D;
		sumXY += temper1->a1*temper2->a2-temper2->a1*temper1->a2;
		sumYZ += temper1->a2*temper2->a3-temper2->a2*temper1->a3;
		sumZX += temper1->a3*temper2->a1-temper2->a3*temper1->a1;
	}
	point3D * normalV = new point3D;
	point3D * upV3D = new point3D;
	point3D * upV2D = new point3D;
	point3D * normalV2D = new point3D;

	normalV->a1 = 0.5*sumYZ;
	normalV->a2 = 0.5*sumZX;
	normalV->a3 = 0.5*sumXY;      // now get the normal vector of basering.


	upV3D->a1 = baseringOut->next3D->a1 - baseringOut->a1;
	upV3D->a2 = baseringOut->next3D->a2 - baseringOut->a2;
	upV3D->a3 = baseringOut->next3D->a3 - baseringOut->a3;    // up vector in 3D.

	upV2D->a1 = localinfo->nextTri->a1 - localinfo->nextTri->c1;
	upV2D->a2 = localinfo->nextTri->a2 - localinfo->nextTri->c2;
	upV2D->a3 = 0.0;     // up vector in 2D.

	normalV2D->a1 = localinfo->b1;
	normalV2D->a2 = localinfo->b2;
	normalV2D->a3 = 0.0;       // normal vector in 2D.


	// cross product between normalV and normalV2D.
	float nv1 = normalV->a2*normalV2D->a3 - normalV2D->a2*normalV->a3;
	float nv2 = normalV->a3*normalV2D->a1 - normalV2D->a3*normalV->a1;
	float nv3 = normalV->a1*normalV2D->a2 - normalV2D->a1*normalV->a2;
	float sum_nv = sqrt(nv1*nv1+nv2*nv2+nv3*nv3);

	    //==========================


	nv1 = nv1/sum_nv;
	nv2 = nv2/sum_nv;
	nv3 = nv3/sum_nv;

	float Angle1 = acos((normalV->a1*normalV2D->a1+normalV->a2*normalV2D->a2)/sqrt(
		normalV->a1*normalV->a1+normalV->a2*normalV->a2+normalV->a3*normalV->a3)/
		sqrt(normalV2D->a1*normalV2D->a1+normalV2D->a2*normalV2D->a2));

	float m[3][3];    // 1st rotation matrix.
	rotatematrix(-Angle1, nv1, nv2, nv3, m);

	float up2Dx = m[0][0]*upV2D->a1 + m[0][1]*upV2D->a2 + m[0][2]*upV2D->a3;
	float up2Dy = m[1][0]*upV2D->a1 + m[1][1]*upV2D->a2 + m[1][2]*upV2D->a3;
	float up2Dz = m[2][0]*upV2D->a1 + m[2][1]*upV2D->a2 + m[2][2]*upV2D->a3;

	float Angle2 = acos((upV3D->a1*up2Dx + upV3D->a2*up2Dy + upV3D->a3*up2Dz)/sqrt(
		upV3D->a1*upV3D->a1+upV3D->a2*upV3D->a2+upV3D->a3*upV3D->a3)/sqrt(
		up2Dx*up2Dx+up2Dy*up2Dy+up2Dz*up2Dz));

	sum_nv = sqrt(normalV->a1*normalV->a1+normalV->a2*normalV->a2+normalV->a3*normalV->a3);
	nv1 = normalV->a1/sum_nv;
	nv2 = normalV->a2/sum_nv;
	nv3 = normalV->a3/sum_nv;

	float n[3][3];    // 2ed rotation matrix.
	rotatematrix(-Angle2, nv1, nv2, nv3, n);

	//=======find the translation distance=======
	float N2Dx = m[0][0]*localinfo->a1 + m[0][1]*localinfo->a2 + m[0][2]*0.0;   // find the center points(relate with gravity) after TWO rotation only.
	float N2Dy = m[1][0]*localinfo->a1 + m[1][1]*localinfo->a2 + m[1][2]*0.0;
	float N2Dz = m[2][0]*localinfo->a1 + m[2][1]*localinfo->a2 + m[2][2]*0.0;
	float Nx = N2Dx, Ny = N2Dy, Nz = N2Dz;


	N2Dx = n[0][0]*Nx + n[0][1]*Ny + n[0][2]*Nz;
	N2Dy = n[1][0]*Nx + n[1][1]*Ny + n[1][2]*Nz;
	N2Dz = n[2][0]*Nx + n[2][1]*Ny + n[2][2]*Nz;

	float normaL = sqrt(normalV->a1*normalV->a1+normalV->a2*normalV->a2+normalV->a3*normalV->a3);  // gravity3D + (-)normalV: give the depth.
	float na1 = -2*normalV->a1/normaL;
	float na2 = -2*normalV->a2/normaL;
	float na3 = -2*normalV->a3/normaL;

	//======================================


	float transX = baseringOut->a1 + na1 - N2Dx;
	float transY = baseringOut->a2 + na2 - N2Dy;
	float transZ = baseringOut->a3 + na3 - N2Dz;     // translation vector.

	trimesh * currentT = triEx;
	while (currentT!=NULL)
	{
		transtriEX(transX,transY,transZ,m,n,currentT);
		currentT = currentT->nextTri;
	}
}

//============ end rotation and translating============



//============ Eraser =============


// rotate one current of tri3D with xRot and yRot.   sign>0 than xRot,yRot>0.
void rotatepoint3D (float Xroat, float Yroat, point3D * currentP, float sign)
{
	if (sign<0.0)
	{
		Xroat = -Xroat;
		Yroat = -Yroat;
	}
	Xroat = Xroat*3.1416/180;
	Yroat = Yroat*3.1416/180;
	float a1,a2,a3,b1,b2,b3,c1,c2,c3;
	
	a1 = currentP->a1;
	a2 = currentP->a2 * cos(Xroat) - currentP->a3 * sin(Xroat);
	a3 = currentP->a2 * sin(Xroat) + currentP->a3 * cos(Xroat);

	currentP->a1 = a3*sin(Yroat)+a1*cos(Yroat);
	currentP->a2 = a2;
	currentP->a3 = a3*cos(Yroat)-a1*sin(Yroat);
}

// rotate  first Y then X.
void rotatepoint3DYX (float Xroat, float Yroat, point3D * currentP, float sign)
{
	if (sign<0.0)
	{
		Xroat = -Xroat;
		Yroat = -Yroat;
	}
	Xroat = Xroat*3.1416/180;
	Yroat = Yroat*3.1416/180;
	float a1,a2,a3,b1,b2,b3,c1,c2,c3;
	
	a1 = currentP->a3 * sin(Yroat) + currentP->a1 * cos(Yroat);
	a2 = currentP->a2;
	a3 = currentP->a3 * cos(Yroat) - currentP->a1 * sin(Yroat);

	currentP->a1 = a1;
	currentP->a2 = a2*cos(Xroat) - a3*sin(Xroat);
	currentP->a3 = a2*sin(Xroat) + a3*cos(Xroat);
}


//  given int X, delete the strokes with index = X.
struct point3D * deletePainting( int X )
{
	pointOut->front3D = NULL;
	point3D * current = pointOut, * next;
	if ( current->next3D==NULL )
	{
		delete current;
		pointOut = NULL;
		return NULL;
	}

	while (current->next3D!=NULL)   // go through pointOut.
	{
		if ( current->X == X )
		{
			next = current->next3D;
			if (current->front3D==NULL)   // if the stroke is the first stroke( need to delete pointOut ).
			{
				delete current;
			    current = next;
				pointOut = current;
				current->front3D=NULL;
			    continue;
			}
			else     // not the first stroke ( no need to handle pointOut ).
			{
				next->front3D = current->front3D;
				delete current;
			    current = next;
				continue;
			}

		}
		current = current->next3D;
	}

	if (current->X == X)
	{
		next = current;
		current = current->front3D;
		delete next;
		if (current == NULL)
		{
			pointOut=NULL;
			return NULL;
		}
		current->next3D = NULL;
		return current;
	}

	current = pointOut;
	while (current!=NULL)   // need to return the current point3D if not the last stroke.
	{
		if (current->X == X+1)
			return current;
		current = current->next3D;
	}
}

// remove all Painting strokes if the strokes in header do intersected with painting.
void antiPainting( void )
{
	mypoint * current = header;
	if (header->nextpoint==NULL)
		return;

	point3D * currentP =pointOut;   // the pointOut must be already rotated.
	mypoint a,b,m,n,Int;
	int ind, type, indOut;

	indOut = current->index;
	while ( current->nextpoint!=NULL )   // for every line seg in header.
	{
		while ( current->nextpoint!=NULL && current->nextpoint->index == indOut )
		{
		  a.x = current->x;
		  a.y = current->y;
		  current = current->nextpoint;
		  b.x = current->x;
		  b.y = current->y;

		  ind = pointOut->X;
		  currentP = pointOut;
		  while ( currentP->next3D!=NULL )   // for every polyline in pointOut.
		  {
			  while ( currentP->next3D!=NULL && currentP->next3D->X == ind )
			  {
				m.x = currentP->a1;
				m.y = currentP->a2;
				currentP = currentP->next3D;
				n.x = currentP->a1;
				n.y = currentP->a2;
				type = Intersection(a,b,m,n,Int);
				if (type==2 || type==3 || type==1)
				  {
					currentP = deletePainting(ind);
					if (pointOut==NULL)
						return;
					if (currentP->X!=ind)
						goto out;
				  }
				
			  }
			
			  if (currentP->next3D==NULL)   // the end of pointOut, break.
				  break;
			  currentP = currentP->next3D;  // the end of one stroke in pointOut, change to next stroke.
			  ind = currentP->X;
		    }
		  out:;
		  }
		  if (current->nextpoint==NULL)
			  break;
		  current = current->nextpoint;
		  indOut = current->index;
		}
}

// eraising those strokes that do not want.
void eraser(float Xroat, float Yroat)
{
	point3D * current = pointOut;
	while (current!=NULL)
	{
		rotatepoint3DYX(Xroat,Yroat,current,1.0);
		current = current->next3D;
	}

	antiPainting();

	current = pointOut;

	while (current!=NULL)
	{
		rotatepoint3D(Xroat,Yroat,current,-1.0);
		current = current->next3D;
	}
}


//============ end Eraser ============



//============ Cutting ============

// rotate one current of tri3D with xRot and yRot.   sign>0 than xRot,yRot>0.
void rotatetrimesh (float Xroat, float Yroat, trimesh * currentP, float sign)
{
	if (sign<0.0)
	{
		Xroat = -Xroat;
		Yroat = -Yroat;
	}
	Xroat = Xroat*3.1416/180;
	Yroat = Yroat*3.1416/180;
	float a1,a2,a3,b1,b2,b3,c1,c2,c3;
	
	a1 = currentP->a1;
	a2 = currentP->a2 * cos(Xroat) - currentP->a3 * sin(Xroat);
	a3 = currentP->a2 * sin(Xroat) + currentP->a3 * cos(Xroat);

	currentP->a1 = a3*sin(Yroat)+a1*cos(Yroat);
	currentP->a2 = a2;
	currentP->a3 = a3*cos(Yroat)-a1*sin(Yroat);

	b1 = currentP->b1;
	b2 = currentP->b2 * cos(Xroat) - currentP->b3 * sin(Xroat);
	b3 = currentP->b2 * sin(Xroat) + currentP->b3 * cos(Xroat);

	currentP->b1 = b3*sin(Yroat)+b1*cos(Yroat);
	currentP->b2 = b2;
	currentP->b3 = b3*cos(Yroat)-b1*sin(Yroat);
	
	c1 = currentP->c1;
	c2 = currentP->c2 * cos(Xroat) - currentP->c3 * sin(Xroat);
	c3 = currentP->c2 * sin(Xroat) + currentP->c3 * cos(Xroat);

	currentP->c1 = c3*sin(Yroat)+c1*cos(Yroat);
	currentP->c2 = c2;
	currentP->c3 = c3*cos(Yroat)-c1*sin(Yroat);

}

// rotate  first Y then X.
void rotatetrimeshYX (float Xroat, float Yroat, trimesh * currentP, float sign)
{
	if (sign<0.0)
	{
		Xroat = -Xroat;
		Yroat = -Yroat;
	}
	Xroat = Xroat*3.1416/180;
	Yroat = Yroat*3.1416/180;
	float a1,a2,a3,b1,b2,b3,c1,c2,c3;
	
	a1 = currentP->a3 * sin(Yroat) + currentP->a1 * cos(Yroat);
	a2 = currentP->a2;
	a3 = currentP->a3 * cos(Yroat) - currentP->a1 * sin(Yroat);

	currentP->a1 = a1;
	currentP->a2 = a2*cos(Xroat) - a3*sin(Xroat);
	currentP->a3 = a2*sin(Xroat) + a3*cos(Xroat);

	b1 = currentP->b3 * sin(Yroat) + currentP->b1 * cos(Yroat);
	b2 = currentP->b2;
	b3 = currentP->b3 * cos(Yroat) - currentP->b1 * sin(Yroat);

	currentP->b1 = b1;
	currentP->b2 = b2*cos(Xroat) - b3*sin(Xroat);
	currentP->b3 = b2*sin(Xroat) + b3*cos(Xroat);
	
	c1 = currentP->c3 * sin(Yroat) + currentP->c1 * cos(Yroat);
	c2 = currentP->c2;
	c3 = currentP->c3 * cos(Yroat) - currentP->c1 * sin(Yroat);

	currentP->c1 = c1;
	currentP->c2 = c2*cos(Xroat) - c3*sin(Xroat);
	currentP->c3 = c2*sin(Xroat) + c3*cos(Xroat);

}


// find the nearest line seg of that given point. 
// return the current(that line seg) of stroke1(need to delete after using).
// if A->B  then output->xy = A,  output->xf,yf = B.
struct mypoint * findnearline ( mypoint * stroke1, float x, float y )
{
	mypoint * current = stroke1;
	float temp1, temp2, c1,c2, Distance, tempD=-1.0;
	int count;


	while (current!=NULL)
	{
		temp1 = current->x;
		temp2 = current->y;
		current = current->nextpoint;
		if (current==NULL)
			break;
		c1 = (temp1+current->x)/2;
		c2 = (temp2+current->y)/2;  // center of the line seg.

		Distance = sqrt((x-c1)*(x-c1)+(y-c2)*(y-c2));
		if (tempD<0)
		{
			tempD = Distance;
			count = current->num;
			continue;
		}
		else if (tempD> Distance)   // see if short than previous one.
		{
			tempD = Distance;     // if so, replace it.
			count = current->num;
		}
	}
	current = stroke1;
	while (current!=NULL)
	{
		if (current->num == count-1)     // find the nearest line seg.
			break;
		current = current->nextpoint;
	}

	mypoint * output = new mypoint;
	output->x = current->x;
	output->y = current->y;
	output->xf = current->nextpoint->x;
	output->yf = current->nextpoint->y;

	return output;
}



// A->B  C:x1,x2
// determine if C on the left side of A->B.
 bool isontheleft(float a1, float a2, float b1, float b2, float x1, float x2)
 {
	 float ax1, ax2, m1, m2;
	 m1 = b2 - a2;   //  Vab2.
	 m2 = a1 - b1;   // -Vab1.
	 ax1 = x1 - a1;
	 ax2 = x2 - a2;

	 if ((m1*ax1+m2*ax2)<0)
		 return true;
	 else
		 return false;
 }


 // find the cut line's correspoinding projective vertex on the object,
 // then store then in trimesh in which it has two pair of end points with an ordered num.
 // always the near surface's point are on top( output->a1~b3).
 struct trimesh * findlineintersection( mypoint * stroke1 )
 {
	float t, u, v, a, b, c, m, n ,p, a1, b1, c1, m1, n1, p1;
	point3D dir;
	dir.a1 = 0.0 ;
	dir.a2 = 0.0 ;
	dir.a3 = -1000.0 ;
	mypoint * current=stroke1;
	point3D orig;
	trimesh * currentri=tri3D;


	bool type = false;
	float d = 0;
	int count = 0;
	int round = 0;
	int num = 2;
	trimesh * currentO, * trans=NULL;
	trimesh * output = new trimesh;
	currentO = output;

	while (current!=NULL)
	{
		orig.a1=current->x;
		orig.a2=current->y;
		orig.a3=1000.0;
		currentri = tri3D;
		while (currentri!=NULL)
		{
			type = TriangleIntersect(currentri, orig, dir, &t, &u, &v);
			if (type==true)
			{
				if ( round==0)
				{
				    if (count==0)
				    {
					   a = orig.a1+dir.a1*t;    // use orig + dir * d to calculate stroke point's location.
		               b = orig.a2+dir.a2*t;
		               c = orig.a3+dir.a3*t;
					   d = t;
					   count = 1;
			   	     }
				    else
					{
					   m = orig.a1+dir.a1*t;    // use orig + dir * d to calculate stroke point's location.
		               n = orig.a2+dir.a2*t;
		               p = orig.a3+dir.a3*t;
					   if (d>t)                 // always make the near vertex on top.
					   {
						  swap(a,m);
						  swap(b,n);
						  swap(c,p);
					   }
					   count = 0;
					   round = 1;
					   break;
				     }
				 }
				if ( round == 1 )
				{
					if (count==0)
				    {
					   a1 = orig.a1+dir.a1*t;    // use orig + dir * d to calculate stroke point's location.
		               b1 = orig.a2+dir.a2*t;
		               c1 = orig.a3+dir.a3*t;
					   d = t;
					   count = 1;
			   	     }
				    else
					{
					   m1 = orig.a1+dir.a1*t;    // use orig + dir * d to calculate stroke point's location.
		               n1 = orig.a2+dir.a2*t;
		               p1 = orig.a3+dir.a3*t;
					   if (d>t)                 // always make the near vertex on top.
					   {
						  swap(a1,m1);
						  swap(b1,n1);
						  swap(c1,p1);
					   }
					   count = 0;
					   round = 2;
					   currentO->a1 = a;
					   currentO->a2 = b;
					   currentO->a3 = c;
					   currentO->b1 = m;
					   currentO->b2 = n;
					   currentO->b3 = p;
					   currentO->c1 = a1;
					   currentO->c2 = b1;
					   currentO->c3 = c1;
					   currentO->d1 = m1;
					   currentO->d2 = n1;
					   currentO->d3 = p1;
					   currentO->X = num;    // for 1st line segment.
					   currentO->link3D = NULL;
					   num++;
					   currentO->nextTri = new trimesh;
					   trans = currentO;
					   currentO = currentO->nextTri;
					   currentO->frontTri = trans;
					   break;
				     }
				}
				if ( round == 2 )
				{
					if (count==0)
				    {
						swap(a1,a);
						swap(b1,b);
						swap(c1,c);
						swap(m1,m);
						swap(n1,n);
						swap(p1,p);

						a1 = orig.a1+dir.a1*t;    // use orig + dir * d to calculate stroke point's location.
		                b1 = orig.a2+dir.a2*t;
		                c1 = orig.a3+dir.a3*t;
					    d = t;
					    count = 1;
			   	     }
				    else
					{
					   m1 = orig.a1+dir.a1*t;    // use orig + dir * d to calculate stroke point's location.
		               n1 = orig.a2+dir.a2*t;
		               p1 = orig.a3+dir.a3*t;
					   if (d>t)                 // always make the near vertex on top.
					   {
						  swap(a1,m1);
						  swap(b1,n1);
						  swap(c1,p1);
					   }
					   count = 0;
					   currentO->a1 = a;
					   currentO->a2 = b;
					   currentO->a3 = c;
					   currentO->b1 = m;
					   currentO->b2 = n;
					   currentO->b3 = p;
					   currentO->c1 = a1;
					   currentO->c2 = b1;
					   currentO->c3 = c1;
					   currentO->d1 = m1;
					   currentO->d2 = n1;
					   currentO->d3 = p1;
					   currentO->X = num;    // for 1st line segment.
					   currentO->link3D = NULL;
					   num++;
					   currentO->nextTri = new trimesh;
					   trans = currentO;
					   currentO = currentO->nextTri;
					   currentO->frontTri = trans;
					   break;
				     }
				}
			}
			currentri = currentri->nextTri ;
		}
		current = current->nextpoint ;
	}

	current = stroke1;
	trimesh * currentO2;
	if (trans!=NULL)
	{
	while ( current!=NULL )
	{
		if (abs(current->x - trans->c1)<0.001 && abs(current->y - trans->c2)<0.001)   // if find the final point in output.
		{
			currentO->a1 = a1;
			currentO->a2 = b1;
			currentO->a3 = c1;
			currentO->b1 = m1;
			currentO->b2 = n1;
			currentO->b3 = p1;
			currentO->c1 = current->nextpoint->x;
			currentO->c2 = current->nextpoint->y;
			currentO->c3 = 0;
			currentO->d1 = current->nextpoint->x;
			currentO->d2 = current->nextpoint->y;
			currentO->d3 = 0;
			currentO->X = num;
			currentO->link3D = NULL;
			currentO->nextTri = NULL;
		}
		if (abs(current->x - output->a1)<0.001 && abs(current->y - output->a2)<0.001)   // if find the start point in output.
		{
			currentO2 = new trimesh;

			currentO2->a1 = current->frontpoint->x;
			currentO2->a2 = current->frontpoint->y;
			currentO2->a3 = 1000.0;
			currentO2->b1 = current->frontpoint->x;
			currentO2->b2 = current->frontpoint->y;
			currentO2->b3 = 1000.0;
			currentO2->c1 = output->a1;
			currentO2->c2 = output->a2;
			currentO2->c3 = output->a3;
			currentO2->d1 = output->b1;
			currentO2->d2 = output->b2;
			currentO2->d3 = output->b3;
			currentO2->X = 1;
			currentO2->link3D = NULL;

			currentO2->nextTri = output;
			output->frontTri = currentO2;
			output = currentO2;
		}
		current = current->nextpoint;
	}
	}
	else
	{
		current = stroke1;
		while ( current!=NULL )
	   {
		if (abs(current->x - a)<0.001 && abs(current->y - b)<0.001)   // if find the first point in output.
		{
			currentO->a1 = current->frontpoint->x;
			currentO->a2 = current->frontpoint->y;
			currentO->a3 = 0;
			currentO->b1 = current->frontpoint->x;
			currentO->b2 = current->frontpoint->y;
			currentO->b3 = 0;
			currentO->c1 = a;
			currentO->c2 = b;
			currentO->c3 = c;
			currentO->d1 = m;
			currentO->d2 = n;
			currentO->d3 = p;
			currentO->X = 1;
			currentO->link3D = NULL;
			currentO->nextTri = new trimesh;
			num++;
			currentO2 = currentO->nextTri;
			currentO2->a1 = a;
			currentO2->a2 = b;
			currentO2->a3 = c;
			currentO2->b1 = m;
			currentO2->b2 = n;
			currentO2->b3 = p;
			currentO2->c1 = current->nextpoint->x;
			currentO2->c2 = current->nextpoint->y;
			currentO2->c3 = 1000;
			currentO2->d1 = current->nextpoint->x;
			currentO2->d2 = current->nextpoint->y;
			currentO2->d3 = 1000;
			currentO2->X = 2;
			currentO2->link3D = NULL;
			currentO2->frontTri = currentO;
			currentO2->nextTri = NULL;
		}
		current = current->nextpoint;
		}
	}

	return output;
}

// for every triangle in tri3D, search if any of the three vertexes are located inside the delete region,
// if so, delete the triangle in tri3D,
// if only 1~2 points inside, store the copy with count number in triEx.
// in triEx, inside points always on the top, that is ->a and/or ->b.
void searchboundarypoint( mypoint * stroke1 )
{
	trimesh * current = tri3D;
	trimesh * temp;
	float a1,a2,a3,b1,b2,b3,c1,c2,c3;
	mypoint * sline;
	bool type;
	int count=0;

	while ( current!= NULL )
	{
		a1 = current->a1;
		a2 = current->a2;
		a3 = current->a3;
		b1 = current->b1;
		b2 = current->b2;
		b3 = current->b3;
		c1 = current->c1;
		c2 = current->c2;
		c3 = current->c3;

		sline = findnearline( stroke1, current->a1, current->a2 );  // find nearest line of a1 a2.
		type = isontheleft(sline->x, sline->y, sline->xf, sline->yf, current->a1, current->a2);

		if ( type==true )   // for first point a.
		{
			count++;
			type = false;
		}
		delete sline;

		sline = findnearline( stroke1, current->b1, current->b2 );
		type = isontheleft(sline->x, sline->y, sline->xf, sline->yf, current->b1, current->b2);

		if ( type==true )    // for second point b.
		{
			type = false;
			if ( count==0 )
			{
				swap(a1,b1);
				swap(a2,b2);
				swap(a3,b3);
			}
			count++;
		}
		delete sline;

		sline = findnearline( stroke1, current->c1, current->c2 );
		type = isontheleft(sline->x, sline->y, sline->xf, sline->yf, current->c1, current->c2);

		if ( type==true )    // for third point c.
		{
			type = false;
			if (count==0)
			{
				swap(a1,c1);
				swap(a2,c2);
				swap(a3,c3);
			}
			if (count==1)
			{
				swap(b1,c1);
				swap(b2,c2);
				swap(b3,c3); 
			}
			count++;      // now, count with varify between 0~3.
		}
		delete sline;

		if ( count==3 ) // all points are in the delete domain.
		{
			if (current==tri3D)
			{
				tri3D = tri3D->nextTri;
				tri3D->frontTri=NULL;
				delete current;
				current = tri3D;
				count = 0;
				continue;
			}
			temp = current->frontTri;
			current = current->nextTri;
			delete temp->nextTri;
			temp->nextTri = current;
			if ( current!=NULL )
				current->frontTri = temp;
		}

		if ( count==1 || count==2 )  // if there are points inside the delete region but not all.
		{
			if (current==tri3D)     // first delete this triangle as above one.
			{
				tri3D = tri3D->nextTri;
				tri3D->frontTri=NULL;
				delete current;
				current = tri3D;
				Save2trimeshEX(a1,a2,a3,b1,b2,b3,c1,c2,c3,count);
				count = 0;
				continue;
			}
			temp = current->frontTri;
			current = current->nextTri;
			delete temp->nextTri;
			temp->nextTri = current;
			if ( current!=NULL )
				current->frontTri = temp;

			Save2trimeshEX(a1,a2,a3,b1,b2,b3,c1,c2,c3,count);
		}

		if ( count == 0 )
			current = current->nextTri;

		count = 0;
	}
}


// input the current of the output of findlineintersection( mypoint * stroke1 ).
// a1 a2 a3 are that intersection between cutline and triangle in tri3D.
// X denotes the index of cutting surface.
void point3Dtrimesh( trimesh * current, float a1, float a2, float a3, int X)
{
	point3D * P = current->link3D;

	if (P!=NULL)
		{
			point3D * trans = new point3D ;
	        point3D * tempOut ;
	        trans->a1 = a1;
			trans->a2 = a2;
			trans->a3 = a3;
			trans->X = X ;

			tempOut = P->next3D ;
			trans->front3D = P ;
	        P->next3D = trans ;
			trans->next3D = tempOut ;
	        if (tempOut!=NULL)
				tempOut->front3D = trans ;
		}

		if (P==NULL)
		{
			P = new point3D ;
			P->a1 = a1;
			P->a2 = a2;
			P->a3 = a3;
			P->X = X ;

			P->next3D = NULL ;
			current->link3D = P;
		}
}


// search through whole line3D to find an intersection with a,b,
// the output is current line3D and the intersected point in x1 and x2.
struct trimesh * searchcutline( trimesh * line3D, float a1, float a2, float b1, float b2, float * x1, float * x2)
{
	trimesh * current = line3D;
	mypoint p1,p2,p3,p4,Int;
	int type;
	while (current!=NULL)
	{
		p1.x = a1;
		p1.y = a2;
		p2.x = b1;
		p2.y = b2;
		p3.x = current->a1;
		p3.y = current->a2;
		p4.x = current->c1;
		p4.y = current->c2;
		type = Intersection(p1,p2,p3,p4,Int);

		if (type == 2 || type == 3 || type == 1)  // do intersected.
		{
			*x1 = Int.x;
			*x2 = Int.y;
			return current;
		}
		current = current->nextTri;
	}
	return NULL;
}



// refine tri3D by introduce new triangles in it.
// refine line3D by adding corresponding intersections.
// point3D->X of line3D will be set to 1 as UnUSEed, if Used, should change point3D->X to 4, and ->X = 2 or 3 repesent top and bottom.
void refineline3D ( trimesh * line3D )
{
	trimesh * currentE = triEx;  // here stores the problem triangles.
	trimesh * temp, *trans;
	float x1, x2, y1, y2;
	float a1,a2,b1,b2,c1,c2,D,Z1,Z2;
	float TYPE;
	float md1, md2, md3;
	int indexa=54;

	while ( currentE!=NULL )
	{
		a1 = currentE->a1;
		a2 = currentE->a2;
		b1 = currentE->b1;
		b2 = currentE->b2;
		c1 = currentE->c1;
		c2 = currentE->c2;
		indexa = 54;
		if (currentE->X==1)   // a problem triangle with only ONE inner point.
		{
			temp = searchcutline(line3D,a1,a2,b1,b2,&x1,&x2);
			D = sqrt((a1-x1)*(a1-x1)+(a2-x2)*(a2-x2))/sqrt((a1-b1)*(a1-b1)+(a2-b2)*(a2-b2));
			Z1 = currentE->a3-(currentE->a3 - currentE->b3)*D;    // find Z-axis of intersection.
			point3Dtrimesh(temp, x1, x2, Z1, 1);
			trans = temp;

			temp = searchcutline(line3D,a1,a2,c1,c2,&y1,&y2);
			if ( temp == trans )
				TYPE=1;    // same line.
			else
				TYPE=2;    // two different cutting lines.

			D = sqrt((a1-y1)*(a1-y1)+(a2-y2)*(a2-y2))/sqrt((a1-c1)*(a1-c1)+(a2-c2)*(a2-c2));
			Z2 = currentE->a3-(currentE->a3 - currentE->c3)*D;    // find Z-axis of intersection.
			point3Dtrimesh(temp, y1, y2, Z2, 1);

			// need to construct new triangles to refine the cutting boundary.
			if (TYPE==1)
			{
				Save2trimesh3D(b1,b2,currentE->b3,x1,x2,Z1,y1,y2,Z2,54);    // num2 54 denotes that the third edge need to be shown up.
				Save2trimesh3D(y1,y2,Z2,b1,b2,currentE->b3,c1,c2,currentE->c3,1);
			}
			if (TYPE==2)
			{
				if (trans->X > temp->X)
				{
					md1 = trans->a1;
					md2 = trans->a2;
					md3 = trans->a3;   // need to comfirm.

					if (abs(md3-currentE->c3)>abs(trans->b3-currentE->c3))
						md3 = trans->b3;
					if (trans->X == 1)
					{
						md1 = trans->c1;
						md2 = trans->c2;
						md3 = trans->c3;
						if (abs(md3-currentE->c3)>abs(trans->b3-currentE->c3))
						md3 = trans->d3;
						indexa = 1;
					}

				}
				else
				{
					md1 = temp->a1;
					md2 = temp->a2;
					md3 = temp->a3;    // need to comfirm.
					if (abs(md3-currentE->c3)>abs(temp->b3-currentE->c3))
						md3 = temp->b3;
					if (temp->X == 1)
					{
						md1 = temp->c1;
						md2 = temp->c2;
						md3 = temp->c3;
						if (abs(md3-currentE->c3)>abs(trans->b3-currentE->c3))
						md3 = trans->d3;
						indexa = 1;
						
					}
				}

				Save2trimesh3D(b1,b2,currentE->b3,md1,md2,md3,x1,x2,Z1,indexa);
				Save2trimesh3D(md1,md2,md3,c1,c2,currentE->c3,b1,b2,currentE->b3,1);
				Save2trimesh3D(c1,c2,currentE->c3,md1,md2,md3,y1,y2,Z2,indexa);
			}
		}

		if (currentE->X==2)   // a problem triangle with TWO inner points.
		{
			temp = searchcutline(line3D,a1,a2,c1,c2,&x1,&x2);
			D = sqrt((a1-x1)*(a1-x1)+(a2-x2)*(a2-x2))/sqrt((a1-c1)*(a1-c1)+(a2-c2)*(a2-c2));
			Z1 = currentE->a3-(currentE->a3 - currentE->c3)*D;    // find Z-axis of intersection.
			point3Dtrimesh(temp, x1, x2, Z1, 1);
			trans = temp;

			temp = searchcutline(line3D,b1,b2,c1,c2,&y1,&y2);
			if ( temp == trans )
				TYPE=2;    // two different cutting lines.
			else
				TYPE=1;    // same line.

			D = sqrt((b1-y1)*(b1-y1)+(b2-y2)*(b2-y2))/sqrt((b1-c1)*(b1-c1)+(b2-c2)*(b2-c2));
			Z2 = currentE->b3-(currentE->b3 - currentE->c3)*D;    // find Z-axis of intersection.
			point3Dtrimesh(temp, y1, y2, Z2, 1);

			// need to construct new triangles to refine the cutting boundary.
			if (TYPE==1)
			{
				Save2trimesh3D(c1,c2,currentE->c3,x1,x2,Z1,y1,y2,Z2,54);
			}
			if (TYPE==2)
			{
				if (trans->X > temp->X)
				{
					md1 = trans->a1;
					md2 = trans->a2;
					md3 = trans->a3;
					if (abs(md3-currentE->c3)>abs(trans->b3-currentE->c3))
						md3 = trans->b3;
					if (trans->X == 1)
					{
						md1 = trans->c1;
						md2 = trans->c2;
						md3 = trans->c3;
						if (abs(md3-currentE->c3)>abs(trans->b3-currentE->c3))
						md3 = trans->d3;
						indexa = 1;
					}
				}
				else
				{
					md1 = temp->a1;
					md2 = temp->a2;
					md3 = temp->a3;
					if (abs(md3-currentE->c3)>abs(temp->b3-currentE->c3))
						md3 = temp->b3;
					if (temp->X == 1)
					{
						md1 = temp->c1;
						md2 = temp->c2;
						md3 = temp->c3;
						if (abs(md3-currentE->c3)>abs(trans->b3-currentE->c3))
						md3 = temp->d3;
						indexa = 1;
					}
				}
				Save2trimesh3D(c1,c2,currentE->c3,x1,x2,Z1,md1,md2,md3,indexa);
				Save2trimesh3D(c1,c2,currentE->c3,y1,y2,Z2,md1,md2,md3,indexa);
			}
		}

		currentE = currentE->nextTri;
	}
}

// delete redundent(same value) point in point3D.
// current will be current->next3D or NULL if at the end.
struct point3D * deletecurrentPoint3D( point3D * current )
{
	point3D * temp, * temp2;
	if (current->next3D==NULL)
	{
		temp = current->front3D;
		delete current;
		temp->next3D = NULL;
		current = NULL;
	}
	else
	{
		temp = current->front3D;
		temp2 = current->next3D;
		delete current;
		current = temp2;
		current->front3D = temp;
		temp->next3D = current;
	}
	return current;
}
		


// delete redundent points in line3D->link3D in that one intersection will appear twice by two adjacent triangle.
void refine2line3D ( trimesh * line3D )
{
	trimesh * currentri = line3D;
	point3D * current1,*current2;
	float a1,a2,a3;

	while ( currentri!=NULL )
	{
		current1 = currentri->link3D;
		current2 = current1->next3D;
		while ( current1!=NULL )
		{
			a1 = current1->a1;
			a2 = current1->a2;
			a3 = current1->a3;
			while ( current2!=NULL )
			{
				if (a1==current2->a1 && a2==current2->a2 && a3==current2->a3)
				{
					current2 = deletecurrentPoint3D( current2 );
					continue;
				}
				current2 = current2->next3D;
			}
			current1 = current1->next3D;
			if (current1==NULL || current1->next3D==NULL)
				break;
			current2 = current1->next3D;
		}
		currentri = currentri->nextTri;
	}
}



// input the current line3D (must on the boundary),
// output the nearest unused point to cutting point top-Z.
// if no further point3D founded, return NULL.
struct point3D * findnearedgeP(trimesh * currentLine3D, bool istart)
{
	float Z, D;
	int type=0;
	if (istart == true)
	{
		Z = currentLine3D->c3;
		D = Z - currentLine3D->d3;
	}
	else
	{
		Z = currentLine3D->a3;
		D = Z - currentLine3D->b3;
	}

	point3D * current = currentLine3D->link3D;
	point3D * temp;
	while ( current!=NULL )
	{
		if ( current->X==1 )
		{
			if ( (Z - current->a3)<D )
			{
				D = Z - current->a3;
				temp = current;
				type = 1;
			}
		}
		current = current->next3D;
	}

	if (type==0)
		return NULL;

	temp->X = 4;
	return temp;
}

// cutting surface triangles ->X = 2;
// store all edge triangles in the tri3D.
// related to the first and the end line in line3D.
void storeedgetri( trimesh * line3D )
{
	trimesh * current = line3D;    //line3D is the start line.
	while (current->nextTri!=NULL)
	{
		current = current->nextTri;
	}   // now, current is the end line.
	// for start line:
	float a1,a2,a3,b1,b2,b3;
	a1 = line3D->c1;
	a2 = line3D->c2;
	a3 = line3D->c3;
	b1 = line3D->d1;
	b2 = line3D->d2;
	b3 = line3D->d3;
	point3D * temp1, *trans;
	int count = 1;
	while (1)
	{
		temp1 = findnearedgeP(line3D, true);
		if (temp1==NULL)
			break;
		if (count==1)
		{
			Save2trimesh3D(a1,a2,a3,temp1->a1,temp1->a2,temp1->a3,a1,a2,temp1->a3,2);
			count++;
			trans = temp1;
			continue;
		}

		Save2trimesh3D(trans->a1,trans->a2,trans->a3,a1,a2,temp1->a3,a1,a2,trans->a3,2);
		Save2trimesh3D(trans->a1,trans->a2,trans->a3,a1,a2,temp1->a3,temp1->a1,temp1->a2,temp1->a3,2);
		trans = temp1;
	}
	Save2trimesh3D(b1,b2,b3,trans->a1,trans->a2,trans->a3,b1,b2,trans->a3,2);

	// for the end line.
	a1 = current->a1;
	a2 = current->a2;
	a3 = current->a3;
	b1 = current->b1;
	b2 = current->b2;
	b3 = current->b3;

	point3D * temp2;
	count = 1;

	while (1)
	{
		temp2 = findnearedgeP(current, false);
		if (temp2==NULL)
			break;
		if (count==1)
		{
			Save2trimesh3D(a1,a2,a3,temp2->a1,temp2->a2,temp2->a3,a1,a2,temp2->a3,2);
			count++;
			trans = temp2;
			continue;
		}

		Save2trimesh3D(trans->a1,trans->a2,trans->a3,a1,a2,temp2->a3,a1,a2,trans->a3,2);
		Save2trimesh3D(trans->a1,trans->a2,trans->a3,a1,a2,temp2->a3,temp2->a1,temp2->a2,temp2->a3,2);
		trans = temp2;
	}
	Save2trimesh3D(b1,b2,b3,trans->a1,trans->a2,trans->a3,b1,b2,trans->a3,2);

}


// using X,Y and Z-axis to determine the closest points to ->a or ->b.
// if isup==ture, then use ->a as base,
// else, use ->b as base.
// return the found point in point3D.
struct point3D * findclosestedgeP(trimesh * currentLine3D, bool isup)
{
	float a1,a2,a3,D,d;
	int type=0;
	int check=0;
	if (isup==true)
	{
		a1= currentLine3D->a1;
	    a2= currentLine3D->a2;
		a3= currentLine3D->a3;
		type = 2;   // related to top vertexes.
		D = sqrt(pow(a1-currentLine3D->c1,2)+pow(a2-currentLine3D->c2,2)+pow(a3-currentLine3D->c3,2));
		D = D*3.0;
	}
	else
	{
		a1= currentLine3D->b1;
	    a2= currentLine3D->b2;
		a3= currentLine3D->b3;
		type = 3;   // related to down vertexes.
		D = sqrt(pow(a1-currentLine3D->d1,2)+pow(a2-currentLine3D->d2,2)+pow(a3-currentLine3D->d3,2));
		D = D*3.0;
	}
	
	point3D * current = currentLine3D->link3D;
	point3D * temp;
	while ( current!=NULL )
	{
		if ( current->X==type )
		{
			d = sqrt(pow(a1-current->a1,2)+pow(a2-current->a2,2)+pow(a3-current->a3,2));
			if ( d < D )
			{
				D = d ;
				temp = current;
				check = 1;
			}
		}
		current = current->next3D;
	}

	if (check==0)
		return NULL;

	temp->X = 4;    //  set to used.
	return temp;
}


//  for every given current in line3D which also inside the object, triangulate this plane, 
//  then store into tri3D.
void storeinedgetri( trimesh * current )
{
	float ma1,ma2,ma3,mb1,mb2,mb3,a1,a2,a3,c1,c2,c3;
	float o1,o2,o3;
	ma1 = (current->a1+current->b1)/2;
	ma2 = (current->a2+current->b2)/2;
	ma3 = (current->a3+current->b3)/2;
	mb1 = (current->c1+current->d1)/2;
	mb2 = (current->c2+current->d2)/2;
	mb3 = (current->c3+current->d3)/2;

	o1 = (ma1+mb1)/2;
	o2 = (ma2+mb2)/2;
	o3 = (ma3+mb3)/2;

	a1 = current->a1;
	a2 = current->a2;
	a3 = current->a3;
	c1 = current->c1;
	c2 = current->c2;
	c3 = current->c3;

	point3D * currentP = current->link3D;
	while( currentP!=NULL )
	{
		if (currentP->a3>o3)
			currentP->X = 2;   // unused points on the top.
		else
			currentP->X = 3;   // unused points at the botton.
		currentP = currentP->next3D;
	}

	Save2trimesh3D(current->a1,current->a2,current->a3,o1,o2,o3,ma1,ma2,ma3,2);
	Save2trimesh3D(current->b1,current->b2,current->b3,o1,o2,o3,ma1,ma2,ma3,2);
	Save2trimesh3D(current->c1,current->c2,current->c3,o1,o2,o3,mb1,mb2,mb3,2);
	Save2trimesh3D(current->d1,current->d2,current->d3,o1,o2,o3,mb1,mb2,mb3,2);

	// for those top vertexes:
	point3D * temp1, *trans=NULL;
	int count = 1;
	int num=1;
	bool isup= true;

	while (num!=3)
	{
	   count = 1;
	   while (1)
	  {
		  temp1 = findclosestedgeP(current, isup);
	  	  if (temp1==NULL)
			  break;
		  if (count==1)
		  {
			  Save2trimesh3D(a1,a2,a3,temp1->a1,temp1->a2,temp1->a3,o1,o2,o3,2);
			  count++;
			  trans = temp1;
			  continue;
		  }

		  Save2trimesh3D(trans->a1,trans->a2,trans->a3,o1,o2,o3,temp1->a1,temp1->a2,temp1->a3,2);
		  trans = temp1;
	  }
	   if (trans!=NULL)
		   Save2trimesh3D(c1,c2,c3,trans->a1,trans->a2,trans->a3,o1,o2,o3,2);
	   else
		   Save2trimesh3D(c1,c2,c3,a1,a2,a3,o1,o2,o3,2);

	  a1 = current->b1;
	  a2 = current->b2;
	  a3 = current->b3;
	  c1 = current->d1;
	  c2 = current->d2;
	  c3 = current->d3;
	  num++;
	  isup = false;
	  trans = NULL;
	}
}


//  make use of all cutting functions to apply it.
void cutting( mypoint * header, float Xroat, float Yroat )
{
	mypoint * stroke1 = refine2stroke( header, 1, 20 );   // refine stroke so that becomes to a polyline with fixed line-seg length.

	trimesh * currentri = tri3D;
	while ( currentri!=NULL )    // first rotate the whole object(tri3D)
	{
		rotatetrimeshYX(Xroat, Yroat, currentri, 1.0);
		currentri = currentri->nextTri;
	}

	trimesh * line3D = findlineintersection( stroke1 );   // project cutting line and find out their location on the object(store in line3D for every pair).

	searchboundarypoint( stroke1 );  // generate the triEx for stroing boundary triangles which intersected with cutting line.

	refineline3D( line3D );   // 1) add new triangles into tri3D;  2) add point3D to line3D.

	refine2line3D( line3D );

	storeedgetri( line3D );   // store first and end cutting planes relate to the first and last nods of line3D.

	trimesh * current = line3D->nextTri;

	while ( current->nextTri!=NULL)
	{
		storeinedgetri( current );   // store inner cutting planes one by one into tri3D.
		current = current->nextTri;
	}
	
	currentri = tri3D;
	while ( currentri!=NULL )    // finally rotate the whole object(tri3D) back to normal.
	{
		rotatetrimesh(Xroat, Yroat, currentri, -1.0);
		currentri = currentri->nextTri;
	}

	Delete_mypoint(header);    // clear up.
	header = NULL;
	index_all=1;
	num_all=1;
	Delete_trimesh(triEx);
	Delete_trimesh(line3D);
	triEx = NULL;
}




//============ For test =============


//============== GL FUNCTION =================

void init (void)
{
	glViewport (0, 0, winWidth, winHeight);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	//gluPerspective(45.0, (float)winWidth/winHeight, 0.001f, 1000.0);
	gluOrtho2D(0.0, 300.0, 0.0, 300.0);

	glClearColor (1.0, 1.0, 1.0, 1.0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(0.0f, 0.0f, 0.0f);
	
}

void displayFcn (void)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 
	glColor3f (0.0, 1.0, 0.0); // Set point color to green.
	glPointSize (1.0);  //Set point size to 1.0.

		 glPushMatrix();

		 if(buttonState == 1)
		 {
			 xRot+=(xRotLength-xRot)*0.1f;
		     yRot+=(yRotLength-yRot)*0.1f;
			 cout<< "Y   :"<<yRot<<endl;
			 cout<< "X   :"<<xRot<<endl;
		 }
		 //glTranslatef(xTrans, yTrans, zTrans);
         glRotatef(xRot, 1.0f, 0.0f, 0.0f);
         glRotatef(yRot, 0.0f, 1.0f, 0.0f); 

		 if (switchmesh==true)
			 drawtri3D();
		 else
		 {
			 glColor3f (1.0, 1.0, 1.0);
			 drawtri3Dreal();
			 drawcuttingedge();
			 drawstrokes_basering();
		 }

		 drawstrokes();

		 if (header!=NULL)
			 drawLineSegment();

		 glPopMatrix();
		 glutSwapBuffers();
}


void winReshapeFcn (GLint newWidth, GLint newHeight)
{
	/* Reset viewport and projection parameters */
	glViewport (0, 0, newWidth, newHeight);
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity ();
	//gluPerspective(45.0, (float)winWidth/winHeight, -600.000f, 1000.0);
    //gluOrtho2D (0.0, GLdouble (newWidth), 0.0, GLdouble(newHeight));
	glOrtho(-newWidth/2,newWidth/2,-newHeight/2, newHeight/2, -1200.0f, 1200.0f);

	/*Reset display window size parameters. */
	winWidth = newWidth;
	winHeight = newHeight;

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(0.0f, 0.0f, 0.0f);
}

void specialkey_func(int key, int x, int y)
{
	if(key == GLUT_KEY_UP)
	{
		zTrans += 2.0f;
	}

	if(key == GLUT_KEY_DOWN)
	{
		zTrans -= 2.0f;
	}

	if(key == GLUT_KEY_LEFT)
	{
		xTrans -= 2.0f;
	}

	if(key == GLUT_KEY_RIGHT)
	{
		xTrans += 2.0f;
	}

	glutPostRedisplay();
}


void drawLineSegment (scrPt endPt1, int x, int y)
{
	glBegin (GL_LINES);
	  glVertex2i (endPt1.x, winHeight-endPt1.y);  // previous point in movement.
	  glVertex2i (x, winHeight-y);  // current point location.
    glEnd();
}


void drawPoints (int x, int y)
{
	if (buttonState ==2)
	{

	num_all++;
	current->nextpoint = new mypoint;
	current = current->nextpoint;
	current->x=x-winWidth/2;
	current->y=winHeight-y-winHeight/2;
	current->num = num_all;
	current->index = index_all;
	current->nextpoint = NULL;

	tempstrokes = rotateOnepoint(xRot, yRot, current, -1.0);

	cout<<"X: "<< (int)current->x<<"\t"<< "Y: "<<current->y<<endl;
	cout<<"Num: "<< (int)current->num<<endl;
	cout<<"Index: "<<(int) current->index<<endl;

	glutPostRedisplay();
	}
	
		float dx, dy;
        dx = (float)(x - ox);
        dy = (float)(y - oy);

	    if (buttonState == 1) 
	    {
		    xRotLength += dy / 5.0f;
		    yRotLength += dx / 5.0f;
	    
            ox = x; oy = y;
            glutPostRedisplay();
		}
}


void drawFreeLines (GLint button, GLint action, GLint xMouse, GLint yMouse)
{
	if (button== GLUT_LEFT_BUTTON && action == GLUT_DOWN)
	{
		//endPt1.x= xMouse;
		//endPt1.y= yMouse;
		buttonState = 2;
		if (index_all==1){
			header = new mypoint;
			header->x =xMouse-winWidth/2;
			header->y =winHeight-yMouse-winHeight/2;
			header->num = 1;
			header->index = 1;
			header->nextpoint = NULL;
			
			tempstrokes = rotateOnepoint(xRot, yRot, header, -1.0);

				cout<<"X: "<< header->x<<"\t"<< "Y: "<<header->y<<endl;
	            cout<<"Num: "<< header->num<<endl;
	            cout<<"Index: "<< header->index<<endl;

				current = header;

		}
		else{
			current->nextpoint = new mypoint;
			current = current->nextpoint;
			current->x=xMouse-winWidth/2;
	        current->y=winHeight-yMouse-winHeight/2;
	        current->num = num_all;
	        current->index = index_all;
			current->nextpoint = NULL;

			tempstrokes = rotateOnepoint(xRot, yRot, current, -1.0);
	        
				cout<<"X: "<< current->x<<"\t"<< "Y: "<<current->y<<endl;
	            cout<<"Num: "<< current->num<<endl;
	            cout<<"Index: "<< current->index<<endl;

		}
		
	}
	else if (button == GLUT_LEFT_BUTTON && action == GLUT_UP)
	{
		buttonState = 0;
		num_all=1;
		index_all++;
	}
	else if (button == GLUT_RIGHT_BUTTON)
	{
		if (action == GLUT_DOWN)
			buttonState = 1;
		else if (action == GLUT_UP)
			buttonState = 0;

		ox = xMouse; oy = yMouse;
		
		glutPostRedisplay();
	}

}

void GL_DrawLine(mypoint * current2)
{
	glClear (GL_COLOR_BUFFER_BIT);
	glBegin(GL_LINE_LOOP);
    while(current2!=NULL)
     {
		glColor3f (0.0, 0.0, 1.0);  
		glVertex2f((GLfloat)current2->x ,(GLfloat)current2->y);
		current2 = current2->nextpoint ;
	 }
     glEnd();
     glFlush ();
}

void mymenu (int value) //build menu to concile the datases of stroke points.
{
	//CDT->generating object.
	if (value==1){      
		mypoint * header2 = changelength (header);
		mypoint * current2 = header2 ; 
		header = header2;
	    GL_DrawLine(current2);

		end_mypoint = Fin_end( header2, end_mypoint );
		output* headerOut = CDT(header);
		output* junctriangle = adjust_chain(headerOut);

		headerOut = reset_X(headerOut);
		junctriangle = reset_X (junctriangle);
		output * headerresult =termianl_chain(headerOut,junctriangle,header);

		spineOut = MakeSpine(headerOut , junctriangle, headerresult );

		Delete_mypoint(header);
		header = NULL;
		index_all=1;
		num_all=1;

		Delete_point3D(tempstrokes);
		tempstrokes = NULL;
		spine2p = NULL;
		spineTri = NULL;
		spineOut = NULL;
	}
	//Painting.
	if (value==2){
		pointOut = locate3Dstrokes(tri3D, xRot, yRot, header, pointOut);
		Delete_mypoint(header);
		header = NULL;
		index_all=1;
		num_all=1;

		Delete_point3D(tempstrokes);
		tempstrokes = NULL;
	}
	//Erasering.
	if (value==3){
		eraser(xRot,yRot);
		Delete_mypoint(header);
		header = NULL;
		index_all=1;
		num_all=1;
		Delete_point3D(tempstrokes);
		tempstrokes = NULL;
	}
	//Extrusion.
	if (value==4){
        triEX2tri3D(header, xRot,yRot);

		Delete_mypoint(header);
		header = NULL;
		index_all=1;
		num_all=1;
		Delete_trimesh(triEx);
		triEx = NULL;
		Delete_point3D(tempstrokes);
		tempstrokes = NULL;
		baseringP = NULL;
	}
	//Cutting.
	if (value==5){
		cutting(header, xRot, yRot);
		Delete_point3D(tempstrokes);
		tempstrokes = NULL;
	}
	//Undo.
	if (value==6){
		Delete_mypoint(header);
		header = NULL;
		index_all=1;
		num_all=1;
		Delete_point3D(tempstrokes);
		tempstrokes = NULL;
	}
	//SWITCH
	if (value==7){
		if (switchmesh==false)
			switchmesh=true;
		else
			switchmesh=false;
		if ( switchmesh==false )
		{
			//glEnable(GL_LIGHTING); 
	        //glEnable(GL_LIGHT1);							// enable light.
	        //glEnable(GL_LIGHT2);
	       // glEnable(GL_LIGHT3);
	   
		}
		else
		{
			//glDisable(GL_LIGHTING);
		}

	}


}



//=================== GL FINISH ==================


int main (int argc, char** argv)
{
	current = header;

	glutInit (&argc, argv);
	//glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB);
	//if (seeround==true)
	glutInitDisplayMode(GLUT_DEPTH  | GLUT_RGB | GLUT_DOUBLE);
	
    //glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB);
	//glutInitDisplayMode(GLUT_DEPTH  | GLUT_RGB | GLUT_DOUBLE);
	glutInitWindowPosition (100, 100);
	glutInitWindowSize (winWidth, winHeight);
	glutCreateWindow ("Free Draw");

	glLightfv(GL_LIGHT1, GL_AMBIENT, LightAmbient);				// setting enviroment light.
	glLightfv(GL_LIGHT1, GL_DIFFUSE, LightDiffuse);				// diffuse light.
	glLightfv(GL_LIGHT1, GL_POSITION,LightPosition1);			// position.
	glLightfv(GL_LIGHT2, GL_POSITION,LightPosition2);
	glLightfv(GL_LIGHT3, GL_POSITION,LightPosition3);
	glLightfv(GL_LIGHT4, GL_POSITION,LightPosition4);
	glLightfv(GL_LIGHT5, GL_POSITION,LightPosition5);
	glLightfv(GL_LIGHT6, GL_POSITION,LightPosition6);
	

	init();

	glutDisplayFunc (displayFcn);
	glutReshapeFunc (winReshapeFcn);
	glutSpecialFunc(specialkey_func);
	glutMouseFunc (drawFreeLines);
	glutMotionFunc (drawPoints);

	int id;
	  id = glutCreateMenu (mymenu);  // use gl_func to create MENU
	  glutAddMenuEntry ("CREATING",1);
	  glutAddMenuEntry ("PAINTING",2);
	  glutAddMenuEntry ("ERASER",3);
	  glutAddMenuEntry ("EXTRUSION",4);
	  glutAddMenuEntry ("CUTTING",5);
	  glutAddMenuEntry ("UNDO",6);
	  glutAddMenuEntry ("SWITCH",7);
	  glutAttachMenu(GLUT_MIDDLE_BUTTON);  

	  glEnable(GL_DEPTH_TEST);


	glutMainLoop();
	//exit(0);
}
