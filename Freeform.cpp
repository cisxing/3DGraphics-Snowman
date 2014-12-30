
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include<string>
#include <vector>
using std::vector;

 
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__)
// Needed on MsWindows
#include <windows.h>
#endif // Win32 platform
 
#include <GL/gl.h>
#include <GL/glu.h>
// Download glut from: http://www.opengl.org/resources/libraries/glut/
#include <GL/glut.h>

class float2
{
public:
	float x;
	float y;

	float2()
	{
		x = 0.0f;
		y = 0.0f;
	}


	float2(float x, float y):x(x),y(y){}

	float2 operator-() const
	{
		return float2(-x, -y);
	}


	float2 operator+(const float2& addOperand) const
	{
		return float2(x + addOperand.x, y + addOperand.y);
	}

	float2 operator-(const float2& operand) const
	{
		return float2(x - operand.x, y - operand.y);
	}

	float2 operator*(const float2& operand) const
	{
		return float2(x * operand.x, y * operand.y);
	}
	
	float2 operator*(float operand) const
	{
		return float2(x * operand, y * operand);
	}

	void operator-=(const float2& a)
	{
		x -= a.x;
		y -= a.y;
	}

	void operator+=(const float2& a)
	{
		x += a.x;
		y += a.y;
	}

	void operator*=(const float2& a)
	{
		x *= a.x;
		y *= a.y;
	}

	void operator*=(float a)
	{
		x *= a;
		y *= a;
	}

	float norm()
	{
		return sqrtf(x*x+y*y);
	}

	float norm2()
	{
		return x*x+y*y;
	}

	float2 normalize()
	{
		float oneOverLength = 1.0f / norm();
		x *= oneOverLength;
		y *= oneOverLength;
		return *this;
	}

};

class Object
{
protected:
    float2 scaleFactor;
    float2 position;
    float orientation;
public:
	Object* translate(float2 offset){
        position += offset; return this;
    }
    Object* scale(float2 factor){
        scaleFactor *= factor; return this;
    }
    Object* rotate(float angle){
        orientation += angle; return this;
    }
    Object():orientation(0.0f), scaleFactor(1.0,1.0){}
    virtual ~Object(){}
    virtual void draw()
    {
        // apply scaling, translation and orientation
        glPushMatrix();
        glTranslatef(position.x, position.y, 0.0f);
        glRotatef(orientation, 0.0f, 0.0f, 1.0f);
        glScalef(scaleFactor.x, scaleFactor.y, 1.0f);
        drawModel();
        glPopMatrix();
    }
    virtual void drawModel()=0;
    virtual void move(double t, double dt)=0;
};

class Asteroid : public Object
{
    void drawModel()
    {
        glColor3d(1.0, 0.5, 0.0);
        glBegin(GL_TRIANGLE_FAN);
             glVertex2d(0.0, 0.0);
             for(double phi=0.0; phi<6.5; phi+=0.32)
               glVertex2d( 
                   cos(phi)*(0.8 + 0.2*sin(phi*7) ),
                   sin(phi)*(0.8 + 0.2*sin(phi*7) ));
        glEnd();
    }
	void move(double t, double dt)
	{
		translate(float2(dt,dt));
		
		scale(float2(0.5,0.5));
		rotate(dt);	
	}
};

class Scene
{
    std::vector<Object*> objects;
public:
    void addObject(Object* object) {
        objects.push_back(object);
    }
    ~Scene() {
        for(unsigned int i=0; i<objects.size(); i++)
            delete objects.at(i); 
    }
    void draw() {
        for(unsigned int i=0; i<objects.size(); i++)
            objects.at(i)->draw();
    }
    void move(double t, double dt) {
        for(unsigned int i=0; i<objects.size(); i++)
            objects.at(i)->move(t, dt);
    }
};

Scene scene;

void onDisplay()
{
	scene.draw();
	//have to call this method or otherwise it will not work
	glutSwapBuffers();

}


void onIdle()
{
	scene.move(0,0.001); 
	glutPostRedisplay();

}


int main(int argc, char *argv[])
{
	//scene.addObject(new Asteroid);
	/*for(int i=0; i<20;i++)
	{
		scene.addObject((new Asteroid())
        ->scale(float2(0.1f+i/20, 0.1f+i/20))
        ->rotate(45.0f+i*5)
        ->translate(float2(0.1+i*0.1, 0.2+i*0.1))
		);
	}*/
	scene.addObject((new Asteroid())
        ->scale(float2(0.1f, 0.1f))
        ->rotate(45.0f)
        ->translate(float2(0.1, 0.2))
		);

	glutInit(&argc, argv);                 		// GLUT initialization
    glutInitWindowSize(700, 700);				// Initial resolution of the MsWindows Window is 600x600 pixels 
    glutInitWindowPosition(200, 0);            // Initial location of the MsWindows window
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);    // Image = 8 bit R,G,B + double buffer + depth buffer
 
    glutCreateWindow("Asteroid");        	// Window is born
 									// Register event handlers
	 glutDisplayFunc(onDisplay);
	 glutIdleFunc(onIdle);
	 //glutMouseFunc(mouseClicks);
	 //glutKeyboardFunc(keyPressed);
	 //glutKeyboardUpFunc(keyReleased);
	//glutMotionFunc(onMouseDrag);
	 glutMainLoop();  
    return 0;


}