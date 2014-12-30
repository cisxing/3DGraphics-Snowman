class BezierCurve : public Freeform
{
 static double bernstein(int i, int n, double t) { 
	 if(n == 1) { if(i == 0) return 1-t; if(i == 1) return t; return 0; } if(i < 0 || i > n) return 0; return (1 - t) * bernstein(i, n-1, t) + t * bernstein(i-1, n-1, t); 
 }


 float2 getPoint(float t) { 
	 float2 r(0.0, 0.0); // for every control point // compute weight using the Bernstein formula // add control point to r, weighted return r
	  }
};

int main(int argc, char *argv[]) {
    glutInit(&argc, argv);                 		// GLUT initialization
    glutInitWindowSize(500, 500);				// Initial resolution of the MsWindows Window is 600x600 pixels 
    glutInitWindowPosition(100, 100);            // Initial location of the MsWindows window
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);    // Image = 8 bit R,G,B + double buffer + depth buffer
 
    glutCreateWindow("Circle");        	// Window is born
 									// Register event handlers
	 glutDisplayFunc(onDisplay);       
	 glutMouseFunc(mouseClicks);
	 glutMainLoop();  
    return 0;
}