/*!
 *  @file main.cpp
 *  @brief Kalman filter demo in C++ using OpenGL GUI
 *  @author Yan Yang
 *
 */

#include "kalman.h"
#include "inc.h"

kalman_filter mykf;
int error_rate;

#define trajectory 0

#if trajectory
	int pointsize = 5;
	int linewidth = 2;
	int crosssize = 5;
#else
	int pointsize = 15;
	int linewidth = 5;
	int crosssize = 10;
#endif

double rand_noise()
/*!
 * @brief generate uniformly distributed random noise in the range [-1, 1]
 * */
{
	return 2*((rand()/(double)RAND_MAX) - 0.5);
}

void test()
/*!
 * @brief a simple example of kalman filtering for a 1D dynamic system
 * */
{
	int n_iter = 50;
	double x = 0.5; //real value

	srand(0);
	kalman_filter kf = kalman_filter();
	kf.init(0.1);
	for(int k = 0; k < n_iter; k++)
	{
		// adding random noise to the real value as measurement
		double measurement = x + rand_noise() * 0.1;
		double t = kf.predict(measurement);
		cout << measurement  << " " << t << endl;
	}
}

void handleKeypress(unsigned char key, int x, int y) {
	switch (key) {
		case 27: //Escape key
			exit(0);
	}
}

void initRendering() {
	glClearColor(0,0,0,0);
	glClear(GL_COLOR_BUFFER_BIT);
	glViewport( 0,0, 500, 500 );
	glMatrixMode( GL_PROJECTION );
	glOrtho( 0.0, 500.0, 0.0, 500.0, 1.0, -1.0 );
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void drawPoints(int x, int y, string colour)
{
	glPointSize(pointsize);
    glBegin( GL_POINTS );
    if(colour == "RED")
    	glColor3f( 1.0,0.0,0.0 );
    else if(colour == "BLUE")
    	glColor3f( 0.0,0.0,1.0 );
    else if(colour == "GREEN")
    	glColor3f( 0.0,1.0,0.0 );
    glVertex2f(x, 500-y);
    glEnd();
    glutSwapBuffers();
#if !trajectory
    glutPostRedisplay();
#endif
}

void drawCross(int x, int y, string colour)
{
	if(colour == "RED")
		glColor3f( 1.0,0.0,0.0 );
	else if(colour == "BLUE")
		glColor3f( 0.0,0.0,1.0 );
	else if(colour == "GREEN")
		glColor3f( 0.0,1.0,0.0 );

	glLineWidth(linewidth);
	glBegin( GL_LINES );
	glVertex2f(x-crosssize, 500-y-crosssize);
	glVertex2f(x+crosssize, 500-y+crosssize);
	glEnd();

	glBegin( GL_LINES );
	glVertex2f(x+crosssize, 500-y-crosssize);
	glVertex2f(x-crosssize, 500-y+crosssize);
	glEnd();
	glutSwapBuffers();
#if !trajectory
	glutPostRedisplay();
#endif
}
void motionPassive(int x, int y)
{
	cout << "Mouse moved at "
		<< "(" << x << "," << y << ")" << endl;
	drawPoints(x, y, "RED");
	usleep(5e+3);
	mykf.predict();

	// adding noise to the measurement
	x = x + rand_noise() * error_rate;
	y = y + rand_noise() * error_rate;
	vec2 myestimated = mykf.correct(x, y);
	drawCross(myestimated.x, myestimated.y, "GREEN");
	cout << "my estimation: "
			<< myestimated.x << " " << myestimated.y << endl;
}

void display()
{
	glClear(GL_COLOR_BUFFER_BIT);
}

void opengltest()
{
	mykf.init(0.0, 0.0);
	glutInitDisplayMode(GLUT_SINGLE|GLUT_RGB);
	glutInitWindowSize(500, 500); //Set the window size
	//Create the window
	glutCreateWindow("Graphic interface for 2D kalman filter tracking");
	glutKeyboardFunc(handleKeypress);
	glutDisplayFunc(display);
	glutPassiveMotionFunc(motionPassive);
	initRendering(); //Initialize rendering
	glutMainLoop();
}

int main (int argc, char** argv)
/* the error_rate introduce measurement error
 * to the mouse position
 * e.g. error_rate = 50, the error of measurement is
 *  in the range of [-50, 50] pixels
 *
 * alternative option to run test();
 * see how Kalman filter works for one-dimensional data
 */
{
	error_rate = 50;
	glutInit(&argc, argv);
	opengltest();
	return 0;
}
