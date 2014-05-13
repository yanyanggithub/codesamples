/*!
 *  @file inc.h
 *  @brief the necessary includes and definition of static const variables
 *  @author Yan Yang
 *  Created on: 09/05/2014
 *
 */

#ifndef INC_H_
#define INC_H_
#include <iostream>
#include <vector>
#include <math.h>
#include <unistd.h>

#ifdef __APPLE__
#include <GLUT/glut.h>          /* Open GL Util    APPLE */
#else
#include <GL/glut.h>            /* Open GL Util    OpenGL*/
#endif
using namespace std;

static const double pi = 3.14159265;

#endif /* INC_H_ */
