/*!
 *  @file kalman.h
 *  @brief 1d and 2d kalman filter, including basic matrix functions
 *
 *  @author Yan Yang
 *
 *  References:
 *  [1] Greg Welch and Gary Bishop (2006), "An Introduction to the Kalman
 *   Filter". http://www.cs.unc.edu/~welch/kalman/kalmanIntro.html
 *  [2] Kalman, R. E. (1960). "A New Approach to Linear Filtering and
 *   Prediction Problems".Journal of Basic Engineering 82 (1)
 *  [3] en.wikipedia.org/wiki/Kalman_filter
 *
 *  Created on: 08/05/2014
 */


#ifndef KALMAN_H_
#define KALMAN_H_
#include "inc.h"

class vec2
/*!
 * @brief 1x2 vector [x, y]
 */
{
public:
	double x ,y;
	vec2 (double x = 0.0, double y = 0.0) : x(x), y(y) {};
	vec2 operator +(vec2 v) { return vec2(x + v.x, y + v.y); }
	vec2 operator -(vec2 v) { return vec2(x - v.x, y - v.y); }
};

class mat2
/*!
 * @brief 2x2 matrix [a11, a12; a21, a22]
 *
 */
{
public:
	double a11, a12, a21, a22;
public:
	mat2 (double a11 = 0.0, double a12 = 0.0, double a21 = 0.0,	double a22 = 0.0) : a11(a11), a12(a12), a21(a21), a22(a22) {};
	mat2 t() { return mat2(a11, a21, a12, a22); }
	mat2 operator +(mat2 m) { return mat2(a11 + m.a11, a12 + m.a12, a21 + m.a21, a22 + m.a22); }
	mat2 operator -(mat2 m) { return mat2(a11 - m.a11, a12 - m.a12, a21 - m.a21, a22 - m.a22); }
	mat2 operator *(mat2 m) { return mat2(a11*m.a11 + a12*m.a21, a11*m.a12 + a12*m.a22, a21*m.a11 + a22*m.a21, a21*m.a12 + a22*m.a22); }
	vec2 operator *(vec2 v) { return vec2(a11*v.x + a12*v.y, a21*v.x + a22*v.y); }
	mat2 inv()
	{
		double d = 1.0 / (a11*a22 - a12*a21);
		return mat2(d*a22, -d*a12, -d*a21, d*a11);
	}
};


void printmat2(mat2 &x)
{
	cout << "mat2: [" << x.a11 << ", " << x.a12 << "; " << x.a21 << ", "
			<< x.a22 << "]" << endl;
}

void printvec2(vec2 &x)
{
	cout << "vec2: [" << x.x << ", " << x.y << "]" << endl;
}

class kalman_filter
/*!
 * @brief kalman filter for a linear dynamic system
 */
{
private:
	/*1d variables*/
	double Q; // process variance
	double R; // estimate of measurement variance
	double P; // posteri error estimate
	double Pminus; // priori error estimate
	double xhat; // postori estimate of x
	double xhatminus; // priori estimate of x
	double K; // kalman gain or blending factor

	/*2d variables*/
	vec2 statePre; // x'(k)
	vec2 statePost; // x(k)
	vec2 measurement;

	mat2 processNoiseCov; //Q
	mat2 measurementMat; //H
	mat2 measurementNoiseCov; // R
	mat2 errorCovPre;  // P'(k)
	mat2 errorCovPost; // P(k)
	mat2 transitionMatrix; // A
	mat2 gain; // K(k)

public:
	kalman_filter ();
	void init (double var);
	void init (double x, double y);

	double predict(double x);
	void update();

	vec2 predict();
	vec2 correct(double x, double y);
};

kalman_filter::kalman_filter()
{
	Q = 1e-5;
	R = 0.0;
	P = 1.0;
	Pminus = 0.0;
	xhat = 0.0;
	xhatminus = 0.0;
	K = 0.0;

	statePre = vec2();
	statePost = vec2();

	measurementMat = mat2(1.0, 0.0, 0.0, 1.0);
	transitionMatrix = mat2(1.0, 0.0, 0.0, 1.0);
	processNoiseCov = mat2(1.0, 1e-4, 1e-4, 1.0);
	measurementNoiseCov = mat2(1.0, 1e-1, 1e-1 , 1.0);
	errorCovPre = mat2(1.0, 1e-1, 1e-1, 1.0);
	errorCovPost = mat2(1.0, 1e-1, 1e-1, 1.0);
}

void kalman_filter::init (double var)
/*!
 * @brief kalman filter with known variance
 **/
{
	R = var;
}


void kalman_filter::init (double x, double y)
/*!
 * @brief kalman filter with 2d measurement
 */
{
	statePre = vec2(x, y);
}

void kalman_filter::update()
/*!
 * @brief update the state status
 */
{
	//time update
	xhatminus = xhat;
	Pminus = P;

	//update the kalman gain
	K = Pminus / (Pminus + R);
}

double kalman_filter::predict(double x)
/*!
 * @brief predict the correct value according to the current measurement
 */
{
	update();
	xhat = xhatminus + K * (x - xhatminus);
	P = (1 - K) * Pminus;
	return xhat;
}

vec2 kalman_filter::predict()
{
	//update the state x'(k) = A*x(k)
	statePre = transitionMatrix * statePost;
	mat2 temp1;

	//update error covariance matrix temp1 = A*P(k)
	temp1 = transitionMatrix * errorCovPost;

	//P'(k) = temp1*A.t() + Q.t()
	errorCovPre =  temp1 * transitionMatrix.t() + processNoiseCov.t();

	statePost = statePre;
	errorCovPost = errorCovPre;

	return statePre;
}

vec2 kalman_filter::correct(double x, double y)
{
	mat2 temp1, temp2, temp3;
	vec2 temp4;

	vec2 measurement(x, y);

	//temp1 = H * P'(k)
	temp1 = measurementMat * errorCovPre;

	//temp2 = temp1*H.t() + R.t()
	temp2 = temp1 * measurementMat.t() + measurementNoiseCov.t();

	//Kt(k) = temp2.inv() * temp1
	temp3 = temp2.inv() * temp1;

	//K(t)
	gain = temp3.t();

	//temp4 = z(k) - H*x'(k)
	temp4 = measurement - measurementMat * statePre;

	//x(k) = x'(k) + K(k)*temp4
	statePost = statePre + gain * temp4;

	//P(k) = P'(k) - K(k)*temp1
	errorCovPost = errorCovPre - gain * temp1;

	return statePost;
}

#endif /* KALMAN_H_ */
