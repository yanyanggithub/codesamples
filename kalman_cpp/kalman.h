/*!
 *  @file kalman.h
 *  @brief 1d and 2d kalman filter, including basic matrix functions
 *
 *  @author Yan Yang
 *
 *  reference:
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
	double x;
	double y;

public:
	vec2 ();
	vec2 (double x, double y);
	vec2 operator = (const vec2 &v);
};

vec2::vec2 ()
{
	this->x = 0.0;
	this->y = 0.0;
}

vec2::vec2 (double x, double y)
{
	this->x = x;
	this->y = y;
}

vec2 vec2::operator = (const vec2 &v)
{
	this->x = v.x;
	this->y = v.y;
	return *this;
}

class mat2
/*!
 * @brief 2x2 matrix [a11, a12; a21, a22]
 *
 */
{
public:
	double a11, a12, a21, a22;
public:
	mat2 ();
	mat2 (double a11, double a12, double a21, double a22);
	mat2 operator = (const mat2 &m);
};

mat2::mat2()
{
	this->a11 = 0.0;
	this->a12 = 0.0;
	this->a21 = 0.0;
	this->a22 = 0.0;
}

mat2::mat2(double a11, double a12, double a21, double a22)
{
	this->a11 = a11;
	this->a12 = a12;
	this->a21 = a21;
	this->a22 = a22;
}

mat2 mat2::operator = (const mat2 &m)
{
	this->a11 = m.a11;
	this->a12 = m.a12;
	this->a21 = m.a21;
	this->a22 = m.a22;
	return *this;
}

void printmat2(mat2 &x)
{
	cout << "mat2: [" << x.a11 << ", " << x.a12 << "; " << x.a21 << ", "
			<< x.a22 << "]" << endl;
}

void printvec2(vec2 &x)
{
	cout << "vec2: [" << x.x << ", " << x.y << "]" << endl;
}

vec2 add(vec2 &a, vec2 &b)
{
	vec2 result;
	result.x = a.x + b.x;
	result.y = a.y + b.y;
	return result;
}

vec2 subtract(vec2 &a, vec2 &b)
{
	vec2 result;
	result.x = a.x - b.x;
	result.y = a.y - b.y;
	return result;
}

mat2 transpose(mat2 &a)
{
	mat2 result;
	result.a11 = a.a11;
	result.a22 = a.a11;
	result.a12 = a.a21;
	result.a21 = a.a12;
	return result;
}

mat2 add(mat2 &a, mat2 &b)
{
	mat2 result;
	result.a11 = a.a11 + b.a11;
	result.a12 = a.a12 + b.a12;
	result.a21 = a.a21 + b.a21;
	result.a22 = a.a22 + b.a22;
	return result;
}

mat2 addT(mat2 &a, mat2 &b)
{
	mat2 result, temp;
	temp = transpose(b);
	result = add(a, temp);
	return result;
}

mat2 subtract(mat2 &a, mat2 &b)
{
	mat2 result;
	result.a11 = a.a11 - b.a11;
	result.a12 = a.a12 - b.a12;
	result.a21 = a.a21 - b.a21;
	result.a22 = a.a22 - b.a22;
	return result;
}

mat2 multiply(mat2 &a, mat2 &b)
{
	mat2 result;
	result.a11 = a.a11*b.a11 + a.a12*b.a21;
	result.a12 = a.a11*b.a12 + a.a12*b.a22;
	result.a21 = a.a21*b.a11 + a.a22*b.a21;
	result.a22 = a.a21*b.a12 + a.a22*b.a22;
	return result;
}

mat2 multiplyT(mat2 &a, mat2 &b)
{
	mat2 result, temp;
	temp = transpose(b);
	result = multiply(a, temp);
	return result;
}

vec2 multiply(mat2 &a, vec2 &b)
{
	vec2 result;
	result.x = a.a11*b.x + a.a12*b.y;
	result.y = a.a21*b.x + a.a22*b.y;
	return result;
}

mat2 invert(mat2 &a)
{
	mat2 result;
	double d = 1.0 / (a.a11*a.a22 - a.a12*a.a21);
	result.a11 = d * a.a22;
	result.a12 = -d * a.a12;
	result.a21 = -d * a.a21;
	result.a22 = d * a.a11;
	return result;
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
	void init (double R);
	void init (double x, double y);

	double predict(double x);
	void update();

	vec2 predict();
	vec2 correct(double x, double y);
};

kalman_filter::kalman_filter()
{
	this->Q = 1e-5;
	this->R = 0.0;
	this->P = 1.0;
	this->Pminus = 0.0;
	this->xhat = 0.0;
	this->xhatminus = 0.0;
	this->K = 0.0;

	this->statePre = vec2();
	this->statePost = vec2();

	this->measurementMat = mat2(1.0, 0.0, 0.0, 1.0);
	this->transitionMatrix = mat2(1.0, 0.0, 0.0, 1.0);
	this->processNoiseCov = mat2(1.0, 1e-4, 1e-4, 1.0);
	this->measurementNoiseCov = mat2(1.0, 1e-1, 1e-1 , 1.0);
	this->errorCovPre = mat2(1.0, 1e-1, 1e-1, 1.0);
	this->errorCovPost = mat2(1.0, 1e-1, 1e-1, 1.0);
}

void kalman_filter::init (double R)
/*!
 * @brief kalman filter with known variance
 **/
{
	this->R = R;
}


void kalman_filter::init (double x, double y)
/*!
 * @brief kalman filter with 2d measurement
 */
{
	this->statePre = vec2(x, y);
}

void kalman_filter::update()
/*!
 * @brief update the state status
 */
{
	//time update
	this->xhatminus = this->xhat;
	this->Pminus = this->P;

	//update the kalman gain
	this->K = this->Pminus / (this->Pminus + R);
}

double kalman_filter::predict(double x)
/*!
 * @brief predict the correct value according to the current measurement
 */
{
	update();
	this->xhat = this->xhatminus + this->K * (x - this->xhatminus);
	this->P = (1 - this->K) * this->Pminus;
	return this->xhat;
}

vec2 kalman_filter::predict()
{
	//update the state x'(k) = A*x(k)
	this->statePre = multiply(this->transitionMatrix, this->statePost);
	mat2 temp1, temp2;

	//update error covariance matrix temp1 = A*P(k)
	temp1 = multiply(this->transitionMatrix, this->errorCovPost);

	//P'(k) = temp1*A.t() + Q.t()
	temp2 = multiplyT(temp1, this->transitionMatrix);
	this->errorCovPre = addT(temp2, this->processNoiseCov);

	this->statePost = this->statePre;
	this->errorCovPost = this->errorCovPre;

	return this->statePre;
}

vec2 kalman_filter::correct(double x, double y)
{
	mat2 temp1, temp2, temp3, temp4, temp5;
	vec2 temp6, temp7, temp8;
	mat2 temp9;

	vec2 measurement(x, y);

	//temp1 = H * P'(k)
	temp1 = multiply(this->measurementMat, this->errorCovPre);

	//temp3 = temp1*H.t() + R.t()
	temp2 = multiplyT(temp1, this->measurementMat);
	temp3 = addT(temp2, this->measurementNoiseCov);

	//temp5 = inv(temp3)*temp1
	temp4 = invert(temp3);
	temp5 = multiply(temp4, temp1);

	//K(t)
	this->gain = transpose(temp5);

	//temp7 = z(k) - H*x'(k)
	temp6 = multiply(this->measurementMat, this->statePre);
	temp7 = subtract(measurement, temp6);

	//x(k) = x'(k) + K(k)*temp7
	temp8 = multiply(this->gain, temp7);
	this->statePost = add(this->statePre, temp8);

	//P(k) = P'(k) - K(k)*temp1
	temp9 = multiply(this->gain, temp1);
	this->errorCovPost = subtract(this->errorCovPre, temp9);

	return this->statePost;
}

#endif /* KALMAN_H_ */
