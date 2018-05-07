#pragma warning( disable : 4996 )
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "Vector3f.h"
#include "bitmap_image.hpp"

#define PI 3.14159265

using namespace std;

typedef struct light {
	Vector3f direction;
	Vector3f intensity;
	Vector3f spotlightPosition;
	float cutoffAngle;
	bool spotlight;
} light;

typedef struct spher {
	Vector3f center;
	float radius;
	Vector3f Ka;
	Vector3f Kd;
	Vector3f Ks;
	float shininess;
} spher;

typedef struct plane {
	Vector3f normal;
	Vector3f center;
	float width;
	float length;
	Vector3f Ka;
	Vector3f Kd;
	Vector3f Ks;
	float shininess;
	bool isFloor;
	bool isSide;
	bool isMirror;
} plane;

typedef struct scene {
	Vector3f Pc;
	Vector3f up;
	Vector3f right;
	float w;
	float h;
	int Rx;
	int Ry;
	Vector3f camera;
	Vector3f ambient;
	vector<light> lights;
	vector<spher> sphers;
	vector<plane>planes;
	int numOfLights;
	int numOfSphers;
	int numOfPlanes;
} scene;

//intersection is either with a spher or a plane
typedef struct intersection {
	float t;
	bool isSpher;
	spher spher;
	plane plane;
	bool whiteTile;
}intersection;

Vector3f constructRayThroughPixel(scene myScene, int i, int j);
float intersectWithSphere(Vector3f ray, Vector3f raySourcePoint, spher mySpher, scene myScene);
float intersectWithPlane(Vector3f ray, Vector3f raySourcePoint, plane myPlane, scene myScene);
intersection findIntersection(Vector3f ray, Vector3f raySourcePoint, scene myScene);
Vector3f calcAmbientColor(scene myScene, intersection hit);
Vector3f calcDiffuseColor(Vector3f ray, Vector3f raySourcePoint, light myLight, scene myScene, intersection hit);
Vector3f calcSpecularColor(Vector3f ray, Vector3f raySourcePoint, light myLight, scene myScene, intersection hit);
rgb_t calcReflectedColor(Vector3f ray, Vector3f raySourcePoint, scene myScene, intersection hit);
bitmap_image rayTrace(scene myScene);
bool intersectWithTriangle(Vector3f ray, Vector3f raySourcePoint, float t, scene myScene, Vector3f c1, Vector3f c2, Vector3f c3);
bool tileColor(plane myPlane, Vector3f hitPoint);

int main(int argc, const char* argv[]) {
	FILE *fd;
	char resultImage[128];
	char objectType[10];
	float vectorValues[3];
	char tmpChar = 'a';		// Default value
	scene sceneA;
	sceneA.numOfLights = 0;
	sceneA.numOfSphers = 0;
	sceneA.numOfPlanes = 0;
	light tmpLight;
	spher tmpSpher;
	plane tmpPlane;
	Vector3f tmpVector;
	bitmap_image img;

	sceneA.camera.makeZero();

	if (argc < 3) {
		cout << "Insufficient amount of parameters!" << endl;
		return 1;
	}

	fd = fopen(argv[1], "r");
	if (fd == NULL) {	// File open failed
		cout << "File open failed!" << endl;
		return 1;
	}

	strcpy(resultImage, argv[2]);	//Saving destination image name from parameters

	while (fgetc(fd) != EOF) {
		fseek(fd, -1, SEEK_CUR);	//Correction, because we use fgetc in the while condition
		fscanf(fd, "%s", objectType);
		fseek(fd, 1, SEEK_CUR);

		if (strcmp(objectType, "scene") == 0) {
			fscanf(fd, "%f", &vectorValues[0]);
			fseek(fd, 1, SEEK_CUR);
			fscanf(fd, "%f", &vectorValues[1]);
			fseek(fd, 1, SEEK_CUR);
			fscanf(fd, "%f", &vectorValues[2]);
			fseek(fd, 1, SEEK_CUR);
			// Make new vector Pc from vectorValues
			sceneA.Pc = Vector3f(vectorValues);

			fscanf(fd, "%f", &vectorValues[0]);
			fseek(fd, 1, SEEK_CUR);
			fscanf(fd, "%f", &vectorValues[1]);
			fseek(fd, 1, SEEK_CUR);
			fscanf(fd, "%f", &vectorValues[2]);
			fseek(fd, 1, SEEK_CUR);
			// Make new vector up from vectorValues
			sceneA.up = Vector3f(vectorValues);

			fscanf(fd, "%f", &sceneA.w);
			fseek(fd, 1, SEEK_CUR);
			fscanf(fd, "%d", &sceneA.Rx);
			fseek(fd, 1, SEEK_CUR);
			fscanf(fd, "%d", &sceneA.Ry);
			fseek(fd, 1, SEEK_CUR);
			fscanf(fd, "%f", &vectorValues[0]);
			fseek(fd, 1, SEEK_CUR);
			fscanf(fd, "%f", &vectorValues[1]);
			fseek(fd, 1, SEEK_CUR);
			fscanf(fd, "%f", &vectorValues[2]);
			// Make new vector ambient from vectorValues
			sceneA.ambient = Vector3f(vectorValues);
			sceneA.right = Vector3f::crossProduct(sceneA.Pc - sceneA.camera, sceneA.up); //Notice, Pc-camera gives us Vtowards vector of the camera\scene
			sceneA.h = (sceneA.Ry / sceneA.Rx)*sceneA.w;
		}

		if (strcmp(objectType, "light") == 0) {
			tmpLight.spotlight = 0;
			fscanf(fd, "%f", &vectorValues[0]);
			fseek(fd, 1, SEEK_CUR);
			fscanf(fd, "%f", &vectorValues[1]);
			fseek(fd, 1, SEEK_CUR);
			fscanf(fd, "%f", &vectorValues[2]);
			fseek(fd, 1, SEEK_CUR);
			// Make new vector lightDirection from vectorValues
			tmpLight.direction = Vector3f(vectorValues);

			fscanf(fd, "%f", &vectorValues[0]);
			fseek(fd, 1, SEEK_CUR);
			fscanf(fd, "%f", &vectorValues[1]);
			fseek(fd, 1, SEEK_CUR);
			fscanf(fd, "%f", &vectorValues[2]);
			// Make new vector lightIntensity from vectorValues
			tmpLight.intensity = Vector3f(vectorValues);
			tmpChar = fgetc(fd);
			if (tmpChar != '\n') {	// Means this is a spotlight, and we need to get more values
				tmpLight.spotlight = 1;
				fscanf(fd, "%f", &vectorValues[0]);
				fseek(fd, 1, SEEK_CUR);
				fscanf(fd, "%f", &vectorValues[1]);
				fseek(fd, 1, SEEK_CUR);
				fscanf(fd, "%f", &vectorValues[2]);
				fseek(fd, 1, SEEK_CUR);
				// Make new vector spotLightPosition from vectorValues 
				tmpLight.spotlightPosition = Vector3f(vectorValues);
				fscanf(fd, "%f", &tmpLight.cutoffAngle);

				sceneA.numOfLights++;
				sceneA.lights.push_back(tmpLight);
			}
			else {
				sceneA.numOfLights++;
				sceneA.lights.push_back(tmpLight);
				continue;	// Already read the \n delimiter, skipping the following fgetc
			}
		}

		if (strcmp(objectType, "spher") == 0) {
			fscanf(fd, "%f", &vectorValues[0]);
			fseek(fd, 1, SEEK_CUR);
			fscanf(fd, "%f", &vectorValues[1]);
			fseek(fd, 1, SEEK_CUR);
			fscanf(fd, "%f", &vectorValues[2]);
			fseek(fd, 1, SEEK_CUR);
			// Make new vector spherCenter from vectorValues
			tmpSpher.center = Vector3f(vectorValues);

			fscanf(fd, "%f", &tmpSpher.radius);
			fseek(fd, 1, SEEK_CUR);

			fscanf(fd, "%f", &vectorValues[0]);
			fseek(fd, 1, SEEK_CUR);
			fscanf(fd, "%f", &vectorValues[1]);
			fseek(fd, 1, SEEK_CUR);
			fscanf(fd, "%f", &vectorValues[2]);
			fseek(fd, 1, SEEK_CUR);
			// Make new vector Ka from vectorValues}  
			tmpSpher.Ka = Vector3f(vectorValues);

			fscanf(fd, "%f", &vectorValues[0]);
			fseek(fd, 1, SEEK_CUR);
			fscanf(fd, "%f", &vectorValues[1]);
			fseek(fd, 1, SEEK_CUR);
			fscanf(fd, "%f", &vectorValues[2]);
			fseek(fd, 1, SEEK_CUR);
			// Make new vector Kd from vectorValues}  
			tmpSpher.Kd = Vector3f(vectorValues);

			fscanf(fd, "%f", &vectorValues[0]);
			fseek(fd, 1, SEEK_CUR);
			fscanf(fd, "%f", &vectorValues[1]);
			fseek(fd, 1, SEEK_CUR);
			fscanf(fd, "%f", &vectorValues[2]);
			fseek(fd, 1, SEEK_CUR);
			// Make new vector Ks from vectorValues}  
			tmpSpher.Ks = Vector3f(vectorValues);

			fscanf(fd, "%f", &tmpSpher.shininess); // float ??

			sceneA.numOfSphers++;
			sceneA.sphers.push_back(tmpSpher);

		}

		if (strcmp(objectType, "plane") == 0) {
			tmpPlane.isMirror = 0;
			tmpPlane.isFloor = 0;
			tmpPlane.isSide = 0;
			fscanf(fd, "%f", &vectorValues[0]);
			fseek(fd, 1, SEEK_CUR);
			fscanf(fd, "%f", &vectorValues[1]);
			fseek(fd, 1, SEEK_CUR);
			fscanf(fd, "%f", &vectorValues[2]);
			fseek(fd, 1, SEEK_CUR);
			// Make new vector planeNormal from vectorValues}
			tmpPlane.normal = Vector3f(vectorValues);

			fscanf(fd, "%f", &vectorValues[0]);
			fseek(fd, 1, SEEK_CUR);
			fscanf(fd, "%f", &vectorValues[1]);
			fseek(fd, 1, SEEK_CUR);
			fscanf(fd, "%f", &vectorValues[2]);
			fseek(fd, 1, SEEK_CUR);
			// Make new vector planeCenter from vectorValues}  
			tmpPlane.center = Vector3f(vectorValues);

			fscanf(fd, "%f", &tmpPlane.width);   
			fseek(fd, 1, SEEK_CUR);
			fscanf(fd, "%f", &tmpPlane.length);
			fseek(fd, 1, SEEK_CUR);

			fscanf(fd, "%f", &vectorValues[0]);
			fseek(fd, 1, SEEK_CUR);
			fscanf(fd, "%f", &vectorValues[1]);
			fseek(fd, 1, SEEK_CUR);
			fscanf(fd, "%f", &vectorValues[2]);
			fseek(fd, 1, SEEK_CUR);
			// Make new vector Ka from vectorValues}  
			tmpPlane.Ka = Vector3f(vectorValues);

			fscanf(fd, "%f", &vectorValues[0]);
			fseek(fd, 1, SEEK_CUR);
			fscanf(fd, "%f", &vectorValues[1]);
			fseek(fd, 1, SEEK_CUR);
			fscanf(fd, "%f", &vectorValues[2]);
			fseek(fd, 1, SEEK_CUR);
			// Make new vector Kd from vectorValues}  
			tmpPlane.Kd = Vector3f(vectorValues);

			fscanf(fd, "%f", &vectorValues[0]);
			fseek(fd, 1, SEEK_CUR);
			fscanf(fd, "%f", &vectorValues[1]);
			fseek(fd, 1, SEEK_CUR);
			fscanf(fd, "%f", &vectorValues[2]);
			fseek(fd, 1, SEEK_CUR);
			// Make new vector Ks from vectorValues}  
			tmpPlane.Ks = Vector3f(vectorValues);

			fscanf(fd, "%f", &tmpPlane.shininess); 

			if (tmpPlane.normal.p[0] == 0 && tmpPlane.normal.p[2] == 0) {
				tmpPlane.isFloor = 1;
			}
			if (tmpPlane.normal.p[1] == 0 && tmpPlane.normal.p[2] == 0) {
				tmpPlane.isSide = 1;
			}

			tmpChar = fgetc(fd);
			if (tmpChar != '\n') {	// Check if this is a Mirror plane
				tmpChar = fgetc(fd);
				if (tmpChar == 'M') {
					tmpPlane.isMirror = 1;
				}
				else {
					fseek(fd, -1, SEEK_CUR);
				}
			}

			sceneA.numOfPlanes++;
			sceneA.planes.push_back(tmpPlane);
			fseek(fd, -1, SEEK_CUR);
		}
		tmpChar = fgetc(fd);	// Reading the \n delimiter
	}

	img = rayTrace(sceneA);
	img.save_image(resultImage);
	return 0;
}

Vector3f constructRayThroughPixel(scene myScene, int i, int j) {
	// Determinig pixel size according to scene width
	float pixelWidthRatio = myScene.w / myScene.Rx;
	float pixelHeightRatio = myScene.h / myScene.Ry;
	Vector3f right = Vector3f(myScene.right);
	Vector3f up = Vector3f(myScene.up);
	Vector3f ray = Vector3f(myScene.Pc);

	//Normalizing right and up vectors
	right.normalize();
	up.normalize();
	//Normalizing i to be between (-Rx/2,Rx/2) and j to be between (-Ry/2,Ry/2)
	i = i - floorf(myScene.Rx / 2);
	j = j - floorf(myScene.Ry / 2);
	ray = ray + right*i*pixelWidthRatio + up*j*pixelHeightRatio;	//Calculating the point p on the scene grid
	ray = ray - myScene.camera;		//Creating the ray from the point found
	ray.normalize();
	return ray;
}

float intersectWithSphere(Vector3f ray, Vector3f raySourcePoint, spher mySpher, scene myScene) {
	float tm, th, d_squared;
	Vector3f L, V;
	float t1, t2;

	L = (mySpher.center - raySourcePoint);
	V = ray;
	V.normalize();

	tm = Vector3f::dotProduct(L, V);
	d_squared = powf(L.getLength(), 2.0) - powf(tm, 2.0);

	if (d_squared > (powf(mySpher.radius, 2.0)))
		return -1.0;

	th = sqrtf((powf(mySpher.radius, 2.0)) - d_squared);
	t1 = tm - th;
	t2 = tm + th;
	if (t1 < 0) {
		if (t2 > 0) {
			return t2;
		}
		else
			return -1.0;
	}
	else
		return t1;
}

float intersectWithPlane(Vector3f ray, Vector3f raySourcePoint, plane myPlane, scene myScene) {
	float t;
	Vector3f dir1, dir2, V, intersectionPoint;
	Vector3f c1, c2, c3, c4;

	myPlane.normal.getTwoOrthogonals(dir1, dir2);

	V = ray;
	V.normalize();

	t = Vector3f::dotProduct(myPlane.normal, ((myPlane.center - raySourcePoint) / (Vector3f::dotProduct(myPlane.normal, V))));
	if (t < 0)
		return t;

	intersectionPoint = raySourcePoint + t*ray;
	intersectionPoint -= myPlane.center;

	if (myPlane.isFloor) {
		if (fabsf(intersectionPoint.p[0]) <= (myPlane.width / 2) && fabsf(intersectionPoint.p[2]) <= (myPlane.length / 2)) {
			return t;
		}
	}
	else {
		if (myPlane.isSide) {
			if (fabsf(intersectionPoint.p[2]) <= (myPlane.width / 2) && fabsf(intersectionPoint.p[1]) <= (myPlane.length / 2)) {
				return t;
			}
		}
		else {
			if (fabsf(intersectionPoint.p[0]) <= (myPlane.width / 2) && fabsf(intersectionPoint.p[1]) <= (myPlane.length / 2)) {
				return t;
			}
		}
	}

	return INFINITY;
}


bool intersectWithTriangle(Vector3f ray, Vector3f raySourcePoint, float t, scene myScene, Vector3f c1, Vector3f c2, Vector3f c3) {
	Vector3f v1, v2, N, V, P;
	P = raySourcePoint + t*ray;
	V = P - raySourcePoint;

	v1 = c1 - raySourcePoint;
	v2 = c2 - raySourcePoint;
	N = Vector3f::crossProduct(v1, v2);
	N.normalize();

	if (Vector3f::dotProduct(V, N) < 0)
		return false;

	v1 = c2 - raySourcePoint;
	v2 = c3 - raySourcePoint;
	N = Vector3f::crossProduct(v1, v2);
	N.normalize();
	
	if (Vector3f::dotProduct(V, N) < 0)
		return false;

	v1 = c3 - raySourcePoint;
	v2 = c1 - raySourcePoint;
	N = Vector3f::crossProduct(v1, v2);
	N.normalize();
	
	if (Vector3f::dotProduct(V, N) < 0)
		return false;
	
	return true;
}

intersection findIntersection(Vector3f ray, Vector3f raySourcePoint, scene myScene) {
	float min_t = INFINITY;		//Also used as default value for non intersection with any primitive for a given ray
	float t;
	spher min_spher;
	plane min_plane;
	bool isMinSpher = 0;
	bool isWhiteTile = 0;
	intersection ans;
	Vector3f newRaySourcePoint;


	float stepSize = (myScene.w / myScene.Rx)*(float)0.95;
	newRaySourcePoint = raySourcePoint + ray*stepSize;

	for (int i = 0; i < myScene.numOfSphers; i++) {
		t = intersectWithSphere(ray, newRaySourcePoint, myScene.sphers.at(i), myScene);
		if (t > 0 && t < min_t) {
			min_spher = myScene.sphers.at(i);
			min_t = t;
			isMinSpher = 1;
		}
	}

	for (int i = 0; i < myScene.numOfPlanes; i++) {
		t = intersectWithPlane(ray, newRaySourcePoint, myScene.planes.at(i), myScene);
		if (t > 0 && t < min_t) {
			min_plane = myScene.planes.at(i);
			min_t = t;
			if (!min_plane.isMirror) {
				isWhiteTile = tileColor(myScene.planes.at(i), ray*t);
			}
			isMinSpher = 0;
		}
	}
	ans.isSpher = isMinSpher;
	ans.plane = min_plane;
	ans.spher = min_spher;
	ans.t = min_t;
	ans.whiteTile = isWhiteTile;
	return ans;
}

bool tileColor(plane myPlane, Vector3f hitPoint) {
	Vector3f dir1, dir2, planeStart_hit, planeStartPoint;
	float i, j;

	myPlane.normal.getTwoOrthogonals(dir1, dir2);

	if (myPlane.isFloor) {
		planeStartPoint = myPlane.center - (myPlane.width / 2)*dir1 - (myPlane.length / 2)*dir2;
		planeStart_hit = hitPoint - planeStartPoint;
		i = Vector3f::projectOntoVector(planeStart_hit, dir1).getLength();
		j = Vector3f::projectOntoVector(planeStart_hit, dir2).getLength();
	}
	else {
		planeStartPoint = myPlane.center - (myPlane.width / 2)*dir2 - (myPlane.length / 2)*dir1;
		planeStart_hit = hitPoint - planeStartPoint;
		i = Vector3f::projectOntoVector(planeStart_hit, dir2).getLength();
		j = Vector3f::projectOntoVector(planeStart_hit, dir1).getLength();
	}

	i = floorf(i / myPlane.width * 32);
	j = floorf(j / myPlane.length * 32);
	return ((int)(i + j) % 2 == 1);
}

rgb_t calcColor(Vector3f ray, Vector3f raySourcePoint, scene myScene, intersection hit) {
	rgb_t color;
	Vector3f colorVector, L;
	intersection block;
	colorVector.makeZero();
	bool shadow;

	if (!hit.isSpher && hit.plane.isMirror) {
		return calcReflectedColor(ray, raySourcePoint, myScene, hit);
	}
	colorVector = (calcAmbientColor(myScene, hit));
	for (int i = 0; i < myScene.numOfLights; i++) {
		if (!myScene.lights.at(i).spotlight) { //directional light 
			L = -1 * myScene.lights.at(i).direction;
			L.normalize();

			//shot ray from hit point thowrds light source and make sure
			//there aren't any blockink objects
			block = findIntersection(L, raySourcePoint + ray*hit.t, myScene);

			//if (block.t <= 0.001 || block.t == INFINITY)
			if (block.t <= 0 || block.t == INFINITY)
				shadow = false;
			else
				shadow = true;
		}
		else { //spotlight 
			L = myScene.lights.at(i).spotlightPosition - (raySourcePoint + ray*hit.t);
			L.normalize();

			//shot ray from hit point thowrds light source and make sure
			//there aren't any blockink objects
			block = findIntersection(L, raySourcePoint + ray*hit.t, myScene);

			//check that it don't hit anything 
			//or that the blocking object is farther than the spotlight position

			//			if (block.t == INFINITY || block.t <= 0.001 || ((L*block.t + (raySourcePoint + ray*hit.t)).getLength() > (myScene.lights.at(i).spotlightPosition - (raySourcePoint + ray*hit.t)).getLength()))
			if (block.t == INFINITY || block.t <= 0 || ((L*block.t + (raySourcePoint + ray*hit.t)).getLength() > (myScene.lights.at(i).spotlightPosition - (raySourcePoint + ray*hit.t)).getLength()))
				shadow = false;
			else
				shadow = true;
		}

		if (!shadow) {
			colorVector += (calcDiffuseColor(ray, raySourcePoint, myScene.lights.at(i), myScene, hit));
			colorVector += (calcSpecularColor(ray, raySourcePoint, myScene.lights.at(i), myScene, hit));
		}
	}
	if (colorVector.p[0] > 1.0)
		colorVector.p[0] = 1.0;
	if (colorVector.p[1] > 1.0)
		colorVector.p[1] = 1.0;
	if (colorVector.p[2] > 1.0)
		colorVector.p[2] = 1.0;

	color.red = colorVector.p[0] * 255;
	color.green = colorVector.p[1] * 255;
	color.blue = colorVector.p[2] * 255;

	return color;
}

Vector3f calcAmbientColor(scene myScene, intersection hit) {
	Vector3f colorVec;

	colorVec.makeZero();
	if (hit.isSpher)
		colorVec = hit.spher.Ka*myScene.ambient;
	else {// with plane
		if (hit.whiteTile)
			colorVec = Vector3f(1, 1, 1)*myScene.ambient;
		else
			colorVec = Vector3f(0, 0, 0)*myScene.ambient;
	}
	if (colorVec.p[0] < 0 || colorVec.p[1] < 0 || colorVec.p[2] < 0)
		colorVec.makeZero();

	return colorVec;

}

Vector3f calcDiffuseColor(Vector3f ray, Vector3f raySourcePoint, light myLight, scene myScene, intersection hit) {
	Vector3f colorVec;
	Vector3f D;
	Vector3f N;
	Vector3f L;
	Vector3f IL;
	Vector3f kd;
	float cosTheta;

	colorVec.makeZero();

	if (hit.isSpher) {
		N = (raySourcePoint + ray*hit.t) - hit.spher.center;
		kd = hit.spher.Kd;
	}

	else {//hit with plane
		N = hit.plane.normal;
		kd = hit.plane.Kd;
	}
	N.normalize();

	if (!myLight.spotlight) { //directional light 
		L = -1 * myLight.direction;
		L.normalize();
		IL = myLight.intensity*Vector3f::dotProduct(L, L);/////??????
	}
	else { //spotlight 
		L = myLight.spotlightPosition - (ray*hit.t + raySourcePoint);
		L.normalize();

		D = myLight.direction;
		D.normalize();
		cosTheta = Vector3f::dotProduct(D, (-1 * L));
		if ((acos(cosTheta) * (float)180.0 / PI) <= myLight.cutoffAngle)
			IL = myLight.intensity*cosTheta;
		else
			IL.makeZero();
	}

	cosTheta = Vector3f::dotProduct(N, L);

	colorVec = kd*cosTheta*IL;

	if (colorVec.p[0] < 0 || colorVec.p[1] < 0 || colorVec.p[2] < 0)
		colorVec.makeZero();
	return colorVec;
}

Vector3f calcSpecularColor(Vector3f ray, Vector3f raySourcePoint, light myLight, scene myScene, intersection hit) {
	Vector3f colorVec;
	Vector3f N, L, V, R, D;
	Vector3f IL;
	Vector3f ks;
	float cosTheta, n;

	colorVec.makeZero();
	IL.makeZero();
	V = -1 * ray;

	if (hit.isSpher) {
		N = (ray*hit.t + raySourcePoint) - hit.spher.center;
		ks = hit.spher.Ks;
		n = hit.spher.shininess;
	}

	else {//hit with plane
		N = hit.plane.normal;
		ks = hit.plane.Ks;
		n = hit.plane.shininess;
	}
	N.normalize();

	if (!myLight.spotlight) { //directional light 
		L = -1 * myLight.direction;//from hit point thordws light source
		L.normalize();
		IL = myLight.intensity*Vector3f::dotProduct(L, L);//????????
		L.makeNegative();
	}
	else { //spotlight 
		L = (ray*hit.t + raySourcePoint) - myLight.spotlightPosition;
		L.normalize();
		D = myLight.direction;
		D.normalize();
		cosTheta = Vector3f::dotProduct(D, L);

		if ((acos(cosTheta) * (float)180.0 / PI) <= myLight.cutoffAngle)
			IL = myLight.intensity*cosTheta;
		else
			IL.makeZero();

	}

	R = L - 2 * N* Vector3f::dotProduct(N, L);
	R.normalize();
	V.normalize();
	cosTheta = Vector3f::dotProduct(R, V);

	cosTheta = max(cosTheta, (float)0.0);
	colorVec = ks*powf(cosTheta, n)*IL;

	if (colorVec.p[0] < 0 || colorVec.p[1] < 0 || colorVec.p[2] < 0)
		colorVec.makeZero();

	return colorVec;
}


rgb_t calcReflectedColor(Vector3f ray, Vector3f raySourcePoint, scene myScene, intersection hit) {
	Vector3f reflectedRay, N, newSourcePoint;
	intersection newIntersectoin;
	rgb_t color;
	color.red = 0;
	color.green = 0;
	color.blue = 0;

	N = hit.plane.normal;
	N.normalize();
	newSourcePoint = raySourcePoint + ray*hit.t;
	reflectedRay = ray - 2 * N* Vector3f::dotProduct(N, ray);
	reflectedRay.normalize();
	newIntersectoin = findIntersection(reflectedRay, newSourcePoint, myScene);
	
	if (hit.t > 0 && hit.t != INFINITY) {
		color = calcColor(reflectedRay, newSourcePoint, myScene, newIntersectoin);
	}
	return color;
}

bitmap_image rayTrace(scene myScene) {
	Vector3f ray;
	bitmap_image img(myScene.Rx, myScene.Ry);
	intersection hit;

	for (int y = 0; y < myScene.Ry; y++) {
		for (int x = 0; x < myScene.Rx; x++) {
			ray = constructRayThroughPixel(myScene, x, y);
			hit = findIntersection(ray, myScene.camera, myScene);

			if (hit.t > 0 && hit.t != INFINITY) {
				img.set_pixel(x, (myScene.Ry - 1) - y, calcColor(ray, myScene.camera, myScene, hit));
			}
		}
	}

	return img;
}