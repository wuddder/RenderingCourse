#include "stdafx.h"
#include "cameras/realistic.h"
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <cstdlib>
#include "cmath"
#include "core/sampler.h"
#include "../core/sampler.h"
#include "../shapes/sphere.h"
#include "../core/montecarlo.h"

using namespace std;

RealisticCamera::RealisticCamera(const AnimatedTransform &cam2world,
				 float hither, float yon, 
				 float sopen, float sclose, 
				 float filmdistance, float aperture_diameter, string specfile, 
				 float filmdiag, Film *f)
	: Camera(cam2world, sopen, sclose, f) // pbrt-v2 doesnot specify hither and yon
{	
	_hither=hither;
	_yon=yon;
	_filmdistance=filmdistance;
	_aperture_diameter=aperture_diameter;
	_filmdiag=filmdiag;
	_distance=ReadSpecFile(specfile);
	// compute transform: raster to camera, film to camera
	float diag = sqrtf(f->xResolution * f->xResolution + f->yResolution * f->yResolution);
	float scale = filmdiag / diag;
	float X = scale * f->xResolution * 0.5f;
	float Y = scale * f->yResolution * 0.5f;
	RasterToCamera = Translate(Vector(0.f, 0.f, -(_filmdistance + _distance))) * 
					 Translate(Vector(X, -Y, 0.f)) *
					 Scale(-scale, scale, 1) ;

	
}
float RealisticCamera::ReadSpecFile(const string filename){
	float radius=0;
	float axpos=0;
	float N=0;
	float aperture=0;
	float distance=0;
	ifstream spec(filename);
	string str;
	while(getline(spec,str)){
		if(str[0]=='#')
			continue;
		else if(str[0]!='#'){
			Len lens;
			sscanf(str.c_str(),"%f%f%f%f",&radius,&axpos,&N,&aperture);
			lens.isStop = (N==0);
			lens.radius = radius;
			lens.z = -distance;
			lens.n = (N==0) ? 1 : N;
			lens.aperture = aperture/2;
			distance += axpos;
			lensgroup.push_back(lens);
		}		
	
	}
	spec.close();
	return distance;
}
float RealisticCamera::GenerateRay(const CameraSample &sample, Ray *ray) const {
	// Generate raster and camera samples
    Point Pras(sample.imageX, sample.imageY, 0);
    Point Pcamera;
    RasterToCamera(Pras, &Pcamera);
    //setup iteration
    int current=lensgroup.size()-1;
    float lensU,lensV,lensZ;
    // Method 1: pbrt-v2 provide
	//ConcentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);
	//lensU *= lensgroup[current].aperture;
	//lensV *= lensgroup[current].aperture;
	//lensZ = -_distance;

	// Method 2: textbook, r = \sqrt{lensU}, theta = 2 \pi lensV
	float r = sqrtf(sample.lensU), theta = 2 * M_PI * sample.lensV;
	lensU = lensgroup[current].aperture * r * cosf(theta);
	lensV = lensgroup[current].aperture * r * sinf(theta);
	lensZ = -_distance;


	Point hit = Point(lensU, lensV, lensZ);
	Vector T = hit - Pcamera;

	// Tracing the ray	
	for(int i = current; i >= 0; --i){
		if (lensgroup[i].isStop) {
			float deltaZ = lensgroup[i].z - hit.z;
			T = Normalize(T);
			float t = deltaZ / T.z;
			hit = hit + t * T;
			if (hit.x * hit.x + hit.y * hit.y > lensgroup[i].aperture * lensgroup[i].aperture) return 0.f;
		}
		else {
			float n2 = (i == 0) ? 1 :lensgroup[i-1].n;
			if (!Refraction(hit, T, lensgroup[i], n2)) return 0.f;
		}
	}
	ray->o = hit;
	ray->d = Normalize(T);
	ray->mint = _hither;
	ray->maxt = (_yon - _hither) / ray->d.z;
	ray->time = Lerp(sample.time, shutterOpen, shutterClose);
	float cosTheta = Dot(Normalize(ray->o - Pcamera), Vector(0, 0, 1));
	float Z = _filmdistance;
	float weight = (lensgroup[current].aperture *lensgroup[current].aperture * M_PI) / ( Z * Z);
	weight = weight * cosTheta * cosTheta * cosTheta * cosTheta;
	CameraToWorld(*ray, ray);
	ray->d = Normalize(ray->d);
	return weight;

}
bool RealisticCamera::Refraction(Point& O, Vector& D, Len surface, float n2) const{
	/*
	 I: incident vector I (normalize), 
	 C: incident circle center
	OC + I * t = OP
	| OP | = radius
	| OC + I* t | = | OP |
	| OC + I * t | = radius
	| OC | ^ 2 + 2 * ( OC • I ) * t + | I * t | ^ 2 - radius * radius = 0
	get t value by solving the equation above
	*/
	Vector I = Normalize(D);
	Point C = Point(0.f, 0.f, surface.z - surface.radius);
	Vector OC = O - C;
	float b = Dot(OC, I);
	float c = OC.LengthSquared() - surface.radius * surface.radius;
	float determine = b * b - c;
	float t = 0;
	if (determine < 0) {
		return false;
	}
	else {
		float root = sqrtf(determine);
		t = (surface.radius > 0) ? (-b + root) : (-b - root);
	}
	O = O + t * I;
	if (surface.aperture * surface.aperture < O.y * O.y + O.x * O.x) return false;
	Vector N = (surface.radius > 0.f) ? Normalize(C - O) : Normalize(O - C);

	// Heckber's Method
	float n_ratio = surface.n / n2;
	float c1 = -Dot(I, N);
	float c2 = 1.f - n_ratio * n_ratio * (1.f - c1 * c1);
	if (c2 <= 0.f) return false;
	else c2 = sqrtf(c2);
	D = n_ratio * I + (n_ratio * c1 - c2) * N;

	return true;
}
RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film) {
	// Extract common camera parameters from \use{ParamSet}
	float hither = params.FindOneFloat("hither", -1);
	float yon = params.FindOneFloat("yon", -1);
	float shutteropen = params.FindOneFloat("shutteropen", -1);
	float shutterclose = params.FindOneFloat("shutterclose", -1);

	// Realistic camera-specific parameters
	string specfile = params.FindOneString("specfile", "");
	float filmdistance = params.FindOneFloat("filmdistance", 70.0); // about 70 mm default to film
 	float fstop = params.FindOneFloat("aperture_diameter", 1.0);	
	float filmdiag = params.FindOneFloat("filmdiag", 35.0);

	Assert(hither != -1 && yon != -1 && shutteropen != -1 &&
		shutterclose != -1 && filmdistance!= -1);
	if (specfile == "") {
	    Severe( "No lens spec file supplied!\n" );
	}
	return new RealisticCamera(cam2world, hither, yon,
				   shutteropen, shutterclose, filmdistance, fstop, 
				   specfile, filmdiag, film);
}