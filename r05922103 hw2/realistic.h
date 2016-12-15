#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CAMERAS_REALISTIC_H
#define PBRT_CAMERAS_REALISTIC_H

// #include "../core/camera.h"
#include <stdlib.h>
#include "camera.h"
#include "paramset.h"
#include "film.h"
#include "vector"
#include <cstdlib>
using namespace std;

struct Len {
	bool isStop;
	float radius, axpos, n, aperture;
	bool end;
	float z;
};

// RealisticCamera Declarations
class RealisticCamera : public Camera {
public:
	// RealisticCamera Public Methods
	RealisticCamera(const AnimatedTransform &cam2world,
						float hither, float yon, float sopen,
						float sclose, float filmdistance, float aperture_diameter, string specfile,
						float filmdiag, Film *film);
	float GenerateRay(const CameraSample &sample, Ray *) const;
  
private:
	// RealisticCamera Private Member
	float _hither, _yon, _aperture_diameter, _filmdistance, _filmdiag;
	float _distance;
	vector<Len> lensgroup;
	Transform RasterToCamera;
	// RealisticCamera Private Methods
	float ReadSpecFile(string fileName);
	bool Refraction(Point& O, Vector& D, Len surface, float n2) const;
};


RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film);


#endif	// PBRT_CAMERAS_REALISTIC_H