
/**
x^T A x
Quadratic surface: 4*4 matrix. f(x,y,z) = 0 implicit equation Inside the surface, positive, outside the surface, negative
so 4*4: x,y,z,1. Why do we need heterogenius? Think about a sphere: x^2+y^2+z^2 -1 = 0
In such case: {1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,-1} we need the constant term that is why we added the last row
What about ellipsoid: Still a diagonal matrix but not with constants
what about cylindar: in z direction it is invariant so we should set that value to 0: {1,0,0,0},{0,1,0,0},{0,0,0,0},{0,0,0,-1}
Just need to define the matrix and then should be easy. 
What is the coefficient of the matrix:
Ex. Sphere: two roots: corresponding to the entry point and exit point. (Choose the closest point for shading?)

How to make a plane:
Find the reference point that we find the line that is perpendicular to the plane. It is just the dot product and exactly the implicit function.

Diffuse surfaces:
ex. diffuse sphere: what is a diffuse material: No reflection.
It receives the most light when the direction of the light source is in the perpendicular direction as the given point 
So we calculate the dot product max(0,N (dot) L) to see where is the illumination point, multiply the diffusion coefficient depending on the material
Iterate all the light sources, then multiply the function above

What is the inverse of the scaling matrix? {s_1,0,0,0},{0,s_2,0,0},{0,0,s_3,0},{0,0,0,1} then it is {1/s_1,0,0,0},{0,1/s_2,0,0},{0,0,1/s_3,0},{0,0,0,1}
how about rotation? if we rotate in z axis: {cos x,sin x,0,0},{-sin x,cos x,0,0},{0,0,1,0},{0,0,0,1} then it is {cos (-x),sin (-x),0,0},{-sin (-x),cos (-x),0,0},{0,0,1,0},{0,0,0,1}
											which is just the transpose of the matrix because sin(-x) = -sin(x) and cos(-x) = cos(x)
how about transformation:{1,0,0,0},{0,1,0,0},{0,0,0,0},{dealta x,delta y,delta z,1} to get x'= x+ (delta x) and so on. 
											Inverse would be {1,0,0,0},{0,1,0,0},{0,0,0,0},{-dealta x,-delta y,-delta z,1}


if we have a shiny surface, we want to find the normal vector and then
need to calculate the half vector then we can measure the variation H (dot) N, {max(0, H dot N)}^n the shininess if really shiny, then fone highlight is small(the n)
The half vector is the viewing direction(?) when the H is in the same direction as N , get the maximum intensity
To figure out the triangle intersection 



**/

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
 
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__)
// Needed on MsWindows
#define NOMINMAX
#include <windows.h>
#endif // Win32 platform
 
#include <GL/gl.h>
#include <GL/glu.h>
// Download glut from: http://www.opengl.org/resources/libraries/glut/
#include <GL/glut.h>

#include "float2.h"
#include "float3.h"
#include "float4.h"
#include "float4x4.h"
#include <vector>
#include <algorithm>
#include "perlin.h"
#include <string>














// For procedural texturing -> what does this mean?
Perlin perlin;
// for quadrics, so that we do not need a float4.cpp
const float4x4 float4x4::identity(
	1, 0, 0, 0,
	0, 1, 0, 0,
	0, 0, 1, 0,
	0, 0, 0, 1);

// Abstract base class for light sources
class LightSource
{
public:
	virtual float3 getPowerDensityAt(float3 x)=0;
	virtual float3 getLightDirAt(float3 x)=0;
	virtual float  getDistanceFrom(float3 x)=0;

};

// TO BE CREATED AT PRACTICAL
class DirectionalLightSource : public LightSource
{
	float3 direction;
	float3 density;

public:
	DirectionalLightSource(float3 dir, float3 den)
	{
		direction = dir;
		density = den;

	}
	float3 getPowerDensityAt(float3 x)
	{
		return density;
	}
	float3 getLightDirAt(float3 x){
		return direction.normalize();
	}
	float getDistanceFrom(float3 x){
		return FLT_MAX;

	}

};

// TO BE CREATED AT PRACTICAL
class PointLightSource : public LightSource
{
protected:
	float3 power;
	float3 density;

public:
	PointLightSource(float3 center, float3 den)
	{
		power = center;
		density = den;

	}
	float3 getPowerDensityAt(float3 x)
	{
		float densityDistance = 4*3.14159*getDistanceFrom(x)*getDistanceFrom(x);
		return density*float3(1/densityDistance,1/densityDistance,
			1/densityDistance);
	}

	//dont really know why we make the direction pointing from the point to the lightsource
	float3 getLightDirAt(float3 x){
		float3 difference = power-(x);
		return difference.normalize();
	}
	float getDistanceFrom(float3 x){
		
		float3 difference = power-(x);
		return difference.norm();
	}
};






















// Skeletal Material class. Feel free to add methods e.g. for illumination computation (shading).
class Material
{
public:
	bool reflective;
	bool refractive;
	bool textured;
	float3 minReflectance;		// Fresnel coefficient
	float refractiveIndex;			// index of refraction
	float3 kd;			// diffuse reflection coefficient
	float3 ks;			// specular reflection coefficient
	float shininess;	// specular exponent
	Material()
	{
		reflective = false;
		refractive = false;
		textured = false;
		minReflectance = float3(0.93, 0.85, 0.4);
		refractiveIndex = 1;
		kd = float3(0.5, 0.5, 0.5) + kd * 0.5;
		ks = float3(1, 1, 1);
		shininess = 15;
	}

	float3 reflect(float3 inDir, float3 normal) {
		return inDir - normal * normal.dot(inDir) * 2;
	}
	float3 refract(float3 inDir, float3 normal) {
		float ri = refractiveIndex;
		float cosa = -normal.dot(inDir);
		if(cosa < 0) { cosa = -cosa; normal = -normal; ri = 1 / ri; }
		float disc = 1 - (1 - cosa * cosa) / ri / ri;
		if(disc < 0) return reflect(inDir, normal);
		return inDir * (1.0 / ri) + normal * (cosa / ri - sqrt(disc));
	}
	float3 getReflectance(float3 inDir, float3 normal) {
		float cosa = fabs(normal.dot(inDir));
		return minReflectance + (float3(1, 1, 1) - minReflectance) * pow(1 - cosa, 5);
	}

	virtual float3 shade(
		float3 position,
		float3 normal,
		float3 viewDir,
		float3 lightDir,
		float3 lightPowerDensity)
		{
			// TO BE IMPLEMENTED AT PRACTICAL
			//this is for diffuse only, copied from the class ppt
			float cosTheta = normal.dot(lightDir);
			if(cosTheta<0) return float3(0,0,0);
			float3 halfway = (viewDir +lightDir).normalize();
			float cosDelta= normal.dot(halfway);
			float3 diffuse = kd*lightPowerDensity*cosTheta;

			if(cosDelta<0) return diffuse;
			return diffuse+lightPowerDensity *ks *pow(cosDelta,shininess);

		}
};

// Skeletal Camera class. Feel free to add custom initialization, set aspect ratio to fit viewport dimensions, or animation.
class Camera
{
	float3 eye;

	float3 lookAt;
	float3 right;
	float3 up;

public:
	float3 getEye()
	{
		return eye;
	}
	Camera()
	{
		eye = float3(0, 0, 3);
		lookAt = float3(0, 0, 2);
		right = float3(1, 0, 0);
		up = float3(0, 1, 0);
	}

	float3 rayDirFromNdc(const float2 ndc) {
		return (lookAt - eye
			+ right * ndc.x
			+ up    * ndc.y
			).normalize();
	}
};

// Ray structure.
class Ray
{
public:
    float3 origin;
    float3 dir;
    Ray(float3 o, float3 d)
    {
        origin = o;
        dir = d;
    }
};

// Hit record structure. Contains all data that describes a ray-object intersection point.
class Hit
{
public:
	Hit()
	{
		t = -1;
	}
	float t;
	float3 position;
	float3 normal;
	Material* material;
};

// Object abstract base class.
class Intersectable
{
protected:
	Material* material;
public:
	Intersectable(Material* material):material(material) {}
    virtual Hit intersect(const Ray& ray)=0;
	virtual void rotate(float3 axis, float angle)=0;
    
    virtual void scale(float3 scaleFactor) =0;
    
    virtual void translate(float3 translation)=0;

	virtual	void translateB(float3 translation)=0;
	virtual void rotateB(float3 axis, float angle)=0;
};

// Object realization.
class Sphere : public Intersectable
{
	float3 center;
	float radius;
public:
    Sphere(const float3& center, float radius, Material* material):
		Intersectable(material),
		center(center),
		radius(radius)
    {
    }
    Hit intersect(const Ray& ray)
    {
        float3 diff = ray.origin - center;
        double a = ray.dir.dot(ray.dir);
        double b = diff.dot(ray.dir) * 2.0;
        double c = diff.dot(diff) - radius * radius;
 
        double discr = b * b - 4.0 * a * c;
        if ( discr < 0 ) 
            return Hit();
        double sqrt_discr = sqrt( discr );
        double t1 = (-b + sqrt_discr)/2.0/a;
        double t2 = (-b - sqrt_discr)/2.0/a;
 
		float t = (t1<t2)?t1:t2;
		if(t < 0)
			t = (t1<t2)?t2:t1;
		if (t < 0)
            return Hit();

		Hit h;
		h.t = t;
		h.material = material;
		h.position = ray.origin + ray.dir * t;
		h.normal = h.position - center;
		h.normal.normalize();

		return h;

    }
	void rotate(float3 axis, float angle) {

    }
    
    void scale(float3 scaleFactor) {
    }
    
    void translate(float3 translation) {
    }

	void translateB(float3 translation){
		
	}
	void rotateB(float3 axis, float angle){
		
	}
};































// TO BE CREATED AT PRACTICAL
class Plane : public Intersectable 
 {
 public:
	 float3 normal;
	 float3 x0;
	 
	 Plane(float3 planeNormal, float3 position, Material* material):
	 Intersectable(material),normal(planeNormal),x0(position)
	 {

	 }

	 Hit intersect(const Ray& ray){
		 float3 rayd= ray.dir;
		 float denom = rayd.dot(normal);
		 if(denom ==0)
		 {
			 return Hit();
		 }
		 else
		 {
			 Hit hit;
			 //do not know if this is correct
			 float d = (x0-ray.origin).dot(normal)/denom;
			 hit.t = d;
			 hit.material = material;
			 hit.normal= normal;
			 hit.normal.normalize();
			 hit.position = ray.dir*d+ray.origin;
			 return hit;
		 }
	 }
	 void rotate(float3 axis, float angle) {

    }
    
    void scale(float3 scaleFactor) {
    }
    
    void translate(float3 translation) {
    }

	void translateB(float3 translation){
		
	}
	void rotateB(float3 axis, float angle){
		
	}

};

// TO BE CREATED AT PRACTICAL
class Quadric : public Intersectable
{
	float4x4 A;

public:
	float radius;
	float3 center;
  Quadric(Material* material,float radius1, float3 position):
    radius(radius1),center(position),Intersectable(material)
  { // ellipsoid hardwired here
    // you should add methods
    // to set up different quadrics


    A = float4x4::identity;
	//this is the original stuff
	A._11 = 2;
	A._33 = -radius;
	A._30 = center.x;
	A._31 = center.y;
	A._32 = center.z;


  }


	Hit intersect(const Ray& ray) {
		// ray in homo coords
		float4 e = float4(ray.origin.x,
		ray.origin.y, ray.origin.z, 1);
		float4 d = float4(ray.dir.x,
		ray.dir.y, ray.dir.z, 0);
		// quadratic coeffs.
		double a = d.dot( A * d );
		double b = e.dot( A * d ) 
           + d.dot( A * e );
		double c = e.dot( A * e );
 
        double discr = b * b - 4.0 * a * c;
        if ( discr < 0 ) 
            return Hit();
        double sqrt_discr = sqrt( discr );
        double t1 = (-b + sqrt_discr)/2.0/a;
        double t2 = (-b - sqrt_discr)/2.0/a;
 
		float t = (t1<t2)?t1:t2;
		if(t < 0)
			t = (t1<t2)?t2:t1;
		if (t < 0)
            return Hit();

		Hit h;
		h.t = t;
		h.material = material;
		h.position = ray.origin + ray.dir * t;

		// homo position
		float4 hPos = float4(h.position.x,
		h.position.y, h.position.z, 1);
		// homo normal per quadric normal formula
		float4 hNormal = A * hPos +  hPos * A;
		// Cartesian normal
		h.normal = float3(hNormal.x, hNormal.y, hNormal.z).normalize();
		h.normal.normalize();
		return h;
 }

	void rotate(float3 axis, float angle) {
        A = float4x4::rotation(axis, -angle) * A * float4x4::rotation(axis, angle);
    }
    
    void scale(float3 scaleFactor) {
        scaleFactor = float3(1 / scaleFactor.x, 1 / scaleFactor.y, 1 / scaleFactor.z);
        A = float4x4::scaling(scaleFactor) * A * float4x4::scaling(scaleFactor);
    }
    
    void translate(float3 translation) {
        A = float4x4::translation(translation) * A * float4x4::translation(translation).transpose();
    }

	void translateB(float3 translation){
		
	}
	void rotateB(float3 axis, float angle){
		
	}

};

class ClippedQuadric : public Intersectable
{
public:
  float4x4 A;
  float4x4 B;
public:
  ClippedQuadric(Material* material, int type)
    :Intersectable(material) {
    // infinite cylinder hardwired
    A = float4x4::identity;
    A._11 = 0;
    A._33 = -0.25;
    // sphere or radius 2 hardwired
    B = float4x4::identity;
    B._33 = -0.35;
	if(type == 0)
	{
		getCone();
	}
  } // add methods to change quadric

  void getCone()
  {
	A = float4x4::identity;
	A._00=36;
	A._11=-1;
	A._22=36;
	A._33=0;
	//A._33= 0;
    // sphere or radius 2 hardwired
    B = float4x4::identity;
	//B._13=5;
    B._33 =-3;
  }

	void rotate(float3 axis, float angle) {
	  //both the quadric and the sphere rotate in the same direction
        A = float4x4::rotation(axis, -angle) * A * float4x4::rotation(axis, angle);
        B = float4x4::rotation(axis, -angle) * B * float4x4::rotation(axis, angle);
    }
    
    void scale(float3 scaleFactor) {
        scaleFactor = float3(1 / scaleFactor.x, 1 / scaleFactor.y, 1 / scaleFactor.z);
        A = float4x4::scaling(scaleFactor) * A * float4x4::scaling(scaleFactor);
        B = float4x4::scaling(scaleFactor) * B * float4x4::scaling(scaleFactor);
    }
    
    void translate(float3 translation) {
        A = float4x4::translation(translation) * A * float4x4::translation(translation).transpose();
        B = float4x4::translation(translation) * B * float4x4::translation(translation).transpose();
    }

	void translateB(float3 translation){
		 B = float4x4::translation(translation) * B * float4x4::translation(translation).transpose();
	}
	void rotateB(float3 axis, float angle){
		B = float4x4::rotation(axis, -angle) * B * float4x4::rotation(axis, angle);
	}


  Hit intersect(const Ray& ray)
  {
	float4 e = float4(ray.origin.x,
		ray.origin.y, ray.origin.z, 1);
		float4 d = float4(ray.dir.x,
		ray.dir.y, ray.dir.z, 0);
		// quadratic coeffs.
		double a = d.dot( A * d );
		double b = e.dot( A * d ) 
           + d.dot( A * e );
		double c = e.dot( A * e );
 
        double discr = b * b - 4.0 * a * c;
        if ( discr < 0 ) 
            return Hit();
        double sqrt_discr = sqrt( discr );
        double t1 = (-b + sqrt_discr)/2.0/a;
        double t2 = (-b - sqrt_discr)/2.0/a;
 
		  // get t1, t2 for quadric A

		float4 hit1 = e + d * t1;
		if(hit1.dot(B * hit1) > 0) // if not in B
			t1 = -1;				 // invalidate
			float4 hit2 = e + d * t2;
		if(hit2.dot(B * hit2) > 0) // if not in B
			t2 = -1; 				 // invalidate
		float t = (t1<t2)?t1:t2;
		
		if(t < 0)
			t = (t1<t2)?t2:t1;
		if (t < 0)
            return Hit();

		Hit h;
		h.t = t;
		h.material = material;
		h.position = ray.origin + ray.dir * t;

		// homo position
		float4 hPos = float4(h.position.x,
		h.position.y, h.position.z, 1);
		// homo normal per quadric normal formula
		float4 hNormal = A * hPos +  hPos * A;
		// Cartesian normal
		h.normal = float3(hNormal.x, hNormal.y, hNormal.z).normalize();
		h.normal.normalize();
		return h;

  }


};

class HoleQuadric : public Intersectable
{
  float4x4 A;
  float4x4 B;
public:
  HoleQuadric(Material* material)
    :Intersectable(material) {
    // infinite cylinder hardwired
		//A is an ellipsoid
	A = float4x4::identity;
	B = float4x4::identity;
  	B._22 = 0;
    B._33 = -0.6;
	//B is cylindar
    
	A._11=2;
    A._33 = -2.44;


  } // add methods to change quadric



	void rotate(float3 axis, float angle) {
	  //both the quadric and the sphere rotate in the same direction
        A = float4x4::rotation(axis, -angle) * A * float4x4::rotation(axis, angle);
        B = float4x4::rotation(axis, -angle) * B * float4x4::rotation(axis, angle);
    }
    
    void scale(float3 scaleFactor) {
        scaleFactor = float3(1 / scaleFactor.x, 1 / scaleFactor.y, 1 / scaleFactor.z);
        A = float4x4::scaling(scaleFactor) * A * float4x4::scaling(scaleFactor);
        B = float4x4::scaling(scaleFactor) * B * float4x4::scaling(scaleFactor);
    }
    
    void translate(float3 translation) {
        A = float4x4::translation(translation) * A * float4x4::translation(translation).transpose();
        B = float4x4::translation(translation) * B * float4x4::translation(translation).transpose();
    }

	void translateB(float3 translation){
		 B = float4x4::translation(translation) * B * float4x4::translation(translation).transpose();
	}
	void rotateB(float3 axis, float angle){
		B = float4x4::rotation(axis, -angle) * B * float4x4::rotation(axis, angle);
	}


  Hit intersect(const Ray& ray)
  {
	float4 e = float4(ray.origin.x,
		ray.origin.y, ray.origin.z, 1);
		float4 d = float4(ray.dir.x,
		ray.dir.y, ray.dir.z, 0);
		// quadratic coeffs.
		double a = d.dot( A * d );
		double b = e.dot( A * d ) 
           + d.dot( A * e );
		double c = e.dot( A * e );
 
        double discr = b * b - 4.0 * a * c;
        if ( discr < 0 ) 
            return Hit();
        double sqrt_discr = sqrt( discr );
        double t1 = (-b + sqrt_discr)/2.0/a;
        double t2 = (-b - sqrt_discr)/2.0/a;
 
		  // get t1, t2 for quadric A

		float4 hit1 = e + d * t1;
		if(hit1.dot(B * hit1) < 0) // if in B
			t1 = -1;				 // invalidate
			float4 hit2 = e + d * t2;
		if(hit2.dot(B * hit2) < 0) // if in B
			t2 = -1; 				 // invalidate
		float t = (t1<t2)?t1:t2;
		
		if(t < 0)
			t = (t1<t2)?t2:t1;
		if (t < 0)
            return Hit();

		Hit h;
		h.t = t;
		h.material = material;
		h.position = ray.origin + ray.dir * t;

		// homo position
		float4 hPos = float4(h.position.x,
		h.position.y, h.position.z, 1);
		// homo normal per quadric normal formula
		float4 hNormal = A * hPos +  hPos * A;
		// Cartesian normal
		h.normal = float3(hNormal.x, hNormal.y, hNormal.z).normalize();
		h.normal.normalize();
		return h;

  }


};



class Scene
{
	Camera camera;
	std::vector<LightSource*> lightSources;
	std::vector<Intersectable*> objects;
	std::vector<Material*> materials;
public:
	Scene()
	{
		// ADD LIGHT SOURCES HERE
		//direction and density
		lightSources.push_back(new DirectionalLightSource(float3(1,2.5,1),float3(1,1,1)));
		//center and density
		//lightSources.push_back(new PointLightSource(float3(3,5,7),float3(30,30,30)));
		// ADD MATERIALS HERE
		materials.push_back(new Material());
		materials[0]->kd= float3(1,1,1);
		materials[0]->shininess = 20;
		
		materials.push_back(new Material());
		materials[1]->kd= float3(0,.25,.5);
		materials[1]->ks = float3(0,.25,.5);
		materials[1]->reflective= true;
		materials[1]->minReflectance=float3(1,1,1);
		//index of refraction of Copper
		materials[1]->refractiveIndex= 0.4609;
		//materials[1] properties

		materials.push_back(new Material());
		materials[2]->kd= float3(1,0.5,0.0);
		materials[2]->ks = float3(1,.5,0.0);
		materials[2]->minReflectance=float3(0.1,0.05,0);

		materials.push_back(new Material());
		materials[3]->kd= float3(1,0.0,0.0);
		materials[3]->ks = float3(1,0.0,0.0);
		materials[3]->minReflectance=float3(0.1,0,0);

		// ADD OBJECTS HERE
		//working code
		
		objects.push_back(new Plane(float3(0,1,0),float3(0,-2.2,0),materials[0]));
		
		//objects.push_back(new Quadric(materials[0],1.2,float3(0,3,0)));
		objects.push_back(new Quadric(materials[0],0.8,float3(0,-2,0)));
		objects.push_back(new Sphere(float3(-0.7,1.8,1),0.3, materials[3]));
		objects.push_back(new ClippedQuadric(materials[1],1));
		objects[3]->translate(float3(0,-1.5,0));
		objects.push_back(new ClippedQuadric(materials[2],0));
		objects[4]->translate(float3(0,2,0));
		//I am moving it 90 degrees here,but I do not like it
		//objects[5]->rotate(float3(1,0,0),3.14159/2);
		objects[4]->rotate(float3(1,0,0),3.5/2);
		objects[4]->translateB(float3(0,0,1.74));
		objects[4]->translate(float3(0,0,-8));
		objects[4]->scale(float3(0.3,0.3,0.3));

		//testing
		objects.push_back(new HoleQuadric(materials[0]));
		objects[5]->translate(float3(0,1,0));

		//a joke\\objects.push_back(new Quadric(materials[2],2.2,float3(0,-5,0)));
		/*objects.push_back(new Quadric(materials[0],1.0,float3(0,-2,0)));*/
		//objects.push_back(new Sphere(float3(0,-0.5,0),1.2, materials[0]));
		//objects.push_back(new Sphere(float3(0,1,0),0.9,materials[0]));
				
		
	}
	~Scene()
	{
		for (std::vector<LightSource*>::iterator iLightSource = lightSources.begin(); iLightSource != lightSources.end(); ++iLightSource)
			delete *iLightSource;
		for (std::vector<Material*>::iterator iMaterial = materials.begin(); iMaterial != materials.end(); ++iMaterial)
			delete *iMaterial;
		for (std::vector<Intersectable*>::iterator iObject = objects.begin(); iObject != objects.end(); ++iObject)
			delete *iObject;		
	}

public:
	Camera& getCamera()
	{
		return camera;
	}

	// IMPLEMENTED FOR YOUR CONVENIENCE, CALL THIS WHEN APPROPRIATE
	Hit firstIntersect(const Ray& ray)
	{
		Hit bestHit;
		bestHit.t = FLT_MAX;
		for(int oi=0; oi < objects.size(); oi++)
		{
			Hit hit = objects[oi]->intersect(ray);
			if(hit.t > 0 && hit.t < bestHit.t)
				bestHit = hit;
		}
		if(bestHit.t == FLT_MAX)
			return Hit();
		return bestHit;
	}

	float3 trace(const Ray& ray, int depth)
	{
		//Sphere s(float3(0, 0, 0), 1, NULL);
		Hit hit = firstIntersect(ray);

		//if(hit.t < 0||depth>5)
		if(hit.t<0)
		{
			Perlin per;
			return float3(0,3*ray.dir.x+ cos(8*ray.dir.y),sin(.3*ray.dir.x)+ 8*ray.dir.y)*per.marble(hit.position);
			//return ray.dir*ray.dir;
		}
			
	//
		float3 contribution = float3(0,0,0);
		for(int i = 0; i<lightSources.size(); i++){
			Ray shadowRay(hit.position+hit.normal*.00001,lightSources[i]->getLightDirAt(hit.position));
			Hit shadowHit = firstIntersect(shadowRay);
			float dist = (shadowHit.position-(hit.position)).norm();
			if(shadowHit.t <0|| dist>lightSources[i]->getDistanceFrom(hit.position)){
				contribution = contribution + (hit.material)->shade(hit.position,hit.normal, -ray.dir,
					lightSources[i]->getLightDirAt(hit.position),lightSources[i]->getPowerDensityAt(hit.position));
			}
		}


		if(hit.material->reflective){
			if(depth==0)
			{
				contribution+= float3(0,0,0);
			}
			else
			{
				float3 reflectionDir = hit.material->reflect(ray.dir, hit.normal);
				Ray reflectedRay(hit.position + hit.normal*0.0001, reflectionDir );
				//k^2 WHAT ON EARTH IS K HERE AWWWWWWWW I am guessing it is 0 in our case
				// WHAT IS GOING ON ??float3 f0=hit.material->minReflectance;
				float f0 = ((hit.material->refractiveIndex-1)*(hit.material->refractiveIndex-1))/((hit.material->refractiveIndex+1)*(hit.material->refractiveIndex+1));
				float3 normalDir = ray.dir;
				float costheta = normalDir.normalize().dot(hit.normal.normalize());
				contribution += trace(reflectedRay,depth-1)
                            * (f0+(1-f0)*(1-costheta)*(1-costheta)*(1-costheta)*(1-costheta)*(1-costheta));
			}
		}

			return contribution;
	
	}
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
// global application data

// screen resolution
const int screenWidth = 600;
const int screenHeight = 600;
// image to be computed by ray tracing
float3 image[screenWidth*screenHeight];

Scene scene;

bool computeImage()
{
	static unsigned int iPart = 0;

	if(iPart >= 64)
		return false;
    for(int j = iPart; j < screenHeight; j+=64)
	{
        for(int i = 0; i < screenWidth; i++)
		{
			float3 pixelColor = float3(0, 0, 0);
			float2 ndcPixelCentre( (2.0 * i - screenWidth) / screenWidth, (2.0 * j - screenHeight) / screenHeight );

			Camera& camera = scene.getCamera();
			Ray ray = Ray(camera.getEye(), camera.rayDirFromNdc(ndcPixelCentre));
			
			image[j*screenWidth + i] = scene.trace(ray,5);
		}
	}
	iPart++;
	return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// OpenGL starts here. In the ray tracing example, OpenGL just outputs the image computed to the array.

// display callback invoked when window needs to be redrawn
void onDisplay( ) {
    glClearColor(0.1f, 0.2f, 0.3f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // clear screen

	if(computeImage())
		glutPostRedisplay();
    glDrawPixels(screenWidth, screenHeight, GL_RGB, GL_FLOAT, image);
 
    glutSwapBuffers(); // drawing finished
}

int main(int argc, char **argv) {
    glutInit(&argc, argv);						// initialize GLUT
    glutInitWindowSize(600, 600);				// startup window size 
    glutInitWindowPosition(100, 100);           // where to put window on screen
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);    // 8 bit R,G,B,A + double buffer + depth buffer
 
    glutCreateWindow("Ray caster");				// application window is created and displayed
 
    glViewport(0, 0, screenWidth, screenHeight);

    glutDisplayFunc(onDisplay);					// register callback
 
    glutMainLoop();								// launch event handling loop
    
    return 0;
}

