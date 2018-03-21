precision highp float;

struct PointLight {
  vec3 position;
  vec3 color;
};

struct Material {
  vec3  diffuse;
  vec3  specular;
  float glossiness;
  //Added 
  float reflection;
  float refraction;
  float refractionIndex;
};


struct Sphere {
  vec3 position;
  float radius;
  Material material;
};

struct Plane {
  vec3 normal;
  float d;
  Material material;
};

struct Cylinder {
  vec3 position;
  vec3 direction;  
  float radius;
  Material material;
};

const int lightCount = 2;
const int sphereCount = 3;
const int planeCount = 1;
const int cylinderCount = 2;

struct Scene {
  vec3 ambient;
  PointLight[lightCount] lights;
  Sphere[sphereCount] spheres;
  Plane[planeCount] planes;
  Cylinder[cylinderCount] cylinders;
};

struct Ray {
  vec3 origin;
  vec3 direction;
};

// Contains all information pertaining to a ray/object intersection
struct HitInfo {
  bool hit;
  float t;
  vec3 position;
  vec3 normal;
  Material material;
};

HitInfo getEmptyHit() {
  return HitInfo(
    false, 
    0.0, 
    vec3(0.0), 
    vec3(0.0), 
  	// Depending on the material definition extension you make, this constructor call might need to be extened as well
    Material(vec3(0.0), vec3(0.0), 0.0, 0.0, 0.0, 0.0)
	);
}

// Sorts the two t values such that t1 is smaller than t2
void sortT(inout float t1, inout float t2) {
  // Make t1 the smaller t
  if(t2 < t1)  {
    float temp = t1;
    t1 = t2;
    t2 = temp;
  }
}

// Tests if t is in an interval
bool isTInInterval(const float t, const float tMin, const float tMax) {
  return t > tMin && t < tMax;
}

// Get the smallest t in an interval
bool getSmallestTInInterval(float t0, float t1, const float tMin, const float tMax, inout float smallestTInInterval) {
  
  sortT(t0, t1);
  
  // As t0 is smaller, test this first
  if(isTInInterval(t0, tMin, tMax)) {
  	smallestTInInterval = t0;
    return true;
  }
  
  // If t0 was not in the interval, still t1 could be
  if(isTInInterval(t1, tMin, tMax)) {
  	smallestTInInterval = t1;
    return true;
  }  
  
  // None was
  return false;
}

HitInfo intersectSphere(const Ray ray, const Sphere sphere, const float tMin, const float tMax) {
              
    vec3 to_sphere = ray.origin - sphere.position;
  
    float a = dot(ray.direction, ray.direction);
    float b = 2.0 * dot(ray.direction, to_sphere);
    float c = dot(to_sphere, to_sphere) - sphere.radius * sphere.radius;
    float D = b * b - 4.0 * a * c;
    if (D > 0.0)
    {
		float t0 = (-b - sqrt(D)) / (2.0 * a);
		float t1 = (-b + sqrt(D)) / (2.0 * a);
      
      	float smallestTInInterval;
      	if(!getSmallestTInInterval(t0, t1, tMin, tMax, smallestTInInterval)) {
          return getEmptyHit();
        }
      
      	vec3 hitPosition = ray.origin + smallestTInInterval * ray.direction;      

      	vec3 normal = 
          	length(ray.origin - sphere.position) < sphere.radius + 0.001? 
          	-normalize(hitPosition - sphere.position) : 
      		normalize(hitPosition - sphere.position);      

        return HitInfo(
          	true,
          	smallestTInInterval,
          	hitPosition,
          	normal,
          	sphere.material);
    }
    return getEmptyHit();
}

HitInfo intersectPlane(const Ray ray,const Plane plane, const float tMin, const float tMax) {
	// First we have to verify the plane is not parallel to our ray. 
    // If the dot product is equal to Zero it means the ray does not hit the plane.   	
    float parallel = dot(plane.normal,ray.direction);
    if (parallel != 0.0)
    {  
  	
      //We find the line from the ray to the plane. 
      //We compute t using the formula to find the intersection between line and a plane   
      float t = -1.0 *dot(plane.normal, ray.origin+plane.d)/ parallel; 
      
      //We check that the t value is between our range. 
      // If t is too high it means the ray is hitting a plane too far in the scene
      // If t is too low it means the ray is not hitting anything in tghe scene.
      // For both of these cases we return that the ray does not intersect with anything
      if(!isTInInterval(t,tMin,tMax)) {
        return getEmptyHit();
      }
      
      // We compute the hit position using the equation of a line.
      vec3 hitPosition = ray.origin + t * ray.direction;  
      
        return HitInfo(
          	true,
          	t,
          	hitPosition,
          	plane.normal,
          	plane.material);     
    } 
   return getEmptyHit();
}



float lengthSquared(vec3 x) {
  return dot(x, x);
}

HitInfo intersectCylinder(const Ray ray, const Cylinder cylinder, const float tMin, const float tMax) {
	//We are intersecting a line with an infinite cylinder.
  	//A Cylinder has equation x^2 + z^2 - r^2 = 0.
  	//http://mrl.nyu.edu/~dzorin/rend05/lecture2.pdf   next I used the information from this page.
  	//Using that information, we know to solve for t we have to solve the quadratic equation At^2+Bt+C = 0
  	//The fact that the equation for t is a quadratic equation makes sense because t will hit the cylinder twice
  	//Our task now is too find the first hit point. This is the point with the smallest t value.
  
  	vec3 p = ray.origin-cylinder.position;
  	float VdotVa = dot(ray.direction,cylinder.direction);
  
  	float A = dot((ray.direction-VdotVa*cylinder.direction),(ray.direction-VdotVa*cylinder.direction));
  	float B = 2.0*dot((ray.direction-VdotVa*cylinder.direction),p-(dot(p,cylinder.direction)*cylinder.direction));
  	float C = dot((p-(dot(p,cylinder.direction)*cylinder.direction)),(p-(dot(p,cylinder.direction)*cylinder.direction)))-cylinder.radius*cylinder.radius;
  
  	//As it is a quadratic formula we compute the discriminant
  	float d = B*B-4.0*A*C;
  	
  	//t will only have solutions when the discrimant is higher or equal to one
  	if(d >= 0.0){
      
      // We compute the two solutions to the quadratic equation
      float t0 = (-B - sqrt(d)) / (2.0 * A);
		float t1 = (-B + sqrt(d)) / (2.0 * A);
        
      	// the t we are looking for is the smallest value between the two t values and the one inside the range
  		// If we do not find a t value which corresponds to our conditions it means we have not hit the cylinder	
      float smallestTInInterval;
      	if(!getSmallestTInInterval(t0, t1, tMin, tMax, smallestTInInterval)) {
          return getEmptyHit();
        }
      	
      //We compute the hit position now that we know the t value of the equation.
      	vec3 hitPosition = ray.origin + smallestTInInterval * ray.direction;      
		
      //This part is to find the normal of the hit point
      // We find the zposition of the hitpoint on the cylinder
   		float zPos = dot((hitPosition-cylinder.position),cylinder.direction);
      // The normal is then equal to the difference between the hitpoint and the point on the cylinder's direction at the same height
      	vec3 normal = normalize(hitPosition-(cylinder.position+zPos*cylinder.direction));

        return HitInfo(
          	true,
          	smallestTInInterval,
          	hitPosition,
          	normal,
          	cylinder.material);
    }
    return getEmptyHit();
}


HitInfo getBetterHitInfo(const HitInfo oldHitInfo, const HitInfo newHitInfo) {
	if(newHitInfo.hit)
  		if(newHitInfo.t < oldHitInfo.t)  // No need to test for the interval, this has to be done per-primitive
          return newHitInfo;
  	return oldHitInfo;
}

HitInfo intersectScene(const Scene scene, const Ray ray, const float tMin, const float tMax) {
  HitInfo bestHitInfo;
  bestHitInfo.t = tMax;
  bestHitInfo.hit = false;
  for (int i = 0; i < cylinderCount; ++i) {
    bestHitInfo = getBetterHitInfo(bestHitInfo, intersectCylinder(ray, scene.cylinders[i], tMin, tMax));
  }
  for (int i = 0; i < sphereCount; ++i) {
    bestHitInfo = getBetterHitInfo(bestHitInfo, intersectSphere(ray, scene.spheres[i], tMin, tMax));
  }
  for (int i = 0; i < planeCount; ++i) {
    bestHitInfo = getBetterHitInfo(bestHitInfo, intersectPlane(ray, scene.planes[i], tMin, tMax));
  }
  
  return bestHitInfo;
}

vec3 shadeFromLight(
  const Scene scene,
  const Ray ray,
  const HitInfo hit_info,
  const PointLight light)
{ 
  vec3 hitToLight = light.position - hit_info.position;
  
  vec3 lightDirection = normalize(hitToLight);
  vec3 viewDirection = normalize(hit_info.position - ray.origin);
  vec3 reflectedDirection = reflect(viewDirection, hit_info.normal);
  float diffuse_term = max(0.0, dot(lightDirection, hit_info.normal));
  float specular_term  = pow(max(0.0, dot(lightDirection, reflectedDirection)), hit_info.material.glossiness);
   
  float visibility = 1.0;
  
  //We cast a ray from the hitpoint to the light
  Ray intersectToLight;
  intersectToLight.origin = hit_info.position;
  intersectToLight.direction = lightDirection;
  
  HitInfo hitsBetweenHitandLight;
  //We then intersect the ray with the scene 
  hitsBetweenHitandLight = intersectScene(scene, intersectToLight, 0.0001, length(hitToLight)); 
  
  //If the ray hits an object in the scene we set the visbility to zero
  if(hitsBetweenHitandLight.hit){
 	visibility = 0.0;   
  }

     
  return 	visibility * 
    		light.color * (
    		specular_term * hit_info.material.specular +
      		diffuse_term * hit_info.material.diffuse);
}

vec3 background(const Ray ray) {
  // A simple implicit sky that can be used for the background
  return vec3(0.2) + vec3(0.8, 0.6, 0.5) * max(0.0, ray.direction.y);
}

// It seems to be a WebGL issue that the third parameter needs to be inout instea dof const on Tobias' machine
vec3 shade(const Scene scene, const Ray ray, inout HitInfo hitInfo) {
  
  	if(!hitInfo.hit) {
  		return background(ray);
  	}
  
    vec3 shading = scene.ambient * hitInfo.material.diffuse;
    for (int i = 0; i < lightCount; ++i) {
        shading += shadeFromLight(scene, ray, hitInfo, scene.lights[i]); 
    }
    return shading;
}


Ray getFragCoordRay(const vec2 frag_coord) {
  	float sensorDistance = 1.0;
  	vec2 sensorMin = vec2(-1, -0.5);
  	vec2 sensorMax = vec2(1, 0.5);
  	vec2 pixelSize = (sensorMax- sensorMin) / vec2(800, 400);
  	vec3 origin = vec3(0, 0, sensorDistance);
    vec3 direction = normalize(vec3(sensorMin + pixelSize * frag_coord, -sensorDistance));  
  
  	return Ray(origin, direction);
}

float fresnel(const vec3 viewDirection, const vec3 normal) {
  	// fresnel computes the intensty of the reflection depending on the angle of the incoming ray
 	// We want the reflection to be the strongest when the angle gets greater
  	//In other words in the middle of the sphere there is little reflection but more on the sides
  	// Thus we take the absolute value of the dot product of the normal and the direction
  	// We know that this value is high low when the angle gets closer to 90 degrees which represents the sides of the sphere
  // Since we want the opposite relationship i.e. high reflectivity on the sides, we do 1.0 - the dot product
	return 1.0 - abs(dot(normalize(viewDirection), normalize(normal)));
}

vec3 colorForFragment(const Scene scene, const vec2 fragCoord) {
      
    Ray initialRay = getFragCoordRay(fragCoord);  
  	HitInfo initialHitInfo = intersectScene(scene, initialRay, 0.0001, 10000.0);  
  	vec3 result = shade(scene, initialRay, initialHitInfo);
	
  	Ray currentRay;
  	HitInfo currentHitInfo;
  	
  	// Compute the reflection
  	currentRay = initialRay;
  	currentHitInfo = initialHitInfo;
  	
  	// The initial strength of the reflection
  	float reflectionWeight = 1.0;
  	
  	const int maxReflectionStepCount = 2;
  	for(int i = 0; i < maxReflectionStepCount; i++) {
      
      if(!currentHitInfo.hit) break;
     
      //We add the fresnel constant here so that the reflectivity at this hit is affected
      reflectionWeight *= fresnel(currentRay.direction, currentHitInfo.normal)*currentHitInfo.material.reflection;
      
      Ray nextRay;
	  //We set the next ray's direction equal to reflection of the current ray coming on the object
      // To compute the reflection we need to know the normal and the direction of the incoming ray.
      // The angle of reflection of the incoming ray and the normal will be equal to the angle between the reflected ray and the normal
      nextRay.origin = currentHitInfo.position;
      nextRay.direction = reflect(currentRay.direction, currentHitInfo.normal);


      currentRay = nextRay;
      
      currentHitInfo = intersectScene(scene, currentRay, 0.0001, 10000.0);      
            
      result += reflectionWeight * shade(scene, currentRay, currentHitInfo);
    }
  
  	// Compute the refraction
  	currentRay = initialRay;  
  	currentHitInfo = initialHitInfo;
   
  	// The initial medium is air
  	float currentIOR = 1.0;

  	// The initial strength of the refraction.
  	float refractionWeight = 1.0;
  
  	const int maxRefractionStepCount = 2;
  	for(int i = 0; i < maxRefractionStepCount; i++) {
      
      if(!currentHitInfo.hit) break;

      // add the fresnel value here. This is the refraction fresnel value so we use the relationship between the two
      // ReflectionConstant = 1 - RefractionConstant
      refractionWeight *= (1.0-fresnel(currentRay.direction, currentHitInfo.normal))*currentHitInfo.material.refraction;
    
      float index = currentHitInfo.material.refractionIndex;
      
      //we have to care about the 2 cases of reflection:
      //- when we are outside the object and going in; we invert the refraction index according to snell's law
      //- when we are inside and coming out; we invert the normal
      if(dot(currentRay.direction,currentHitInfo.normal) < 0.0){
       index=1.0/index;
      } 
      else {
       currentHitInfo.normal = -currentHitInfo.normal;
      }
      
      //the next ray is then equal to refraction of the current direction and the normal and the index of refraction 
      Ray nextRay;
      nextRay.origin = currentHitInfo.position;
      nextRay.direction = refract(currentRay.direction, currentHitInfo.normal, index);
      
      currentRay = nextRay;
      
      currentHitInfo = intersectScene(scene, currentRay, 0.001, 10000.0);
            
      result += refractionWeight * shade(scene, currentRay, currentHitInfo);
      
      
    }
  return result;
}

Material getDefaultMaterial() {
  // Will need to update this to match the new Material definition
  return Material(vec3(0.3), vec3(0), 1.0, 0.0, 0.0, 0.0);
}

Material getPaperMaterial() {
	//Paper is white and compact 
   // diffuse = it is bright to relatively high value
  //  specular = not strong
  // glossiness = low paper is not reflective  
  // reflection = nothing
  // refraction = nothing in this case because it is very low in the real world
  // refractionIndex = none

  return Material(vec3(0.6), vec3(0.2), 2.0, 0.0, 0.0, 0.0);
}

Material getPlasticMaterial() {
  // Plastic is yellow and glossy
   //  diffuse = diffuse towards a yellow color
  //  specular = pretty high 
  // glossiness = plastic is very high
  // reflection = very low because the plastic is not clear here. It is yellow
  // refraction = none
  // refractionIndex = none not translucid
  return Material(vec3(1,1,0), vec3(1.0), 20.0, 0.1, 0.0, 0.0);
}

Material getGlassMaterial() {
  // Glass is transparent and not moderatily reflective.
  //  diffuse = nothing because it is transparent
  //  specular = none
  // glossiness = none
  // reflection = average, glass is more see through than reflective
  // refraction = high
  // refractionIndex = the official value is 1.55. However, I reduced it to make it appear more natural
  return Material(vec3(0.05), vec3(0.01), 0.01, 0.79, 1.0, 1.2);
}

Material getSteelMirrorMaterial() {
 // Steel is compact and reflective
  // Replace by your definition of a steel mirror material
  //  diffuse = very low
  //  specular = none
  // glossiness = none
  // reflection = high
  // refraction = none 
  // refractionIndex = none
  return Material(vec3(0.01), vec3(0.0), 0.01, 0.6, 0.0, 0.0);
}

vec3 tonemap(const vec3 radiance) {
  const float monitorGamma = 2.0;
  return pow(radiance, vec3(1.0 / monitorGamma));
}

void main()
{
    // Setup scene
    Scene scene;
  	scene.ambient = vec3(0.12, 0.15, 0.2);
  
    // Lights
    scene.lights[0].position = vec3(5, 15, -5);
    scene.lights[0].color    = 0.5 * vec3(0.8, 0.6, 0.5);
    
  	scene.lights[1].position = vec3(-15, 10, 2);
    scene.lights[1].color    = 0.5 * vec3(0.5, 0.7, 1.0);
  
    // Primitives
    scene.spheres[0].position            	= vec3(8, -2, -13);
    scene.spheres[0].radius              	= 4.0;
    scene.spheres[0].material 				= getPaperMaterial();
    
  	scene.spheres[1].position            	= vec3(-7, -1, -13);
    scene.spheres[1].radius             	= 4.0;
    scene.spheres[1].material				= getPlasticMaterial();
  
    scene.spheres[2].position            	= vec3(0, 0.5, -5);
    scene.spheres[2].radius              	= 2.0;
    scene.spheres[2].material   			= getGlassMaterial();

  	scene.planes[0].normal            		= vec3(0, 1, 0);
  	scene.planes[0].d              			= 4.5;
    scene.planes[0].material				= getSteelMirrorMaterial();
  
  	scene.cylinders[0].position            	= vec3(-1, 1, -18);
  	scene.cylinders[0].direction            = normalize(vec3(-1, 2, -1));
  	scene.cylinders[0].radius         		= 1.5;
    scene.cylinders[0].material				= getPaperMaterial();
  
  	scene.cylinders[1].position            	= vec3(3, 1, -5);
  	scene.cylinders[1].direction            = normalize(vec3(1, 4, 1));
  	scene.cylinders[1].radius         		= 0.25;
    scene.cylinders[1].material				= getPlasticMaterial();

  // compute color for fragment
  gl_FragColor.rgb = tonemap(colorForFragment(scene, gl_FragCoord.xy));
  gl_FragColor.a = 1.0;
}
