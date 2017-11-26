#define PROJECTION
#define RASTERIZATION
#define CLIPPING
#define INTERPOLATION
//#define ZBUFFERING
//#define ANIMATION

precision highp float;
uniform float time;

// Polygon / vertex functionality
const int MAX_VERTEX_COUNT = 8;

uniform ivec2 viewport;

struct Vertex {
    vec3 position;
    vec3 color;
};

struct Polygon {
    // Numbers of vertices, i.e., points in the polygon
    int vertexCount;
    // The vertices themselves
    Vertex vertices[MAX_VERTEX_COUNT];
};

// Appends a vertex to a polygon
void appendVertexToPolygon(inout Polygon polygon, Vertex element) {
    for (int i = 0; i < MAX_VERTEX_COUNT; ++i) {
        if (i == polygon.vertexCount) {
            polygon.vertices[i] = element;
        }
    }
    polygon.vertexCount++;
}

// Copy Polygon source to Polygon destination
void copyPolygon(inout Polygon destination, Polygon source) {
    for (int i = 0; i < MAX_VERTEX_COUNT; ++i) {
        destination.vertices[i] = source.vertices[i];
    }
    destination.vertexCount = source.vertexCount;
}

// Get the i-th vertex from a polygon, but when asking for the one behind the last, get the first again
Vertex getWrappedPolygonVertex(Polygon polygon, int index) {
    if (index >= polygon.vertexCount) index -= polygon.vertexCount;
    for (int i = 0; i < MAX_VERTEX_COUNT; ++i) {
        if (i == index) return polygon.vertices[i];
    }
}

// Creates an empty polygon
void makeEmptyPolygon(out Polygon polygon) {
  polygon.vertexCount = 0;
}

// Clipping part

#define ENTERING 0
#define LEAVING 1
#define OUTSIDE 2
#define INSIDE 3

int getCrossType(Vertex poli1, Vertex poli2, Vertex wind1, Vertex wind2) {
#ifdef CLIPPING

  	//calculate half space equation between wind1 and wind2 y= slope x + b
  	//to do that we have to find the slope and the c constant
    float slope =  (wind2.position[1]-wind1.position[1]) / (wind2.position[0]-wind1.position[0]);  
  	float c = wind2.position[1]- (slope*wind2.position[0]);
  	
  	//find value for each of the coordinates given
  	float resultPoli1 = (slope*poli1.position[0] + c - poli1.position[1])*(wind2.position[0]-wind1.position[0]);
  	float resultPoli2 = (slope*poli2.position[0] + c - poli2.position[1])*(wind2.position[0]-wind1.position[0]);  	
  	
  	//For each of the different cases
	if(resultPoli1 > 0.0 && resultPoli2 > 0.0){
     		return INSIDE;
    }
  	if(resultPoli1 < 0.0 && resultPoli2 < 0.0){
         	return OUTSIDE;
    } 
  	if(resultPoli1 > 0.0 && resultPoli2 < 0.0){
      		return LEAVING;
    }
  	if(resultPoli1 < 0.0 && resultPoli2 > 0.0){
      		return ENTERING;
    }     
  	
#else
    return INSIDE;
#endif
}

// This function assumes that the segments are not parallel or collinear.
Vertex intersect2D(Vertex a, Vertex b, Vertex c, Vertex d) {
#ifdef CLIPPING
	
  	//we have two parametric equations intersecting
  	// p1 = c + s*CD
  	// p2 = a + t*AB
  	// we have to find the values for which s = t

  	float AB1 = b.position[0]-a.position[0];
  	float AB2 = b.position[1]-a.position[1];
  	float AB3 = b.position[2]-a.position[2];
 	float a1 = a.position[0];
    float a2 = a.position[1];
  	float a3 = a.position[2];
  	
  	float CD1 = d.position[0]-c.position[0];
  	float CD2 = d.position[1]-c.position[1];
  	float CD3 = d.position[2]-c.position[2];
    float c1 = c.position[0];
    float c2 = c.position[1];
  	float c3 = c.position[2];
  
  	float s = (AB1*(a2-c2) + AB2*(c1 - a1)) / (CD2*AB1-CD1*AB2);
  	
  	//now that we have s plug first parametric equation to get x y z values
	Vertex intersection;
	intersection.position[0] = c1 + s*CD1;
    intersection.position[1] = c2 + s*CD2;
  	 
    // z adjusted for interpolation
  	float d3 = d.position[2];
  	intersection.position[2] = 1.0 / ((1.0/c3) + s*((1.0/d3)-(1.0/c3)));
  	return intersection;
  	
#else
    return a;
#endif
}

void sutherlandHodgmanClip(Polygon unclipped, Polygon clipWindow, out Polygon result) {
    Polygon clipped;
    copyPolygon(clipped, unclipped);

    // Loop over the clip window
    for (int i = 0; i < MAX_VERTEX_COUNT; ++i) {
        if (i >= clipWindow.vertexCount) break;

        // Make a temporary copy of the current clipped polygon
        Polygon oldClipped;
        copyPolygon(oldClipped, clipped);

        // Set the clipped polygon to be empty
        makeEmptyPolygon(clipped);

        // Loop over the current clipped polygon
        for (int j = 0; j < MAX_VERTEX_COUNT; ++j) {
            if (j >= oldClipped.vertexCount) break;
            
            // Handle the j-th vertex of the clipped polygon. This should make use of the function 
            // intersect() to be implemented above.
#ifdef CLIPPING
			Vertex p0 = getWrappedPolygonVertex(oldClipped , j);
            Vertex p1 = getWrappedPolygonVertex(oldClipped , j+1);
          	
          	Vertex clipWindow1 = getWrappedPolygonVertex(clipWindow , i);;
          	Vertex clipWindow2 = getWrappedPolygonVertex(clipWindow , i+1); ;
          	
          	if(getCrossType(p0,p1,clipWindow1, clipWindow2) == 0) {
              //entering
              // add intersection and p1
              appendVertexToPolygon(clipped, intersect2D(p0, p1, clipWindow1, clipWindow2));
              appendVertexToPolygon(clipped, p1);
              
            }
          	if(getCrossType(p0,p1,clipWindow1, clipWindow2) == 1) {
              //leaving
              appendVertexToPolygon(clipped, intersect2D(p0, p1, clipWindow1, clipWindow2));
            } 
          	if(getCrossType(p0,p1,clipWindow1, clipWindow2) == 2) {
              //outside
              //do not do anything
            } 
          	if(getCrossType(p0,p1,clipWindow1, clipWindow2) == 3){
          	  //inside
              appendVertexToPolygon(clipped, p1);
        	}
          	
#else
            appendVertexToPolygon(clipped, getWrappedPolygonVertex(oldClipped, j));
#endif
        }
    }

    // Copy the last version to the output
    copyPolygon(result, clipped);
}

// Rasterization and culling part

#define INNER_SIDE 0
#define OUTER_SIDE 1

// Assuming a clockwise (vertex-wise) polygon, returns whether the input point 
// is on the inner or outer side of the edge (ab)
int edge(vec2 point, Vertex a, Vertex b) {
#ifdef RASTERIZATION
  	//calculate half space equation involves finding y = ax+c
    float slope =  (b.position[1]-a.position[1]) / (b.position[0]-a.position[0]);
    
  	//b = y/(-ax)
  	float c = a.position[1]- (slope*a.position[0]);
  	
  	//Equation = ax + b - y
  	float result = (slope*point[0] + c - point[1])*(b.position[0]-a.position[0]);
  	
  	if(result > 0.0){
      return INNER_SIDE;
    }
  
#endif
    return OUTER_SIDE;
}

// Returns if a point is inside a polygon or not
bool isPointInPolygon(vec2 point, Polygon polygon) {
    // Don't evaluate empty polygons
    if (polygon.vertexCount == 0) return false;
    // Check against each edge of the polygon
    bool rasterise = true;
    for (int i = 0; i < MAX_VERTEX_COUNT; ++i) {
        if (i < polygon.vertexCount) {
#ifdef RASTERIZATION
          	Vertex a = getWrappedPolygonVertex(polygon,i);
            Vertex b = getWrappedPolygonVertex(polygon,i+1);
  			int halfTest = edge(point, a, b);
          	if(halfTest == OUTER_SIDE){
              rasterise = false;
            }
#else
      	 rasterise = false;
#endif
        }
    }
    return rasterise;
}

bool isPointOnPolygonVertex(vec2 point, Polygon polygon) {
    for (int i = 0; i < MAX_VERTEX_COUNT; ++i) {
        if (i < polygon.vertexCount) {
          	ivec2 pixelDifference = ivec2(abs(polygon.vertices[i].position.xy - point) * vec2(viewport));
          	int pointSize = viewport.x / 200;
            if( pixelDifference.x <= pointSize && pixelDifference.y <= pointSize) {
              return true;
            }
        }
    }
    return false;
}

float triangleArea(vec2 a, vec2 b, vec2 c) {
    // https://en.wikipedia.org/wiki/Heron%27s_formula
    float ab = length(a - b);
    float bc = length(b - c);
    float ca = length(c - a);
    float s = (ab + bc + ca) / 2.0;
    return sqrt(max(0.0, s * (s - ab) * (s - bc) * (s - ca)));
}

Vertex interpolateVertex(vec2 point, Polygon polygon) {
    float weightSum = 0.0;
    vec3 colorSum = vec3(0.0);
    vec3 positionSum = vec3(0.0);
    float depthSum = 0.0;
    for (int i = 0; i < MAX_VERTEX_COUNT; ++i) {
        if (i < polygon.vertexCount) {
#if defined(INTERPOLATION) || defined(ZBUFFERING)
    		
          //For current vertex we will draw a triangle between 1 verteces before and after the current vertex and the point given.
          //we then find the area of this triangle which is equal to the weight of how much the color of point A affects the point given. 
          //We add the weighted color of A to the color of the point given
          //we repeat this for all verteces in the polygon
      
          //A is the current point we are weighting
          vec2 A = vec2(getWrappedPolygonVertex(polygon,i).position[0],getWrappedPolygonVertex(polygon,i).position[1]);
          vec3 Acolor = getWrappedPolygonVertex(polygon,i).color;
          
          //First point in the triangle
          vec2 B = vec2(getWrappedPolygonVertex(polygon,i+1).position[0],getWrappedPolygonVertex(polygon,i+1).position[1]);
          
          //Second point in the triangle
          vec2 C; 
          if(i == 0){
            C = vec2(getWrappedPolygonVertex(polygon,polygon.vertexCount - 1).position[0],getWrappedPolygonVertex(polygon,polygon.vertexCount - 1).position[1]);
          }else {
          	C = vec2(getWrappedPolygonVertex(polygon,i-1).position[0],getWrappedPolygonVertex(polygon,i-1).position[1]); 
          }
         
          // find the area of triangle
          float weightOfA = triangleArea(B,C,point);
          
          
          
#else
#endif
#ifdef ZBUFFERING
    // Put your code here
#endif
#ifdef INTERPOLATION
   		  //add the weighted color
          weightSum += weightOfA;
          colorSum += Acolor*weightOfA;
#endif
        }
    }
    
    Vertex result = polygon.vertices[0];
  
#ifdef INTERPOLATION
   	result.color = colorSum / weightSum;
#endif
#ifdef ZBUFFERING
    // Put your code here
#endif
#if !defined(INTERPOLATION) && !defined(ZBUFFERING)
    // Put your code here
#endif

  return result;
}

// Projection part

// Used to generate a projection matrix.
mat4 computeProjectionMatrix() {
    mat4 projectionMatrix = mat4(1);
  
  float aspect = float(viewport.x) / float(viewport.y);  
  float imageDistance = 0.5;

#ifdef PROJECTION
    //using standard perspective projection seen in https://en.wikibooks.org/wiki/GLSL_Programming/Vertex_Transformations
  	// aspect ratio a is already given and n is imageDistance.
    // d is found with angle of 0.69 radians obtained through trial and error to look like given picture
 	float d = 1.0/tan(0.69/2.0);    
    //assume f is infinite
    float f = 10000.0;
      
  	projectionMatrix[0] = vec4(d/aspect,0,0,0);
  	projectionMatrix[1] = vec4(0,d,0,0);
  	projectionMatrix[2] = vec4(0,0,(imageDistance+f)/(imageDistance-f),-2.0*imageDistance*f/(imageDistance-f));
  	projectionMatrix[3] = vec4(0,0,-1,0);  
  	
#endif
  
    return projectionMatrix;
}

// Used to generate a simple "look-at" camera. 
mat4 computeViewMatrix(vec3 VRP, vec3 TP, vec3 VUV) {
    mat4 viewMatrix = mat4(1);

#ifdef PROJECTION
  
 	//This was done following the slides on creating a general camera
     
  	// Define VPN as vector pointing away from the camera
  	vec3 VPN = TP - VRP;
  	
  	// Generate the camera axes.
    vec3 n = normalize(VPN);
    vec3 u = normalize(cross(VUV, n));
    vec3 v = normalize(cross(n, u));
  	//reference point
  	vec4 q = vec4(0,0,0,1.0);
  	vec3 t = vec3(- dot(VRP, u), - dot(VRP, v), - dot(VRP, n));
  	
 	// add the camera axes to the view matrix including 
  	viewMatrix[0]  = vec4(u[0], v[0], n[0], q[0]);
    viewMatrix[1]  = vec4(u[1], v[1], n[1], q[1]);
  	viewMatrix[2]  = vec4(u[2], v[2], n[2], q[2]);
  	viewMatrix[3]  = vec4(t, q[3]);
#endif
    return viewMatrix;
}

vec3 getCameraPosition() {  
#ifdef ANIMATION
    // Put your code here
#else
    return vec3(0, 0, 10);
#endif
}

// Takes a single input vertex and projects it using the input view and projection matrices
vec3 projectVertexPosition(vec3 position) {

  // Set the parameters for the look-at camera.
    vec3 TP = vec3(0, 0, 0);
  	vec3 VRP = getCameraPosition();
    vec3 VUV = vec3(0, 1, 0);
  
    // Compute the view matrix.
    mat4 viewMatrix = computeViewMatrix(VRP, TP, VUV);
  
  // Compute the projection matrix.
    mat4 projectionMatrix = computeProjectionMatrix();
  
#ifdef PROJECTION
  	mat4 transformation = projectionMatrix*viewMatrix;
  	vec4 newpos = transformation*vec4(position,1.0);
 	return vec3(newpos)/newpos[3];
#endif
}

// Projects all the vertices of a polygon
void projectPolygon(inout Polygon projectedPolygon, Polygon polygon) {
    copyPolygon(projectedPolygon, polygon);
    for (int i = 0; i < MAX_VERTEX_COUNT; ++i) {
        if (i < polygon.vertexCount) {
            projectedPolygon.vertices[i].position = projectVertexPosition(polygon.vertices[i].position);
        }
    }
}

// Draws a polygon by projecting, clipping, ratserizing and interpolating it
void drawPolygon(
  vec2 point, 
  Polygon clipWindow, 
  Polygon oldPolygon, 
  inout vec3 color, 
  inout float depth)
{
    Polygon projectedPolygon;
    projectPolygon(projectedPolygon, oldPolygon);  
  
    Polygon clippedPolygon;
    sutherlandHodgmanClip(projectedPolygon, clipWindow, clippedPolygon);

    if (isPointInPolygon(point, clippedPolygon)) {
      
        Vertex interpolatedVertex = 
          interpolateVertex(point, projectedPolygon);
#if defined(ZBUFFERING)    
    // Put your code here
#else
      // Put your code to handle z buffering here
      color = interpolatedVertex.color;
      depth = interpolatedVertex.position.z;      
#endif
   }
  
   if (isPointOnPolygonVertex(point, clippedPolygon)) {
        color = vec3(1);
   }
}

// Main function calls

void drawScene(vec2 pixelCoord, inout vec3 color) {
    color = vec3(0.3, 0.3, 0.3);
  
  	// Convert from GL pixel coordinates 0..N-1 to our screen coordinates -1..1
    vec2 point = 2.0 * pixelCoord / vec2(viewport) - vec2(1.0);

    Polygon clipWindow;
    clipWindow.vertices[0].position = vec3(-0.65,  0.95, 1.0);
    clipWindow.vertices[1].position = vec3( 0.65,  0.75, 1.0);
    clipWindow.vertices[2].position = vec3( 0.75, -0.65, 1.0);
    clipWindow.vertices[3].position = vec3(-0.75, -0.85, 1.0);
    clipWindow.vertexCount = 4;
  
  	// Draw the area outside the clip region to be dark
    color = isPointInPolygon(point, clipWindow) ? vec3(0.5) : color;

    const int triangleCount = 2;
    Polygon triangles[triangleCount];
  
    triangles[0].vertices[0].position = vec3(-2, -2, 0.0);
    triangles[0].vertices[1].position = vec3(4, 0, 3.0);
    triangles[0].vertices[2].position = vec3(-1, 2, 0.0);
    triangles[0].vertices[0].color = vec3(1.0, 0.5, 0.2);
    triangles[0].vertices[1].color = vec3(0.8, 0.8, 0.8);
    triangles[0].vertices[2].color = vec3(0.2, 0.5, 1.0);
    triangles[0].vertexCount = 3;
  
    triangles[1].vertices[0].position = vec3(3.0, 2.0, -2.0);
  	triangles[1].vertices[2].position = vec3(0.0, -2.0, 3.0);
    triangles[1].vertices[1].position = vec3(-1.0, 2.0, 4.0);
    triangles[1].vertices[1].color = vec3(0.2, 1.0, 0.1);
    triangles[1].vertices[2].color = vec3(1.0, 1.0, 1.0);
    triangles[1].vertices[0].color = vec3(0.1, 0.2, 1.0);
    triangles[1].vertexCount = 3;

    float depth = 10000.0;
    // Project and draw all the triangles
    for (int i = 0; i < triangleCount; i++) {
        drawPolygon(point, clipWindow, triangles[i], color, depth);
    }   
}

void main() {
    drawScene(gl_FragCoord.xy, gl_FragColor.rgb);
    gl_FragColor.a = 1.0;
}
