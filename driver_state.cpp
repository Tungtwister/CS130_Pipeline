#include "driver_state.h"
#include <cstring>
#include <cstdlib>
#include <algorithm>

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_color=0;
    state.image_depth=0;
    state.image_width=width;
    state.image_height=height;
    
    state.image_depth = new float[width*height];
    state.image_color = new pixel[width*height];
   	for(int i = 0; i < width*height; i++)
    {
      state.image_color[i] = make_pixel(0,0,0);
		}
    //delete state.image_color;
    //std::cout<<"TODO: allocate and initialize state.image_color and state.image_depth."<<std::endl;
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{
    //std::cout<<"TODO: implement rendering."<<std::endl;
    switch(type)
    {
      case render_type::triangle:
      {
        const data_geometry** g = new const data_geometry*[3];
        const data_geometry** h = new const data_geometry*[3];
        
        g[0] = new data_geometry;
  			g[1] = new data_geometry;
  			g[2] = new data_geometry;
  			h[0] = new data_geometry;
  			h[1] = new data_geometry;
  			h[2] = new data_geometry;
        
   			const_cast<data_geometry*>(g[0])->data = new float[MAX_FLOATS_PER_VERTEX];
  			const_cast<data_geometry*>(g[1])->data = new float[MAX_FLOATS_PER_VERTEX];
  			const_cast<data_geometry*>(g[2])->data = new float[MAX_FLOATS_PER_VERTEX];
  			const_cast<data_geometry*>(h[0])->data = new float[MAX_FLOATS_PER_VERTEX];
  			const_cast<data_geometry*>(h[1])->data = new float[MAX_FLOATS_PER_VERTEX];
  			const_cast<data_geometry*>(h[2])->data = new float[MAX_FLOATS_PER_VERTEX];
    		for(int i = 0; i < state.num_vertices*state.floats_per_vertex; i += 3*state.floats_per_vertex)
        {
    			for(int j = 0; j < 3; j++)
          {
    				g[j]->data[0] = state.vertex_data[(j *state.floats_per_vertex)];
    				g[j]->data[1] = state.vertex_data[(j *state.floats_per_vertex)+ 1];
    				g[j]->data[2] = state.vertex_data[(j *state.floats_per_vertex)+ 2];
    				h[j]->data[0] = state.vertex_data[i+(j *state.floats_per_vertex)];
    				h[j]->data[1] = state.vertex_data[i+(j *state.floats_per_vertex)+ 1];
    				h[j]->data[2] = state.vertex_data[i+(j *state.floats_per_vertex)+ 2];
          }
          //std::cout<< g[0]->data[0] <<std::endl;
          //std::cout<< h[0]->data[0] <<std::endl;
        }
        rasterize_triangle(state,g);
        rasterize_triangle(state,h);
        
        			delete [] g[0]->data;
        			delete [] g[1]->data;
        			delete [] g[2]->data;
        			delete g[0];
        			delete g[1];
        			delete g[2];
        			delete [] g;
        			delete [] h[0]->data;
        			delete [] h[1]->data;
        			delete [] h[2]->data;
        			delete h[0];
        			delete h[1];
        			delete h[2];
        			delete [] h;
      }
            break;
      case render_type::indexed:
            break;
      case render_type::fan:
            break;
      case render_type::strip:
            break;
      default: 
            break;
    }
}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry* in[3],int face)
{
    if(face==6)
    {
        rasterize_triangle(state, in);
        return;
    }
    std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    clip_triangle(state,in,face+1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
//    //std::cout<<"TODO: implement rasterization"<<std::endl;
//    //std::cout<< in[0]->gl_Position[0] <<std::endl;
    data_geometry* out = new data_geometry[3];  
    data_vertex v;
  
    for(int i = 0; i < 3; i++)
    {
    	v.data = in[i]->data;
    
    	state.vertex_shader(v, out[i], state.uniform_data);
    
    	out[i].gl_Position[0] /= out[i].gl_Position[3];
    	out[i].gl_Position[1] /= out[i].gl_Position[3];
      out[i].gl_Position[2] /= out[i].gl_Position[3];
      
      int x = state.image_width/2.0 * out[i].gl_Position[0] + state.image_width/2.0 - (0.5);
      int y = state.image_height/2.0 * out[i].gl_Position[1] + state.image_height/2.0 - (0.5);
      
      state.image_color[x+y*state.image_width] = make_pixel(255, 255, 255);
    }
//    std::cout<< out[0].gl_Position[0] <<std::endl;
//    std::cout<< out[0].gl_Position[1] <<std::endl;
//    std::cout<< out[1].gl_Position[0] <<std::endl;
//    std::cout<< out[1].gl_Position[1] <<std::endl;
//    std::cout<< out[2].gl_Position[0] <<std::endl;
//    std::cout<< out[2].gl_Position[1] <<std::endl;
//    std::cout<< "\n";
    int Ax = (state.image_width/2.0)*out[0].gl_Position[0] + (state.image_width/2.0) - (0.5);
    int Ay = (state.image_height/2.0)*out[0].gl_Position[1] + (state.image_height/2.0) - (0.5);
    
    int Bx = (state.image_width/2.0)*out[1].gl_Position[0] + (state.image_width/2.0) - (0.5);
    int By = (state.image_height/2.0)*out[1].gl_Position[1] + (state.image_height/2.0) - (0.5);
    
    int Cx = (state.image_width/2.0)*out[2].gl_Position[0] + (state.image_width/2.0) - (0.5);
    int Cy = (state.image_height/2.0)*out[2].gl_Position[1] + (state.image_height/2.0) - (0.5);
    
    int min_x = 0;
    int max_x = 0;
    int min_y = 0;
    int max_y = 0;
    
    min_x = std::min(Ax,Bx);
    min_x = std::min(min_x,Cx);
    
    max_x = std::max(Ax,Bx);
    max_x = std::max(max_x,Cx);
    
    min_y = std::min(Ay,By);
    min_y = std::min(min_y,Cy);
    
    max_y = std::max(Ay,By);
    max_y = std::max(max_y,Cy);
    
 
    float abc = 0.5 * ((Bx*Cy - Cx*By) - (Ax*Cy - Cx*Ay) + (Ax*By - Bx*Ay));
    //std::cout << abc <<std::endl;
    
    for(int i = min_x; i < max_x; i++)
    {
      for(int j = min_y; j < max_y; j++)
      {
        float pbc = 0.5 * ((Bx*Cy - Cx*By) - (i*Cy - Cx*j) + (i*By - Bx*j));
        float apc = 0.5 * ((i*Cy - Cx*j) - (Ax*Cy - Cx*Ay) + (Ax*j - i*Ay));
        float abp = 0.5 * ((Bx*j - i*By) - (Ax*j - i*Ay) + (Ax*By - Bx*Ay));		
    
    		float alpha = pbc / abc;
    		float beta = apc / abc;
    		float gamma = abp / abc;		
        //barycentric
    		if(alpha >= 0 && beta >= 0 && gamma >= 0)
        {
          //std::cout << "color" <<std::endl;
    			state.image_color[i+j*state.image_width] = make_pixel(255, 255, 255);
    		}
      }
    }
}

