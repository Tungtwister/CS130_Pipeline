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
    for(int i = 0; i < width*height; i++)
    {
      state.image_depth[i] = 2.0;
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
        for(int i = 0; i < state.num_vertices * state.floats_per_vertex; i += 3*state.floats_per_vertex)
        {
            const data_geometry * geo[3];
            data_geometry g0,g1,g2;
            data_vertex v0,v1,v2;
            
            g0.data = &state.vertex_data[i];
            g1.data = &state.vertex_data[i + state.floats_per_vertex];
            g2.data = &state.vertex_data[i + (state.floats_per_vertex * 2)];
  
            v0.data = &state.vertex_data[i];
            v1.data = &state.vertex_data[i + state.floats_per_vertex];
            v2.data = &state.vertex_data[i + (state.floats_per_vertex * 2)];
            
            state.vertex_shader(v0,g0,state.uniform_data);
            state.vertex_shader(v1,g1,state.uniform_data);
            state.vertex_shader(v2,g2,state.uniform_data);
            
            for(int j = 0; j < 3; j++)
            {
              g0.gl_Position[j] /= g0.gl_Position[3];
              g1.gl_Position[j] /= g1.gl_Position[3];
              g2.gl_Position[j] /= g2.gl_Position[3];
            }
            
            geo[0] = &g0;
            geo[1] = &g1;
            geo[2] = &g2;
  
            rasterize_triangle(state,geo);
        }
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


//    std::cout<< out[0].gl_Position[0] <<std::endl;
//    std::cout<< out[0].gl_Position[1] <<std::endl;
//    std::cout<< out[1].gl_Position[0] <<std::endl;
//    std::cout<< out[1].gl_Position[1] <<std::endl;
//    std::cout<< out[2].gl_Position[0] <<std::endl;
//    std::cout<< out[2].gl_Position[1] <<std::endl;
//    std::cout<< "\n";
    int Ax = state.image_width/2.0 * in[0]->gl_Position[0] + (state.image_width/2.0) - (0.5);
    int Ay = state.image_height/2.0 * in[0]->gl_Position[1] + (state.image_height/2.0) - (0.5);
    
    int Bx = state.image_width/2.0 * in[1]->gl_Position[0] + (state.image_width/2.0) - (0.5);
    int By = state.image_height/2.0 * in[1]->gl_Position[1] + (state.image_height/2.0) - (0.5);
    
    int Cx = state.image_width/2.0 * in[2]->gl_Position[0] + (state.image_width/2.0) - (0.5);
    int Cy = state.image_height/2.0 * in[2]->gl_Position[1] + (state.image_height/2.0) - (0.5);
    
    int min_x = 0;
    int max_x = 0;
    int min_y = 0;
    int max_y = 0;
    
    min_x = std::min(0,Ax);
    min_x = std::min(min_x,Bx);
    min_x = std::min(min_x,Cx);
    
    max_x = std::max(0,Ax);
    max_x = std::max(max_x,Bx);
    max_x = std::max(max_x,Cx);
    
    min_y = std::min(0,Ay);
    min_y = std::min(min_y,By);
    min_y = std::min(min_y,Cy);
    
    max_y = std::max(0,Ay);
    max_y = std::max(max_y,By);
    max_y = std::max(max_y,Cy);
    
 
    float abc = 0.5 * ((Bx*Cy - Cx*By) - (Ax*Cy - Cx*Ay) + (Ax*By - Bx*Ay));
    //std::cout << abc <<std::endl;
    
    for(int i = 0; i < state.image_width; i++)
    {
      for(int j = 0; j < state.image_height; j++)
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
          //do interpolation
          
          data_fragment frag;
          float frag_data[MAX_FLOATS_PER_VERTEX];
          frag.data = frag_data;
          
          float zBuffer = alpha * in[0]->gl_Position[2]  + beta * in[1]->gl_Position[2] + gamma * in[2]->gl_Position[2];
          if(zBuffer < state.image_depth[i + j * state.image_width])
          {
            for(int k = 0; k < state.floats_per_vertex; k++)
            {
              float p = alpha / in[0]->gl_Position[3] + beta / in[1]->gl_Position[3] + gamma / in[2]->gl_Position[3];
              float a = alpha / (in[0]->gl_Position[3] * p);
              float b = beta / (in[1]->gl_Position[3] * p);
              float c = gamma / (in[2]->gl_Position[3] * p);
              
              switch(state.interp_rules[k])
              {
                case interp_type::flat:
                  frag_data[k] = in[0]->data[k];
                  //std::cout<< frag_data[k] <<std::endl; 
                  break;
                case interp_type::smooth:
                  frag_data[k] = (a * in[0]->data[k] + b * in[1]->data[k] + c * in[2]->data[k]);
                  break;
                case interp_type::noperspective:
                  frag_data[k] = (alpha * in[0]->data[k] + beta * in[1]->data[k] + gamma * in[2]->data[k]);
                  break;
                case interp_type::invalid:
                  break;
              }
            }
            
            
            
            //call fragment shader
            data_output dOut;
            state.fragment_shader(frag, dOut,state.uniform_data); 
            //std::cout<< color.output_color <<std::endl;
            int red = 255 * dOut.output_color[0];
            int green = 255 * dOut.output_color[1];
            int blue = 255 * dOut.output_color[2];
            
            //std::cout<< red <<std::endl;
            state.image_depth[i + j * state.image_width] = zBuffer;
      			state.image_color[i+j*state.image_width] = make_pixel(red, green, blue);
            
          }
    		}
      }
    }
}

