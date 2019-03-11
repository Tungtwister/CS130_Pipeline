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
            
//            for(int j = 0; j < 3; j++)
//            {
//              g0.gl_Position[j] /= g0.gl_Position[3];
//              g1.gl_Position[j] /= g1.gl_Position[3];
//              g2.gl_Position[j] /= g2.gl_Position[3];
//            }
            
            geo[0] = &g0;
            geo[1] = &g1;
            geo[2] = &g2;
            
            //std::cout<<  " zero" <<geo[2]->gl_Position[0] << " one" <<geo[2]->gl_Position[1] << " two" << geo[2]->gl_Position[2] << " w" <<geo[2]->gl_Position[3] <<std::endl;
            clip_triangle(state,geo,0);
            //rasterize_triangle(state,geo);
        }
      }
            break;
      case render_type::indexed:
      {
        for(int i = 0; i < 3*state.num_triangles; i += 3)
        {
            const data_geometry * geo[3];
            data_geometry g0,g1,g2;
            data_vertex v0,v1,v2;
            
            g0.data = &state.vertex_data[state.index_data[i] * state.floats_per_vertex];
            g1.data = &state.vertex_data[state.index_data[i+1] * state.floats_per_vertex];
            g2.data = &state.vertex_data[state.index_data[i+2] * state.floats_per_vertex];
  
            v0.data = &state.vertex_data[state.index_data[i] * state.floats_per_vertex];
            v1.data = &state.vertex_data[state.index_data[i+1] * state.floats_per_vertex];
            v2.data = &state.vertex_data[state.index_data[i+2] * state.floats_per_vertex];
            
            state.vertex_shader(v0,g0,state.uniform_data);
            state.vertex_shader(v1,g1,state.uniform_data);
            state.vertex_shader(v2,g2,state.uniform_data);
            
//            std::cout<< g0.gl_Position[0] << " " << g0.gl_Position[1] << " " << g0.gl_Position[2] <<std::endl;
//            std::cout<< g1.gl_Position[0] << " " << g1.gl_Position[1] << " " << g1.gl_Position[2] <<std::endl;
//            std::cout<< g2.gl_Position[0] << " " << g2.gl_Position[1] << " " << g2.gl_Position[2] <<std::endl;
            
//            for(int j = 0; j < 3; j++)
//            {
//              g0.gl_Position[j] /= g0.gl_Position[3];
//              g1.gl_Position[j] /= g1.gl_Position[3];
//              g2.gl_Position[j] /= g2.gl_Position[3];
//            }
            
            geo[0] = &g0;
            geo[1] = &g1;
            geo[2] = &g2;
  
            rasterize_triangle(state,geo);
        }
      }
            break;
      case render_type::fan:
      {
        for(int i = 0; i < state.num_vertices * state.floats_per_vertex - 2*state.floats_per_vertex; i += state.floats_per_vertex)
        {
            const data_geometry * geo[3];
            data_geometry g0,g1,g2;
            data_vertex v0,v1,v2;
            
            g0.data = &state.vertex_data[0];
            g1.data = &state.vertex_data[i+ state.floats_per_vertex];
            g2.data = &state.vertex_data[i + (state.floats_per_vertex * 2)];
  
            v0.data = &state.vertex_data[0];
            v1.data = &state.vertex_data[i + state.floats_per_vertex];
            v2.data = &state.vertex_data[i + (state.floats_per_vertex * 2)];
            
            state.vertex_shader(v0,g0,state.uniform_data);
            state.vertex_shader(v1,g1,state.uniform_data);
            state.vertex_shader(v2,g2,state.uniform_data);
            
//            for(int j = 0; j < 3; j++)
//            {
//              g0.gl_Position[j] /= g0.gl_Position[3];
//              g1.gl_Position[j] /= g1.gl_Position[3];
//              g2.gl_Position[j] /= g2.gl_Position[3];
//            }
            
            geo[0] = &g0;
            geo[1] = &g1;
            geo[2] = &g2;
  
            rasterize_triangle(state,geo);
        }
      }
            break;
      case render_type::strip:
      {
        for(int i = 0; i < state.num_vertices * state.floats_per_vertex  - 2*state.floats_per_vertex; i += state.floats_per_vertex)
        {
            const data_geometry * geo[3];
            data_geometry g0,g1,g2;
            data_vertex v0,v1,v2;
            int odd = i % 2;
      
            if(odd == 0)
            {
              g0.data = &state.vertex_data[i];
              g1.data = &state.vertex_data[i+ state.floats_per_vertex];
              g2.data = &state.vertex_data[i + (state.floats_per_vertex * 2)];
    
              v0.data = &state.vertex_data[i];
              v1.data = &state.vertex_data[i + state.floats_per_vertex];
              v2.data = &state.vertex_data[i + (state.floats_per_vertex * 2)];
            }
            else
            {
              g1.data = &state.vertex_data[i];
              g0.data = &state.vertex_data[i+ state.floats_per_vertex];
              g2.data = &state.vertex_data[i + (state.floats_per_vertex * 2)];
    
              v1.data = &state.vertex_data[i];
              v0.data = &state.vertex_data[i + state.floats_per_vertex];
              v2.data = &state.vertex_data[i + (state.floats_per_vertex * 2)];
            }
            
            
            state.vertex_shader(v0,g0,state.uniform_data);
            state.vertex_shader(v1,g1,state.uniform_data);
            state.vertex_shader(v2,g2,state.uniform_data);
            
            //std::cout<< g0.gl_Position[0] << " " << g0.gl_Position[1] << " " << g0.gl_Position[2] <<std::endl;
            
//            for(int j = 0; j < 3; j++)
//            {
//              g0.gl_Position[j] /= g0.gl_Position[3];
//              g1.gl_Position[j] /= g1.gl_Position[3];
//              g2.gl_Position[j] /= g2.gl_Position[3];
//            }
            
            geo[0] = &g0;
            geo[1] = &g1;
            geo[2] = &g2;
  
            rasterize_triangle(state,geo);
        }
      }
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
    //couldnt get clip triangle to work 
    if(face==6)
    {
        rasterize_triangle(state, in);
        return;
    }
    else
    {
      vec4 A = in[0]->gl_Position;
      vec4 B = in[1]->gl_Position;
      vec4 C = in[2]->gl_Position;
      //std::cout<< "C " << C[0] << " " << C[1] << " " << C[2] << " " << C[3] <<std::endl;
      //std::cout<< "in[0] " << in[0]->gl_Position[0] << in[0]->gl_Position[1] << in[0]->gl_Position[2] <<std::endl;
      //determine which plane we are looking at
      float axis;
      float sign;
      sign =  2 * (face % 2) - 1;
      axis = face % 3;
      //std::cout<< sign << axis <<std::endl;
      bool inA, inB, inC = false;

      //std::cout<< sign <<std::endl;
      if(sign == 1)
      {
        if(A[axis] <= A[3])
        {
          inA = true;
        }          
        if(B[axis] <= B[3])
        {
          inB = true;
        }          
        if(C[axis] <= C[3])
        {
          inC = true;
        }
      }
      else if(sign == -1)
      {
        if(A[axis] >= -A[3])
        {
          inA = true;
        }          
        if(B[axis] >= -B[3])
        {
          inB = true;
        }    
        if(C[axis] >= -C[3])
        {
          inC = true;
        }
      }
      //std::cout<< face << inA << inB << inC << std::endl;
      //if A,B,C inside triangle
      if(inA == true && inB == true && inC == true)
      {
         clip_triangle(state,in,face+1); 
      }
      //if A,B,C outside triangle
      if(inA == false && inB == false && inC == false)
      {
        return;
      }
      //if A inside, but B, C outside
      if(inA == true && inB == false && inC == false)
      {

      }
      //if B inside, but A, C outside
      if(inA == true && inB == false && inC == false)
      {
      }
      //if C inside, but A,B outside
      if(inA == false && inB == false && inC == true)
      {
      }
      
      //if A,B inside, but C outside
      if(inA == true && inB == true && inC == false)
      {
//        data_geometry inter1;
//        data_geometry inter2;
//        
//        inter1.data = new float[MAX_FLOATS_PER_VERTEX];
//        inter2.data = new float[MAX_FLOATS_PER_VERTEX];
//        
//        float alphaAC = (sign*C[3] - C[axis]) / (A[axis] - sign*A[3] + sign*C[3] - C[axis]);
//        float alphaCB = (sign*B[3] - B[axis]) / (C[axis] - sign*C[3] + sign*B[3] - B[axis]);
//        
//        vec4 AC = alphaAC * A + (1 - alphaAC) * C;
//        vec4 CB = alphaCB * C + (1 - alphaCB) * B;
//        
//        inter1.gl_Position = AC;
//        inter2.gl_Position = CB;
//        //std::cout<< inter1.gl_Position <<std::endl;
//        for(int i = 0; i < state.floats_per_vertex; i++)
//        {
//          switch(state.interp_rules[i])
//          {
//            case interp_type::flat:
//              inter1.data[i] = in[0]->data[i];
//              inter2.data[i] = in[0]->data[i];
//              break;
//            case interp_type::smooth:
//              break;
//            case interp_type::noperspective:
//              break;
//            case interp_type::invalid:
//              break;
//          }
//        }
//        in[2] = &inter1;
//        clip_triangle(state,in,face+1);
//        
//        in[0] = &inter1;
//        in[2] = &inter2;
//        clip_triangle(state,in,face+1);
      }
      //if A,C inside, but B outside
      if(inA == true && inB == false && inC == true)
      {
      }
      //if B,C inside, but A outside
      if(inA == false && inB == true && inC == true)
      {
      }
       
      
      
    }
    //std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
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
    data_geometry* out = new data_geometry[3];
    for(int i = 0; i < 3; i++)
    {
      out[i] = *in[i];
      
      out[i].gl_Position[0] /= out[i].gl_Position[3];
      out[i].gl_Position[1] /= out[i].gl_Position[3];
      out[i].gl_Position[2] /= out[i].gl_Position[3];
      
    }
    
    float Ax = state.image_width/2.0 * out[0].gl_Position[0] + (state.image_width/2.0) - (0.5);
    float Ay = state.image_height/2.0 * out[0].gl_Position[1] + (state.image_height/2.0) - (0.5);
    
    float Bx = state.image_width/2.0 * out[1].gl_Position[0] + (state.image_width/2.0) - (0.5);
    float By = state.image_height/2.0 * out[1].gl_Position[1] + (state.image_height/2.0) - (0.5);
    
    float Cx = state.image_width/2.0 * out[2].gl_Position[0] + (state.image_width/2.0) - (0.5);
    float Cy = state.image_height/2.0 * out[2].gl_Position[1] + (state.image_height/2.0) - (0.5);
    
    float min_x = 0.0;
    float max_x = 0.0;
    float min_y = 0.0;
    float max_y = 0.0;
    
    min_x = std::min(min_x,Ax);
    min_x = std::min(min_x,Bx);
    min_x = std::min(min_x,Cx);
    
    max_x = std::max(max_x,Ax);
    max_x = std::max(max_x,Bx);
    max_x = std::max(max_x,Cx);
    
    min_y = std::min(min_y,Ay);
    min_y = std::min(min_y,By);
    min_y = std::min(min_y,Cy);
    
    max_y = std::max(max_y,Ay);
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
          
          float zBuffer = alpha * out[0].gl_Position[2]  + beta * out[1].gl_Position[2] + gamma * out[2].gl_Position[2];
          if(zBuffer < state.image_depth[i + j * state.image_width])
          {
            for(int k = 0; k < state.floats_per_vertex; k++)
            {
              float p = alpha / out[0].gl_Position[3] + beta / out[1].gl_Position[3] + gamma / out[2].gl_Position[3];
              float a = alpha / (out[0].gl_Position[3] * p);
              float b = beta / (out[1].gl_Position[3] * p);
              float c = gamma / (out[2].gl_Position[3] * p);
              
              switch(state.interp_rules[k])
              {
                case interp_type::flat:
                  frag_data[k] = out[0].data[k];
                  //std::cout<< frag_data[k] <<std::endl; 
                  break;
                case interp_type::smooth:
                  frag_data[k] = (a * out[0].data[k] + b * out[1].data[k] + c * out[2].data[k]);
                  break;
                case interp_type::noperspective:
                  frag_data[k] = (alpha * out[0].data[k] + beta * out[1].data[k] + gamma * out[2].data[k]);
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

