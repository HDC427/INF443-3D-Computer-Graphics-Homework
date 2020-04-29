#pragma once

#include "main/scene_base/base.hpp"

#ifdef SCENE_3D_GRAPHICS

// Stores some parameters that can be set from the GUI
struct gui_scene_structure
{
    bool wireframe;
};

struct vec3t{
    vcl::vec3 p; // position
    float t;     // time
};

struct scene_model : scene_base
{

    /** A part must define two functions that are called from the main function:
     * setup_data: called once to setup data before starting the animation loop
     * frame_draw: called at every displayed frame within the animation loop
     *
     * These two functions receive the following parameters
     * - shaders: A set of shaders.
     * - scene: Contains general common object to define the 3D scene. Contains in particular the camera.
     * - data: The part-specific data structure defined previously
     * - gui: The GUI structure allowing to create/display buttons to interact with the scene.
    */

    void setup_data(std::map<std::string,GLuint>& shaders, scene_structure& scene, gui_structure& gui);
    void frame_draw(std::map<std::string,GLuint>& shaders, scene_structure& scene, gui_structure& gui);

    void set_gui();

    // Exercise 2.3, function to initialize tree positions
    void update_tree_position();

    // Exercise 3.3, function to initialize grass positions
    void update_grass_position();

    // Exercise 4.2, inner bird animation
    void update_bird_inner(float t);

    // Exercise 4.3, initialize bird position and timer
    void initialize_bird();

     // Exercise 4.3, outer bird animation
    void update_global_position(vcl::hierarchy_mesh_drawable_node &obj, float t, vcl::buffer<vec3t> const& v);
    float K = 0.5;

    // visual representation of a surface
    vcl::mesh_drawable terrain;

    std::vector<vcl::vec3> tree_position;
    std::vector<vcl::vec3> grass_position;
    vcl::mesh_drawable trunk;
    vcl::mesh_drawable foliage;
    vcl::mesh_drawable grass;
    GLuint texture_id;
    GLuint billboard_id;

    vcl::hierarchy_mesh_drawable bird;
    // Data (p_i,t_i)
    vcl::buffer<vec3t> keyframes; // Given (position,time)

    gui_scene_structure gui_scene;
    vcl::timer_interval timer;
};

#endif


