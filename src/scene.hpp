#pragma once

/**
	Introduction to the use of CGP library
*/

#include "cgp/cgp.hpp"



struct gui_parameters {
	bool display_frame     = true;
	bool display_wireframe = false;
};



// The structure of the custom scene
struct scene_structure {
	
	// ****************************** //
	// Elements and shapes of the scene
	// ****************************** //

	
	cgp::mesh_drawable global_frame;  // The standard global frame
	
	// Different primitives displayed as mesh_drawable
	cgp::mesh_drawable cube;     
	cgp::mesh_drawable ground;
	cgp::mesh_drawable cylinder;
	cgp::mesh_drawable sphere;

	cgp::curve_drawable curve;    // A set of vertices displayed as a curve

	cgp::mesh shape;
	cgp::mesh floor;
	cgp::buffer<cgp::vec3> initial_position;
	float initial_time;
	cgp::mesh_drawable shape_visual;
	cgp::mesh_drawable floor_visual;
	float wind_str = 1.0;
	float wind_angle = 1.6;

	cgp::timer_basic timer; // A timer to have access to the elapsed time
	cgp::scene_environment_basic_camera_spherical_coords environment; // Standard environment controler
	gui_parameters gui;                       // Standard GUI element storage

	// ****************************** //
	// Functions
	// ****************************** //

	void initialize();  // Standard initialization to be called before the animation loop
	void display();     // The frame display to be called within the animation loop
	void display_gui(); // The display of the GUI, also called within the animation loop

	
	void evolve_shape();

};
float K_integration(float K ,float x0, float (*h)(float, float, float));

float sloped_floor(float x, float limit, float slope);
float atan_floor(float x);
float valley_floor(float x, float position, float width);
float valley_wall(float x, float position, float width);



