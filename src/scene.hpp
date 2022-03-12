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
	cgp::buffer<cgp::vec3> initial_position;
	float initial_time;
	cgp::mesh_drawable shape_visual;

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





