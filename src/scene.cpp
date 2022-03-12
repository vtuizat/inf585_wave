#include "scene.hpp"


using namespace cgp;





void scene_structure::initialize()
{
	// Initialize the shapes of the scene
	// ***************************************** //

	// Set the behavior of the camera and its initial position
	environment.camera.axis = camera_spherical_coordinates_axis::z;
	environment.camera.look_at({ 5.0f,-4.0f,2.0f }, { 0,0,0 });

	// Create a visual frame representing the coordinate system
	global_frame.initialize(mesh_primitive_frame(), "Frame");

	// Create the ground plane
	ground.initialize(mesh_primitive_quadrangle({-2,-2,-1}, {2,-2,-1}, {2, 2,-1}, {-2, 2,-1}), "Ground");

	// send data to GPU and store it into a curve_drawable structure


	int N = 100;
	shape = mesh_primitive_grid({ 0,0,0 }, { 1,0,0 }, { 1,1,0 }, { 0,1,0 }, N, 5* N);
	for (int i = 0; i < N * N * 5; i++){
		shape.position.at(i)=vec3(2 * shape.position.at(i).x, 5 * shape.position.at(i).y, shape.position.at(i).z);
	}
	initial_position = shape.position;
	initial_time = timer.t;
	shape_visual.initialize(shape, "Deforming shape");
	shape_visual.shading.color = { 0.6f, 0.6f, 0.9f };

		// Reset the color of the shape to white (only the texture image will be seen)
	shape_visual.shading.color = {1,1,1};

	// Load the image and associate the texture id to the structure
	//shape_visual.texture = opengl_load_texture_image("assets/squirrel.jpg");

}


void scene_structure::display()
{
	// Set the light to the current position of the camera
	environment.light = environment.camera.position();
	// Update the current elapsed time
	timer.update();


	// the general syntax to display a mesh is:
	//   draw(mesh_drawableName, environment);
	//     Note: scene is used to set the uniform parameters associated to the camera, light, etc. to the shader
	draw(ground, environment);

	
	// conditional display of the global frame (set via the GUI)
	if (gui.display_frame)
		draw(global_frame, environment);


	draw(shape_visual, environment);
	if (gui.display_wireframe)
		draw_wireframe(shape_visual, environment, { 0,0,0 });
	
	evolve_shape();
	shape_visual.update_position(shape.position);

	// Recompute normals on the CPU (given the position and the connectivity currently in the mesh structure)
	shape.compute_normal();
	// Send updated normals on the GPU
	shape_visual.update_normal(shape.normal);

}


void scene_structure::display_gui()
{
	ImGui::Checkbox("Frame", &gui.display_frame);
	ImGui::SliderFloat("Time Scale", &timer.scale, 0.0f, 2.0f, "%.1f");
	ImGui::Checkbox("Wireframe", &gui.display_wireframe);
}

void scene_structure::evolve_shape()
{
    size_t const N = initial_position.size();

	//wind
	vec3 wind_dir = normalize(vec3(1.0, 1.0, 0.0));
	float wind_str = 3.f;

	float R = 0.007065f * std::pow(wind_str, 2.5) / 2; //0.06;
	float w = 9.81f * std::sqrt(2.f/3.f) / wind_str ; //4.0;

	float K = w * w ;/// 9.81f;   //10.0;
	float lambda = 1.5;

	//std::cout << "R = " << R << " w = " << w << " K = " << K <<"\n";
    for(size_t k=0; k<N; ++k)
    {
        
		vec3 const& p0 = initial_position[k];
        vec3& p        = shape.position[k];
		
        //p.x = p0.x + R * std::sin(K * dot(wind_dir, vec3(1,0,0)) * p0.x - w * timer.t - lambda * (p.z - p0.z) * (timer.t - initial_time));
		p.y = p0.y + R * std::sin(K * dot(wind_dir, vec3(0,1,0)) * p0.y - w * timer.t - lambda * (p.z - p0.z) * (timer.t - initial_time));
		p.z = p0.z - R* std::cos(K * p0.y - w * timer.t - lambda * (p.z - p0.z) * (timer.t - initial_time));
		
    }
}


