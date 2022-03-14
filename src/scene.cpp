#include "scene.hpp"


using namespace cgp;





void scene_structure::initialize()
{
	// Initialize the shapes of the scene
	// ***************************************** //
	int N = 100;
	// Set the behavior of the camera and its initial position
	environment.camera.axis = camera_spherical_coordinates_axis::z;
	//environment.camera.look_at({ 5.0f,-4.0f,2.0f }, { 0,0,0 });

	// Create a visual frame representing the coordinate system
	global_frame.initialize(mesh_primitive_frame(), "Frame");

	// Create the ground plane
	ground.initialize(mesh_primitive_quadrangle({-2,-2,-1}, {2,-2,-1}, {2, 2,-1}, {-2, 2,-1}), "Ground");

	// send data to GPU and store it into a curve_drawable structure


	
	shape = mesh_primitive_grid({ 0,0,0 }, { 1,0,0 }, { 1,1,0 }, { 0,1,0 }, N, 5* N);
	for (int i = 0; i < N * N * 5; i++){
		shape.position.at(i)=vec3(2 * shape.position.at(i).x, 5 * shape.position.at(i).y, shape.position.at(i).z);
	}
	floor = mesh_primitive_grid({ 0,0,0 }, { 1,0,0 }, { 1,1,0 }, { 0,1,0 }, N, 5* N);
	for (int i = 0; i < N * N * 5; i++){
		floor.position.at(i)=vec3(2 * floor.position.at(i).x, 5 * floor.position.at(i).y, -1 * sloped_floor(5 * floor.position.at(i).y, 2.0, 0.13) - 0.1);
	}
	environment.camera.look_at({ 5.0f,-4.0f,2.0f }, shape.position.at(N * N * 5 - N/2 ));
	//shape.position.at(N/2, 5*N/2, n/2)
	initial_position = shape.position;
	initial_time = timer.t;
	shape_visual.initialize(shape, "Deforming shape");
	shape_visual.shading.color = { 0.6f, 0.6f, 0.9f };
	floor_visual.initialize(floor, "Deforming floor");
	floor_visual.shading.color = { 0.6f, 0.6f, 0.0f };
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
	//draw(ground, environment);

	
	// conditional display of the global frame (set via the GUI)
	if (gui.display_frame)
		draw(global_frame, environment);


	draw(floor_visual, environment);
	if (gui.display_wireframe)
		draw_wireframe(floor_visual, environment, { 0,0,0 });

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
	ImGui::SliderFloat("Wind Strenght", &wind_str, 0.0f, 10.0f, "%.1f");
	ImGui::Checkbox("Wireframe", &gui.display_wireframe);
	bool const restart = ImGui::Button("Restart");

	if (restart)
		initialize();
	if (restart)
		evolve_shape();
}

void scene_structure::evolve_shape()
{
    size_t const N = initial_position.size();

	//wind
	vec3 wind_dir = normalize(vec3(1.0, 1.0, 0.0));
	//wind_str = 3.f;

	float R = 0.007065f * std::pow(wind_str, 2.5) / 2; //0.06;
	float w = 9.81f * std::sqrt(2.f/3.f) / wind_str ; //4.0;

	float K = w * w ;/// 9.81f;   //10.0;
	float lambda = 1.5;

	//std::cout << "R = " << R << " w = " << w << " K = " << K <<"\n";
    for(size_t k=0; k<N; ++k)
    {
        
		vec3 const& p0 = initial_position[k];
        vec3& p        = shape.position[k];
		float depth = 0.1 + sloped_floor( p.y, 2.0, 0.13);//std::log(p.y);
		float K_depth = K / std::sqrt( std::tanh(K * depth));

        //p.x = p0.x + R * std::sin(K * dot(wind_dir, vec3(1,0,0)) * p0.x - w * timer.t - lambda * (p.z - p0.z) * (timer.t - initial_time));
		p.y = p0.y + R * std::sin( K_depth * p0.y + w * timer.t + lambda * (p.z - p0.z) * (timer.t - initial_time));
		p.z = p0.z - R* std::cos(K_depth * p0.y + w * timer.t + lambda * (p.z - p0.z) * (timer.t - initial_time));
		
    }
}

float sloped_floor(float x, float limit, float slope){
	if (x < limit){
		return slope * x;
	}
	else {
		return slope * limit;
	}
}

