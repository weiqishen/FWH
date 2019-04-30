// Reading all the input files

#include "Reader.h"
#include "param_reader.h"


using namespace std;

Reader::Reader(char *input_fnameC)
{
	input_fnameS.assign(input_fnameC);
}

void Reader::read()
{
	hid_t fid, attr_id, dataspace_id;
	hid_t str_ftype,str_type;
	char **temp_field;
	hsize_t n_fields;
	ndarray<string> fields;

	//read input parameters
	read_input();

	//open hdf5 data file
	fid = H5Fopen(input.fwh_surf_fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	if (fid < 0)
		Fatal_Error("Can't open HDF5 data file");

	//read geometry
	read_face_ctr(fid);
	read_face_normal(fid);
	read_face_area(fid);
	faces.n_eles=faces.A.get_len();

	//read fields
	attr_id = H5Aopen(fid, "fields", H5P_DEFAULT);
	if (attr_id < 0)
		Fatal_Error("Can't open fields");
	dataspace_id = H5Aget_space(attr_id);
	H5Sget_simple_extent_dims(dataspace_id, &n_fields, NULL);
	fields.setup(n_fields); 
	temp_field = new char *[n_fields];
    str_ftype = H5Aget_type(attr_id);
    str_type = H5Tget_native_type(str_ftype, H5T_DIR_ASCEND);
	H5Aread(attr_id, str_type, temp_field);
	for (size_t i = 0; i < (size_t)n_fields; i++)
	{
		fields(i).assign(temp_field[i]);
		delete[] temp_field[i];
	}
	delete[] temp_field;
	H5Tclose(str_ftype);
    H5Tclose(str_type);
	H5Sclose(dataspace_id);
	H5Aclose(attr_id);
	//test if fields are 'rho u v w pressure'
	if (fields(0) != "rho" || fields(1) != "u" || fields(2) != "v" || fields(3) != "w" || fields(4) != "pressure")
		Fatal_Error("Fields in the fwh surface file are not in the order of \"rho u v w pressure\"");

	//read dt and set up source time
	attr_id = H5Aopen(fid, "dt", H5P_DEFAULT);
	if (attr_id < 0)
		Fatal_Error("Can't open dt");
	H5Aread(attr_id, H5T_NATIVE_DOUBLE, &faces.dt);
	H5Aclose(attr_id);

	//read field variables
	read_data(fid);

	//close handles
	H5Fclose(fid);
}

void Reader::read_face_ctr(hid_t fid)
{
	hid_t ctr_id, dataspace_id;
	hsize_t dim[2];
	cout << " -> Reading the face Centers data from the file... " << flush;

	//open dataset
	ctr_id = H5Dopen2(fid, "coord", H5P_DEFAULT);
	if (ctr_id < 0)
		Fatal_Error("Can't open dataset \"coord\"");
	dataspace_id = H5Dget_space(ctr_id);
	H5Sget_simple_extent_dims(dataspace_id, dim, NULL);

	// Allocating the memory for face center coordinates
	faces.center.setup({(size_t)dim[1], (size_t)dim[0]});

	// Reading the face Centers file and storing the data in global variables.
	H5Dread(ctr_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, faces.center.get_ptr());

	H5Sclose(dataspace_id);
	H5Dclose(ctr_id);

	cout << "Done." << endl;
}

void Reader::read_face_normal(hid_t fid)
{
	hid_t normal_id, dataspace_id;
	hsize_t dim[2];
	cout << " -> Reading the face Normals data from the file... " << flush;

	//open dataset
	normal_id = H5Dopen2(fid, "normal", H5P_DEFAULT);
	if (normal_id < 0)
		Fatal_Error("Can't open dataset \"normal\"");
	dataspace_id = H5Dget_space(normal_id);
	H5Sget_simple_extent_dims(dataspace_id, dim, NULL);

	// Allocating the memory for face normals
	faces.normal.setup({(size_t)dim[1], (size_t)dim[0]});

	// Reading the face Centers file and storing the data in global variables.
	H5Dread(normal_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, faces.normal.get_ptr());

	H5Sclose(dataspace_id);
	H5Dclose(normal_id);

	cout << "Done." << endl;
}

void Reader::read_face_area(hid_t fid)
{
	hid_t area_id, dataspace_id;
	hsize_t dim;
	cout << " -> Reading the face Areas data from the file... " << flush;

	//open dataset
	area_id = H5Dopen2(fid, "area", H5P_DEFAULT);
	if (area_id < 0)
		Fatal_Error("Can't open dataset \"area\"");
	dataspace_id = H5Dget_space(area_id);
	H5Sget_simple_extent_dims(dataspace_id, &dim, NULL);

	// Allocating the memory for face normals
	faces.A.setup((size_t)dim);

	// Reading the face Centers file and storing the data in global variables.
	H5Dread(area_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, faces.A.get_ptr());

	H5Sclose(dataspace_id);
	H5Dclose(area_id);

	cout << "Done." << endl;
}

void Reader::read_data(hid_t fid)
{
	hid_t dataset_id, dataspace_id;
	hsize_t dim[3];

	cout << " -> Reading field variables on the surface file... " << flush;

	dataset_id = H5Dopen2(fid, "data", H5P_DEFAULT);
	if (dataset_id < 0)
		Fatal_Error("Can't open dataset");
	dataspace_id = H5Dget_space(dataset_id);
	H5Sget_simple_extent_dims(dataspace_id, dim, NULL);

	//setup source time
	faces.tau.setup((size_t)dim[2]);
	for (size_t i = 0; i < faces.tau.get_len(); i++)
		faces.tau(i) = i * faces.dt;

	//initialize data array
	faces.data.setup({(size_t)dim[2], (size_t)dim[1], (size_t)dim[0]});//time*face*field

	//read data
	H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, faces.data.get_ptr());

	//close handles
	H5Sclose(dataspace_id);
	H5Dclose(dataset_id);

	cout << "Done." << endl;
}

void Reader::read_input()
{
	cout << " -> Reading input parameters... " << flush;

	param_reader pr(input_fnameS);
	pr.openFile();

	//read ambient parameters
	pr.getScalarValue("T_static", input.T_static, 300.);
	pr.getScalarValue("gamma", input.gamma, 1.4);
	pr.getScalarValue("R_gas", input.R_gas, 286.9);
	pr.getScalarValue("fwh_surf_fname", input.fwh_surf_fname);
	pr.getScalarValue("output_fname", input.output_fname);
	//read observer positions
	pr.getVectorValue("observer_x", microphone.x);
	pr.getVectorValue("observer_y", microphone.y);
	pr.getVectorValue("observer_z", microphone.z);
	microphone.n_oberver = microphone.x.get_len();

	if (microphone.y.get_len() != microphone.n_oberver || microphone.z.get_len() != microphone.n_oberver)
		Fatal_Error("Number of observers not agree in each dimension");

	//endcap averaging
	pr.getScalarValue("endcap_avg", input.endcap_avg, 0);
	if (input.endcap_avg)
	{
		pr.getVectorValue("endcap_x", input.endcap_x);
	}
	pr.closeFile();
	cout << "Done.\n";
}

Reader::~Reader()
{
	
}