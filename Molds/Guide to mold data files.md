Each mold must be defined as a MATLAB function using the following system:

[X,Y,Z,z,F,DT,MoldEdge] = <Mold name>(Inp)

where:
- X,Y,Z are grid matrices than when passed to the surf function plots the inside of the net boundary of the mold.
- DT, z are respectively a 2D Delayney triangulation object (in XY-plane) and corresponding z-coordinates of each vertex.
- F: interpolation object created with the built-in scatteredInterpolant function.
- MoldEdge: cell array with three cells where each cell defines a net boundary segment by means of points:
	---> MoldEdge{1}: left net boundary (in ascending y-coordinates)
	---> MoldEdge{2}: right net boundary (in ascending y-coordinates)
	---> MoldEdge{1}: bottom/root net boundary (in ascending x-coordinates)

With multiple layers, the above information must be stored in a cell array with one cell for each layer, e.g.
Mold{3}.MoldEdge{2}: Right mold net boundary line of 3rd layer.

The MATLAB function must be located in a folder of the same name, i.e. from the root directory: 
./Molds/<Mold name>/<Mold name>.m

The idea is that for each mold subfolder, there is a similar subfolder with input files under ./Input files/. It is possible to have different variants of the mold that all use the same input files, e.g. with different layer offsets or different discretizations. This can be done by adding three underscores after the main mold name and afterwards the extension, e.g. <Mold name>___<Extension>.m. This MATLAB function must still be located in the subfolder <Mold name>.

X,Y,Z and MoldEdge will subsequently be stored in the struct Mold. The remaining variables were not stored due to performance considerations because they are called very many times.
