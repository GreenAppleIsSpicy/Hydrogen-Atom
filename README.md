# Hydrogen Atom Orbital Viewer

Allows you to visualize the orbitals in the Hydrogen Atom in 2D or in 3D as a pointcloud which can be exported as a .PLY file and opened in programs like Blender to be rendered.
You will need to have an extra file, '\Electron_Orbitals\' in the same directory as the python file is running, it will be the saving location for your PLYs.

You will need the following Python libraries:
- Numpy
- Math
- Matplotlib
- Open3D (Only if you plan to export as a .ply)

## Rendering your files in Blender

### For experienced Blender users
- Import your file.
- Use "Mesh to Points" Geometry node and set a material.
- The color of the Pointcloud can be found in the attributes and can be assigned to the material in the shading tab.

### Plotting
Ensure you are using the function `Plot_Hydrogen_3D()` and have the argument `plysave` set to `True`. 
Once run, the file `'Electron_Positions_(n, l, m).ply'` will be created in the folder `'Electron_Positions'` with your specified values for `n`, `l` and `m`.

I reccomend running with the argument `plot` also set to `True` to ensure that your plot came out as you intended.

### Importing to Blender
In Blender navigate to File > Import > Stanford (.ply), then in the window navigate to `'Electron_Positions_(n, l, m).ply'`. Once selected you can hit import and your PLY will enter Blender.

### Coloring
From here navigate to the "Geometry Nodes" tab and ensure your pointcloud is selected. Select "New" in the lower window then in the same window Add > Mesh > Operations > Mesh to Points. Select and drag in between the line connecting "Group Input" and "Group Output." Then, Add > Material > Set Material, and drag between "Mesh to Points" and "Group Output." In "Set Material" click on 'Material' and select the option 'Material.'

#### If you had the heatmap enabled and want to use those colors:
- Head to the "Shading" tab, ensure your pointcloud is selected. Select "New" then Add > Input > Attribute, in "Attribute" click 'Name' and enter "Col". Drag from the dot labeled 'Color' on "Attribute" to the dot labeled 'Base Color' on "Principled BSDF." Then in  Return to the Geometry Nodes tab and in the "Set Material" Node select your newly created material (likely called Material.001).

Return to the "Layout" tab and in the top right of the viewport window, select 'Viewport Shading.' The points should become colored (if you enabled the heatmap colors) and become spheres. 
- If the points do not turn into spheres:
  - On the right hand side lower window.
  - Go to the "Render" menu.
  - Set 'Render Engine' to "Cycles."
- If the points do not become colored:
  - On the right hand side lower window.
  - Go to the "Data" menu.
  - Go to the "Attribute" dropdown.
  - Ensure that the name for the attribute Vertex>Color matches the name in the Attribute node in the Shading tab. (Should be "Col")

### Now you are free to roam Blender and all the cool stuff you can do with lighting and nodes!
