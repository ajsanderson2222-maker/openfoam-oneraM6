"""
Render 3D wing surface mesh wireframe using paraview offscreen rendering.
"""
from paraview.simple import *
import os

case_dir = "/home/andrew/openfoam_jobs/openfoam-oneraM6/simulation"
out_dir  = "/home/andrew/openfoam_jobs/openfoam-oneraM6"

foam = OpenFOAMReader(registrationName="case.foam",
                      FileName=os.path.join(case_dir, "case.foam"))
foam.MeshRegions = ["patch/wing"]
foam.CellArrays  = []
foam.UpdatePipeline(time=2000)

edges = ExtractEdges(Input=foam)
edges.UpdatePipeline()

view = CreateRenderView()
view.ViewSize               = [2400, 1600]
view.Background             = [1, 1, 1]
view.UseLight               = 0
view.OrientationAxesVisibility = 0

disp = Show(edges, view)
disp.Representation = "Wireframe"
disp.AmbientColor   = [0, 0, 0]
disp.DiffuseColor   = [0, 0, 0]
disp.LineWidth      = 0.8

Render(view)

# Top-down parallel view matching wing_cp_topdown.png
# x=chord (left→right), z=span (bottom→top), camera above in y
# Wing centre: x~0.5, z~0.6; parallel_scale drives half-height in world (z-dir)
# z span=1.22 → half=0.61; window aspect=1.5 → x_shown=0.61*1.5*2=1.83 > chord(1.14) ✓
# Match wing_cp_topdown.png: x/c left→right, span (z) bottom→top
# Camera above in +y looking down, ViewUp = +x so x goes right, z goes up
view.CameraParallelProjection = 1
view.CameraPosition      = [0.5, 10.0, 0.61]
view.CameraFocalPoint    = [0.5,  0.0, 0.61]
view.CameraViewUp        = [0,    0,   1   ]   # z (span) points up in image
view.CameraParallelScale = 0.68   # half-height in z-world = ~half span

Render(view)
SaveScreenshot(os.path.join(out_dir, "wing_mesh_3d.png"), view,
               ImageResolution=[2400, 1600])

# Mirror left/right to match wing_cp_topdown.png orientation
from PIL import Image
import numpy as np
img = Image.open(os.path.join(out_dir, "wing_mesh_3d.png"))
img = img.transpose(Image.FLIP_LEFT_RIGHT)
img.save(os.path.join(out_dir, "wing_mesh_3d.png"))
print("Saved wing_mesh_3d.png")
