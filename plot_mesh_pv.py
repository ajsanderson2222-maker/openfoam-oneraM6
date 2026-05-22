"""
Render mesh wireframe at z=0.034 root section using paraview offscreen rendering.
Zoomed to airfoil near-field: x=[-0.2,1.4], y=[-0.7,0.7].
"""
from paraview.simple import *
import os

case_dir = "/home/andrew/openfoam_jobs/openfoam-oneraM6/simulation"
out_dir  = "/home/andrew/openfoam_jobs/openfoam-oneraM6"

foam = OpenFOAMReader(registrationName="case.foam",
                      FileName=os.path.join(case_dir, "case.foam"))
foam.MeshRegions = ["internalMesh"]
foam.CellArrays  = []
foam.UpdatePipeline(time=2000)

slc = Slice(Input=foam)
slc.SliceType = "Plane"
slc.SliceType.Origin = [0, 0, 0.034]
slc.SliceType.Normal = [0, 0, 1]
slc.UpdatePipeline()

edges = ExtractEdges(Input=slc)
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
disp.LineWidth      = 0.5

# First render to let paraview do its auto-reset
Render(view)

# Override camera after auto-reset
# Target: x=[-0.2, 1.4] width=1.6, y=[-0.7, 0.7] height=1.4
# Window aspect = 2400/1600 = 1.5
# parallel_scale = half_height = 0.7; x_shown = 0.7*1.5*2 = 2.1 — wider than 1.6, fine
view.CameraParallelProjection = 1
view.CameraPosition      = [0.6, 0.0, 10]   # centre of x=[-0.2,1.4]
view.CameraFocalPoint    = [0.6, 0.0, 0.034]
view.CameraViewUp        = [0, 1, 0]
view.CameraParallelScale = 0.7

Render(view)
SaveScreenshot(os.path.join(out_dir, "mesh_root_section.png"), view,
               ImageResolution=[2400, 1600])
print("Saved mesh_root_section.png")
