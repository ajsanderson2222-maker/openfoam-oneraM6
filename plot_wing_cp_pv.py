"""
Render top-down view of wing surface colored by Cp using paraview offscreen rendering.
"""
from paraview.simple import *
import os

case_dir = "/home/andrew/openfoam_jobs/openfoam-oneraM6/simulation"
out_dir  = "/home/andrew/openfoam_jobs/openfoam-oneraM6"

p_inf = 80510.081
rho_inf = 80510.081 / (287.0 * 255.556)
U_inf_mag = (268.654**2 + 14.362**2)**0.5
q_inf = 0.5 * rho_inf * U_inf_mag**2

foam = OpenFOAMReader(registrationName="case.foam",
                      FileName=os.path.join(case_dir, "case.foam"))
foam.MeshRegions = ["patch/wing"]
foam.CellArrays  = ["p"]
foam.UpdatePipeline(time=2000)

# Extract the wing block by name
extract = ExtractBlock(Input=foam)
extract.Selectors = ["/Root/boundary/wing"]
extract.UpdatePipeline()

view = CreateRenderView()
view.ViewSize               = [2400, 1600]
view.Background             = [1, 1, 1]
view.UseLight               = 0
view.OrientationAxesVisibility = 0

disp = Show(extract, view)
disp.Representation = "Surface"
ColorBy(disp, ("CELLS", "p"))

lut = GetColorTransferFunction("p")
lut.RescaleTransferFunction(28920, 121154)
lut.ApplyPreset("Cool to Warm", True)

disp.SetScalarBarVisibility(view, True)
sb = GetScalarBar(lut, view)
sb.Title = "p (Pa)"
sb.ComponentTitle = ""
sb.Orientation = "Vertical"
sb.Position = [0.88, 0.05]
sb.ScalarBarLength = 0.7

Render(view)

# Top-down: looking down -y at upper surface
# Wing x=0..1.14, z=0..1.22; camera above in y
view.CameraParallelProjection = 1
view.CameraPosition      = [0.57, 5.0, 0.61]
view.CameraFocalPoint    = [0.57, 0.0, 0.61]
view.CameraViewUp        = [-1,   0,   0   ]   # chord runs left→right
view.CameraParallelScale = 0.7   # half-height in world units

Render(view)
SaveScreenshot(os.path.join(out_dir, "wing_cp_topdown.png"), view,
               ImageResolution=[2400, 1600])
print("Saved wing_cp_topdown.png")
