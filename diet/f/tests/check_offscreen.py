#!/usr/bin/env python3
"""Check the headless VTK path used by MRI figures without cohort inputs."""
try:
    import resource
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
except ImportError:
    pass
import vtk

sphere = vtk.vtkSphereSource()
mapper = vtk.vtkPolyDataMapper()
mapper.SetInputConnection(sphere.GetOutputPort())
actor = vtk.vtkActor()
actor.SetMapper(mapper)
renderer = vtk.vtkRenderer()
renderer.AddActor(actor)
renderer.ResetCamera()
window = vtk.vtkRenderWindow()
window.SetOffScreenRendering(1)
window.SetSize(64, 64)
window.AddRenderer(renderer)
window.Render()
image = vtk.vtkWindowToImageFilter()
image.SetInput(window)
image.Update()
pixels = image.GetOutput().GetPointData().GetScalars()
if pixels is None or pixels.GetRange()[0] == pixels.GetRange()[1]:
    raise RuntimeError("VTK produced an empty image; install a working Mesa/EGL or OSMesa runtime")
window.Finalize()
print("VTK offscreen rendering: PASS", flush=True)
