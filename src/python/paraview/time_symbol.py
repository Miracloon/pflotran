# Time Symbol Markers for ParaView Line Chart View
# ==================================================
# Draws a moving symbol on each plotted series at the current animation
# timestep. The symbol tracks the interpolated value of each series as the
# animation plays, providing a dot-on-line indicator of the current time.
#
# REQUIREMENTS
# ------------
# Three data sources must be loaded in ParaView:
#   1. The PFLOTRAN HDF5 or VTK output file (transient mesh data). This
#      provides the TIME_STEPS that drive animation.
#   2. The .tsv file produced by dat_to_tsv.py (the mass-balance table).
#      This is displayed in the Line Chart View as the plotted series.
#   3. This Programmable Filter, connected to BOTH sources above, which
#      outputs the interpolated symbol positions at the current timestep.
#
# NOTE: The .tsv file loaded via the CSV reader is a static dataset with no
# TIME_STEPS metadata. The filter must be connected to the PFLOTRAN transient
# source to inherit timesteps — connecting to the .tsv alone will produce a
# symbol that never moves.
#
# SETUP
# -----
# 1. File → Open: load the PFLOTRAN HDF5/VTK transient output.
#
# 2. File → Open: load the .tsv file (produced by dat_to_tsv.py).
#    ParaView will read it with the CSV reader as a static table.
#
# 3. In the Pipeline Browser, Ctrl+click to select BOTH the PFLOTRAN
#    transient source AND the .tsv source simultaneously.
#    Apply: Filters → Programmable Filter
#    Set the output data set type to "vtkTable".
#    The filter will receive:
#      Input 0: PFLOTRAN transient source  (timesteps)
#      Input 1: .tsv table source          (series data)
#
# 4. Load this script into the "Script" box of the Programmable Filter using
#    one of two methods:
#
#    a) Paste: copy and paste the contents of this file directly into the
#       Script box.
#
#    b) Exec (recommended): put only this single line in the Script box:
#           exec(open('/full/path/to/time_symbol.py').read())
#       ParaView will re-read the file each time the filter is evaluated,
#       so edits to the file take effect after clicking Apply — no need to
#       reopen the Script dialog.
#
# 5. Set time_col_name (below) to match the time column in your .tsv file
#    (default: "Time [y]").
#
# 6. Click Apply.
#
# 7. In the Pipeline Browser, HIDE the Programmable Filter in the 3D Render
#    View by clicking its eye icon — a vtkTable cannot be rendered in 3D
#    and will produce errors if left visible there.
#
# 8. Open (or activate) a Line Chart View.
#    - Show the .tsv source there and configure the series to plot.
#    - Also show the Programmable Filter there. In its Display tab:
#        - Set X Array Name to: Time
#        - For each series, set Line Style to "None" and Marker Style to a
#          symbol (circle, square, diamond, etc.) with a larger marker size.
#        - Match the symbol color to the corresponding line from the .tsv.
#
# 9. To hide symbols from the chart legend:
#    In the Display tab, find each series row in the series table and uncheck
#    the legend-visibility column for that row.
#
# NOTES
# -----
# - Values are linearly interpolated between timesteps in the .tsv table.
# - All numeric columns in the .tsv are included as symbol series. To restrict
#   to specific columns, replace the loop below with a list of column names.
# - The output is a single-row vtkTable — one interpolated value per series
#   at the current time t.

import vtk
import numpy as np

executive = self.GetExecutive()
outInfo = executive.GetOutputInformation(0)

t = 0.0
if outInfo.Has(vtk.vtkStreamingDemandDrivenPipeline.UPDATE_TIME_STEP()):
    t = outInfo.Get(vtk.vtkStreamingDemandDrivenPipeline.UPDATE_TIME_STEP())

# Input 0: PFLOTRAN transient source (provides timesteps — not used for data)
# Input 1: .tsv table source (provides series data for interpolation)
table = self.GetInputDataObject(0, 1)
if table is None or not table.IsA('vtkTable'):
    raise RuntimeError(
        "Second input must be the .tsv table source. "
        "Select both the PFLOTRAN transient source and the .tsv source "
        "before applying the Programmable Filter."
    )

# Name of the time column in the .tsv table
time_col_name = "Time [y]"

time_col_idx = None
for i in range(table.GetNumberOfColumns()):
    if table.GetColumn(i).GetName() == time_col_name:
        time_col_idx = i
        break
if time_col_idx is None:
    raise RuntimeError(f"Column '{time_col_name}' not found in the table input.")

n = table.GetNumberOfRows()
times = np.array([table.GetColumn(time_col_idx).GetValue(r) for r in range(n)])

output = self.GetOutput()
output.Initialize()

# Time column for the x-axis
tc = vtk.vtkDoubleArray()
tc.SetName("Time")
tc.InsertNextValue(t)
output.AddColumn(tc)

# One interpolated value per series at time t
for i in range(table.GetNumberOfColumns()):
    col = table.GetColumn(i)
    name = col.GetName()
    if name == time_col_name:
        continue
    values = np.array([col.GetValue(r) for r in range(n)])
    vc = vtk.vtkDoubleArray()
    vc.SetName(name)
    vc.InsertNextValue(float(np.interp(t, times, values)))
    output.AddColumn(vc)

output.GetInformation().Set(vtk.vtkDataObject.DATA_TIME_STEP(), t)
