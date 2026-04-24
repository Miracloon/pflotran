import vtk

executive = self.GetExecutive()
outInfo = executive.GetOutputInformation(0)

t = 0.0
if outInfo.Has(vtk.vtkStreamingDemandDrivenPipeline.UPDATE_TIME_STEP()):
    t = outInfo.Get(vtk.vtkStreamingDemandDrivenPipeline.UPDATE_TIME_STEP())

# Get the input — handle both vtkTable and vtkMultiBlockDataSet
raw_input = self.GetInputDataObject(0, 0)
if raw_input.IsA('vtkTable'):
    input_table = raw_input
elif raw_input.IsA('vtkMultiBlockDataSet'):
    # Walk blocks to find the first vtkTable
    input_table = None
    for i in range(raw_input.GetNumberOfBlocks()):
        blk = raw_input.GetBlock(i)
        if blk is not None and blk.IsA('vtkTable'):
            input_table = blk
            break
    if input_table is None:
        raise RuntimeError("No vtkTable found in MultiBlockDataSet input. "
                           "Connect this filter to a table source (e.g. the .tsv file), "
                           "not the 3D mesh data.")
else:
    raise RuntimeError(f"Unsupported input type: {raw_input.GetClassName()}. "
                       "Connect this filter to a table source (e.g. the .tsv file).")

time_col_name = "Time [y]"
global_min =  float('inf')
global_max = -float('inf')

for i in range(input_table.GetNumberOfColumns()):
    col = input_table.GetColumn(i)
    if col.GetName() == time_col_name:
        continue
    lo, hi = col.GetRange()
    if lo < global_min:
        global_min = lo
    if hi > global_max:
        global_max = hi

output = self.GetOutput()
output.Initialize()

tc = vtk.vtkDoubleArray()
tc.SetName("Time")
tc.InsertNextValue(t)
tc.InsertNextValue(t)

mc = vtk.vtkDoubleArray()
mc.SetName("Marker")
mc.InsertNextValue(global_min)
mc.InsertNextValue(global_max)

output.AddColumn(tc)
output.AddColumn(mc)
output.GetInformation().Set(vtk.vtkDataObject.DATA_TIME_STEP(), t)