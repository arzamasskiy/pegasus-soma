"""
Read Pegasus .vtk output files

Last update  2016-02-26  leva@astro.princetone.edu
"""

# Python modules
import numpy as np
from vtk import vtkStructuredPointsReader
from vtk.util import numpy_support as VN

#=======================================================================================

def vtk_read(filename):
  reader = vtkStructuredPointsReader()
  reader.SetFileName(filename)
  reader.ReadAllVectorsOn()
  reader.ReadAllScalarsOn()
  numS = reader.GetNumberOfScalarsInFile()
  numV = reader.GetNumberOfVectorsInFile()
  reader.Update()
	
  data = reader.GetOutput()
  dim = data.GetDimensions()
  vec = list(dim)
  for i in range(len(dim)):
    if vec[i]>1:
      vec[i] = vec[i]-1

  result = dict()

  for i in range(0, numS):
    name = reader.GetScalarsNameInFile(i)
    array = VN.vtk_to_numpy(data.GetCellData().GetArray(name))
    result.update([(name, array.reshape(vec, order='F'))])
  vec.append(3)
  for i in range(0, numV):
    name = reader.GetVectorsNameInFile(i)
    array = VN.vtk_to_numpy(data.GetCellData().GetArray(name))
    result.update([(name, array.reshape(vec, order='F'))])
	
  x = np.zeros(data.GetNumberOfPoints())
  y = np.zeros(data.GetNumberOfPoints())
  z = np.zeros(data.GetNumberOfPoints())

  for i in range(data.GetNumberOfPoints()):
    x[i],y[i],z[i] = data.GetPoint(i)

  x = x.reshape(dim,order='F')
  y = y.reshape(dim,order='F') 
  z = z.reshape(dim,order='F')

  result.update([('x', x)])
  result.update([('y', y)])
  result.update([('z', z)])

  return result

class PegasusError(RuntimeError):
  """General exception class for Pegasus read functions."""
  pass
