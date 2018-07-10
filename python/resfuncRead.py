import numpy as np

"""read resolution functions from text files"""
def getRFunc(infile):

  #open the file
  f = open(infile)

  #setup output
  out = {}

  for i,line in enumerate(f.readlines()):
    #print(line)
    det = np.uint32(line.split()[0])
    out[det] = {}
    sqrt_vec = np.ones((6,),dtype=np.float64)
    sqrt_vec[0] = np.float64(line.split()[1])
    sqrt_vec[1] = np.float64(line.split()[2])
    sqrt_vec[2] = np.float64(line.split()[3])
    sqrt_vec[3] = np.float64(line.split()[4])
    sqrt_vec[4] = np.float64(line.split()[5])
    sqrt_vec[5] = np.float64(line.split()[6])
    out[det]['sqrt'] = sqrt_vec 
    lin_vec = np.ones((6,),dtype=np.float64)
    lin_vec[0] = np.float64(line.split()[7])
    lin_vec[1] = np.float64(line.split()[8])
    lin_vec[2] = np.float64(line.split()[9])
    lin_vec[3] = np.float64(line.split()[10])
    out[det]['lin'] = lin_vec 

  return out

def makeFunc(vec,islin = False):

  if not islin:
    f = lambda x: np.sqrt(vec[0] + x*vec[2] + x**2*vec[4]) 
  else:
    f = lambda x: vec[0] + x*vec[2]

  return f
