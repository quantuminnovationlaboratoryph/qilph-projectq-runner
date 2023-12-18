import math
import sys
import numpy as np
from projectq          import MainEngine
from projectq.backends import Simulator
from projectq.ops      import *
from qiskit            import *

#=======================================#

debug = False

gQasmFilename = ""

for parameter in sys.argv:
  if parameter.endswith(".qasm"):
    gQasmFilename = parameter

#=======================================#

gRegisters      = {}
simulator1      = Simulator()
engine1         = MainEngine(backend=simulator1)

gQuantumCircuit = QuantumCircuit.from_qasm_file(gQasmFilename)
gateCounts      = 0
gCircuitQubits  = gQuantumCircuit.num_qubits
gCircuitDepth   = gQuantumCircuit.depth()
gCircuitGates   = gQuantumCircuit.count_ops()

print("")
print("======================================================================")
print("INFO: QUBITS: " + str(gCircuitQubits))
print("INFO: DEPTH : " + str(gCircuitDepth))


print("")
print("======================================================================")
print(gQuantumCircuit)


print("")
print("======================================================================")
print("INFO: QREGS: ", gQuantumCircuit.qregs)

for qreg in gQuantumCircuit.qregs:
  qregSize = qreg.size
  qregName = qreg.name
  gRegisters[qregName] = engine1.allocate_qureg(qregSize)
  print(qreg)

for reg in gRegisters:
  print(reg, gRegisters[reg])


print("")
print("======================================================================")

qasmOneQubitGateLabels = [ 
"u3"  ,
"u"   ,
"u2"  ,
"rx"  ,
"ry"  ,
"u1"  ,
"p"   ,
"rz"  ,
"u0"  ,
"id"  ,
"t"   , 
"tdg" ,
"s"   ,
"sdg" ,
"z"   ,
"x"   ,
"y"   ,
"h"   ,
"sx"  ,
"sxdg",
]

qasmMultiQubitGateLabels = [ 
"cx"  ,
"cu"  ,
]

# This function takes a small threshold number and a complex number cNum = p + qi, where a and b are 
# the two real components of cNum, and return new complex number cNum' = r + si where r is equal to
# p if the absolute value of p's fractional component is greater than the threshold value otherwise
# r is the integral part of p and s is equal to q if the absolute value of q's fractional component
# is greater than the threshold value otherwise s is the integral part of q.

def zeroDownFrac (cNum, threshold):
  realComp = np.modf(np.real(cNum))
  realInt  = realComp[1]
  realFrac = realComp[0]
  imagComp = np.modf(np.imag(cNum))
  imagInt  = imagComp[1]
  imagFrac = imagComp[0]
  realFrac = 0.0 if np.fabs(realFrac) < threshold else realFrac
  imagFrac = 0.0 if np.fabs(imagFrac) < threshold else imagFrac
  return complex(realInt + realFrac, imagInt + imagFrac)

def qasmU3Matrix ( _t, _p, _l, thres):
  _k = _p + _l
  p  = np.cos(_t/2.0)
  q  = np.sin(_t/2.0)
  e0 = 1
  e1 = -(np.cos(_l) + np.sin(_l)*0j)
  e2 =  (np.cos(_p) + np.sin(_p)*0j)
  e3 =  (np.cos(_k) + np.sin(_k)*0j)
  w  = p * e0
  x  = q * e1
  y  = q * e2
  z  = p * e3
  ww = zeroDownFrac(w, thres)
  xx = zeroDownFrac(x, thres)
  yy = zeroDownFrac(y, thres)
  zz = zeroDownFrac(z, thres)
  if debug == True:
    print("theta   =", _t)
    print("phi     =", _p)
    print("lambda  =", _l)
    print("p (cos) = ", p)
    print("q (sin) = ", q)
    print("e0      = ", e0)
    print("e1      = ", e1)
    print("e2      = ", e2)
    print("e3      = ", e3)
    print("w       = ", w)
    print("x       = ", x)
    print("y       = ", y)
    print("z       = ", z)
    print("ww      = ", ww)
    print("xx      = ", xx)
    print("yy      = ", yy)
    print("zz      = ", zz)
  return np.array([[ww, xx],[yy, zz]])


thres = 0.0000000001
#print(H.matrix)
#print(myGate.matrix)
#print("Hadmard Matrix from U3:\n", qasmU3Matrix(math.pi/2.0, 0, math.pi, thres), "\n")
#print("NOT Matrix from U3:\n", qasmU3Matrix(math.pi, 0, math.pi, thres), "\n")

idGate = MatrixGate([[1,0],
                     [0,1]])

SqrtXdag = get_inverse(SqrtX)

# Mapping of QASM gate labels to ProjectQ gate functions
qasmGateToProjectQ = {
# qubits: 1, parameters: 3

"u3"      :   MatrixGate,  # CASE: A
"u"       :   MatrixGate,  # CASE: A

# qubits: 1, parameters: 2

"u2"      :   MatrixGate,  # CASE: A

# qubits: 1, parameters: 1

"rx"      :           Rx,
"ry"      :           Ry,
"u1"      :   MatrixGate,
"p"       :           Ph,
"rz"      :           Rz,

# qubits: 1, parameters: 0

"u0"      :       idGate,
"id"      :       idGate,
"t"       :            T, 
"tdg"     :         Tdag, 
"s"       :            S, 
"sdg"     :         Sdag,
"z"       :            Z, 
"x"       :            X, 
"y"       :            Y, 
"h"       :            H, 
"sx"      :        SqrtX,
"sxdg"    :     SqrtXdag,
"barrier" :      Barrier,
"measure" :      Measure,

# qubits: 2, parameters: 0

"cx"      :         C(X),
"cz"      :         C(Z),
"cy"      :         C(Y),
"ch"      :         C(H),
"csx"     :     C(SqrtX), 
"swap"    :         Swap,

# qubits: 2, parameters: 1

"crx"     :           Rx,
"cry"     :           Ry,
"crz"     :           Rz,
"cp"      :           Ph, 
"cu1"     :           Ph,
"rzz"     :          Rzz, 
"rxx"     :          Rxx,

# qubits: 2, parameters: 3

"cu3"     :   MatrixGate,  # CASE: A

# qubits: 2, parameters: 4

"cu"      :   MatrixGate,  # CASE: A

# qubits: 3, parameters: 0

"cswap"   :      C(Swap),
"ccx"     :      Toffoli,
"rccx"    :  C(idGate,2),

# qubits: 4, parameters: 0

"rc3x"    :  C(idGate,3), 
"rcccx"   :  C(idGate,3), 
"c3x"     :       C(X,3), 
"c3sx"    :   C(SqrtX,3),

# qubits: 5, parameters: 0

"c4x"     :       C(X,4),

}



#gate = MatrixGate([[0, 1], [1, 0]])

for gate in gQuantumCircuit.data:
  gateName      = gate.operation.name.lower()
  gateParams    = gate.operation.params
  gateNumQubits = gate.operation.num_qubits
  gateNumClbits = gate.operation.num_clbits
  gateQubits    = gate.qubits

  # For Debugging
  if gateName in ["barrier", "measure"]:
    continue

  print("\n")
  print(gate)
  print("Gate Name      : " + gateName)
  print("Gate Parameters: " + str(gateParams))
  print("No. of Qubits  : " + str(gateNumQubits))
  print("No. of Clbits  : " + str(gateNumClbits))

 
  engine1.flush()
   
  #================================================#
  # Covert the list of qubits to a tuple of qubits #
  #================================================#
  qubitsTuple = ()
  for aQubit in gateQubits:
    qubitRegName = gQuantumCircuit.find_bit(aQubit)[1][0][0].name
    qubitIndex   = gQuantumCircuit.find_bit(aQubit).registers[0][1]
    theQubit     = gRegisters[qubitRegName][qubitIndex]
    qubitsTuple  = qubitsTuple + (theQubit,)
    if debug==True:
      print("aQubit         :", aQubit)
      print("Qubit Reg Name : " + qubitRegName)
      print("Qubit Index    : " + str(qubitIndex))
      print("theQubit       :", theQubit)
      print("theQubit(type) :", type(theQubit))

  qubitScope  = len(qubitsTuple)
  paramsCount = len(gateParams)

  print("No. of Params  : " + str(paramsCount))



  #================================================#
  # CASE A: Gates that need a matrix as parameter  #
  #================================================#

  if gateName in ["u3", "u2", "u1", "u", "cu3", "cu"]:
    print("CASE           : A")
    if   len(gateParams) >= 3:
       gateParams2 = (gateParams[0], gateParams[1], gateParams[2], thres)
    elif len(gateParams) == 2:
       gateParams2 = (  math.pi/2.0, gateParams[0], gateParams[1], thres)
    else:
       gateParams2 = ( 0000000000.0, 00000000000.0, gateParams[0], thres)
    matrixParam = qasmU3Matrix(*gateParams2)
    
    if "c" in gateName:
      projQGate   = C(qasmGateToProjectQ[gateName](matrixParam))
    else:
      projQGate   = qasmGateToProjectQ[gateName](matrixParam)
    
    projQGate  | qubitsTuple
    continue

  #================================================#
  # CASE B: Gates with 'no' parameters             #
  #================================================#
  if   paramsCount == 0 or gateName == "u0":
    print("CASE           : B")
    print("gateName       : ", "["+gateName+"]")
    if gateName == "mcx":
      gateName = "c3x" if len(qubitsTuple) == 4  else "c4x"
    projQGate  = qasmGateToProjectQ[gateName]
    projQGate  | qubitsTuple
    continue

  
  #================================================#
  # CASE C: Gates with a single 'angle' parameter  #
  #================================================#
  if paramsCount == 1:
    print("CASE           : C")
    angleParam = gateParams[0] 
    if gateName in ["rxx", "rzz"] + ["rx", "ry", "rz", "p"]:
      projQGate  = qasmGateToProjectQ[gateName](angleParam)
    else:
      projQGate  = C(qasmGateToProjectQ[gateName](angleParam))
    projQGate  | qubitsTuple
    continue

'''
  #==================================================#
  # 1-Qubit Gate Translation                         #
  #==================================================#  
  if gateName in qasmOneQubitGateLabels:

    qubit        = gateQubits[0] 
    qubitRegName = gQuantumCircuit.find_bit(qubit)[1][0][0].name
    qubitIndex   = gQuantumCircuit.find_bit(qubit).registers[0][1]
    targetQubit  = gRegisters[qubitRegName][qubitIndex]

    print("Qubit Reg Name : " + qubitRegName)
    print("Qubit Index    : " + str(qubitIndex))

    # Case A: You need to define a matrix to accomodate the u-gate parameters.
    if gateName in ["u3", "u", "u2", "u1"]:

      if   len(gateParams) == 3:
         gateParams2 = (gateParams[0], gateParams[1], gateParams[2], thres)
      elif len(gateParams) == 2:
         gateParams2 = (  math.pi/2.0, gateParams[0], gateParams[1], thres)
      else:
         gateParams2 = ( 0000000000.0, 00000000000.0, gateParams[0], thres)

      matrixParam = qasmU3Matrix(*gateParams2)
      projQGate   = qasmGateToProjectQ[gateName](matrixParam)

    # CASE B: You need to provide an angle parameter to the gate.
    elif gateName in ["rx", "ry", "rz", "p"]:
      projQGate   = qasmGateToProjectQ[gateName](gateParams[0])

    # CASE C: No parameters. Apply the gate as is to the qubit.
    else:
      projQGate   = qasmGateToProjectQ[gateName]

    projQGate | targetQubit

  #==================================================#
  # Multi-Qubit Gate Translation                     #
  #==================================================#  
  else:
    print("MultiQubit Gate!")
    print("Qubits          :", gateQubits)

    qubitsTuple = ()
    for aQubit in gateQubits:
      qubitRegName = gQuantumCircuit.find_bit(aQubit)[1][0][0].name
      qubitIndex   = gQuantumCircuit.find_bit(aQubit).registers[0][1]
      theQubit     = gRegisters[qubitRegName][qubitIndex]
      qubitsTuple  = qubitsTuple + (theQubit,)
      if debug==True:
        print("aQubit         :", aQubit)
        print("Qubit Reg Name : " + qubitRegName)
        print("Qubit Index    : " + str(qubitIndex))
        print("theQubit       :", theQubit)
        print("theQubit(type) :", type(theQubit))

    if debug==True:
      print("qubitsTuple    :", qubitsTuple)
      for pqQubit in qubitsTuple:
        print("pqQubit        :", pqQubit)
        print("pqQubit(type)  :", type(pqQubit))

    paramsCount = len(gateParams)

    # CASE: Gates with no parameters
    if   paramsCount == 0:
      print("paramsCount =", paramsCount)
      if gateName == "mcx":
        gateName = "c3x" if len(qubitsTuple) == 4  else "c4x"
      projQGate  = qasmGateToProjectQ[gateName]
      projQGate  | qubitsTuple
    
    # CASE: Gates with 1 parameter (angle) 
    elif paramsCount == 1:
      print("paramsCount =", paramsCount)
      angleParam = gateParams[0] 
      if gateName in ["rxx", "rzz"]:
        projQGate  = qasmGateToProjectQ[gateName](angleParam)
      else:
        projQGate  = C(qasmGateToProjectQ[gateName](angleParam))
      projQGate  | qubitsTuple
    else:
      continue
'''    
     
      


print("")
print("======================================================================")

for regGroup in gRegisters:
  print("regGroup       : ", regGroup)
  print("register       : ", gRegisters[regGroup])
  All(Measure) | gRegisters[regGroup]
engine1.flush()


observedString = ""                                                                                     
for regGroup in gRegisters:                                                                             
  print("")                                                                                             
  print("REGISTER GROUP: ", regGroup)                                                                   
  for index in range(0, len(gRegisters[regGroup])):                                                     
    observedString += str(int(gRegisters[regGroup][index]))                                             
    print(f"Index = {index}, value = {int(gRegisters[regGroup][index])}")                            
print("Observed String: ", observedString)                                                           









